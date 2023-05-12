import configparser
from datetime import datetime
import os

import hoomd, hoomd.md
import gsd, gsd.hoomd, gsd.pygsd
import numpy as np

from params import BOND_PARAMS, MD_PARAMS_TYPES, OUTPUT_PARAMS_TYPES, MDParams
import molsys
import timefmt


class MDSystem(molsys.MolSystem):
    def __init__(self):
        super().__init__()
        self.mdpara = MDParams("main", MD_PARAMS_TYPES)
        self.output = MDParams("output", OUTPUT_PARAMS_TYPES)

    def load_mdparams(self, inipath: str):
        if not os.path.isfile(inipath):
            raise OSError(f"Initial configuration file {inipath} can not be found.")
        CONFIG = configparser.ConfigParser()
        CONFIG.read(inipath, encoding='utf-8')
        print("-*-*- SETTINGS FOR SIMULATION -*-*-")
        self.mdpara.digest(CONFIG)
        print(self.mdpara.params)
        self.output.digest(CONFIG)
        print(self.output.params)
    
    def configure_box(self):
        # No negative coordinates
        sxyz = self.P.position
        sxyz -= (sxyz.max(axis=0) + sxyz.min(axis=0)) * .5
        # No particles out of box
        xmax, ymax, zmax = np.max([np.max(sxyz, axis=0), np.abs(np.min(sxyz, axis=0))], axis=0)
        if (self.mdpara.autobox):
            bx = xmax * 2.1
            by = ymax * 2.1
            bz = zmax * 2.1
        else:
            bx = self.mdpara.xbox
            by = self.mdpara.ybox
            bz = self.mdpara.zbox
        if (xmax > bx) | (ymax > by) | (zmax > bz):
            raise ValueError(f"Particles out of box. Box length should be larger than ({xmax}, {ymax}, {zmax}).")

        C = self.S.configuration
        C.dimensions = 3
        C.box = [bx, by, bz, 0., 0., 0.]
        C.step = 0
        # Write initial configuration to file
        with gsd.hoomd.open(name=self.mdpara.initgsd, mode='wb') as f:
            f.append(self.S)
    
    def setup(self):
        # 1 Init
        # Run a simulation using prepared initial configuration file
        device = hoomd.device.GPU() if self.mdpara.gpu else hoomd.device.CPU()
        if self.mdpara.seed < 0:
            self.mdpara.seed = np.random.randint(0, 1 << 16)
        self.system = hoomd.Simulation(device=device, seed=self.mdpara.seed)
        if not os.path.isfile(self.mdpara.initgsd):
            raise OSError(f"Can not find initial frame file {self.mdpara.initgsd}.")
        self.system.create_state_from_gsd(self.mdpara.initgsd)

        # 2 Bonds
        # Set bonds
        harmonic = hoomd.md.bond.Harmonic()
        for btype in self.B.types:
            if btype in BOND_PARAMS:
                bparams = BOND_PARAMS[btype]
                harmonic.params[btype] = dict(k=bparams.k, r0=bparams.r0)
            else:
                raise ValueError("Unexpected bond type name. "
                                 f"Please check whether {btype} is listed in bondparams.")
        
        # Set nonbonded
        nl = hoomd.md.nlist.Cell(buffer=1.5, exclusions=["bond", "body"])   # Neighbor list
        ashbaugh = hoomd.md.pair.Ashbaugh(default_r_cut=3.5, nlist=nl)   # Ash-baugh Van der Waals potential
        yukawa = hoomd.md.pair.Yukawa(default_r_cut=3.5, nlist=nl)
        
        # Calculate pairwise coefficient matrix
        N = len(self.rtypes)
        lambda_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rlambda, dtype=float)
        lambda_mat = (lambda_mat.transpose() + lambda_mat) * .5
        sigma_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rsigma, dtype=float)
        sigma_mat = (sigma_mat.transpose() + sigma_mat) * .05
        epsilon_mat = np.ones((N, N), dtype=float) * np.asarray(self.rcharges, dtype=float)
        epsilon_mat = (epsilon_mat.transpose() * epsilon_mat) * 1.73136

        for i, ri in enumerate(self.rtypes):
            for j, rj in enumerate(self.rtypes):
                ashbaugh.params[(ri, rj)] = dict(sigma=sigma_mat[i, j], epsilon=.8368, lam=lambda_mat[i, j], alpha=1.0)
                yukawa.params[(ri, rj)] = dict(epsilon=epsilon_mat[i, j], kappa=1.)
        
        # 3 Setup MD
        # Setup integrator
        methods = []
        kT = self.mdpara.T * .00831446
        methodkw = self.mdpara.integrator.upper()
        if methodkw == "LD":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v3.10.0/module-md-methods.html#hoomd.md.methods.Langevin
            ld = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT)
            for ri in self.rtypes:
                ld.gamma[ri] = self.rparams[ri].mass * .001
            methods.append(ld)
        elif methodkw == "NVE":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v3.10.0/module-md-methods.html#hoomd.md.methods.NVE
            nve = hoomd.md.methods.NVE(filter=hoomd.filter.All())
            methods.append(nve)
        elif methodkw == "NVT":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v3.10.0/module-md-methods.html#hoomd.md.methods.NVT
            nvt = hoomd.md.methods.NVT(filter=hoomd.filter.All(),
                                       kT=kT, tau=self.mdpara.tau)
            methods.append(nvt)            
        elif methodkw == "NPT":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v3.10.0/module-md-methods.html#hoomd.md.methods.NPT
            npt = hoomd.md.methods.NPT(filter=hoomd.filter.All(),
                                       kT=kT, tau=self.mdpara.tau, S=self.mdpara.p, tauS=self.mdpara.taup, couple="none")
            methods.append(npt)
        else:
            raise ValueError("Unrecognized md integrator type. Only LD/NVE/NVT/NPT are allowed.")
        
        if self.mdpara.em:
            integrator = hoomd.md.minimize.FIRE(
                dt=self.mdpara.dt,
                force_tol=self.mdpara.forcetol,
                angmom_tol=self.mdpara.angmomtol,
                energy_tol=self.mdpara.energytol,
                methods=methods,
                forces=[harmonic, ashbaugh, yukawa])
        else:
            integrator = hoomd.md.Integrator(
                dt=self.mdpara.dt,
                methods=methods,
                forces=[harmonic, ashbaugh, yukawa])
        self.system.operations.integrator = integrator

        # 4 Logging
        # Output system trajectory data and properties recorded in the GSD format
        logger = hoomd.logging.Logger()
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
        self.system.operations.computes.append(thermodynamic_properties)
        logger.add(thermodynamic_properties)
        logger.add(self.system, quantities=["timestep", "walltime"])
        traj_writer = hoomd.write.GSD(filename=self.output.traj,
                                     trigger=hoomd.trigger.Periodic(self.output.period),
                                     mode='xb')
        traj_writer.writer = logger
        self.system.operations.writers.append(traj_writer)

        # Output system particle positions and the box parameters in the DCD format
        hoomd.write.DCD(filename=self.output.checkpoint,
                        trigger=hoomd.trigger.Periodic(self.output.cptperiod))

    def run(self):
        # Run simulation
        stime = datetime.now()
        print("Simulation starts at " + datetime.strftime(stime, '%c'))
        self.system.run(steps=self.mdpara.steps, write_at_start=self.mdpara.newrun)
        etime = datetime.now()
        print("Simulation ends at " + datetime.strftime(etime, '%c'))
        print("Total usage time:", timefmt.strfdelta(etime - stime, "%D days %H hours %M minutes %S seconds"))
