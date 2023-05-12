import configparser
from datetime import datetime
import os

import hoomd, hoomd.md
import gsd, gsd.hoomd, gsd.pygsd
import numpy as np

from params import BOND_PARAMS, MD_PARAMS_TYPES
import molsys
import timefmt


class MDSystem(molsys.MolSystem):
    def __init__(self):
        super().__init__()

    def load_mdparams(self, inipath: str):
        if not os.path.isfile(inipath):
            raise OSError(f"Initial configuration file {inipath} can not be found.")
        CONFIG = configparser.ConfigParser()
        CONFIG.read(inipath, encoding='utf-8')
        print("-*-*- SETTINGS FOR SIMULATION -*-*-")
        for section, args in MD_PARAMS_TYPES.items():
            for nm, dtype in args.items():
                fullnm = nm if section == "main" else section + '_' + nm
                if dtype == int:
                    v = CONFIG.getint(section, nm)
                elif dtype == float:
                    v = CONFIG.getfloat(section, nm)
                elif dtype == bool:
                    v = CONFIG.getboolean(section, nm)
                else:
                    v = CONFIG.get(section, nm)
                setattr(self, fullnm, v)
                print(f"Option {fullnm} = {v}")
    
    def configure_box(self):
        # No negative coordinates
        sxyz = self.P.position
        sxyz -= (sxyz.max(axis=0) + sxyz.min(axis=0)) * .5
        # No particles out of box
        xmax, ymax, zmax = np.max([np.max(sxyz, axis=0), np.abs(np.min(sxyz, axis=0))], axis=0)
        if (self.autobox):
            self.xbox = xmax * 2.1
            self.ybox = ymax * 2.1
            self.zbox = zmax * 2.1
        if (xmax > self.xbox) | (ymax > self.ybox) | (zmax > self.zbox):
            raise ValueError(f"Particles out of box. Box length should be larger than ({xmax}, {ymax}, {zmax}).")

        C = self.S.configuration
        C.dimensions = 3
        C.box = [self.xbox, self.ybox, self.zbox, 0., 0., 0.]
        C.step = 0
        # Write initial configuration to file
        with gsd.hoomd.open(name=self.initgsd, mode='wb') as f:
            f.append(self.S)
    
    def setup(self):
        # 1 Init
        # Run a simulation using prepared initial configuration file
        device = hoomd.device.GPU() if self.gpu else hoomd.device.CPU()
        if self.seed < 0:
            self.seed = np.random.randint(0, 1 << 16)
        self.system = hoomd.Simulation(device=device, seed=self.seed)
        if not os.path.isfile(self.initgsd):
            raise OSError(f"Initial frame file {self.initgsd} can not be found.")
        self.system.create_state_from_gsd(self.initgsd)

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
        kT = self.T * .00831446
        methodkw = self.integrator.upper()
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
                                       kT=kT, tau=self.tau)
            methods.append(nvt)            
        elif methodkw == "NPT":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v3.10.0/module-md-methods.html#hoomd.md.methods.NPT
            npt = hoomd.md.methods.NPT(filter=hoomd.filter.All(),
                                       kT=kT, tau=self.tau, S=self.p, tauS=self.taup, couple="none")
            methods.append(npt)
        else:
            raise ValueError("Unrecognized md integrator type. Only LD/NVE/NVT/NPT are allowed.")
        integrator = hoomd.md.Integrator(dt=self.dt, methods=methods, forces=[harmonic, ashbaugh, yukawa])
        
        # 4 Logging
        # Output system trajectory data and properties recorded in the GSD format
        logger = hoomd.logging.Logger()
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
        self.system.operations.computes.append(thermodynamic_properties)
        logger.add(thermodynamic_properties)
        logger.add(self.system, quantities=["timestep", "walltime"])
        traj_writer = hoomd.write.GSD(filename=self.output_traj,
                                     trigger=hoomd.trigger.Periodic(self.output_period),
                                     mode='xb')
        traj_writer.writer = logger
        self.system.operations.writers.append(traj_writer)

        # Output system particle positions and the box parameters in the DCD format
        hoomd.write.DCD(filename=self.output_checkpoint,
                        trigger=hoomd.trigger.Periodic(self.output_cptperiod))

    def run(self):
        # Run simulation
        stime = datetime.now()
        print("Simulation starts at " + datetime.strftime(stime, '%c'))
        self.system.run(steps=self.steps, write_at_start=self.newrun)
        etime = datetime.now()
        print("Simulation ends at " + datetime.strftime(etime, '%c'))
        print("Total usage time:", timefmt.strfdelta(etime - stime, "%D days %H hours %M minutes %S seconds"))