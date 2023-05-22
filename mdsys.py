import configparser
from datetime import datetime
import os

import hoomd, hoomd.md
import gsd, gsd.hoomd, gsd.pygsd
import numpy as np

from params import BOND_PARAMS, MDParams
import molsys
import timefmt


class MDSystem(molsys.MolSystem):
    def __init__(self, version: str):
        super().__init__()
        # Fetch HOOMD-Blue version
        self.version_major = int(version.split('.', maxsplit=1)[0])
        self.version_minor = int(version.split('.', maxsplit=2)[1])
        # Setup parameter types dictionary
        MD_PARAMS_TYPES = dict(
            gpu=dict(vtype=bool, defv=True),
            seed=dict(vtype=int),
            initgsd=dict(vtype=str, defv="md_init.gsd"),
            autobox=dict(vtype=bool, defv=True),
            automargin=dict(vtype=float, defv=0.1, required="autobox", expect=True),
            xbox=dict(vtype=float, required="autobox", expect=False),
            ybox=dict(vtype=float, required="autobox", expect=False),
            zbox=dict(vtype=float, required="autobox", expect=False),
            cutoff=dict(vtype=float, defv=3.5),
            em=dict(vtype=bool, defv=False),
            forcetol=dict(vtype=float, required="em", expect=True),
            angmomtol=dict(vtype=float, required="em", expect=True),
            energytol=dict(vtype=float, required="em", expect=True),
            integrator=dict(vtype=str),
            steps=dict(vtype=int),
            newrun=dict(vtype=bool, defv=True),
            dt=dict(vtype=float),
            T=dict(vtype=float, required="integrator", expect=["ld", "nvt", "npt"]),
            tau=dict(vtype=float, required="integrator", expect=["nvt", "npt"]),
            p=dict(vtype=float, required="integrator", expect="npt"),
            taup=dict(vtype=float, required="integrator", expect="npt"))
        OUTPUT_PARAMS_TYPES = dict(
            period=dict(vtype=int),
            traj=dict(vtype=str, defv="md_traj.gsd"),
            writecpt=dict(vtype=bool, defv=True),
            cptperiod=dict(vtype=int, required="writecpt", expect=True),
            checkpoint=dict(vtype=str, defv="md_cpt.dcd", required="writecpt", expect=True))
        MD_PARAMS_TYPES
        if self.version_major == 2:
            for kw in ["em", "forcetol", "angmomtol", "energytol", "newrun"]:
                MD_PARAMS_TYPES.pop(kw)
            MD_PARAMS_TYPES.update(tlimit=dict(vtype=int, defv=240))
            OUTPUT_PARAMS_TYPES.update(log=dict(vtype=str, defv="md_run.log"))
        self.mdpara = MDParams("main", MD_PARAMS_TYPES)
        self.output = MDParams("output", OUTPUT_PARAMS_TYPES)

    def load_mdparams(self, inipath: str):
        if not os.path.isfile(inipath):
            raise OSError(f"Initial configuration file {inipath} can not be found.")
        CONFIG = configparser.ConfigParser()
        CONFIG.read(inipath, encoding='utf-8')
        print("-*-*- SETTINGS FOR SIMULATION -*-*-")
        self.mdpara.digest(CONFIG)
        self.output.digest(CONFIG)
    
    def configure_box(self):
        # No negative coordinates
        sxyz = self.P.position
        sxyz -= (sxyz.max(axis=0) + sxyz.min(axis=0)) * .5
        # No particles out of box
        xmax, ymax, zmax = np.max([np.max(sxyz, axis=0), np.abs(np.min(sxyz, axis=0))], axis=0)
        if (self.mdpara.autobox):
            bx = 2. * (1. + self.mdpara.automargin) * (xmax + self.mdpara.cutoff)
            by = 2. * (1. + self.mdpara.automargin) * (ymax + self.mdpara.cutoff)
            bz = 2. * (1. + self.mdpara.automargin) * (zmax + self.mdpara.cutoff)
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
        if self.version_major == 2:
            self.setup_hoomd_v2()
        elif self.version_major == 3:
            self.setup_hoomd_v3()
        else:
            raise ValueError("Unsupported HOOMD-Blue version. Failed to setup simulation details.")

    def setup_hoomd_v2(self):
        import azplugins
        # 1 Init
        inargs = '--mode=cpu' if self.mdpara.gpu else '--mode=cpu'
        hoomd.context.initialize(args=inargs)
        # Run a simulation using prepared initial configuration file
        if not os.path.isfile(self.mdpara.initgsd):
            raise OSError(f"Can not find initial frame file {self.mdpara.initgsd}.")
        hoomd.init.read_gsd(self.mdpara.initgsd)

        # 2 Bonds
        # Set bonds
        harmonic = hoomd.md.bond.harmonic()
        for btype in self.B.types:
            if btype in BOND_PARAMS:
                bparams = BOND_PARAMS[btype]
                harmonic.bond_coeff.set(btype, k=bparams.k, r0=bparams.r0)
            else:
                raise ValueError("Unexpected bond type name. "
                                 f"Please check whether {btype} is listed in bondparams.")
        
        # Set nonbonded
        nl = hoomd.md.nlist.cell()   # Neighbor list
        nb = azplugins.pair.ashbaugh(r_cut=self.mdpara.cutoff, nlist=nl)
        yukawa = hoomd.md.pair.yukawa(r_cut=self.mdpara.cutoff, nlist=nl)
        
        # Calculate pairwise coefficient matrix
        N = len(self.rtypes)
        lambda_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rlambda, dtype=float)
        lambda_mat = (lambda_mat.transpose() + lambda_mat) * .5
        sigma_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rsigma, dtype=float)
        sigma_mat = (sigma_mat.transpose() + sigma_mat) * .5
        epsilon_mat = np.ones((N, N), dtype=float) * np.asarray(self.rcharges, dtype=float)
        epsilon_mat = (epsilon_mat.transpose() * epsilon_mat) * 1.73136

        for i, ri in enumerate(self.rtypes):
            for j, rj in enumerate(self.rtypes):
                nb.pair_coeff.set(ri, rj, lam=lambda_mat[i, j], epsilon=0.8368, sigma=sigma_mat[i, j])
                yukawa.pair_coeff.set(ri, rj, epsilon=epsilon_mat[i, j], kappa=1.0)
        
        # 3 Setup MD
        # Group particles
        grpall = hoomd.group.all()
        # Setup integrator
        hoomd.md.integrate.mode_standard(dt=self.mdpara.dt)
        methodkw = self.mdpara.integrator.upper()
        if methodkw in ["LD", "NVT", "NPT"]:
            kT = self.mdpara.T * .00831446
        if methodkw == "LD":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-integrate.html#hoomd.md.integrate.langevin
            integrator = hoomd.md.integrate.langevin(group=grpall, kT=kT, seed=self.mdpara.seed)
            for ri in self.rtypes:
                integrator.set_gamma(ri, gamma=self.rparams[ri].mass * .001)
        elif methodkw == "NVE":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-integrate.html#hoomd.md.integrate.nve
            integrator = hoomd.md.methods.NVE(group=grpall)
        elif methodkw == "NVT":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-integrate.html#hoomd.md.integrate.nvt
            integrator = hoomd.md.integrate.nvt(group=grpall, kT=kT, tau=self.mdpara.tau)
        elif methodkw == "NPT":
            # For more details, see https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-integrate.html#hoomd.md.integrate.npt
            integrator = hoomd.md.integrate.npt(group=grpall, kT=kT, tau=self.mdpara.tau, tauP=self.mdpara.taup, P=self.mdpara.p)
        else:
            raise ValueError("Unrecognized md integrator type. Only LD/NVE/NVT/NPT are allowed.")

        # 4 Logging
        hoomd.analyze.log(
            filename=self.output.log,
            quantities=['potential_energy', 'pressure_xx', 'pressure_yy', 'pressure_zz', 'temperature','lx','ly','lz', 'pressure_xy', 'pressure_xz', 'pressure_yz'],
            period=self.output.period,
            overwrite=False,
            header_prefix='#')
        hoomd.dump.gsd(
            filename=self.output.traj,
            period=self.output.period,
            group=grpall)
        if self.output.writecpt:
            hoomd.dump.dcd(
                filename=self.output.checkpoint,
                period=self.output.cptperiod,
                group=grpall,
                overwrite=False)
    
    def setup_hoomd_v3(self):
        # 1 Init
        device = hoomd.device.GPU() if self.mdpara.gpu else hoomd.device.CPU()
        if self.mdpara.seed < 0:
            self.mdpara.seed = np.random.randint(0, 1 << 16)
        self.system = hoomd.Simulation(device=device, seed=self.mdpara.seed)
        # Run a simulation using prepared initial configuration file
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
        ashbaugh = hoomd.md.pair.Ashbaugh(default_r_cut=self.mdpara.cutoff, nlist=nl)   # Ash-baugh Van der Waals potential
        yukawa = hoomd.md.pair.Yukawa(default_r_cut=self.mdpara.cutoff, nlist=nl)
        
        # Calculate pairwise coefficient matrix
        N = len(self.rtypes)
        lambda_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rlambda, dtype=float)
        lambda_mat = (lambda_mat.transpose() + lambda_mat) * .5
        sigma_mat = np.zeros((N, N), dtype=float) + np.asarray(self.rsigma, dtype=float)
        sigma_mat = (sigma_mat.transpose() + sigma_mat) * .5
        epsilon_mat = np.ones((N, N), dtype=float) * np.asarray(self.rcharges, dtype=float)
        epsilon_mat = (epsilon_mat.transpose() * epsilon_mat) * 1.73136

        for i, ri in enumerate(self.rtypes):
            for j, rj in enumerate(self.rtypes):
                ashbaugh.params[(ri, rj)] = dict(sigma=sigma_mat[i, j], epsilon=0.8368, lam=lambda_mat[i, j], alpha=1.0)
                yukawa.params[(ri, rj)] = dict(epsilon=epsilon_mat[i, j], kappa=1.0)
        
        # 3 Setup MD
        # Setup integrator
        methods = []
        methodkw = self.mdpara.integrator.upper()
        if methodkw in ["LD", "NVT", "NPT"]:
            kT = self.mdpara.T * .00831446
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
        if self.version_major == 2:
            hoomd.run_upto(step=self.mdpara.steps, limit_hours=self.mdpara.tlimit)
        elif self.version_major == 3:
            self.system.run(steps=self.mdpara.steps, write_at_start=self.mdpara.newrun)
        else:
            raise ValueError("Unsupported HOOMD-Blue version. Failed to setup simulation run task.")
        etime = datetime.now()
        print("Simulation ends at " + datetime.strftime(etime, '%c'))
        print("Total usage time:", timefmt.strfdelta(etime - stime, "%D days %H hours %M minutes %S seconds"))
