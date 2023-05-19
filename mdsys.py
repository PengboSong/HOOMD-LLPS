import configparser
from datetime import datetime
import os

import hoomd, hoomd.md
import azplugins
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
        self.output.digest(CONFIG)
    
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
        nb = azplugins.pair.ashbaugh(r_cut=3.5, nlist=nl)
        yukawa = hoomd.md.pair.yukawa(r_cut=3.5, nlist=nl)
        
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
                nb.pair_coeff.set(ri, rj, lam=lambda_mat[i, j], epsilon=.8368, sigma=sigma_mat[i, j])
                yukawa.pair_coeff.set(ri, rj, epsilon=epsilon_mat[i, j], kappa=1.)
        
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

    def run(self):
        # Run simulation
        stime = datetime.now()
        print("Simulation starts at " + datetime.strftime(stime, '%c'))
        hoomd.run_upto(step=self.mdpara.steps, limit_hours=self.mdpara.tlimit)
        etime = datetime.now()
        print("Simulation ends at " + datetime.strftime(etime, '%c'))
        print("Total usage time:", timefmt.strfdelta(etime - stime, "%D days %H hours %M minutes %S seconds"))
