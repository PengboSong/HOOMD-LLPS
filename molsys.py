import os.path
from typing import Tuple

import gsd.hoomd
import numpy as np

from params import AA_PARAMS, DNA_PARAMS, RNA_PARAMS, MOL_PARAMS


class MolSystem:
    def __init__(self):
        self.S = gsd.hoomd.Snapshot()   # Snapshot is replaced by Frame since gsd 2.8.0
        self.P = self.S.particles
        self.P.N = 0
        self.P.types = []
        self.P.typeid = []
        self.P.mass = []
        self.P.charge = []
        self.P.position = np.array([], dtype=np.float32).reshape(-1, 3)
        self.B = self.S.bonds
        self.B.N = 0
        self.B.types = []
        self.B.typeid = []
        self.B.group = np.array([], dtype=np.uint32).reshape(-1, 2)

        self.rindex = {}
        self.rparams = {}
        self.rtypes = []
        self.rcharges = []
        self.rlambda = []
        self.rsigma = []
    
    def load_oneletterseq(self, seq: str, seqtype: str, xyzoffset: Tuple[float, float, float]):
        # Choose parameter set according to seqtype
        if seqtype.upper() == "DNA":
            PARAMS_SET = DNA_PARAMS
            rtype = "deoxynucleotide"
            bond = .50
            btype = "NT_bond"
        elif seqtype.upper() == "RNA":
            PARAMS_SET = RNA_PARAMS
            rtype = "nucleotide"
            bond = .50
            btype = "NT_bond"
        else:
            PARAMS_SET = AA_PARAMS
            rtype = "amino acid"
            bond = .38
            btype = "AA_bond"

        # One-letter residue set & count
        rset = set()
        rseq = []
        rcount = 0
        for i, letter in enumerate(seq):
            if letter in PARAMS_SET:
                rset.add(letter)
                rseq.append(PARAMS_SET[letter].residue)
                # Add system particles count by 1
                rcount += 1
            else:
                raise ValueError(f"Unexpected {rtype} residue letter. "
                                 f"Please check letter {letter} at seq {i + 1}.")
        
        sid = self.P.N
        tid = len(self.rindex)
        for i, letter in enumerate(rset):
            params = PARAMS_SET[letter]
            res = params.residue
            # Add system particles type / ID / mass / charge
            self.P.types.append(res)
            self.rindex[res] = tid + i
            self.rparams[res] = params

            self.rtypes.append(params.residue)
            self.rcharges.append(params.charge)
            self.rlambda.append(params._lambda)
            self.rsigma.append(params.sigma)
        
        self.P.N += rcount
        for res in rseq:
            self.P.typeid.append(self.rindex[res])
            self.P.mass.append(self.rparams[res].mass)
            self.P.charge.append(self.rparams[res].charge)
        
        # Initialize N x 3 particles position matrix
        sxyz = self.P.position
        xyz = np.zeros((rcount, 3), dtype=float) + np.asarray(xyzoffset)
        zcenter = -.5 * rcount * bond
        '''
        zmax = sxyz.max(axis=0)[2]
        zcenter = zmax - .5 * rcount * bond
        '''
        for i in range(rcount):
            xyz[i, :] += np.array([0., 0., zcenter + i * bond])
        # Merge sub position matrix to main matrix
        assert (sxyz.ndim == 2) & (sxyz.shape[1] == 3)
        self.P.position = np.vstack([sxyz, xyz])

        # Set bonds
        if btype not in self.B.types:
            self.B.types.append(btype)
        bid = self.B.types.index(btype)
        bcount = rcount - 1
        self.B.N += bcount
        # Initialize N - 1 x 2 bonds pair matrix
        bpairs = np.zeros((bcount, 2), dtype=int)
        for i in range(bcount):
            self.B.typeid.append(bid)
            bpairs[i, :] = [sid + i, sid + i + 1]
        # Merge sub pair matrix to main matrix
        bgroup = self.B.group
        assert (bgroup.ndim == 2) & (bgroup.shape[1] == 2)
        self.B.group = np.vstack([bgroup, bpairs])
    

    def load_pdb(self, PDB: str, btype: str="AA_bond"):
        if not os.path.isfile(PDB):
            raise FileNotFoundError(f"Can not find PDB file {PDB}")
        
        linen  = 0       # Currently reading line row number 
        rset = set()     # Unique residue names
        rseq = []        # List of residue names in order
        rcoord = []      # 2-D list of residue C-alpha coordinates
        rcount = 0       # Number of residues
        bpairs = []      # 2-D list of residue connections
        resids = set()   # Unique residue IDs
        with open(PDB, 'r') as f:
            for line in f.readlines():
                linen += 1
                mark = line[:6].strip().upper()
                if mark in ["ATOM", "HETATM"]:
                    res = line[17:20].strip().upper()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if res in MOL_PARAMS:
                        rset.add(res)
                        rseq.append(res)
                        rcoord.append([x, y, z])
                        rcount += 1
                    else:
                        raise ValueError(f"Unexpected residue name. "
                                         f"Please check 3-letter residue name {res} at line {linen}.")
                elif mark == "CONECT":
                    content = line[6:]
                    Warning(f"Non-standard connect line format at line {linen}.")
                    aids = [int(content[5*i:5*i+5]) for i in range(len(content) // 5)]
                    resids |= set(aids)
                    if len(aids) < 2:
                        raise ValueError(f"Wrong connect line format with less than 2 atoms in one line at {linen}.")
                    ai = aids[0]
                    for aj in aids[1:]:
                        if ai < aj:
                            bpairs.append((ai, aj))

        sid = self.P.N
        tid = len(self.rindex)
        for i, res in enumerate(rset):
            params = MOL_PARAMS[res]
            # Add system particles type / ID / mass / charge
            self.P.types.append(res)
            self.rindex[res] = tid + i
            self.rparams[res] = params

            self.rtypes.append(params.residue)
            self.rcharges.append(params.charge)
            self.rlambda.append(params._lambda)
            self.rsigma.append(params.sigma)
            
        self.P.N += rcount
        for res in rseq:
            self.P.typeid.append(self.rindex[res])
            self.P.mass.append(self.rparams[res].mass)
            self.P.charge.append(self.rparams[res].charge)
        
        # Set N x 3 particles position matrix
        # PDB coordinates unit in A -> System coordinates unit in nm
        xyz = np.array(rcoord, dtype=float) * .1
        # Merge sub position matrix to main matrix
        sxyz = self.P.position
        assert (sxyz.ndim == 2) & (sxyz.shape[1] == 3)
        self.P.position = np.vstack([sxyz, xyz])
        
        # Set bonds
        if btype not in self.B.types:
            self.B.types.append(btype)
        bid = self.B.types.index(btype)
        # Remapping index of bpairs
        # New index should start from tid (included)
        atomidmap = dict(zip(sorted(resids), range(sid, sid + len(resids))))
        bpairs = list(map(lambda bij: [atomidmap[bij[0]], atomidmap[bij[1]]], bpairs))
        bpairs = np.asarray(bpairs, dtype=int)
        # Bond pairs matrix should have `bcount` rows and 2 columns
        assert (bpairs.ndim == 2) & (bpairs.shape[1] == 2)
        bcount = bpairs.shape[0]
        self.B.N += bcount
        self.B.typeid.extend([bid] * bcount)
        # Merge sub pair matrix to main matrix
        bgroup = self.B.group
        assert (bgroup.ndim == 2) & (bgroup.shape[1] == 2)
        self.B.group = np.vstack([bgroup, bpairs])

