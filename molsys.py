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
        self.P.typeid = np.array([], dtype=np.uint32)
        self.P.mass = np.array([], dtype=np.float32)
        self.P.charge = np.array([], dtype=np.float32)
        self.P.position = np.array([], dtype=np.float32).reshape(-1, 3)
        self.B = self.S.bonds
        self.B.N = 0
        self.B.types = []
        self.B.typeid = np.array([], dtype=np.uint32)
        self.B.group = np.array([], dtype=np.uint32).reshape(-1, 2)

        self.rindex = {}
        self.rparams = {}
        self.rtypes = []
        self.rcharges = []
        self.rlambda = []
        self.rsigma = []
    
    def load_oneletterseq(self, seq: str, seqtype: str, xyzoffset: Tuple[float, float, float], zigzagdeg: float = 45.):
        # Choose parameter set according to seqtype
        tag = seqtype[:3].upper()
        placeholder = seqtype[3]
        if tag in ["DNA", "RNA", "PRO"] and placeholder == '|':
            if tag == "DNA":
                PARAMS_SET = DNA_PARAMS
                rtype = "deoxynucleotide"
                bond = .500
                btype = "NT_bond"
            elif tag == "RNA":
                PARAMS_SET = RNA_PARAMS
                rtype = "nucleotide"
                bond = .500
                btype = "NT_bond"
            else:
                PARAMS_SET = AA_PARAMS
                rtype = "amino acid"
                bond = .381
                btype = "AA_bond"
        else:
            raise ValueError("Expect the first 4 letters of sequence to be identifier "
                             "with the format of `TAG|`, where TAG should be one of "
                             "DNA, RNA, PRO")

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
        self.P.typeid = np.append(self.P.typeid, [self.rindex[res] for res in rseq])
        self.P.mass = np.append(self.P.mass, [self.rparams[res].mass for res in rseq])
        self.P.charge = np.append(self.P.charge, [self.rparams[res].charge for res in rseq])
        
        # Initialize N x 3 particles position matrix
        sxyz = self.P.position
        xyz = np.zeros((rcount, 3), dtype=np.float32) + np.asarray(xyzoffset)
        # zigzagdeg is the degree between the first bond vector and the positive z-axis
        # 0 ~ 180 deg refers to the positive y-axis side, 180 ~ 360 def refers to the negative y-axis side
        zigzagrad = min(max(zigzagdeg, 0.), 360.) / 180. * np.pi
        dy = bond * np.sin(zigzagrad)
        dz = bond * np.cos(zigzagrad)
        for i in range(rcount):
            xyz[i, :] += np.array([0., (i & 1) * dy, i * dz])
        # Merge sub position matrix to main matrix
        assert (sxyz.ndim == 2) & (sxyz.shape[1] == 3)
        self.P.position = np.vstack([sxyz, xyz])

        # Set bonds
        if btype not in self.B.types:
            self.B.types.append(btype)
        bid = self.B.types.index(btype)
        bcount = rcount - 1
        self.B.N += bcount
        self.B.typeid = np.append(self.B.typeid, [bid] * bcount)
        # Initialize N - 1 x 2 bonds pair matrix
        bpairs = np.zeros((bcount, 2), dtype=np.uint32)
        bpairs[:, 0] = np.arange(bcount, dtype=np.uint32) + sid
        bpairs[:, 1] = np.arange(bcount, dtype=np.uint32) + sid + 1
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
                    atom = line[12:16].strip().upper()
                    res = line[17:20].strip().upper()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if res in MOL_PARAMS:
                        mol = res
                    elif atom in MOL_PARAMS:
                        mol = atom
                    else:
                        raise ValueError(f"Unexpected residue name. Please check atom name {atom} "
                                         f"or residue name {res} at line {linen}.")
                    rset.add(mol)
                    rseq.append(mol)
                    rcoord.append([x, y, z])
                    rcount += 1
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
        self.P.typeid = np.append(self.P.typeid, [self.rindex[res] for res in rseq])
        self.P.mass = np.append(self.P.mass, [self.rparams[res].mass for res in rseq])
        self.P.charge = np.append(self.P.charge, [self.rparams[res].charge for res in rseq])
        
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
        if len(bpairs) == 0:   # NO CONECT lines in PDB file
            # Initialize N - 1 x 2 bonds pair matrix
            bpairs = np.zeros((rcount - 1, 2), dtype=np.uint32)
            # Connect all adjacent particles
            bpairs[:, 0] = np.arange(rcount - 1, dtype=np.uint32) + sid
            bpairs[:, 1] = np.arange(rcount - 1, dtype=np.uint32) + sid + 1
        else:
            bpairs = list(map(lambda bij: [atomidmap[bij[0]], atomidmap[bij[1]]], bpairs))
            bpairs = np.asarray(bpairs, dtype=np.uint32)
        # Bond pairs matrix should have `bcount` rows and 2 columns
        assert (bpairs.ndim == 2) & (bpairs.shape[1] == 2)
        bcount = bpairs.shape[0]
        self.B.N += bcount
        self.B.typeid = np.append(self.B.typeid, [bid] * bcount)
        # Merge sub pair matrix to main matrix
        bgroup = self.B.group
        assert (bgroup.ndim == 2) & (bgroup.shape[1] == 2)
        self.B.group = np.vstack([bgroup, bpairs])
