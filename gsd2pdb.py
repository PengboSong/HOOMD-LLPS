#!/usr/bin/env python3
#coding=utf-8


import argparse
import os.path
from typing import Tuple

import gsd.hoomd
import numpy as np

from gmx.structure.atom_matrix import AtomMatrix


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--file',
        type=str,
        help="GSD format HOOMD trajectory file."
    )
    flags, _ = parser.parse_known_args()
    return flags


def fix_pbc(cm: np.ndarray, bondgrp: np.ndarray, boxpara: Tuple[float, float, float]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fix cooridnate matrix with Periodic Boundary Condition and make molecules whole.
    
    Args:
        cm: Coordinate matrix.
        bondgrp: Bond i, j collections.
        boxpara: Box lengths in X, Y, Z dimensions.
    
    Returns:
        cm: Fixed coordinate matrix that molecules are whole again.
        bondlen: Bond length collections.
    """
    bondlen = []
    boxpara = np.asarray(boxpara, dtype=float)
    halfbox = .5 * boxpara

    for i, j in bondgrp:
        x1 = cm[i]
        x2 = cm[j]
        sign = 1 * (x2 < (x1 - halfbox)) - 1 * (x2 > (x1 + halfbox))
        cm[j] += sign * boxpara
        bondlen.append(np.linalg.norm(cm[j] - cm[i]))
    bondlen = np.asarray(bondlen, dtype=float)
    return cm, bondlen

    
def write_pdb(GSD: str):
    fname = os.path.splitext(GSD)[0]
    with gsd.hoomd.open(name=GSD, mode='rb') as traj:
        for i, snap in enumerate(traj):
            N = snap.particles.N
            R = snap.particles.position
            B = snap.bonds.group
            R, _ = fix_pbc(R, B, boxpara=snap.configuration.box[:3])
            fullname = f"{fname}_{i}.pdb"
            frame = AtomMatrix()
            frame.atomn = N
            frame.atoms = ['CA'] * N
            frame.elements = ['C'] * N
            frame.molids = [m + 1 for m in range(N)]
            frame.molnms = [snap.particles.types[r] for r in snap.particles.typeid]
            # System coordinates unit in nm -> PDB coordinates unit in A
            frame.xyzs = R * 10.
            frame.velos = np.zeros((N, 3), dtype=float)
            frame.bonds = np.asarray(B + 1, dtype=int)
            frame.clean()
            frame.to_pdb(fullname)


def main():
    FLAGS = parseArgs()
    write_pdb(FLAGS.file)


if __name__ == '__main__':
    main()
