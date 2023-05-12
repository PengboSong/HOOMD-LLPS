#!/usr/bin/env python3
# coding=utf-8


import argparse
import os.path

import numpy as np

from gmx.structure.mol_matrix import MolMatrix


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--file',
        type=str,
        help="PDB format file for initial molecules template."
    )
    parser.add_argument(
        '-o', '--out',
        type=str,
        help="Output molecules array conformation file."
    )
    parser.add_argument(
        '-n', '--array',
        type=int,
        nargs='+',
        help="Number of molecule copies in X, Y, Z dimensions (unit: nm)."
    )
    parser.add_argument(
        '-l', '--spacing',
        type=float,
        nargs='+',
        help="Spacing between molecule copies in X, Y, Z dimensions."
    )
    flags, _ = parser.parse_known_args()
    return flags


def gen_mol_array(mols: MolMatrix, nx: int, ny: int, nz: int,
                  lx: float, ly: float, lz: float, auto_center: bool = True) -> MolMatrix:
    """Generate a 3D-array with copies of given molecules.

    Args:
        mols: Input mol matrix with coordinates and connect information.
        nx: Number of molecule copies along X axis.
        ny: Number of molecule copies along Y axis.
        nz: Number of molecule copies along Z axis.
        lx: Spacing between molecule copies along X axis (unit: nm).
        ly: Spacing between molecule copies along Y axis (unit: nm).
        lz: Spacing between molecule copies along Z axis (unit: nm).
        auto_center: Whether to center generated array to origin (Center is calculated by
                     coordinates only, ignoring mass).
    
    Returns:
        Modified mol matrix with residue index column renumbered to array index. Output
        mol matrix should have nx X ny X nz input molecule copies.
    """
    # Number of molecule copies should be positive integers
    nx = max(1, nx)
    ny = max(1, ny)
    nz = max(1, nz)

    # Make copies along X axis
    for i in range(len(mols)):
        for xi in range(1, nx):
            newmol = mols[i].copy()
            newmol.move([xi * lx, 0., 0.])
            mols.append(newmol)
    # Make copies along Y axis
    for i in range(len(mols)):
        for yi in range(1, ny):
            newmol = mols[i].copy()
            newmol.move([0., yi * ly, 0.])
            mols.append(newmol)
    # Make copies along Z axis
    for i in range(len(mols)):
        for zi in range(1, nz):
            newmol = mols[i].copy()
            newmol.move([0., 0., zi * lz])
            mols.append(newmol)
    
    if (auto_center):
        mols.genbox()
        mols.do_shift()
    
    return mols


def main():
    FLAGS = parseArgs()
    assert len(FLAGS.array) == 3, "Expect 3 array size arguments."
    assert len(FLAGS.spacing) == 3, "Expect 3 spacing arguments."
    assert os.path.isfile(FLAGS.file) and os.path.splitext(FLAGS.file)[-1].lower() == ".pdb", "Expect PDB format file."
    mols = MolMatrix()
    mols.from_pdb(FLAGS.file)
    mols = gen_mol_array(mols, *FLAGS.array, *FLAGS.spacing, auto_center=True)
    mols.to_file(FLAGS.out)


if __name__ == '__main__':
    main()
