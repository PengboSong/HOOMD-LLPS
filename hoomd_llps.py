#!/usr/bin/env python3
# coding = utf-8


import argparse
import os.path

import mdsys


def hoomd_version() -> str:
    import hoomd
    try:
        return hoomd.version.version
    except AttributeError:
        pass
    try:
        return hoomd.__version__
    except AttributeError:
        pass


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--file',
        type=str,
        nargs='+',
        help="File path lists of one letter AA/NA sequence / PDB format file "
             "with residue particles only"
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        default='mc.ini',
        help="Molecular dynamics configuration file"
    )
    flags, _ = parser.parse_known_args()
    return flags


def main():
    FLAGS = parseArgs()
    MOLS = mdsys.MDSystem(version=hoomd_version())
    if FLAGS.config and os.path.isfile(FLAGS.config):
        MOLS.load_mdparams(FLAGS.config)
    else:
        raise IOError("Molecular dynamics configuration file is required.")
    if FLAGS.file and len(FLAGS.file) != 0:
        ftype = ""
        for i, f in enumerate(FLAGS.file):
            if not os.path.isfile(f):
                raise IOError(f"Can not find target file {f}.")
            if ftype and os.path.splitext(f)[-1].lower() != ftype:
                raise ValueError("Mixed input file types of seq, pdb and gsd are not supported.")
            ftype = os.path.splitext(f)[-1].lower()
            if ftype not in ['.seq', '.pdb', '.gsd']:
                raise ValueError("Only input types of seq, pdb and gsd are supported.")
            if ftype == ".pdb":
                MOLS.load_pdb(f)
                MOLS.configure_box()
            elif ftype == ".seq":
                with open(f, 'r') as fseq:
                    content = fseq.read().strip()
                    MOLS.load_oneletterseq(
                        seq=content[4:],
                        seqtype=content[:4],
                        xyzoffset=(i * 5., 0., 0.))
                    MOLS.configure_box()
            elif f != MOLS.mdpara.initgsd:
                raise ValueError(f"Only accept initial configuration filename set in {FLAGS.config} by key initgsd.")
    else:
        raise IOError("Input file lists of sequences or conformations are required.")    
    MOLS.setup()
    MOLS.run()


if __name__ == '__main__':
    main()