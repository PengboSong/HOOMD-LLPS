#!/usr/bin/env python3

import argparse

import mdsys


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--file',
        nargs='+',
        help="One letter amino acid sequence file"
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
    MOLS = mdsys.MDSystem()
    MOLS.load_mdparams(FLAGS.config)
    if not FLAGS.file:
        raise IOError("Input sequence file required.")
    '''
    for i, f in enumerate(FLAGS.file):
        with open(f, 'r') as fobj:
            content = fobj.read().strip()
            MOLS.load_oneletterseq(
                seq=content[4:],
                seqtype=content[:3],
                xyzoffset=(i * 5., 0., 0.))
    '''
    for f in FLAGS.file:
        MOLS.load_pdb(f)
    MOLS.configure_box()
    MOLS.setup()
    MOLS.run()


if __name__ == '__main__':
    main()