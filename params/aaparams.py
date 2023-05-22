from .paraset import hoomd_particle as aaparams


KR_AA_PARAMS = {
    'A': aaparams("ALA",  71.08,  0.0, .504, 0.7297),
    'R': aaparams("ARG", 156.20,  1.0, .656, 0.0000),
    'N': aaparams("ASN", 114.10,  0.0, .568, 0.4324),
    'D': aaparams("ASP", 115.10, -1.0, .558, 0.3784),
    'C': aaparams("CYS", 103.10,  0.0, .548, 0.6757),
    'Q': aaparams("GLN", 128.10,  0.0, .602, 0.5135),
    'E': aaparams("GLU", 129.10, -1.0, .592, 0.4595),
    'G': aaparams("GLY",  57.05,  0.0, .450, 0.6486),
    'H': aaparams("HIS", 137.10,  0.5, .608, 0.5946),
    'I': aaparams("ILE", 113.20,  0.0, .618, 0.9730),
    'L': aaparams("LEU", 113.20,  0.0, .618, 0.9730),
    'K': aaparams("LYS", 128.20,  1.0, .636, 0.5135),
    'M': aaparams("MET", 131.20,  0.0, .618, 0.8378),
    'F': aaparams("PHE", 147.20,  0.0, .636, 1.0000),
    'P': aaparams("PRO",  97.12,  0.0, .556, 1.0000),
    'S': aaparams("SER",  87.08,  0.0, .518, 0.5946),
    'T': aaparams("THR", 101.10,  0.0, .562, 0.6757),
    'W': aaparams("TRP", 186.20,  0.0, .678, 0.9459),
    'Y': aaparams("TYR", 163.20,  0.0, .646, 0.8649),
    'V': aaparams("VAL",  99.07,  0.0, .586, 0.8919),
}

URRY_AA_PARAMS = {
    'A': aaparams("ALA",  71.08,  0.0, .504, 0.6029),
    'R': aaparams("ARG", 156.20,  1.0, .656, 0.5588),
    'N': aaparams("ASN", 114.10,  0.0, .568, 0.5882),
    'D': aaparams("ASP", 115.10, -1.0, .558, 0.2942),
    'C': aaparams("CYS", 103.10,  0.0, .548, 0.6471),
    'Q': aaparams("GLN", 128.10,  0.0, .602, 0.5588),
    'E': aaparams("GLU", 129.10, -1.0, .592, 0.0000),
    'G': aaparams("GLY",  57.05,  0.0, .450, 0.5735),
    'H': aaparams("HIS", 137.10,  0.5, .608, 0.7647),
    'I': aaparams("ILE", 113.20,  0.0, .618, 0.7059),
    'L': aaparams("LEU", 113.20,  0.0, .618, 0.7206),
    'K': aaparams("LYS", 128.20,  1.0, .636, 0.3824),
    'M': aaparams("MET", 131.20,  0.0, .618, 0.6765),
    'F': aaparams("PHE", 147.20,  0.0, .636, 0.8235),
    'P': aaparams("PRO",  97.12,  0.0, .556, 0.7588),
    'S': aaparams("SER",  87.08,  0.0, .518, 0.5882),
    'T': aaparams("THR", 101.10,  0.0, .562, 0.5882),
    'W': aaparams("TRP", 186.20,  0.0, .678, 1.0000),
    'Y': aaparams("TYR", 163.20,  0.0, .646, 0.8971),
    'V': aaparams("VAL",  99.07,  0.0, .586, 0.6647),
}

AA_PARAMS = URRY_AA_PARAMS
