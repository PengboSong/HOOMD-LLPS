from .paraset import hoomd_particle as molparams


KR_MOL_PARAMS = {
    # Amino Acid Residues
    "ALA": molparams("ALA",  71.08,  0.0, .504, 0.7297),
    "ARG": molparams("ARG", 156.20,  1.0, .656, 0.0000),
    "ASN": molparams("ASN", 114.10,  0.0, .568, 0.4324),
    "ASP": molparams("ASP", 115.10, -1.0, .558, 0.3784),
    "CYS": molparams("CYS", 103.10,  0.0, .548, 0.6757),
    "GLN": molparams("GLN", 128.10,  0.0, .602, 0.5135),
    "GLU": molparams("GLU", 129.10, -1.0, .592, 0.4595),
    "GLY": molparams("GLY",  57.05,  0.0, .450, 0.6486),
    "HIS": molparams("HIS", 137.10,  0.5, .608, 0.5946),
    "ILE": molparams("ILE", 113.20,  0.0, .618, 0.9730),
    "LEU": molparams("LEU", 113.20,  0.0, .618, 0.9730),
    "LYS": molparams("LYS", 128.20,  1.0, .636, 0.5135),
    "MET": molparams("MET", 131.20,  0.0, .618, 0.8378),
    "PHE": molparams("PHE", 147.20,  0.0, .636, 1.0000),
    "PRO": molparams("PRO",  97.12,  0.0, .556, 1.0000),
    "SER": molparams("SER",  87.08,  0.0, .518, 0.5946),
    "THR": molparams("THR", 101.10,  0.0, .562, 0.6757),
    "TRP": molparams("TRP", 186.20,  0.0, .678, 0.9459),
    "TYR": molparams("TYR", 163.20,  0.0, .646, 0.8649),
    "VAL": molparams("VAL",  99.07,  0.0, .586, 0.8919),
}

URRY_MOL_PARAMS = {
    # Amino Acid Residues
    "ALA": molparams("ALA",  71.08,  0.0, .504, 0.6029),
    "ARG": molparams("ARG", 156.20,  1.0, .656, 0.5588),
    "ASN": molparams("ASN", 114.10,  0.0, .568, 0.5882),
    "ASP": molparams("ASP", 115.10, -1.0, .558, 0.2942),
    "CYS": molparams("CYS", 103.10,  0.0, .548, 0.6471),
    "GLN": molparams("GLN", 128.10,  0.0, .602, 0.5588),
    "GLU": molparams("GLU", 129.10, -1.0, .592, 0.0000),
    "GLY": molparams("GLY",  57.05,  0.0, .450, 0.5735),
    "HIS": molparams("HIS", 137.10,  0.5, .608, 0.7647),
    "ILE": molparams("ILE", 113.20,  0.0, .618, 0.7059),
    "LEU": molparams("LEU", 113.20,  0.0, .618, 0.7206),
    "LYS": molparams("LYS", 128.20,  1.0, .636, 0.3824),
    "MET": molparams("MET", 131.20,  0.0, .618, 0.6765),
    "PHE": molparams("PHE", 147.20,  0.0, .636, 0.8235),
    "PRO": molparams("PRO",  97.12,  0.0, .556, 0.7588),
    "SER": molparams("SER",  87.08,  0.0, .518, 0.5882),
    "THR": molparams("THR", 101.10,  0.0, .562, 0.5882),
    "TRP": molparams("TRP", 186.20,  0.0, .678, 1.0000),
    "TYR": molparams("TYR", 163.20,  0.0, .646, 0.8971),
    "VAL": molparams("VAL",  99.07,  0.0, .586, 0.6647),
}

MOL_PARAMS = URRY_MOL_PARAMS
