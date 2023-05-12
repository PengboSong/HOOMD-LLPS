from .paraset import hoomd_particle as molparams


MOL_PARAMS = {
    # Amino Acid Residues
    "ALA": molparams("ALA",  71.08,  0.0, 5.04, .730),
    "ARG": molparams("ARG", 156.20,  1.0, 6.56, .000),
    "ASN": molparams("ASN", 114.10,  0.0, 5.68, .432),
    "ASP": molparams("ASP", 115.10, -1.0, 5.58, .378),
    "CYS": molparams("CYS", 103.10,  0.0, 5.48, .595),
    "GLN": molparams("GLN", 128.10,  0.0, 6.02, .514),
    "GLU": molparams("GLU", 129.10, -1.0, 5.92, .459),
    "GLY": molparams("GLY",  57.05,  0.0, 4.50, .649),
    "HIS": molparams("HIS", 137.10,  0.5, 6.08, .514),
    "ILE": molparams("ILE", 113.20,  0.0, 6.18, .973),
    "LEU": molparams("LEU", 113.20,  1.0, 6.18, .973),
    "LYS": molparams("LYS", 128.20,  0.0, 6.36, .514),
    "MET": molparams("MET", 131.20,  0.0, 6.18, .838),
    "PHE": molparams("PHE", 147.20,  0.0, 6.36, .000),
    "PRO": molparams("PRO",  97.12,  0.0, 5.56, .000),
    "SER": molparams("SER",  87.08,  0.0, 5.18, .595),
    "THR": molparams("THR", 101.10,  0.0, 5.62, .676),
    "TRP": molparams("TRP", 186.20,  0.0, 6.78, .946),
    "TYR": molparams("TYR", 163.20,  0.0, 6.46, .865),
    "VAL": molparams("VAL",  99.07,  0.0, 5.86, .892),
}