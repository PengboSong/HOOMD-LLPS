from .paraset import hoomd_bond


BOND_PARAMS = {
    "AA_bond": hoomd_bond("AA_bond", 8360.0, .381),
    "NT_bond": hoomd_bond("NT_bond", 8360.0, .500),
}