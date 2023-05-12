class ParaSet(dict):
    def __getattr__(self, __k):
        return self.get(__k)
    
    def __setattr__(self, __k, __v):
        self.update(__k = __v)


def hoomd_particle(residue: str, mass: float, charge: float, sigma: float, lambda_: float):
    return ParaSet({
        "residue": residue,
        "mass": mass,
        "charge": charge,
        "sigma": sigma,
        "_lambda": lambda_})


def hoomd_bond(bond: str, k: float, r0: float):
    return ParaSet({
        "bond": bond,
        "k": k,
        "r0": r0})
