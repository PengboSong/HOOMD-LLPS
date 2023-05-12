MD_PARAMS_TYPES = {
    "main": {
        "gpu": bool,
        "seed": int,
        "initgsd": str,
        "autobox": bool,
        "xbox": float,
        "ybox": float,
        "zbox": float,
        "integrator": str,
        "dt": float,
        "steps": int,
        "newrun": bool,
        "T": float,
        "tau": float,
        "taup": float,
        "p": float,
    },
    "output": {
        "period": int,
        "traj": str,
        "cptperiod": int,
        "checkpoint": str,
    }
}