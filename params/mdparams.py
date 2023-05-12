import configparser
from collections.abc import Iterable
from typing import Any, Dict, List, Optional, Tuple


class TreeNode(object):
    def __init__(self, key: str = '', vtype: type = Any, defv: Any = None):
        self.key = key
        self.vtype = vtype
        self.defv = defv
        self.expect = None
        self.nodes = []
    
    def isleaf(self) -> bool:
        return len(self.nodes) == 0

    def addchild(self, node: 'TreeNode'):
        self.nodes.append(node)


class MDParams(object):
    def __init__(self, section: str, importdict: Dict[str, Tuple[type, Any, str, Any]]):
        self.section = section
        self.params = {}
        self.root = TreeNode()

        self.fromdict(importdict)
    
    def __getitem__(self, __k: str) -> Any:
        return self.params.get(__k)
    
    def __setitem__(self, __k: str, __v: Any):
        self.params[__k] == __v

    def __getattr__(self, __k: str) -> Any:
        return self.params.get(__k)
    
    def __setattr__(self, __k: str, __v: Any):
        self.params[__k] == __v
    
    def setkey(self, key: str, vtype: type, defv: Any = None, required: str = '', expect: Optional[List[Any]] = None):
        newleaf = TreeNode(key=key, vtype=vtype, defv=defv)
        if required and (expect is not None):
            newleaf.expect = expect
            queue = [self.root]
            while (queue):
                node = queue.pop(0)
                if node.val == required:
                    node.addchild(newleaf)
                    break
                if not node.isleaf():
                    queue.extend(node.nodes)
        else:
            self.root.addchild(newleaf)
    
    def fromdict(self, importdict: Dict[str, Tuple[type, Any, str, Any]]):
        for key, params in importdict.items():
            self.setkey(key=key, **params)

    def digest(self, config: configparser.ConfigParser, printoptions: bool = True):
        queue = [(self.root, None)]
        sect = config[self.section]
        vtypemap = {int: sect.getint, float: sect.getfloat, bool: sect.getboolean, str: sect.get}
        while (queue):
            node, val = queue.pop(0)
            if (node.expect is not None) and (val not in node.val):
                continue
            if node.key:
                getmethod = vtypemap[node.vtype] if node.vtype in vtypemap else sect.get
                if node.defv:
                    v = getmethod(node.key, fallback=node.defv)
                else:
                    v = getmethod(node.key)
                self.params[node.key] = v
                if printoptions:
                    print(f"Option {node.key} = {v}")
                if not node.isleaf():
                    for child in node.nodes:
                        queue.append((child, v))


MD_PARAMS_TYPES = {
    "gpu": {"vtype": bool, "defv": True},
    "seed": {"vtype": int},
    "initgsd": {"vtype": str, "defv": "init.gsd"},
    "autobox": {"vtype": bool, "defv": True},
    "xbox": {"vtype": float, "required": "autobox", "expect": [False]},
    "ybox": {"vtype": float, "required": "autobox", "expect": [False]},
    "zbox": {"vtype": float, "required": "autobox", "expect": [False]},
    "em": {"vtype": bool, "defv": False},
    "forcetol": {"vtype": float, "required": "em", "expect": [True]},
    "angmomtol": {"vtype": float, "required": "em", "expect": [True]},
    "energytol": {"vtype": float, "required": "em", "expect": [True]},
    "integrator": {"vtype": str},
    "steps": {"vtype": int, "defv": True},
    "newrun": {"vtype": bool, "defv": True},
    "dt": {"vtype": float},
    "T": {"vtype": float, "required": "integrator", "expect": ["ld", "nve", "nvt", "npt"]},
    "tau": {"vtype": float, "required": "integrator", "expect": ["nvt", "npt"]},
    "p": {"vtype": float, "required": "integrator", "expect": ["npt"]},
    "taup": {"vtype": float, "required": "integrator", "expect": ["npt"]},
}

OUTPUT_PARAMS_TYPES = {
    "period": {"vtype": int},
    "traj": {"vtype": str, "defv": "md_traj.gsd"},
    "writecpt": {"vtype": bool, "defv": True},
    "cptperiod": {"vtype": int, "required": "writecpt", "except": [True]},
    "checkpoint": {"vtype": str, "defv": "md_cpt.dcd", "required": "writecpt", "except": [True]},
}
