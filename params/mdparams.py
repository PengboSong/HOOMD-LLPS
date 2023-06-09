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
    
    def setkey(self, key: str, vtype: type, defv: Any = None, required: str = '', expect: Optional[List[Any]] = None):
        newleaf = TreeNode(key=key, vtype=vtype, defv=defv)
        if required and (expect is not None):
            newleaf.expect = expect
            queue = [self.root]
            while (queue):
                node = queue.pop(0)
                if node.key == required:
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
            if (node.expect is not None) and ((isinstance(node.expect, List) and val not in node.expect) or (not isinstance(node.expect, List) and val != node.expect)):
                continue
            v = ''
            if node.key:
                getmethod = vtypemap[node.vtype] if node.vtype in vtypemap else sect.get
                if node.defv:
                    v = getmethod(node.key, fallback=node.defv)
                else:
                    v = getmethod(node.key)
                self.params[node.key] = v
                setattr(self, node.key, v)
                if printoptions:
                    print(f"Option {node.key} = {v}")
            if not node.isleaf():
                for child in node.nodes:
                    queue.append((child, v))
