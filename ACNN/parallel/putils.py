
class RecArgs(object):
    def __init__(self, args):
        largs = args.split(";")
        self.dict = {}
        for l in largs:
            if "(" in l:
                a, b = l.split("(")[0], "(".join(l.split("(")[1:])
                if a[-1] == "=":
                    a = a[:-1]
                self.dict[a.lower()] = b[:-1]
            elif "=" in l:
                a, b = l.split("=")[0], "=".join(l.split("=")[1:])
                self.dict[a.lower()] = b
            else:
                self.dict[l.lower()] = ""
    
    def __repr__(self):
        largs = []
        for k, v in self.dict.items():
            if v == "":
                largs.append(k)
            elif "(" in v or "=" in v or ("," in v and k != 'fix' and k != 'spfixxy'):
                largs.append("%s(%s)" % (k, v))
            else:
                largs.append("%s=%s" % (k, v))
        return ";".join(largs)

class RefStruct(object):
    def __init__(self, lnumber, xname, xsta, xmul, xidx, pdir):
        self.lnumber = lnumber
        self.xname, self.xsta, self.xmul, self.xidx = xname, xsta, int(xmul), int(xidx)
        self.xdir = "%s/%s.%d/par_%s.xyz.%d" % (pdir, xsta, self.xmul, xname, self.xidx)
        self.xpre = "%s.%d" % (xsta, self.xmul)
    
    @staticmethod
    def solve(l, stage="0", m=0, idx=0, pre_dir="."):
        # format: [stage[.multi]-]name[.idx][:max_number]
        if ":" in l:
            lnumber = int(l.split(":")[1])
            l = l.split(":")[0]
        else:
            lnumber = None
        if "-" not in l:
            lsta = stage
            lmulti = str(m)
        else:
            la, lb = l.split('-')
            if "." in la:
                lsta, lmulti = la.split(".")
            else:
                lsta, lmulti = la, str(m)
            l = lb
        if "." in l:
            lx, ly = l.split('.')
        else:
            lx = l
            ly = str(idx)
        return RefStruct(lnumber, lx, lsta, lmulti, ly, pre_dir)
