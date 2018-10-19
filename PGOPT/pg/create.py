
import os, re
from pg.utils import optcopy

# re for resolving cluster name
rxa = r'^\s*([A-Za-z][a-z]*)\s*([0-9]+)(.*)$'
rxb = r'^\s*\(\s*([^\)]+)\s*\)\s*([0-9]+)(.*)$'
rxc = r'^\s*([A-Za-z][a-z]*)(.*)$'
rxd = r'^\s*$'

# resolve cluster name
def elem_char(name):
    dct = []
    while len(re.findall(rxd, name)) == 0:
        ra = re.findall(rxa, name)
        if len(ra) == 0:
            ra = re.findall(rxb, name)
        if len(ra) != 0:
            ne = ra[0][0][:1].upper() + ra[0][0][1:]
            dct += [ne, ra[0][1]]
            name = ra[0][2]
        else:
            ra = re.findall(rxc, name)
            ne = ra[0][0][:1].upper() + ra[0][0][1:]
            dct += [ne]
            name = ra[0][1]
    return ' '.join(dct)

class PGCreate(object):
    def __init__(self, pre, rseed, dir_name, scr_name):
        self.number = pre["number"]
        self.method = pre["method"]
        self.name = pre["name"]
        if self.method == "ck":
            self.drel = pre["drel"]
        self.rseed = rseed
        self.dir_name = dir_name
        self.scr_name = scr_name
    
    def run(self):
        if self.method == "ck": return self.runck()
        elif self.method == "blda": return self.runblda()
        elif self.method == "bc": return self.runbc()
        elif self.method == "rv": return self.runrv()
        else:
            raise RuntimeError("Unknown method: %s" % self.method)
    
    def runblda(self):
        return [ "create" ]
    
    def runck(self):
        ckpath = './CK/%s' % self.dir_name
        ckopts = { "@GNAME": self.name, "@GDREL": str(self.drel), 
            "@GNUMBER": str(self.number), "@GRSEED": str(abs(self.rseed % (2**30 - 1))) }
        if "AFFCKHOME" not in os.environ:
            raise RuntimeError('must have AFFCKHOME environment variable!')
        if not os.path.exists(ckpath):
            os.makedirs(ckpath)
        if os.path.exists(ckpath + '/ck_structs'):
            os.popen('rm -rf ' + ckpath + '/ck_structs').read()
        optcopy(self.scr_name + '/runck.sh', ckpath + '/runck.sh', ckopts)
        os.chdir(ckpath)
        os.popen('./runck.sh').read()
        os.chdir('./ck_structs')
        os.popen('zip ../ck_structs.zip *.xyz').read()
        os.chdir('../../..')
        os.popen('rm -rf ' + ckpath + '/ck_structs').read()
        return [ "read", "../../CK/%s/ck_structs.zip" % self.dir_name ]
    
    def runbc(self):
        bcpath = './BC/%s' % self.dir_name
        bcopts = { "@GNAME": self.name, "@GNUMBER": str(self.number), 
            "@GRSEED": str(abs(self.rseed % (2**30 - 1))) }
        import imp
        try: imp.find_module('ase')
        except ImportError:
            raise RuntimeError('must have ASE installed!')
        if not os.path.exists(bcpath):
            os.makedirs(bcpath)
        optcopy(self.scr_name + '/runbc.sh', bcpath + '/runbc.sh', bcopts)
        os.chdir(bcpath)
        os.popen('./runbc.sh').read()
        os.chdir('../..')
        return [ "read", "../../BC/%s/bc_structs.xyz" % self.dir_name ]
    
    def runrv(self):
        rvpath = './RV/%s' % self.dir_name
        rvopts = { "@GNAME": self.name, "@GNAMEX": elem_char(self.name), 
            "@GNUMBER": str(self.number) }
        if "CLUSTER_KANTERS_HOME" not in os.environ:
            raise RuntimeError('must have CLUSTER_KANTERS_HOME environment variable!')
        if not os.path.exists(rvpath):
            os.makedirs(rvpath)
        optcopy(self.scr_name + '/runrv.sh', rvpath + '/runrv.sh', rvopts)
        os.chdir(rvpath)
        os.popen('./runrv.sh').read()
        os.chdir('../..')
        return [ "read", "../../RV/%s/rv_structs.xyz" % self.dir_name ]
