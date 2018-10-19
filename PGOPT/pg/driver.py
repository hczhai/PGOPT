
from __future__ import print_function
from shutil import copyfile
import os, re, time, copy
from pg.utils import read_opts, write_opts, read_json, write_json, optcopy
from pg.hosts import get_host, get_model, host_cores, get_job_cmd
from pg.help import PGHelp
from pg.create import PGCreate
from pg.render import ScriptsRender
from subprocess import Popen, PIPE, STDOUT
import shlex

def trans_prog(prog):
    prog_dict = {"tm": "TURBOMOLE 6.6", "vasp": "VASP 5.4.1"}
    if prog in prog_dict:
        return prog_dict[prog]
    else:
        return prog

def test_xran(xran, i):
    if xran[0] == -1 and xran[1] == -1:
        return True
    elif xran[0] == -1 and xran[1] != -1:
        return i < xran[1]
    elif xran[0] != -1 and xran[1] == -1:
        return i >= xran[0]
    else:
        return i in range(xran[0], xran[1])

def add_zero(i, n=3):
    i = str(i)
    while len(i) < n: i = "0" + i
    return i

def read_bool(bs):
    return bs.lower().startswith("t") or bs == "1"

def time_span_short(time):
    time = int(time)
    xmin = (time / 60) % 60
    xhours = time / 60 / 60
    return "%d:%02d" % (xhours, xmin)

def cast_value(orig, x):
    if isinstance(orig, bool): return x in ["True", "true", "T", "t", "1"]
    elif isinstance(orig, float): return float(x)
    elif isinstance(orig, int): return int(x)
    elif isinstance(orig, str) or isinstance(orig, unicode): return x
    elif isinstance(orig, list):
        if len(orig) < len(x.split(",")):
            orig = [''] * len(x.split(","))
        return [cast_value(orig[ix], x) for ix, x in enumerate(x.split(","))]
    else: raise RuntimeError("Unsupport data type: %s" % orig.__class__)

# solve miltiplicity expr: 1-5,9 => [ 1, 3, 5, 9 ]
def solve_multi(ms):
    mk = []
    for l in ms.split(","):
        if '-' not in l: mk.append(int(l))
        else:
            ll = l.split('-')
            mk += range(int(ll[0]), int(ll[1]) + 1, 2)
    return mk

# read current stage and m from json
def read_stage_m(jdict):
    od = jdict["output_dir"]
    odl = od.split("/")[-1].split(".")
    return odl[0], int(old[1])

# solve source configurations expression
def solve_ref(l, stage, m):
    # format: [stage[.multi]-]name.idx[:max_number]
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
    lx, ly = l.split('.')
    if lx == "filter":
        xdir = "../../master/%s.%s/fil_structs.xyz.%s" % (lsta, lmulti, ly)
    else:
        xdir = "../../master/%s.%s/par_%s.xyz.%s" % (lsta, lmulti, lx, ly)
    if lnumber is None: return [ xdir ]
    else: return [ xdir, lnumber ]

def run_cmd(cmd):
    shell = None
    try:
        shell = Popen(shlex.split(cmd), shell=False, universal_newlines=True, \
            bufsize=1, stdout=PIPE, stderr=STDOUT)
        while shell.poll() is None:
            line = shell.stdout.readline().rstrip()
            if line == "":
                continue
            print (line)
        shell = None
    finally:
        if shell is not None:
            shell.kill()

class BaseDriver(object):
    def __init__(self, argv):
        self.this = argv[0]
        self.args = argv[1:]
        if len(self.args) == 0:
            self.args = ["help"]
        self.task, self.args = self.args[0], self.args[1:]
        if "PGOPTHOME" not in os.environ:
            raise RuntimeError("must set PGOPTHOME environ variable!")
        self.scripts_dir = os.environ["PGOPTHOME"] + "/scripts"
        self.scripts_templates_dir = os.environ["PGOPTHOME"] + "/scripts.templates"
        self.scripts_spec_dir = os.environ["PGOPTHOME"] + "/scripts.spec"
        self.surface_dir = os.environ["PGOPTHOME"] + "/surfaces"
        self.tasks = []
        if "PROJECT_NAME" in os.environ:
            self.project_name = os.environ["PROJECT_NAME"]
        else:
            self.project_name = ""
    
    def run(self):
        if "help" in self.task or self.task == "-h":
            self.help(self.args)
        else:
            found = False
            for task in self.tasks:
                if self.task == task:
                    getattr(self, task)(self.args)
                    found = True
                    break
            if not found:
                raise RuntimeError("Unknown task name: %s" % self.task)

    def help(self, args):
        pass

class PGDriver(BaseDriver):
    def __init__(self, argv):
        super(PGDriver, self).__init__(argv)
        self.tasks = ["init", "master", "freq", "relax", "energy", "submit", "torun", "sync",
                      "show", "depend", "report", "log", "check", "final", "create",
                      "get", "set", "filter", "para", "surfgen", "mclog", "neb",
                      "connect", "select", "align", "fdlog", "transfer", "background", "debug",
                      "clean", "enlarge", "draw"]

    def run(self):
        super(PGDriver, self).run()
        if "help" in self.task or self.task == "-h": return
        if self.task in [ "sync", "show", "check", "clean" ]: return
        # pre = self.pre_info()
        self.to_dir(dox="local")
        if os.path.isfile("./CMD-HISTORY"):
            ftime = os.path.getmtime("./CMD-HISTORY")
        else:
            ftime = None
        with open("./CMD-HISTORY", "a") as f:
            if ftime is not None and time.time() - ftime > 3600:
                f.write("# AFTER %d HOURS\n" % int((time.time() - ftime) / 3600))
            f.write("pgopt %s %s\n" % (self.task, " ".join([
                a if ";" not in a and "(" not in a and "=" not in a and "|" not in a else
                ("\"%s\"" % a) for a in self.args])))

    # get local and remote dirs
    # if cur, also return current root dir
    def lr_dirs(self, cur=False):
        ndir = None
        jdir = None
        if os.path.isfile("../../DIRECTORIES"):
            ndir = "../../DIRECTORIES"
            if cur: jdir = "../.."
        elif os.path.isfile("./DIRECTORIES"):
            ndir = "./DIRECTORIES"
            if cur: jdir = "."
        elif os.path.isfile("../DIRECTORIES"):
            ndir = "../DIRECTORIES"
            if cur: jdir = ".."
        else:
            return None
        d = [x for x in open(ndir, 'r').read().split('\n') if len(x) != 0][:2]
        if cur:
            jdir = os.path.abspath(jdir)
            return d + [jdir]
        else:
            return d

    def to_dir(self, dox=None):
        lr = self.lr_dirs()
        if lr == None:
            return False
        else:
            ldir, rdir = lr
            if dox == "local": os.chdir(ldir)
            elif dox == "remote": os.chdir(rdir)
            return True

    def read_arg_opts(self, orig, expr):
        lppre = orig.split(";")
        laarg = {x.split('=')[0].split("(")[0]: x for x in expr.split(";")}
        lpxre = []
        for il, l in enumerate(lppre):
            xx = l.split('=')[0].split("(")[0]
            if xx in laarg.keys():
                lppre[il] = laarg[xx]
                del laarg[xx]
        return ";".join(lppre + laarg.values())

    def get(self, args):
        pre = self.pre_info()
        arg_idx = 0
        if args[0].isdigit():
            arg_idx = int(args[0])
            args = args[1:]
        if args[0] == "args":
            print (pre["parallel"]["arguments"][arg_idx][args[1]])
        else:
            ppre = pre
            for arg in args:
                ppre = ppre[arg]
            print (ppre)

    def set(self, args):
        pre = self.pre_info()
        self.to_dir(dox="local")
        arg_idx = 0
        if args[0].isdigit():
            arg_idx = int(args[0])
            args = args[1:]
        args_temp = read_json(self.scripts_dir + "/args-template.json")
        if args[0] == "args":
            odv = pre["parallel"]["arguments"][arg_idx][args[1]]
            carg = cast_value(odv, args[2])
            pre["parallel"]["arguments"][arg_idx][args[1]] = carg
            print ("[args-%d] %s = %s => %s" % (arg_idx, args[1], str(odv), str(carg)))
        elif args[0] == "molecules+":
            if "molecules" not in pre:
                pre["molecules"] = {}
            pre["molecules"][args[1]] = "../../%s.xyz" % args[1]
            print ("[molecules+] = %s" % pre["molecules"].keys())
            lr = self.lr_dirs()
            if lr is not None and os.path.isfile("./%s.xyz" % args[1]):
                copyfile("./%s.xyz" % args[1], "%s/%s.xyz" % (lr[1], args[1]))
        elif args[0] == "molecules-":
            if "molecules" in pre and args[1] in pre["molecules"]:
                del pre["molecules"][args[1]]
            if len(pre["molecules"]) == 0: del pre["molecules"]
            print ("[molecules-] = %s" % pre.get("molecules", {}).keys())
        elif args[0].startswith("do-"):
            ppre = pre["parallel"]["arguments"][arg_idx]
            ppre[args[0]] = read_bool(args[1])
            print ("[%s %d] = %s" % (args[0], arg_idx, str(ppre[args[0]])))
        elif args[0] == "opts+":
            ppre = pre["parallel"]["arguments"][arg_idx]
            ppre["others"] += ";" + args[1]
            print ("[opts+%d] = %s" % (arg_idx, ppre["others"]))
        elif args[0] == "opts-":
            ppre = pre["parallel"]["arguments"][arg_idx]
            lppre = ppre["others"].split(";")
            ppre["others"] = ";".join([l for l in lppre if l not in args[1].split(";")])
            print ("[opts-%d] = %s" % (arg_idx, ppre["others"]))
        elif args[0] == "opts^":
            ppre = pre["parallel"]["arguments"][arg_idx]
            ppre["others"] = self.read_arg_opts(ppre["others"], args[1])
            print ("[opts^%d] = %s" % (arg_idx, ppre["others"]))
        elif args[0] == "sources+" or args[0] == "refs+":
            spname = args[0][:-1]
            ppre = pre["parallel"]["arguments"][arg_idx]
            spre = args[1].split(",")
            if spre[-1].isdigit(): spre[-1] = int(spre[-1])
            if spre[0] == "inherit": spre[1:] = [ spre[1:] ]
            if spre[0] == "read":
                if '/' not in spre[1]:
                    spre[1] = "../../" + spre[1]
                self.read_structs(spre[1])
            ppre[spname].append(spre)
            print ("[%s+%d] = %s" % (spname, arg_idx, ppre[spname]))
        elif args[0] == "sources-" or args[0] == "refs-":
            spname = args[0][:-1]
            ppre = pre["parallel"]["arguments"][arg_idx]
            spre = args[1].split(",")
            if spre[-1].isdigit(): spre[-1] = int(spre[-1])
            if spre[0] == "inherit": spre[1:] = [ spre[1:] ]
            if spre in ppre[spname]: ppre[spname].remove(spre)
            print ("[%s-%d] = %s" % (spname, arg_idx, ppre[spname]))
        elif args[0] == "surf":
            ppre = pre["parallel"]["arguments"][arg_idx]
            if args[1] == "none":
                if "surface" in ppre:
                    del ppre["surface"]
            else:
                ppre["surface"] = {}
                self.read_surface(args[1], ppre["surface"])
                ppre["others"] = self.read_arg_opts(ppre["others"], ppre["surface"]["others"])
            print ("[surf%d] = %s" % (arg_idx, ppre.get("surface", {"surface": ""})["surface"]))
            print ("[opts^%d] = %s" % (arg_idx, ppre["others"]))
        elif args[0] == "sweep":
            ppre = pre["parallel"]["arguments"][arg_idx]
            if len(args) == 2 and args[1] == "none":
                if "sweep" in ppre:
                    del ppre["sweep"]
            elif len(args) == 1:
                if "sweep" not in ppre:
                    ppre["sweep"] = {"n": 11, "start": 0.0, "end": 1.0,
                        "name": "efield", "copy": True}
            else:
                odv = ppre["sweep"][args[1]]
                carg = cast_value(odv, args[2])
                ppre["sweep"][args[1]] = carg
            print ("[sweep%d] = %s" % (arg_idx, str(ppre.get("sweep", None))))
        elif args[0] == "mc":
            if len(args) == 2 and args[1] == "none":
                if "monte-carlo" in pre:
                    del pre["monte-carlo"]
            elif len(args) == 1:
                pre["monte-carlo"] = { 
                    "max-iter": 200, "step-length": 0.4, "temperature": 300, "swap-site-rate": 1.0, 
                    "short-distance-factor": 0.7, "refuse-separate": True, "light-shell": False,
                    "detailed-balance": False, "light-step-length": 0.6, "final-factor": 1.2, 
                    "swap-site": False, "keep-ch3": False, "swap-site-make-space": [],
                    "keep-co": False, "solid-move": []
                }
            elif len(args) == 3:
                odv = pre["monte-carlo"][args[1]]
                carg = cast_value(odv, args[2])
                pre["monte-carlo"][args[1]] = carg
                print ("[mc] %s = %s => %s" % (args[1], str(odv), str(carg)))
            else:
                print ("[mc] = %s" % pre.get("monte-carlo", {}))
        elif args[0] == "fd":
            if len(args) == 2 and args[1] == "none":
                if "finite-diff" in pre:
                    del pre["finite-diff"]
            elif len(args) == 1:
                pre["finite-diff"] = { "step-length": 0.01 }
            elif len(args) == 3:
                odv = pre["finite-diff"][args[1]]
                carg = cast_value(odv, args[2])
                pre["finite-diff"][args[1]] = carg
                print ("[fd] %s = %s => %s" % (args[1], str(odv), str(carg)))
            else:
                print ("[fd] = %s" % pre.get("finite-diff", {}))
        elif args[0] == "sources^" or args[0] == "refs^":
            spname = args[0][:-1]
            ppre = pre["parallel"]["arguments"][arg_idx]
            spre = args[1].split(",")
            if spre[-1].isdigit(): spre[-1] = int(spre[-1])
            if spre[0] == "inherit": spre[1:] = [ spre[1:] ]
            if spre[0] == "read":
                if '/' not in spre[1]:
                    spre[1] = "../../" + spre[1]
                self.read_structs(spre[1])
            found = False
            for src in ppre[spname]:
                if src[0] == spre[0]:
                    src[1:] = spre[1:]
                    found = True
                    break
            if not found:
                ppre[spname].append(spre)
            print ("[%s^%d] = %s" % (spname, arg_idx, ppre[spname]))
        elif args[0] in args_temp.keys():
            ppre = pre["parallel"]["arguments"]
            aargs = args_temp[args[0]]
            if "creation-surface" in pre:
                for aa in aargs:
                    aa["others"] = self.read_arg_opts(aa["others"],
                                                      pre["creation-surface"]["others"])
            pre["parallel"]["arguments"] = ppre[:arg_idx] + aargs + ppre[arg_idx:]
            print ("[args-%d] created %s" % (arg_idx, args[0]))
        elif args[0] == "none":
            del pre["parallel"]["arguments"][arg_idx]
            print ("[args-%d] %s" % (arg_idx, args[0]))
        else:
            ppre = pre
            for arg in args[:-2]:
                ppre = ppre[arg]
            if isinstance(ppre, list):
                odv = ppre[int(args[-2])]
                ppre[int(args[-2])] = cast_value(ppre[int(args[-2])], args[-1])
            else:
                odv = ppre[args[-2]]
                ppre[args[-2]] = cast_value(ppre[args[-2]], args[-1])
            print ("[%s] %s = %s => %s" % ("/".join(args[:-2]), args[-2], str(odv), args[-1]))
        write_json(pre, "./para-template.json")

    def create(self, args):
        pre = self.pre_info()
        def_pos = { "0": "stage", "1": "multi", "2": "number" }
        opts = { "number": "5000" }
        optl = [ "rseed" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "stage", "multi" ]:
            if k not in opts:
                raise RuntimeError("no %s argument found!" % k)
        
        # multi and random seed list solve
        mk = solve_multi(opts["multi"])
        rs = None
        if "rseed" in opts:
            rs = [int(x) for x in opts["rseed"].split(",")]
            while len(rs) < len(mk): rs += [ rs[-1] ]
        
        msdir = "./master"
        if not os.path.exists(msdir): 
            print ("create %s" % msdir)
            os.makedirs(msdir)
        
        pre["creation"]["number"] = int(opts["number"])
        pre["tasks"] = [ "create" ]
        del pre["parallel"]
        # for each multiplicity
        for im, m in enumerate(mk):
            if rs is None:
                pre["random_seed"] = int(time.time() * 99991) % (2**32 - 1)
            else:
                pre["random_seed"] = rs[im]

            # sources and refs
            xdir_name = "C%s.%d" % (opts["stage"], m)
            pre["output_dir"] = "../../master/%s" % xdir_name

            pgc = PGCreate(pre["creation"], pre["random_seed"], xdir_name, 
                self.scripts_dir)
            pgc.run()

            tmdir = "./ctomaster/C-%s.%d" % (opts["stage"], m)
            print ("create %s" % tmdir)
            if not os.path.exists(tmdir): os.makedirs(tmdir)
            os.chdir(tmdir)
            # para.json
            write_json(pre, "./create.json")
            os.chdir("../..")
    
    def read_structs(self, name):
        lr = self.lr_dirs()
        if lr is None:
            path_pwd, path_full = None, None
        else:
            path_pwd, path_full = lr
        xname = name.replace("../../", "./")
        if os.path.isfile(xname):
            if path_pwd is not None:
                yname = name.replace("../../", "%s/" % path_full)
                copyfile(xname, yname)

    def read_surface(self, name, d):
        lr = self.lr_dirs()
        if lr is None:
            path_pwd, path_full = None, None
        else:
            path_pwd, path_full = lr
        d["surface"] = "../../%s.xyz" % name
        if "fix" in d: del d["fix"]
        if os.path.isfile(('./%s.json') % name):
            ccs = read_json(('./%s.json') % name)
            if path_pwd is not None:
                copyfile(('./%s.xyz') % name, '%s/%s.xyz' % (path_full, name))
        elif os.path.isfile((self.surface_dir + '/%s.json') % name):
            ccs = read_json((self.surface_dir + '/%s.json') % name)
            if path_pwd is not None:
                copyfile((self.surface_dir + '/%s.xyz') % name, \
                    '%s/%s.xyz' % (path_full, name))
            else:
                copyfile((self.surface_dir + '/%s.xyz') % name, './%s.xyz' % name)
        self.to_dir(dox="remote")
        d.update(ccs)
        # para_temp["parallel"]["filter_type"] = "none"
        sff = open('./%s.xyz' % name, 'r').readlines()
        surfx = sff[1].split(' ')
        if "fix" in ccs and ccs["fix"] != "":
            selm = ccs["fix"]
        else:
            selm = list(set([x.strip().split(' ')[0] for x in sff[2:]]))
            selm = ",".join(selm)
        if "(" in surfx[-1]:
            surfx = surfx[:-1]
        scell = ["%.4f" % float(y) for y in [x for x in surfx if len(x) != 0][2:]]
        d.update({"others": "cell=%s;scf(iter=100);fix(%s);nocenter" % \
            (":".join(scell), selm)})
        self.to_dir(dox="local")

    def init(self, args):
        self.to_dir(dox="local")
        opts = {"number": 200, "method": "blda", "charge": 0, "nodes": "1",
                "func": "tpssh", "basis": "def2-SV(P)", "program": "tm",
                "surface": ""}
        optl = ["no-scratch", "model", "cores", "encut", "periodic"] + opts.keys()
        if os.path.isfile("./para-template.json"):
            para_temp = read_json("./para-template.json")
            opts = {"name": para_temp["creation"]["name"]}
        else:
            para_temp = read_json(self.scripts_dir + "/para-template.json")
        def_pos = {"0": "name", "1": "number", "2": "method"}
        optst = read_opts(args, def_pos, optl)
        if "program" in optst and optst["program"] == "vasp":
            opts["basis"] = ""
            opts["func"] = "PBE"
            opts["others"] = "cell=15;scf(iter=300)"
        opts.update(optst)
        for k in [ "number", "charge", "nodes" ]:
            if k in opts: opts[k] = int(opts[k])
        if "name" not in opts:
            raise RuntimeError("name argument must be set!")
        path_pwd = None
        if "no-scratch" not in opts and not os.path.isfile("./DIRECTORIES"):
            path_pwd = os.path.abspath(os.curdir)
            xname = opts["name"].replace("(", "").replace(")", "").replace("|", "")
            path_remote = "%s-%d-%s" % (xname, opts["number"], opts["method"])
            path_id = 0
            if "SCRATCH" in os.environ:
                path_scr = os.environ["SCRATCH"]
            elif "WORKDIR" in os.environ:
                path_scr = os.environ["WORKDIR"]
            else:
                raise RuntimeError("SCRATCH/WORKDIR directory not found in environ!")
            while os.path.exists("%s/%s.%d" % (path_scr, path_remote, path_id)):
                path_id += 1
            path_full = "%s/%s.%d" % (path_scr, path_remote, path_id)
            os.makedirs(path_full)
            f = open('./DIRECTORIES', 'w')
            f.write('%s\n%s\n' % (path_pwd, path_full))
            f.close()
            copyfile('./DIRECTORIES', path_full + "/DIRECTORIES")
        cre_opts = { k: opts[k] for k in ["name", "number", "method"] if k in opts }
        if "method" in cre_opts:
            if cre_opts["method"] == "blda":
                cre_opts["default_sigma"] = 0.1
                cre_opts["2d"] = 0.0
                cre_opts["order"] = 2
                cre_opts["morder"] = 1
                cre_opts["~lowpos"] = -999.0
                cre_opts["~hydro"] = False
                cre_opts["~rmhydro"] = 3
                cre_opts["~nocluselems"] = []
                cre_opts["~nomoleelems"] = []
                cre_opts["~loworderelems"] = []
                # ~sampling = random/all/grad
                # all: for absorbates. generate <given number> of init structs for each isomers in the mole file.
                #      in total <given number> * <number of isomers> generated
                # random: for absorbates. in total <given number> generated. each time select random isomer
                # grid.<m>.<n>: for PES sampling. must use xyfix in DFT parameter. no mole.
                cre_opts["sampling"] = "random"
                if "drel" in para_temp["creation"]:
                    del para_temp["creation"]["drel"]
            elif cre_opts["method"] == "ck":
                cre_opts["drel"] = 0.02
                for p in ["default_sigma", "order", "2d", "morder"]:
                    if p in para_temp["creation"]:
                        del para_temp["creation"][p]
        para_temp["creation"].update(cre_opts)
        if "creation-surface" not in para_temp:
            if "surface" in opts and opts["surface"] != "":
                para_temp["creation-surface"] = {}
                para_temp["creation-surface"]["mix_rate"] = 0.5
                para_temp["creation-surface"]["ncore"] = 4
        else:
            if "surface" in opts and opts["surface"] == "":
                del para_temp["creation-surface"]
        if "creation-periodic" not in para_temp:
            if "periodic" in opts:
                para_temp["creation-periodic"] = {}
                para_temp["creation-periodic"]["max_height"] = 3.5
                para_temp["creation-periodic"]["start_height"] = 1.0
                para_temp["creation-periodic"]["layer_height"] = 0.1
        else:
            if "periodic" not in opts:
                del para_temp["creation-periodic"]
        if "surface" in opts and opts["surface"] != "":
            self.read_surface(opts["surface"], para_temp["creation-surface"])
        xhost = get_host()
        if xhost != "hoffman":
            if "PROJECT_NAME" not in os.environ:
                raise RuntimeError("must set PROJECT_NAME environ variable!")
            self.project_name = os.environ["PROJECT_NAME"]
        if "model" in opts:
            xmodel = opts["model"]
        elif "hosts" in para_temp:
            xmodel = para_temp["hosts"]["model"]
        else:
            xmodel = get_model(xhost)
        if "cores" in opts:
            xcores = int(opts["cores"])
        elif "hosts" in para_temp:
            xcores = para_temp["hosts"]["cores"]
        else:
            xcores = host_cores[xhost]
        if "nodes" in opts:
            xnodes = int(opts["nodes"])
        elif "nodes" in para_temp:
            xnodes = para_temp["hosts"]["nodes"]
        else:
            xnodes = 1
        para_temp["hosts"] = { "name": xhost, "cores": xcores, "model": xmodel,
            "nodes": xnodes }
        write_json(para_temp, "./para-template.json")
    
    def pre_info(self):
        self.to_dir(dox="local")
        if os.path.isfile("./para-template.json"):
            para_temp = read_json("./para-template.json")
        else:
            raise RuntimeError("please run 'pgopt init ...' first!")
        xmodel = para_temp["hosts"]["model"]
        xhost = para_temp["hosts"]["name"]
        self.scripts_render = ScriptsRender(self.scripts_dir, self.scripts_templates_dir, 
                                            self.scripts_spec_dir, xhost, xmodel)
        return para_temp
    
    def submit(self, args):
        return self.submit_or_depend(args, "submit")
    
    def depend(self, args):
        return self.submit_or_depend(args, "depend")
    
    def transfer(self, args):
        return self.submit_or_depend(args, "transfer")
    
    def debug(self, args):
        return self.submit_or_depend(args, "debug")

    def background(self, args):
        return self.submit_or_depend(args, "background")

    def clean(self, args):
        cmd = "mv CMD-HISTORY CMD.bak"
        os.popen(cmd).read()
        cmd = "rm -r CMD-HISTORY DIRECTORIES para-template.json master tomaster"
        os.popen(cmd).read()

    # submit torun/run-batch.sh/run-master.sh
    def submit_or_depend(self, args, task):
        pre = self.pre_info()
        xmodel = pre["hosts"]["model"]
        jcmd = get_job_cmd(xmodel)
        self.to_dir(dox="remote")
        if len(args) == 0:
            raise RuntimeError("Need argument torun/relax/freq/energy/etc!")
        if args[0] == "torun":
            if task == "depend":
                if xmodel not in [ "idataplex", "cray", "bridges", "kalk" ]:
                    raise RuntimeError("depend not implemented for model %s" % xmodel)
            elif task == "transfer" or task == "background" or task == "debug":
                if xmodel not in [ "idataplex", "cray", "kalk" ]:
                    raise RuntimeError("transfer/background/debug not implemented for model %s" % xmodel)
            def_pos = { "0": "runt", "1": "from", "2": "to" }
            opts = { "from": "-1", "to": "-1", "email": "" }
            if task == "depend":
                opts.update({ "jobid": "" })
            optl = opts.keys() + [ 'autodep' ]
            opts.update(read_opts(args[1:], def_pos, optl))
            for k in [ "from", "to" ]:
                if k in opts: opts[k] = int(opts[k])
            if "runt" not in opts or opts["runt"] not in [ "para", "seq" ]:
                raise RuntimeError("no runt=para/seq argument found!")
            os.chdir('./torun/%s' % opts["runt"])
            xran = (opts["from"], opts["to"])
            for l in sorted(os.listdir('.')):
                if l.endswith('.sh'):
                    idx = int(re.findall(r'\-([0-9]+)\.', l)[0])
                    hx = False
                    if idx == 0 and xran[0] != -1 and xran[1] != -1:
                        cmd = "grep '#$ -t ' %s | wc -l" % l
                        if os.popen(cmd).read().strip() == "1":
                            cmd = "sed -i '/#$ -t /c\\#$ -t %d-%d:1' %s" % (xran[0] + 1, xran[1], l)
                            os.popen(cmd).read()
                            hx = True
                    if hx or test_xran(xran, idx):
                        print ("%s %s/%s" % (task, opts["runt"], l))
                        if opts["email"] != "":
                            cmd = "sed -i '/#PBS -m/d' %s" % l
                            os.popen(cmd).read()
                            cmd = "sed -i '/#PBS -M/d' %s" % l
                            os.popen(cmd).read()
                            if opts["email"] != "none":
                                cmd = "sed -i '5 i\\#PBS -M %s' %s" % (opts["email"], l)
                                os.popen(cmd).read()
                                cmd = "sed -i '5 i\\#PBS -m e' %s" % l
                                os.popen(cmd).read()
                        if task == "submit":
                            cmd = '%s \"%s\"' % (jcmd[0], l)
                            print (os.popen(cmd).read().strip())
                        elif task == "background" or task == "debug":
                            if xmodel in [ "idataplex", "cray", "kalk" ]:
                                cmd = "sed -i '/#PBS -q/d' %s" % l
                            else:
                                raise RuntimeError("depend not implemented for model %s" % xmodel)
                            os.popen(cmd).read()
                            cmd = "sed -i '5 i\\#PBS -q %s' %s" % (task, l)
                            os.popen(cmd).read()
                            print ("Q-%s " % task.upper() + l + " -> %s" % task)
                        elif task == "transfer":
                            if xmodel in [ "idataplex", "cray", "kalk" ]:
                                cmd = "sed -i '/#PBS -q/d' %s" % l
                            else:
                                raise RuntimeError("depend not implemented for model %s" % xmodel)
                            os.popen(cmd).read()
                            cmd = "sed -i '5 i\\#PBS -q transfer' %s" % l
                            os.popen(cmd).read()
                            if xmodel in [ "cray", "idataplex" ]:
                                cmd = "grep 'PBS -l select=' %s" % l
                                dd = os.popen(cmd).read().strip().split("=")
                                dd = "%s=%s=1" % (dd[0], dd[1])
                                cmd = "sed -i '/#PBS -l select=/d' %s" % l
                                os.popen(cmd).read()
                                cmd = "sed -i '4 i\\%s' %s" % (dd, l)
                                os.popen(cmd).read()
                                cmd = "sed -i '/module load ccm/d' %s" % l
                                os.popen(cmd).read()
                                cmd = "sed -i 's/ccmrun //g' %s" % l
                                os.popen(cmd).read()
                                cmd = "sed -i '/#PBS -l ccm/d' %s" % l
                                os.popen(cmd).read()
                            print ("Q-TRANSFER " + l + " -> transfer")
                        else:
                            if xmodel in [ "idataplex", "cray", "kalk" ]:
                                cmd = "sed -i '/#PBS -W/d' %s" % l
                            elif xmodel in [ "bridges" ]:
                                cmd = "sed -i '/#SBATCH -d/d' %s" % l
                            else:
                                raise RuntimeError("depend not implemented for model %s" % xmodel)
                            os.popen(cmd).read()
                            if "autodep" in opts:
                                cmd = "qstat -u `whoami` | grep '" + l[:-3] + "' | awk '{print $1}' | tail -1 | cut -d . -f 1"
                                dpid = os.popen(cmd).read().strip()
                                if dpid != "":
                                    if xmodel in [ "idataplex", "cray", "kalk" ]:
                                        cmd = "sed -i '3 i\\#PBS -W depend=afterany:%s' %s" % (dpid, l)
                                    elif xmodel in [ "bridges" ]:
                                        cmd = "sed -i '3 i\\#SBATCH -d afterany:%s' %s" % (dpid, l)
                                    print ("AUTO DEP " + l + " -> " + dpid)
                                else:
                                    print ("AUTO DEP " + l + " -> failed")
                                os.popen(cmd).read()
                            elif opts["jobid"] != "":
                                if xmodel in [ "idataplex", "cray", "kalk" ]:
                                    cmd = "sed -i '3 i\\#PBS -W depend=after:%s' %s" % (opts["jobid"], l)
                                elif xmodel in [ "bridges" ]:
                                    cmd = "sed -i '3 i\\#SBATCH -d after:%s' %s" % (opts["jobid"], l)
                                os.popen(cmd).read()
            os.chdir("../..")
        else:
            def_pos = { "0": "runt", "1": "stage", "2": "multi" }
            if task == "submit": 
                opts = { "email": "" }
                optl = []
            else:
                opts = { "jobid": "", "email": "" }
                optl = opts.keys() + [ 'autodep' ]
            opts.update(read_opts(args, def_pos, optl))
            if "stage" not in opts:
                raise RuntimeError("no stage argument found!")
            # find the most recent batch
            if "multi" not in opts:
                mf, mt, mc = None, None, None
                mfx = "run-batch.sh"
                for l in sorted(os.listdir('./tomaster')):
                    if l.startswith('B-%s-%s' % (opts["runt"], opts["stage"])):
                        if mt is None or mt < os.path.getmtime("./tomaster/%s" % l):
                            mf = "./tomaster/%s" % l
                            mt = os.path.getmtime("./tomaster/%s" % l)
                    elif l == ('%s-%s.0' % (opts["runt"], opts["stage"])):
                        mc = "./tomaster/%s" % l
                if mf is None and mc is None:
                    raise RuntimeError("No run-batch.sh script for %s-%s to submit!" % 
                        (opts["runt"], opts["stage"]))
                if mf is None and mc is not None:
                    mf = mc
                    mfx = "run-master.sh"
            else:
                mk = solve_multi(opts["multi"])
                if len(mk) == 1:
                    mf = "./tomaster/%s-%s.%d" % (opts["runt"], opts["stage"], mk[0])
                    mfx = "run-master.sh"
                else:
                    mf = "./tomaster/B-%s-%s-%s" % (opts["runt"], opts["stage"], opts["multi"])
                    mfx = "run-batch.sh"
            os.chdir(mf)
            if opts["email"] != "":
                l = mfx
                cmd = "sed -i '/#PBS -m/d' %s" % l
                os.popen(cmd).read()
                cmd = "sed -i '/#PBS -M/d' %s" % l
                os.popen(cmd).read()
                if opts["email"] != "none":
                    cmd = "sed -i '5 i\\#PBS -M %s' %s" % (opts["email"], l)
                    os.popen(cmd).read()
                    cmd = "sed -i '5 i\\#PBS -m e' %s" % l
                    os.popen(cmd).read()
            if task == "submit":
                cmd = '%s %s' % (jcmd[0], mfx)
                print ("%s %s/%s" % (task, mf, mfx))
                print (os.popen(cmd).read().strip())
            elif task == "background" or task == "debug":
                l = mfx
                if xmodel in [ "idataplex", "cray", "kalk" ]:
                    cmd = "sed -i '/#PBS -q/d' %s" % l
                else:
                    raise RuntimeError("depend not implemented for model %s" % xmodel)
                os.popen(cmd).read()
                cmd = "sed -i '5 i\\#PBS -q %s' %s" % (task, l)
                os.popen(cmd).read()
                print ("Q-%s " % task.upper() + l + " -> %s" % task)
            elif task == "transfer":
                l = mfx
                if xmodel in [ "idataplex", "cray", "kalk" ]:
                    cmd = "sed -i '/#PBS -q/d' %s" % l
                else:
                    raise RuntimeError("depend not implemented for model %s" % xmodel)
                os.popen(cmd).read()
                cmd = "sed -i '5 i\\#PBS -q transfer' %s" % l
                os.popen(cmd).read()
                if xmodel in [ "cray" ]:
                    cmd = "grep 'PBS -l select=' %s" % l
                    dd = os.popen(cmd).read().strip().split("=")
                    dd = "%s=%s=1" % (dd[0], dd[1])
                    cmd = "sed -i '/#PBS -l select=/d' %s" % l
                    os.popen(cmd).read()
                    cmd = "sed -i '4 i\\%s' %s" % (dd, l)
                    os.popen(cmd).read()
                    cmd = "sed -i '/module load ccm/d' %s" % l
                    os.popen(cmd).read()
                    cmd = "sed -i 's/ccmrun //g' %s" % l
                    os.popen(cmd).read()
                    cmd = "sed -i '/#PBS -l ccm/d' %s" % l
                    os.popen(cmd).read()
                elif xmodel in["idataplex"]:
                    cmd = "grep 'PBS -l select=' %s" % l
                    dd = os.popen(cmd).read().strip().split("=")
                    dd = "%s=%s=1" % (dd[0], dd[1])
                    cmd = "sed -i '/#PBS -l select=/d' %s" % l
                    os.popen(cmd).read()
                    cmd = "sed -i '4 i\\%s' %s" % (dd, l)
                    os.popen(cmd).read()
                print ("Q-TRANSFER " + l + " -> transfer")
            else:
                if xmodel in [ "idataplex", "cray", "kalk" ]:
                    cmd = "sed -i '/#PBS -W/d' %s" % mfx
                elif xmodel in [ "bridges" ]:
                    cmd = "sed -i '/#SBATCH -d/d' %s" % mfx
                elif xmodel in [ "cascade" ]:
                    cmd = "sed -i '/#MSUB -l depend=/d' %s" % mfx
                elif xmodel in [ "mira" ]:
                    cmd = "sed -i '/#COBALT --dependencies/d' %s" % mfx
                else:
                    cmd = "sed -i '/#$ -hold_jid/d' %s" % mfx
                print ("%s %s/%s" % (task, mf, mfx))
                os.popen(cmd).read()
                if "autodep" in opts:
                    xname = pre["creation"]["name"].upper().replace("(", "").replace(")", "").replace("|", "")
                    xmodel = pre["hosts"]["model"]
                    if xmodel in [ "idataplex", "cray", "kalk" ] and len(xname) > 8:
                        xname = xname[:8]
                    for i in range(7, 1, -1):
                        cmd = "qstat -u `whoami` | grep '%s' | awk '{print $1}' | tail -1 | cut -d . -f 1" % \
                            (xname + '-MASTER'[:i])
                        dpid = os.popen(cmd).read().strip()
                        if dpid != "":
                            print("AUTO MASTER DEP -> " + dpid)
                            break
                    if dpid == "":
                        print("AUTO MASTER DEP -> failed")
                    opts["jobid"] = dpid
                if opts["jobid"] != "":
                    if xmodel in [ "idataplex", "cray", "kalk" ]:
                        cmd = "sed -i '3 i\\#PBS -W depend=afterany:%s' %s" % (opts["jobid"], mfx)
                    elif xmodel in [ "hoffman" ]:
                        cmd = "sed -i '3 i\\#$ -hold_jid %s' %s" % (opts["jobid"], mfx)
                    elif xmodel in [ "bridges" ]:
                        cmd = "sed -i '3 i\\#SBATCH -d afterany:%s' %s" % (opts["jobid"], mfx)
                    elif xmodel in [ "cascade" ]:
                        cmd = "sed -i '3 i\\#MSUB -l depend=%s' %s" % (opts["jobid"], mfx)
                    elif xmodel in [ "mira" ]:
                        cmd = "sed -i '3 i\\#COBALT --dependencies %s' %s" % (opts["jobid"], mfx)
                    else:
                        raise RuntimeError("depend not implemented for model %s" % xmodel)
                    os.popen(cmd).read()
            os.chdir("../..")
    
    def align(self, args):
        lrc = self.lr_dirs(True)
        self.to_dir(dox="local")
        if os.path.isfile("./para-template.json"):
            pre = self.pre_info()
            if "molecules" in pre:
                del pre["molecules"]
        else:
            pre = {}
        def_pos = { "0": "sour" }
        opts = { "sour": "" }
        optl = opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        if opts["sour"] == "":
            raise RuntimeError("no %s argument found!" % "sour")
        if not opts["sour"].startswith("/") and lrc is not None:
            opts["sour"] = "%s/%s" % (lrc[2], opts["sour"])
        assert os.path.isfile("%s" % opts["sour"])
        pre["tasks"] = [ "align" ]
        pre["output_dir"] = "./OUT"
        pre["random_seed"] = 0
        pre["align"] = { "input_file": opts["sour"] }
        if "creation-surface" in pre and lrc is not None:
            if pre["creation-surface"]["surface"].startswith("../../"):
                pre["creation-surface"]["surface"] = \
                    pre["creation-surface"]["surface"].replace("../..", lrc[0])
        if lrc is not None:
            os.chdir(lrc[2])
        write_json(pre, "./align.json")
        cmd = os.environ["ACNNHOME"] + '/acnnmain align.json'
        for l in os.popen(cmd):
            print (l.rstrip())
        os.rename('./OUT/ali_structs.xyz.0', './aligned.xyz')
        os.rmdir("./OUT")
        os.remove("./align.json")

    def surfgen(self, args):
        self.to_dir(dox="local")
        if os.path.isfile("./para-template.json"):
            pre = self.pre_info()
            if "molecules" in pre:
                del pre["molecules"]
        else:
            pre = {}
        def_pos = { "0": "sour" }
        opts = { "sname": "", "name": "", "sour": "" }
        optl = [ "fix", "fixz", "gapz", "removez", "cell", "clus", "alldet", \
            "supercell", "dx", "dy", "nospg" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        if opts["sour"] == "":
            raise RuntimeError("no %s argument found!" % "sour")
        assert os.path.isfile("./%s" % opts["sour"])
        pre["tasks"] = [ "surfgen" ]
        pre["output_dir"] = "./OUT"
        pre["random_seed"] = 0
        pre["surfgen"] = { "input_file": opts["sour"] }
        if "fix" in opts:
            pre["surfgen"]["fix"] = opts["fix"]
        if "alldet" in opts:
            pre["surfgen"]["alldet"] = True
        if "nospg" in opts:
            pre["surfgen"]["no_space_group"] = True
        if "cell" in opts:
            pre["surfgen"]["cell"] = [float(x) for x in opts["cell"].split(":")]
        if "fixz" in opts:
            pre["surfgen"]["fixz"] = float(opts["fixz"])
        if "gapz" in opts:
            pre["surfgen"]["gapz"] = float(opts["gapz"])
        if "removez" in opts:
            pre["surfgen"]["removez"] = float(opts["removez"])
        if "supercell" in opts:
            pre["surfgen"]["supercell"] = opts["supercell"]
        if "clus" in opts:
            pre["creation"] = {"name": "(%s)|H" % opts["clus"]}
        if "dx" in opts:
            pre["surfgen"]["dx"] = float(opts["dx"])
        if "dy" in opts:
            pre["surfgen"]["dy"] = float(opts["dy"])
        write_json(pre, "./surfgen.json")
        cmd = os.environ["ACNNHOME"] + '/acnnmain surfgen.json'
        run_cmd(cmd)
        if opts["sname"] != "":
            os.rename('./OUT/surf.json', './%s.json' % opts["sname"])
            os.rename('./OUT/surf.xyz', './%s.xyz' % opts["sname"])
        if opts["name"] != "":
            os.rename('./OUT/orig.xyz', './%s.xyz' % opts["name"])
        if opts["sname"] != "" and opts["name"] != "":
            os.rmdir("./OUT")
        os.remove("./surfgen.json")
    
    def master(self, args, def_runt="relax"):
        pre = self.pre_info()
        self.to_dir(dox="remote")
        lr = self.lr_dirs()
        def_pos = { "0": "stage", "1": "multi" }
        opts = { "time": "24:00:00", "multi": "0", "source": "create", "ref": "", "nprev": "0", 
            "max-config": "-1", "idx": "0", "runt": def_runt, "sec": "", "blocks": "0" }
        optl = [ "sp", "batch", "copy", "rseed" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "nprev", "max-config", "idx", "blocks" ]:
            if k in opts: opts[k] = int(opts[k])
        for k in [ "stage" ]:
            if k not in opts:
                raise RuntimeError("no %s argument found!" % k)
        
        # multi and random seed list solve
        mk = solve_multi(opts["multi"])
        rs = None
        if "rseed" in opts:
            rs = [int(x) for x in opts["rseed"].split(",")]
            while len(rs) < len(mk): rs += [ rs[-1] ]
        
        # for all multiplicities
        pre["parallel"]["secondary"] = opts["sec"]
        pre["parallel"]["arguments"][0]["max_config"] = opts["max-config"]
        pre["parallel"]["arguments"][0]["run_type"] = opts["runt"]
        # pre["parallel"]["arguments"][0]["output_index"] = opts["idx"]
        if "sp" in opts:
            if len(opts["sp"]) == 0:
                pass
            else:
                spx = int(opts["sp"])
                pre["parallel"]["arguments"][0]["step"] = 2
                pre["parallel"]["arguments"][0]["max_step"] = spx
        
        msdir = "./master"
        if not os.path.exists(msdir): 
            print ("create %s" % msdir)
            os.makedirs(msdir)

        # for each multiplicity
        mpr = [ pre["output_dir"], pre["parallel"]["proc_dir"], 
            pre["parallel"]["restart_dir"] ]
        xname = pre["creation"]["name"].upper().replace("(", "").replace(")", "").replace("|", "")
        xmodel = pre["hosts"]["model"]
        if xmodel in [ "idataplex", "cray", "kalk" ] and len(xname) > 8:
            xname = xname[:8]
        if "batch" not in opts:
            for im, m in enumerate(mk):
                if rs is None:
                    pre["random_seed"] = int(time.time() * 99991) % (2**32 - 1)
                else:
                    pre["random_seed"] = rs[im]

                # sources and refs
                xdir_name = "%s.%d" % (opts["stage"], m)
                xsour = []
                for l in opts["source"].split(","):
                    if l == "":
                        continue
                    elif l == "create":
                        pgc = PGCreate(pre["creation"], pre["random_seed"], xdir_name, 
                            self.scripts_dir)
                        xsour.append(pgc.run())
                    else:
                        lx = solve_ref(l, opts["stage"], m)
                        xsour.append([ "read" ] + lx)
                if opts["idx"] == 0 and opts["nprev"] != 0 and im != 0:
                    xsour.insert(0, [ "read", "../../master/%s.%d/par_local.xyz.0" % 
                        (opts["stage"], mk[im - 1]), opts["nprev"]])
                xref = []
                for l in opts["ref"].split(","):
                    if l != "":
                        lx = solve_ref(l, opts["stage"], m)
                        xref.append(lx[0])
                if len(pre["parallel"]["arguments"][0]["sources"]) == 0:
                    pre["parallel"]["arguments"][0]["sources"] = xsour
                if len(pre["parallel"]["arguments"][0].get("refs", [])) == 0:
                    pre["parallel"]["arguments"][0]["refs"] = xref

                # directories
                sopts = { "@RUNT": opts["runt"], "@RUNS": opts["stage"], "@MULT": str(m) }
                mpx = []
                for x in mpr:
                    for k, v in sopts.items():
                        x = x.replace(k, v)
                    mpx.append(x)
                [pre["output_dir"], pre["parallel"]["proc_dir"], 
                    pre["parallel"]["restart_dir"]] = mpx
                
                for xarg in pre["parallel"]["arguments"]:
                    if "multiplicity" not in xarg or xarg["multiplicity"] == 0:
                        xarg["multiplicity"] = m
                    elif isinstance(xarg["multiplicity"], str):
                        if xarg["multiplicity"].startswith("+"):
                            xarg["multiplicity"] = m + int(xarg["multiplicity"][1:])
                        elif xarg["multiplicity"].startswith("-"):
                            xarg["multiplicity"] = m - int(xarg["multiplicity"][1:])
                            if xarg["multiplicity"] <= 0:
                                # if impossible to decrease the multiplicity, ignore that stage
                                xarg["sources"] = []

                # tomaster and restarts dir
                tmdir = "./tomaster/%s-%s.%d" % (opts["runt"], opts["stage"], m)
                rsdir = "./restarts/%s.%d" % (opts["stage"], m)
                print ("create %s" % tmdir)
                if not os.path.exists(tmdir): os.makedirs(tmdir)
                if not os.path.exists(rsdir): os.makedirs(rsdir)
                os.chdir(tmdir)
                # copy sources
                if lr is not None and 'copy' in opts:
                    for arg in pre["parallel"]["arguments"]:
                        for s in arg["sources"]:
                            if s[0] == "read" and not os.path.isfile("%s/%s" % (lr[1], s[1][6:])) and \
                                os.path.isfile("%s/%s" % (lr[0], s[1][6:])):
                                print ("copy source %s" % s[1][6:])
                                rfn = "%s/%s" % (lr[1], s[1][6:])
                                if not os.path.exists(os.path.dirname(rfn)):
                                    os.makedirs(os.path.dirname(rfn))
                                copyfile("%s/%s" % (lr[0], s[1][6:]), rfn)
                        for r in arg["refs"]:
                            if not os.path.isfile("%s/%s" % (lr[1], r[6:])) and \
                                os.path.isfile("%s/%s" % (lr[0], r[6:])):
                                print ("copy reference %s" % r[6:])
                                rfn = "%s/%s" % (lr[1], r[6:])
                                if not os.path.exists(os.path.dirname(rfn)):
                                    os.makedirs(os.path.dirname(rfn))
                                copyfile("%s/%s" % (lr[0], r[6:]), rfn)
                # filter.json
                tidx = 0
                presurf = pre.get("creation-surface", None)
                for iarg, arg in enumerate(pre["parallel"]["arguments"]):
                    if arg["task_type"] == "relax":
                        tidx = iarg
                        if "surface" in arg:
                            presurf = arg["surface"]
                fpre = { "tasks": [ "filter" ], 
                    "output_dir": "../../master/%s.%d" % (opts["stage"], m), 
                    "write_summary": True, 
                    "random_seed": 0, 
                    "filtering": {
                        "input_file": "../../master/%s.%d/par_local.xyz.%d" % (opts["stage"], m, tidx),
                        "max_diff": pre["filtering-report"]["max_diff"], 
                        "max_diff_report": pre["filtering-report"]["max_diff_report"], 
                        "sort": True, 
                        "pre_sort": True, 
                        "align": True, 
                        "para": True
                    }
                }
                if presurf is not None:
                    fpre["creation-surface"] = presurf
                write_json(fpre, "./filter.json")
                # connect.json
                cpre = { "tasks": [ "connect" ], 
                    "output_dir": "../../master/%s.%d" % (opts["stage"], m), 
                    "write_summary": True, 
                    "random_seed": 0, 
                    "reaction-path": copy.deepcopy(pre["reaction-path"])
                }
                cpre["reaction-path"]["input_file"] = "../../master/%s.%d/par_filter.xyz.%d" % (opts["stage"], m, tidx)
                if presurf is not None:
                    cpre["creation-surface"] = presurf
                write_json(cpre, "./connect.json")
                # run-master.sh
                if xmodel == "mira":
                    assert opts["blocks"] != 0
                    ta, tb, _ = opts["time"].split(":")
                    tnc = pre["hosts"]["cores"] * pre["hosts"]["nodes"]
                    assert tnc * opts["blocks"] % 16 == 0
                    assert tnc * opts["blocks"] / 16 % 512 == 0
                    print ("total number of nodes %s" % (tnc * opts["blocks"] / 16))
                    ropts = {
                        "@TIMEMIN": str(int(ta) * 60 + int(tb)), 
                        "@NPROCSTOT": str(tnc * opts["blocks"] / 16), 
                        "@RUNT": pre["parallel"]["proc_type"], 
                        "@MPIPR": str(tnc), 
                        "@BLKPR": str(opts["blocks"]),
                        "@NAME": xname + "-MASTER" }
                    optcopy(self.scripts_render.get("run-block.sh"), "./run-master.sh", ropts)
                    pre["parallel"]["max_time"] = ((int(ta) * 60 + int(tb)) - 6) * 60
                else:
                    ropts = {
                        "@NNODES": "1", 
                        "@NPROCS": str(pre["hosts"]["cores"]), 
                        "@QUEUE": "standard", 
                        "@TIME": opts["time"], 
                        "@PROJ": self.project_name, 
                        "@RUNT": pre["parallel"]["proc_type"], 
                        "@NAME": xname + "-MASTER" }
                    optcopy(self.scripts_render.get("run-master.sh"), "./run-master.sh", ropts)
                # para.json
                write_json(pre, "./para.json")
                # meta_info.txt
                f = open("./meta_info.txt", "w")
                f.write("hostname = %s\n" % pre["hosts"]["name"])
                f.write("cores = %d\n" % pre["hosts"]["cores"])
                xosname = os.popen("uname -s").read().strip()
                if xmodel == "mira":
                    cmd = "cat /proc/cpuinfo | grep 'MHz' | tail -1 | awk '{print $3}'"
                    f.write("MHz = %s\n" % os.popen(cmd).read().strip())
                    cmd = "lscpu | grep 'Model' | awk '{print $2}'"
                    f.write("CPU = %s\n" % os.popen(cmd).read().strip())
                    cmd = "lsb_release -d | awk -F : '{print $2}'"
                    f.write("os = %s\n" % os.popen(cmd).read().strip())
                elif xosname == "Darwin":
                    cmd = "sysctl -n hw.cpufrequency | xargs -I % expr % / 1000000"
                    f.write("MHz = %s\n" % os.popen(cmd).read().strip())
                    cmd = "sysctl -n machdep.cpu.brand_string"
                    f.write("CPU = %s\n" % os.popen(cmd).read().strip())
                    cmd = "echo '-productName' '-productVersion' | xargs -n 1 sw_vers | xargs"
                    f.write("os = %s\n" % os.popen(cmd).read().strip())
                else:
                    cmd = "lscpu | grep 'MHz' | awk '{print $3}'"
                    f.write("MHz = %s\n" % os.popen(cmd).read().split()[0])
                    cmd = "lscpu | grep 'Vendor' | awk '{print $3}'"
                    f.write("CPU = %s\n" % os.popen(cmd).read().strip())
                    if pre["hosts"]["model"] == "cascade":
                        cmd = "cat /etc/issue | head -n 1"
                    else:
                        cmd = "lsb_release -d | awk -F : '{print $2}'"
                    f.write("os = %s\n" % os.popen(cmd).read().strip())
                f.close()
                os.chdir("../..")
        if len(mk) != 1:
            badir = "./tomaster/B-%s-%s-%s" % (opts["runt"], opts["stage"], opts["multi"])
            print ("create %s" % badir)
            if not os.path.exists(badir):
                os.makedirs(badir)
            # run-master.sh
            ropts = { 
                "@NNODES": "1", 
                "@NPROCS": str(pre["hosts"]["cores"]), 
                "@QUEUE": "standard", 
                "@TIME": opts["time"], 
                "@RUNT": opts["runt"], 
                "@RUNS": opts["stage"],
                "@PROT": pre["parallel"]["proc_type"], 
                "@NAME": xname + "-BATCH",
                "@MULTS": ' '.join([str(m) for m in mk]), 
                "@PROJ": self.project_name
            }
            optcopy(self.scripts_render.get("run-batch.sh"), "%s/run-batch.sh" % badir, ropts)
    
    def relax(self, args):
        self.master(args, def_runt="relax")
    
    def neb(self, args):
        self.master(args, def_runt="neb")
    
    def energy(self, args):
        self.master(args, def_runt="energy")
    
    def freq(self, args):
        self.master(args, def_runt="freq")

    # create or delete processors (workers)
    # if no corresponding procs dir, will create it
    def torun(self, args):
        pre = self.pre_info()
        self.to_dir(dox="remote")
        def_pos = { "0": "runt", "1": "from", "2": "to" }
        opts = { "time": "24:00:00", "from": "-1", "to": "-1" }
        optl = [ "clean", "x", "once" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "from", "to" ]:
            if k in opts: opts[k] = int(opts[k])
        if "runt" not in opts or opts["runt"] not in [ "para", "seq" ]:
            raise RuntimeError("no runt=para/seq argument found!")
        if not os.path.exists('./torun/%s' % opts["runt"]):
            os.makedirs('./torun/%s' % opts["runt"])
        if "clean" in opts:
            os.chdir('./torun/%s' % opts["runt"])
            xran = (opts["from"], opts["to"])
            for l in sorted(os.listdir('.')):
                if l.endswith('.sh'):
                    idx = int(re.findall(r'\-([0-9]+)\.', l)[0])
                    if test_xran(xran, idx):
                        os.remove(l)
                        print ('delete torun %s' % l)
            os.chdir("../..")
        else:
            os.chdir('./torun/%s' % opts["runt"])
            if opts["from"] == -1: opts["from"] = 0
            if opts["to"] == -1: opts["to"] = 50
            xran = range(opts["from"], opts["to"])
            xname = pre["creation"]["name"].upper().replace("(", "").replace(")", "").replace("|", "")
            xmodel = pre["hosts"]["model"]
            if xmodel in [ "idataplex", "cray", "kalk" ] and len(xname) > 6:
                xname = xname[:6]
            xopts = {
                "@NNODES": str(pre["hosts"]["nodes"]), 
                "@NPROCS": str(pre["hosts"]["cores"]), 
                "@QUEUE": "standard", 
                "@TIME": opts["time"], 
                "@NAME": xname, 
                "@RUNT": opts["runt"], 
                "@PROJ": self.project_name
            }
            if "once" in opts:
                xopts["@ONCE"] = "TRUE"
            if xmodel == "hoffman":
                xopts["@NPROCS"] = str(pre["hosts"]["cores"] * pre["hosts"]["nodes"])
            if xmodel == "mira":
                pass
            elif opts["runt"] == "para":
                if "x" in opts and (xmodel == "hoffman" or xmodel == "cascade"):
                    xopts.update({
                        "@LOWER": str(opts["from"] + 1),
                        "@UPPER": str(opts["to"]), 
                        "@NAME": xname + "-P",
                        "@NNODES": str(pre["hosts"]["nodes"] * (opts["to"] - opts["from"]))
                    })
                    l = "%s-%s.sh" % (xname, add_zero(0))
                    optcopy(self.scripts_render.get("worker-x.sh"), "./%s" % l, xopts)
                    print ('create torun %s' % l)
                else:
                    for i in xran:
                        xopts.update({
                            "@X": str(i),
                            "@NAME": "%s-%s" % (xname, add_zero(i))
                        })
                        l = "%s-%s.sh" % (xname, add_zero(i))
                        optcopy(self.scripts_render.get("worker-para.sh"), "./%s" % l, xopts)
                        print ('create torun %s' % l)
            elif opts["runt"] == "seq":
                if xmodel == "hoffman":
                    xopts.update({
                        "@LOWER": str(opts["from"] + 1),
                        "@UPPER": str(opts["to"]),
                        "@NAME": xname + "-S",
                    })
                    l = "%s-%s.sh" % (xname, add_zero(0))
                    optcopy(self.scripts_render.get("worker-seq.sh"), "./%s" % l, xopts)
                    print ('create torun %s' % l)
                else:
                    nc = pre["hosts"]["cores"] * pre["hosts"]["nodes"]
                    if opts["to"] % nc != 0:
                        opts["to"] = (opts["to"] / nc + 1) * nc
                    xran = range(opts["from"] / nc, opts["to"] / nc)
                    for i in xran:
                        xopts.update({
                            "@X": str(i),
                            "@NAME": "%s-%s" % (xname, add_zero(i))
                        })
                        l = "S-%s-%s.sh" % (xname, add_zero(i))
                        optcopy(self.scripts_render.get("worker-seq.sh"), "./%s" % l, xopts)
                        print ('create torun %s' % l)
            os.chdir('../..')
            if not os.path.exists('./procs/%s' % opts["runt"]):
                os.makedirs('./procs/%s' % opts["runt"])
            os.chdir('./procs/%s' % opts["runt"])
            xran = range(opts["from"], opts["to"])
            for i in xran:
                l = "./PROC%s" % add_zero(i)
                if not os.path.exists(l):
                    os.mkdir(l)
                    os.chdir(l)
                    optcopy(self.scripts_render.get("template.in"), "./template.in", 
                        { "@NAME": xname })
                    if xmodel == 'cray':
                        optcopy(self.scripts_render.get("worker-main.ccm.sh"), "./run-all.sh", 
                            { "@NPX": "1" if opts["runt"] == "seq" else 
                            str(pre["hosts"]["cores"] * pre["hosts"]["nodes"])})
                        optcopy(self.scripts_render.get("worker-main.pre.sh"), "./worker-main.pre.sh", {})
                        optcopy(self.scripts_render.get("worker-main.post.sh"), "./worker-main.post.sh", {})
                    else:
                        optcopy(self.scripts_render.get("worker-main.sh"), "./run-all.sh", 
                            { "@NPX": "1" if opts["runt"] == "seq" else 
                            str(pre["hosts"]["cores"] * pre["hosts"]["nodes"]) })
                    print ('create proc %s' % l)
                    os.chdir("..")
            os.chdir("../..")
    
    def sync(self, args):
        lr = self.lr_dirs()
        if lr is None:
            raise RuntimeError("DIRECTORIES not found!")
        ldir, rdir = lr
        optl = [ "restarts" ]
        opts = read_opts(args, [], optl)
        sydirs = ['master', 'tomaster']
        if 'restarts' in opts: sydirs += [ 'restarts' ]
        for sd in sydirs:
            for k in os.listdir(rdir + '/' + sd):
                if not os.path.exists(ldir + '/' + sd + '/' + k):
                    os.makedirs(ldir + '/' + sd + '/' + k)
                for j in os.listdir(rdir + '/' + sd + '/' + k):
                    # if j.endswith(".sh"): continue
                    rfn = rdir + '/' + sd + '/' + k + '/' + j
                    lfn = ldir + '/' + sd + '/' + k + '/' + j
                    if "nebstructs" in j and os.path.getsize(rfn) > 2 * 1024**3:
                        continue
                    if not os.path.isfile(lfn) or os.path.getmtime(lfn) < os.path.getmtime(rfn):
                        copyfile(rfn, lfn)
        for k in os.listdir(rdir):
            if k.endswith(".xyz"):
                rfn = rdir + '/' + k
                lfn = ldir + '/' + k
                if not os.path.isfile(lfn) or os.path.getmtime(lfn) < os.path.getmtime(rfn):
                    copyfile(rfn, lfn)
    
    def show(self, args):
        lr = self.lr_dirs(cur=True)
        if lr is None:
            raise RuntimeError("DIRECTORIES not found!")
        ldir, rdir, cdir = lr
        def_pos = { "0": "dir" }
        opts = read_opts(args, def_pos, [])
        if "dir" not in opts:
            if ldir == cdir: opts["dir"] = "remote"
            else: opts["dir"] = "local"
        if opts["dir"] == "local": print (ldir)
        elif opts["dir"] == "remote": print (rdir)
    
    def final(self, args):
        pre = self.pre_info()
        k = self.to_dir(dox="local")
        if k: os.popen(os.environ["PGOPTHOME"] + '/pgopt sync').read()
        if not os.path.exists('./report'):
            os.mkdir('./report')
        base = os.path.basename(os.path.abspath(os.curdir))
        cmd = 'zip -r ./report/%s-%s.zip master tomaster *.xyz' % (pre["hosts"]["name"], base)
        print (os.popen(cmd).read().strip())
    
    def filter(self, args):
        self.filter_or_para(args, "filter")
    
    def connect(self, args):
        self.filter_or_para(args, "connect")
    
    def para(self, args):
        self.filter_or_para(args, "para")

    def filter_or_para(self, args, task="filter"):
        pre = self.pre_info()
        if self.to_dir(dox="local"):
            os.popen(os.environ["PGOPTHOME"] + '/pgopt sync').read()
        self.to_dir(dox="remote")
        lr = self.lr_dirs()
        def_pos = { "0": "stage", "1": "multi" }
        opts = { "multi": "0", "runt": "relax", "path-ref-local": "", "path-ref-max": "" }
        opts.update(read_opts(args, def_pos, opts.keys()))
        for k in [ "stage" ]:
            if k not in opts:
                raise RuntimeError("no %s argument found!" % k)
        mk = solve_multi(opts["multi"])
        if "ACNNHOME" not in os.environ:
            raise RuntimeError('must have ACNNHOME environment variable!')
        for m in mk:
            tmdir = "./tomaster/%s-%s.%d" % (opts["runt"], opts["stage"], m)
            print ("running %s %s" % (task, tmdir))
            os.chdir(tmdir)
            if opts["path-ref-local"] != "" or opts["path-ref-max"] != "":
                if task != "connect":
                    raise RuntimeError('path-ref-* options can only be used for task=connect!')
                cnjn = read_json('%s.json' % task)
                for kop in ["path-ref-local", "path-ref-max"]:
                    if opts[kop] != "":
                        xsour = []
                        for l in opts[kop].split(","):
                            xsour.append("../../" + l)
                        cnjn["reaction-path"][kop] = xsour
                write_json(cnjn, '%s.json' % task)
            cmd = os.environ["ACNNHOME"] + '/acnnmain %s.json' % task
            # for line in os.popen(cmd):
            #     print (line.rstrip())
            run_cmd(cmd)
            os.chdir('../..')
        if self.to_dir(dox="local"):
            os.popen(os.environ["PGOPTHOME"] + '/pgopt sync').read()

    def report(self, args):
        pre = self.pre_info()
        k = self.to_dir(dox="local")
        if k: os.popen(os.environ["PGOPTHOME"] + '/pgopt sync').read()
        if not os.path.exists('./report'):
            os.mkdir('./report')
        base = os.path.basename(os.path.abspath(os.curdir))
        cmd = 'zip -r ./report/final.zip master tomaster *.xyz CMD-*'
        if not os.path.isfile('./report/final.zip'):
            print (os.popen(cmd).read().strip())
        dopts = { "tasks": [ "report" ], "output_dir": "./OUT", 
            "random_seed": 0, "report": { "input_file": "./final.zip" }, 
            "filtering-report": pre["filtering-report"] }
        dopts["report"].update(pre["report"])
        write_json(dopts, './report/report.json')
        if "ACNNHOME" not in os.environ:
            raise RuntimeError('must have ACNNHOME environment variable!')
        os.chdir('./report')
        cmd = os.environ["ACNNHOME"] + '/acnnmain report.json'
        # for l in os.popen(cmd):
        #     print (l.strip())
        run_cmd(cmd)
        if os.path.isfile('./OUT/report.pdf'):
            os.rename('./OUT/report.pdf', './%s-report.pdf' % base)
        else:
            print ("ERROR!!! NO PDF GENERATED!!!")
        # os.rename('./OUT/report_structs.xyz', './%s-structs.xyz' % base)
    
    def fdlog(self, args):
        pre = self.pre_info()
        fdsl = pre["finite-diff"]["step-length"]
        k = self.to_dir(dox="remote")
        jds = []
        for l in os.listdir("./master"):
            ll = "./master/%s" % l
            for k in os.listdir(ll):
                if k.startswith("par_runtime.json."):
                    ja, jb = l.split(".")
                    jc = k.split(".")[-1]
                    jrt = read_json("%s/%s.%s" % (ll, "par_runtime.json", jc))
                    if "fd" not in jrt:
                        continue
                    jrtfd = jrt["fd"]
                    jrtcov = [x if isinstance(x, int) else x[0] for x in jrt["corr_list"].values()]
                    jrtst = jrt["states"]
                    jkeys = sorted([int(m) for m in jrtfd.keys()])
                    for i in jkeys:
                        jf = len(jrtfd[str(i)][1]) * 6
                        je = jf + 1 - len(jrtfd[str(i)][0])
                        if je == 0 and i not in jrtcov:
                            ji = "init"
                        elif je == jf + 1:
                            ji = "finished"
                        elif je != 0 and i not in jrtcov:
                            ji = "paused"
                        else:
                            ji = "running"
                        jlist = [ja, jb, jc, str(i), je - 1, jf, ji]
                        jds.append(jlist)
        jds.sort(key=lambda x: map(int, [x[0], x[1], x[2], x[3]]))
        print ("# FD :: %s - %s - %s" % (pre["creation"]["name"], pre["creation"]["number"], 
            pre["creation"]["method"]))
        print ("# FD :: StepLength = %.5f Angstrom" % fdsl)
        print ("# %s :: %s / %s" % (pre["parallel"]["arguments"][0]["program"], 
            pre["parallel"]["arguments"][0]["functional"], 
            pre["parallel"]["arguments"][0].get("basis", "")))
        tit = "run-id ci fdid  cur/  tot    state"
        nhold = len(tit)
        print ("%s" % ("=" * nhold))
        print (tit)
        print ("%s" % ("=" * nhold))
        for ij, j in enumerate(jds):
            print (("%2s.%2s %3s%5s%5d/%5d %8s") % tuple(j))
            if ij != len(jds) - 1 and ij != 0 and jds[ij + 1][0] != j[0] and jds[ij - 1][0] == j[0]:
                print ("%s" % ("-" * nhold))
        print ("%s" % ("=" * nhold))

    def mclog(self, args):
        pre = self.pre_info()
        mcmi = pre["monte-carlo"]["max-iter"]
        mctp = pre["monte-carlo"]["temperature"]
        mcrs = pre["monte-carlo"].get("refuse-separate", True)
        mcdb = pre["monte-carlo"].get("detailed-balance", False)
        k = self.to_dir(dox="remote")
        jds = []
        for l in os.listdir("./master"):
            ll = "./master/%s" % l
            for k in os.listdir(ll):
                if k.startswith("par_runtime.json."):
                    ja, jb = l.split(".")
                    jc = k.split(".")[-1]
                    jrt = read_json("%s/%s.%s" % (ll, "par_runtime.json", jc))
                    if "mc" not in jrt:
                        continue
                    jrtmc = jrt["mc"]
                    jrtcov = jrt["corr_list"].values()
                    jrtst = jrt["states"]
                    jkeys = sorted([int(m) for m in jrtmc.keys()])
                    for i in jkeys:
                        jd = mcmi
                        je = jrtmc[str(i)][2]
                        jf = jrtmc[str(i)][0]
                        if jrtmc[str(i)][2] != 0:
                            jg = float(jrtmc[str(i)][0]) / jrtmc[str(i)][2]
                        else:
                            jg = 0.0
                        jh = jrtmc[str(i)][3]
                        if je == 0 and i not in jrtcov:
                            ji = "init"
                        elif jrtst[str(i)] != 0:
                            ji = "finished"
                        elif je != 0 and i not in jrtcov:
                            ji = "paused"
                        else:
                            ji = "running"
                        jcc = [0] * 5
                        if len(jrtmc[str(i)]) == 13:
                            jcc = [jrtmc[str(i)][5 + 1], jrtmc[str(i)][5 + 2] + jrtmc[str(i)][5 + 3], 
                                jrtmc[str(i)][5 + 4], jrtmc[str(i)][5 + 5] + jrtmc[str(i)][5 + 6], 
                                jrtmc[str(i)][5 + 7]]
                        jlist = [ja, jb, jc, str(i), jd, je] + jcc + [jf, jg, jh, ji]
                        jds.append(jlist)
        jds.sort(key=lambda x: map(int, [x[0], x[1], x[2], x[3]]))
        jtit = [ "%d K" % int(mctp), "detailed balance", "separation refused" ]
        if not mcdb:
            jtit[1] = "no " + jtit[1]
        if not mcrs:
            jtit[2] = "no " + jtit[2]
        print ("# MC :: %s - %s - %s" % (pre["creation"]["name"], pre["creation"]["number"], 
            pre["creation"]["method"]))
        print ("# MC :: " + " / ".join(jtit))
        print ("# %s :: %s / %s" % (pre["parallel"]["arguments"][0]["program"], 
            pre["parallel"]["arguments"][0]["functional"], 
            pre["parallel"]["arguments"][0].get("basis", "")))
        tit = "run-id ci mcid  tot  cur cnv fai max dup sep  rej/rate stepl    state"
        nhold = len(tit)
        print ("%s" % ("=" * nhold))
        print (tit)
        print ("%s" % ("=" * nhold))
        for ij, j in enumerate(jds):
            print (("%2s.%2s %3s%5s%5d%5d%4d%4d%4d%4d%4d%5d/%.2f %.3f %8s") % tuple(j))
            if ij != len(jds) - 1 and ij != 0 and jds[ij + 1][0] != j[0] and jds[ij - 1][0] == j[0]:
                print ("%s" % ("-" * nhold))
        print ("%s" % ("=" * nhold))

    def log(self, args):
        pre = self.pre_info()
        k = self.to_dir(dox="remote")
        jds = []
        for l in os.listdir("./master"):
            ll = "./master/%s" % l
            for k in os.listdir(ll):
                if k.startswith("par_shlog.txt."):
                    ja, jb = l.split(".")
                    jdl = open("%s/%s" % (ll, k), "r").read().strip()
                    jde, jdf = "     - ", "     - "
                    jd = [ja, jb] + jdl.split(" ")
                    # calculate running time and idle time
                    if os.path.isfile("%s/%s.%s" % (ll, "par_runtime.json", jd[8])):
                        js = read_json("%s/%s.%s" % (ll, "par_runtime.json", jd[8]))["times"].values()
                        jska, jskb = None, None
                        for j in js:
                            for jj in j:
                                if len(jj) >= 1 and (jska is None or jj[0] < jska):
                                    jska = jj[0]
                                if len(jj) >= 3 and (jskb is None or jj[2] > jskb):
                                    jskb = jj[2]
                    if jd[-1] == "waiting" or jd[-1] == "running":
                        if jska is not None:
                            jde = time_span_short(time.time() - jska)
                        if jskb is not None:
                            jdf = time_span_short(time.time() - jskb)
                    elif jd[-1] == "finished":
                        if jska is not None and jskb is not None:
                            jde = time_span_short(jskb - jska)
                    if len(jd) == 10: jd.insert(7, "-")
                    jd += [jde, jdf]
                    jds.append(jd)
        jds.sort(key=lambda x: map(int, [x[0], x[1]]))
        print ("# %s - %s - %s" % (pre["creation"]["name"], pre["creation"]["number"], 
            pre["creation"]["method"]))
        print ("# %s :: %s / %s" % (pre["parallel"]["arguments"][0]["program"], 
            pre["parallel"]["arguments"][0]["functional"], 
            pre["parallel"]["arguments"][0].get("basis", "")))
        tit = "run-id tot rem cnv fai max dup    ci/ti    cp/tp    state   time   idle"
        nhold = len(tit)
        print ("%s" % ("=" * nhold))
        print (tit)
        print ("%s" % ("=" * nhold))
        for ij, j in enumerate(jds):
            print (("%2s.%2s %4s%4s%4s%4s%4s%4s    %2s/%2s    %2s/%2s %8s %6s %6s") % tuple(j))
            if ij != len(jds) - 1 and ij != 0 and jds[ij + 1][0] != j[0] and jds[ij - 1][0] == j[0]:
                print ("%s" % ("-" * nhold))
        print ("%s" % ("=" * nhold))
    
    def draw(self, args):
        def_pos = { "0": "input", "1": "output" }
        opts = {"output": "./draw.pdf", "number": "-1", "ratio": "1.0"}
        optl = ["plot-charge", "plot-force", "surface-depth", "zdepth",
            "title", "rotmat", "prop", "force-factor", "force-image",
            "ortho"] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        if "input" not in opts:
            raise RuntimeError("Must have input file name!")
        assert os.path.isfile("./%s" % opts["input"])
        pre = {
            "tasks": ["draw"],
            "output_dir": ".",
            "random_seed": 1,
            "report": {
                "input_file": opts["input"],
                "output_file": opts["output"],
                "ratio": float(opts["ratio"]),
                "number": int(opts["number"])
            }
        }
        if "rotmat" in opts:
            pre["report"]["rotmat"] = [float(x) for x in opts["rotmat"].split(":")]
        if "title" in opts:
            pre["title"] = opts["title"]
        if "zdepth" in opts:
            pre["report"]["zdepth"] = float(opts["zdepth"])
        if "surface-depth" in opts:
            pre["report"]["surface_depth"] = float(opts["surface-depth"])
        if "plot-charge" in opts:
            pre["report"]["plot_charge"] = True
        if "plot-force" in opts:
            pre["report"]["plot_force"] = True
        if "force-factor" in opts:
            pre["report"]["force_factor"] = float(opts["force-factor"])
        if "force-image" in opts:
            pre["report"]["force_image"] = int(opts["force-image"])
        if "prop" in opts:
            pre["report"]["input_props"] = opts["prop"]
        if "ortho" in opts:
            pre["report"]["perspective"] = False
        write_json(pre, "./draw.json")
        cmd = os.environ["ACNNHOME"] + '/acnnmain draw.json'
        run_cmd(cmd)
        os.remove("./draw.json")

    def enlarge(self, args):
        def_pos = { "0": "input", "1": "cell", "2": "size" }
        opts = {"output": "."}
        opts.update(read_opts(args, def_pos, opts.keys()))
        if "cell" not in opts or "size" not in opts or "input" not in opts:
            raise RuntimeError("Must have input/cell/size parameters!")
        assert os.path.isfile("./%s" % opts["input"])
        pre = {
            "tasks": ["enlarge"],
            "output_dir": opts["output"],
            "random_seed": 1,
            "enlarge": {
                "input_file": opts["input"],
                "size": [int(x) for x in opts["size"].split(':')],
                "cell": [float(x) for x in opts["cell"].split(':')]
            }
        }
        write_json(pre, "./enlarge.json")
        cmd = os.environ["ACNNHOME"] + '/acnnmain enlarge.json'
        run_cmd(cmd)
        os.remove("./enlarge.json")

    def select(self, args):
        def_pos = { "0": "input", "1": "tags", "2": "output" }
        opts = { "output": "selected.xyz", }
        optl = [ "remove", "steps" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        if "tags" not in opts:
            tags = []
        else:
            tags = opts["tags"].split(",")
        if "steps" not in opts:
            steps = []
        else:
            steps = opts["steps"].split(",")
        if "remove" in opts:
            rems = opts["remove"].split(",")
        else:
            rems = []
        structs = []
        ntotal = 0
        nselected = 0
        with open(opts["input"], "r") as f:
            lines = f.readlines()
            n = int(lines[0].strip())
            for i in range(0, len(lines), n + 2):
                cmt = lines[i + 1].split("=")[0].strip().split(":")[1]
                cmts = lines[i + 1].split("=")[0].strip().split(":")[-1]
                sel = True
                if len(tags) == 0 and len(rems) == 0 and len(steps) == 0:
                    sel = True
                else:
                    if len(tags) != 0:
                        if '.' in cmt and '.' not in ''.join(tags):
                            sel = sel and (cmt.split(".")[0] in tags)
                        else:
                            sel = sel and (cmt in tags)
                    if len(rems) != 0:
                        sel = sel and (cmt not in tags)
                    if len(steps) != 0:
                        sel = sel and (cmts in steps)
                if sel:
                    structs.extend(lines[i:i + n + 2])
                    nselected += 1
                ntotal += 1
        with open(opts["output"], "w") as f:
            f.writelines(structs)
        print ("%5d/%5d structs selected." % (nselected, ntotal))

    def check(self, args):
        pre = self.pre_info()
        xmodel = pre["hosts"]["model"]
        jcmd = get_job_cmd(xmodel)
        self.to_dir(dox="remote")
        def_pos = { "0": "runt" }
        opts = { "runt": "para", "maxt": "60" }
        optl = [ "delete" ] + opts.keys()
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "maxt" ]:
            if k in opts: opts[k] = int(opts[k])
        for l in os.listdir("./master"):
            ll = "./master/%s" % l
            for k in os.listdir(ll):
                if k.startswith("par_shlog.txt."):
                    jc = k.split(".")[2]
                    ja, jb = l.split(".")
                    jdl = open("%s/%s" % (ll, k), "r").read().strip().split(" ")
                    fk = "par_log.txt." + jc
                    fxn = "%s/%s" % (ll, fk)
                    if os.path.isfile(fxn) and jdl[-1] != "finished":
                        tt = os.path.getmtime(fxn)
                        if time.time() - tt > opts["maxt"]:
                            print ("master %s.%s not running" % (ll, jc))
        if not os.path.exists("./procs/%s" % opts["runt"]):
            return
        os.chdir("./procs/%s" % opts["runt"])
        if pre["parallel"]["arguments"][0]["program"] in [ trans_prog("tm"), trans_prog("vasp") ] \
            and opts["runt"] == "para":
            for l in sorted(os.listdir(".")):
                idx = re.findall(r'PROC([0-9]+)$', l)
                if len(idx) != 0:
                    idx = int(idx[0])
                    tt = None
                    cc = None
                    ss = None
                    if os.path.isfile("./%s/tmpdir.txt" % l) and os.path.isfile("./%s/RUNNING" % l):
                        tmdir = os.popen("tail -n 1 ./%s/tmpdir.txt" % l).read().strip()
                        cc = os.path.basename(tmdir).split(".")
                        if len(cc) == 2: cc = cc[0]
                        else: cc = '.'.join(cc[0:2])
                        if os.path.isfile("./%s/SLEEPING" % l):
                            ss = True
                        else:
                            if os.path.exists(tmdir):
                                chkneb = False
                                for m in os.listdir(tmdir):
                                    if pre["parallel"]["arguments"][0]["program"] == trans_prog("tm"):
                                        if m.endswith(".chk"):
                                            cfn = "%s/%s/slave1.output" % (tmdir, m)
                                            if os.path.isfile(cfn):
                                                tt = os.path.getmtime(cfn)
                                            break
                                    elif pre["parallel"]["arguments"][0]["program"] == trans_prog("vasp"):
                                        if m.startswith("relax.out."):
                                            cfn = "%s/%s" % (tmdir, m)
                                            tt = os.path.getmtime(cfn)
                                        elif m.startswith("neb.out."):
                                            cfn = "%s/%s" % (tmdir, m)
                                            tt = os.path.getmtime(cfn)
                                            chkneb = True
                                        elif m.startswith("freq.out."):
                                            cfn = "%s/%s" % (tmdir, m)
                                            tt = os.path.getmtime(cfn)
                                if chkneb:
                                    for m in os.listdir(tmdir):
                                        if m.endswith(".chk"):
                                            for mm in os.listdir("%s/%s" % (tmdir, m)):
                                                mmc = "%s/%s/%s/OUTCAR" % (tmdir, m, mm)
                                                if os.path.isdir("%s/%s" % (tmdir, m)) and os.path.isfile(mmc):
                                                    tt = max(tt, os.path.getmtime(mmc))
                                            break
                                        
                    pstr = None
                    if ss is not None:
                        pstr = "%5d %15s sleeping" % (idx, cc)
                    elif tt is not None and time.time() - tt > opts["maxt"]:
                        pstr = "%5d %15s %7s" % (idx, cc, time_span_short(time.time() - tt))
                    if pstr is not None:
                        if "delete" not in opts:
                            print (pstr)
                        else:
                            if '.' in cc and jcmd[1] == "qdel":
                                cc = cc.split('.')
                                os.popen("%s %s -t %s" % (jcmd[1], cc[0], cc[1])).read()
                            elif '.' in cc and jcmd[1] == "canceljob":
                                cc = cc.split('.')
                                os.popen("%s %s" % (jcmd[1], cc[0])).read()
                            else:
                                os.popen("%s %s" % (jcmd[1], cc)).read()
                            os.remove("./%s/RUNNING" % l)
                            print (pstr + " deleted!")
        else:
            raise RuntimeError("Not implemented!")

    def help(self, args):
        PGHelp().help(args)
