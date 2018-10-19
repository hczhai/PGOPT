
from __future__ import print_function
import re, os, json, time, numpy as np, dill
import zipfile, copy
from utils.io import read_json
from formod.at_comp import at_comp
from cluster.base import read_zip_clusters, read_clusters, elem_num
from cluster.base import understand_name, update_stat, Cluster
from cluster.coval import AtomicWeight
from cluster.opt import Filtering
from cluster.align import Align
from surface.base import read_surface, ClusterAtSurface
from surface.surf_comp import surface_compare
from parallel.putils import RefStruct
from parallel.freqs import imag, imag_surface
from parallel.main import MonteCarloParams
from meta.ana_utils.base import FileItem, TimeItem, DictItem, StepAndTime, RawTimeRecord
from meta.ana_utils.thermal import moments_of_inertia, rot_order
from meta.ana_utils.plot_rpath import RPath, RGraph, REdge

httoev = 27.21138505

def formatted_name(name, charge, rise=3, surfn=None):
    name = re.sub(r'([0-9]+)', r'<sub rise=%d>\1</sub>' % rise, name)
    if charge == 1:
        name = name + '<sup rise=%d>+</sup>' % (rise)
    elif charge > 1:
        name = name + '<sup rise=%d>%d+</sup>' % (rise, charge)
    elif charge == -1:
        name = name + '<sup rise=%d>-</sup>' % (rise)
    elif charge < -1:
        name = name + '<sup rise=%d>%d-</sup>' % (rise, -charge)
    if surfn is not None:
        surfn = re.sub(r'([0-9]+)', r'<sub rise=%d>\1</sub>' % rise, surfn)
        name = name + '@' + surfn
    return name

# for comparison
def tint(x):
    g = re.findall(r'([0-9]+)', x)
    if len(g) == 1: return [int(g[0]), x.replace(g[0], '')]
    else: return x

# number to roman
def roman(ip):
    ints = (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
    nums = ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    result = ""
    for i in range(len(ints)):
        count = int(ip / ints[i])
        result += nums[i] * count
        ip -= ints[i] * count
    return result

# structs include par_local and par_max
# must be of ClusterAtSurface type
def read_nebpaths(structs, nlocal):
    lmx = structs
    if len(lmx) != 0 and isinstance(lmx[0], ClusterAtSurface):
        from surface.surf_comp import surface_align
        for cc in lmx:
            surface_align(cc, deep=False)
    nebpaths = {}
    for il, l in enumerate(lmx):
        stid = str(l.tid)
        if stid not in nebpaths:
            nebpaths[stid] = RPath(l.tid, ismax=(il >= nlocal))
        nebpaths[stid].add(l)
    for np in nebpaths.values():
        np.sort()
    # focus the lowest barrier
    lwp = {}
    for np in nebpaths.values():
        ft = "%d.%d" % tuple(np.from_to)
        if ft not in lwp:
            lwp[ft] = [np.barriers[0], np]
        elif lwp[ft][0] > np.barriers[0]:
            lwp[ft] = [np.barriers[0], np]
    for np in nebpaths.values():
        ft = "%d.%d" % tuple(np.from_to)
        if lwp[ft][1] is np:
            np.focus = True
    nebmins = {}
    for np in nebpaths.values():
        if str(np.from_to[0]) not in nebmins:
            nebmins[str(np.from_to[0])] = np.structs[0]
        if str(np.from_to[1]) not in nebmins:
            nebmins[str(np.from_to[1])] = np.structs[-1]
    return nebpaths, nebmins

class StructData(object):
    def __init__(self, tid):
        self.tid = tid
        self.has_freq = False

class SourceData(object):
    def __init__(self, name, start, end, ref_run=None, ref_structs=None):
        self.name = name
        self.ref_run = ref_run
        self.ref_structs = ref_structs
        self.start = start
        self.end = end
        self.number = end - start
        self.range = range(start, end)

# level 1
class RunData(object):
    def __init__(self, name):
        self.name = name
        self.stage, self.m, self.runt = name.split('.')
        self.pname = "%s.%s.%s" % (self.stage, self.m, self.runt)
        self.m = int(self.m)
        self.runt = int(self.runt)
        self.xdata = DictItem()
        for key in ["local", 'max', 'dupl', 'filter']:
            self.xdata[key] = []
        for key in ["props", "runtime"]:
            self.xdata[key] = {}
        self.xdata.runtime = {"corr_list": {}, "times": {}, "states": {}, \
            "steps": {}, "restart": {}}

        self.param = {}
        self.surf = None
        self.focus = None
        self.energy_focus = None
        self.locator = None
        self.mc = None

    def analyze(self, order):
        if order == 0:
            print ("Analyzing Run %s ..." % (self.pname))
            self.times()
            self.procs()
            self.params()
            self.structs()
            self.mcdata()
        elif order == 1:
            self.post_props()
        elif order == 2:
            self.ppost_props()
            self.thermal()
        elif order == 3:
            self.conv()

    def post_props(self):

        # determine ref relationship
        xxsour = [x[0] for x in self.args["sources"]]
        self.sour_name = "unknown"
        if "inherit" in xxsour:
            self.sour_name = "inherit (%d)" % (self.runt - 1)
            xname = "%s.%d.%d" % (self.stage, self.m, self.runt - 1)
            for c in self.optima.values():
                ppt = c.props["prev_tidx"]
                assert self.runt - 1 == ppt[0]
                c.ref_tidx = [self.locator[xname], ppt[1]]
        elif "restart" in xxsour:
            # self.sour_name = "restart"
            restrs = []
            nadd = []
            xhas = True
            for src in self.args["sources"]:
                xx = RefStruct.solve(src[1], idx=self.runt)
                xname = "%s.%d.%d" % (xx.xsta, xx.xmul, xx.xidx)
                nadd.append("%s (%s)" % (xx.xname, xname))
                if xname not in self.locator:
                    xhas = False
                    continue
                xxstr = self.locator[xname].xdata[xx.xname]
                if len(src) >= 3:
                    xxstr = xxstr[:src[2]]
                for x in xxstr:
                    restrs.append([self.locator[xname], x.tidx, x.energy])
            self.sour_name = ",".join(nadd)
            if not xhas:
                for c in self.optima.values():
                    c.ref_tidx = None
            if self.args.get("max_config", -1) != -1:
                restrs.sort(key=lambda x: x[2])
                restrs = restrs[:self.args["max_config"]]
            for c in self.optima.values():
                c.ref_tidx = None
            for ic, c in enumerate(restrs):
                if (str(ic) + ".0") in self.optima:
                    self.optima[str(ic) + ".0"].ref_tidx = c[:2]
        elif len(xxsour) != 0:
            self.sour_name = xxsour[0]
            for c in self.optima.values():
                c.ref_tidx = None

        # initial number of duplicates
        for c in self.optima.values():
            if c.ref_tidx is None:
                c.props["nrep"] = 1
            else:
                targ = c.ref_tidx[0].optima["%d.%d" % tuple(c.ref_tidx[1])]
                c.props["nrep"] = targ.props["nrep"]

        # add number of duplicates
        for c in self.optima.values():
            if c.ttype == "dupl":
                xdupl = c.props["dupl"]
                if isinstance(xdupl[0], int):
                    self.optima["%d.%d" % tuple(xdupl)].props["nrep"] += c.props["nrep"]
                else:
                    xx = RefStruct(None, xdupl[2], xdupl[0], xdupl[1], xdupl[3], '.')
                    xname = "%s.%d.%d" % (xx.xsta, xx.xmul, xx.xidx)
                    self.locator[xname].optima["%d.%d" % tuple(xdupl[-1])].props["nrep"] += c.props["nrep"]

        # add energy and frequency info to focused runs
        if self.focus is not self:
            if self.args["task_type"] == "energy":
                for c in self.optima.values():
                    if c.ref_tidx is not None and "energy" in c.props:
                        targ = c.ref_tidx[0].optima["%d.%d" % tuple(c.ref_tidx[1])]
                        targ.props["energy.%s.%d" % (self.stage, self.runt)] = c.props["energy"]
                        targ.props["time.energy.%s.%d" % (self.stage, self.runt)] = c.props["time"]
                        if self.args["charge"] == c.ref_tidx[0].args["charge"] + 1:
                            targ.props["vip.energy.%s.%d" % (self.stage, self.runt)] = \
                                targ.props["vip.energy"] = c.props["energy"] - targ.props["energy"]
                            c.props["vip.energy"] = targ.props["vip.energy"]
                        if "bader_charges" in c.props:
                            targ.props["charges"] = c.props["bader_charges"]
            if self.args["task_type"] == "freq" or (self.args["task_type"] == "energy" and 
                self.args.get("do-finite-diff", False)):
                for c in self.optima.values():
                    if c.ref_tidx is not None:
                        targ = c.ref_tidx[0].optima["%d.%d" % tuple(c.ref_tidx[1])]
                        for key in ["freqs", "displacements"]:
                            if key in c.props:
                                targ.props["%s.%s.%d" % (key, self.stage, self.runt)] = targ.props[key] = c.props[key]
                        targ.props["time.freqs.%s.%d" % (self.stage, self.runt)] = targ.props["time.freqs"] = c.props["time"]
                        if "freqs" in c.props and c.props["freqs"] is not None:
                            if isinstance(c, Cluster):
                                targ.minimum_type = 2 if imag(c.props["freqs"]) > 0 else 1
                            else:
                                targ.minimum_type = 2 if imag_surface(c.props["freqs"]) > 0 else 1
    
    def ppost_props(self):
        if self.focus is not self:
            for c in self.optima.values():
                if c.ref_tidx is not None and "charges" in c.props:
                    targ = c.ref_tidx[0].optima["%d.%d" % tuple(c.ref_tidx[1])]
                    targ.props["charges"] = c.props["charges"]
                if c.ref_tidx is not None and "freqs" in c.props:
                    targ = c.ref_tidx[0].optima["%d.%d" % tuple(c.ref_tidx[1])]
                    for key in ["freqs", "displacements", "time.freqs"]:
                        if key in c.props:
                            targ.props[key] = c.props[key]
                    targ.minimum_type = c.minimum_type

    def thermal(self):
        for c in self.optima.values():
            if isinstance(c, ClusterAtSurface):
                continue
            c.props["moments.of.inertia"] = moments_of_inertia(c)
            c.props["rot.order"] = rot_order(c.pointgroup())
            c.props["spatial.degeneracy"] = 1 # include point group repr in the future

    def structs(self):
        self.optima = DictItem()
        for key in ["filter"]:
            for c in self.xdata[key]:
                labels = c.label.split(":")
                stid = labels[1]
                c.tid = int(stid.split(".")[0])
                c.tidx = [c.tid, int(stid.split(".")[1])]
                c.stidx = stid
                c.ttype = key
        nf = []
        for key in ["sources", "local", "max", "dupl"]:
            for c in self.xdata.get(key, {}):
                labels = c.label.split(":")
                stid = labels[1]
                c.tid = int(stid.split(".")[0])
                c.tidx = [c.tid, int(stid.split(".")[1])]
                c.stidx = stid
                c.ttype = key
                if stid in self.xdata.props:
                    c.props = self.xdata.props[stid]
                    if "bader_charges" in c.props:
                        c.props["charges"] = c.props["bader_charges"]
                else:
                    c.props = {}
                    if self.args["task_type"] in ["energy", "relax"] or \
                        self.args["task_type"] == "neb" and c.tidx[1] == 0:
                        nf.append(c.tidx)
                if stid in self.time_step:
                    c.props["time"] = self.time_step[stid].total_time()
                    c.props["step"] = self.time_step[stid].total_step()
                    c.mcac = self.time_step[stid].mac
                else:
                    c.props["time"] = 0.0
                    c.props["step"] = 0
                if hasattr(c, 'mag') and c.mag is not None:
                    c.multiplicity = int(np.round(c.mag)) + 1
                else:
                    c.multiplicity = self.m
                c.minimum_type = 0 # 0 UNK, 1 MIN, 2 TS
                c.tname = "%s.%d.%d" % ((self.stage, ) + tuple(c.tidx))
                c.aname = "%s.%s.%d.%d.%d" % (self.parent.parent.parent.name,
                    self.stage, c.multiplicity, c.tidx[0], c.tidx[1])
                self.optima[stid] = c
        if len(nf) != 0:
            g = ", ".join(["%d.%d" % tuple(x) for x in nf])
            print ("WARNING: %s.%d.%d - NOT FINISHED = [%s]" % (self.stage, self.m, self.runt, g))

    # after structs
    def mcdata(self):
        if "monte-carlo" in self.xdata.input:
            ipmc = self.xdata.input["monte-carlo"]
            self.mc = MonteCarloParams.read_from(ipmc)
            if not self.args.get("do-monte-carlo", True):
                self.mc = None
        if self.mc is not None:
            self.mce = {}
            for key in ["local", "max", "dupl"]:
                for c in self.xdata[key]:
                    if not c.mcac:
                        continue
                    stid = str(c.tid)
                    if stid not in self.mce:
                        self.mce[stid] = []
                    self.mce[stid].append([c.tidx[1], c.energy])
            for k, v in self.mce.items():
                v.sort(key=lambda x: x[0])

    def conv(self):
        if self.focus is self and self.args["task_type"] == "relax":
            conv_cut = self.parent.parent.parent.parent.conv_cut
            lall = [x.tid for x in self.xdata.local]
            if len(lall) == 0:
                self.conv_data = [[0, [0, 0]]]
                return
            en_cut = self.xdata.local[0].energy + conv_cut / httoev
            lcut = [x.tid for x in self.xdata.local if x.energy < en_cut]
            self.conv_data = []
            ncut = 0
            nall = 0
            for i in range(0, max(lall) + 1):
                if i in lcut:
                    ncut += 1
                if i in lall:
                    nall += 1
                self.conv_data.append([i, [ncut, nall]])

    # after times
    def procs(self):
        self.xdata.log = [l.strip() for l in self.xdata.log]
        lastl = self.xdata.log[-1]
        shtime = ' '.join(lastl.split(' ')[0:2])
        if len(self.time_item.intervals) == 0:
            time_align = 0
        else:
            sttime = self.time_item.end()
            time_align = sttime - time.mktime(time.strptime(shtime, "%y-%m-%d %H:%M:%S"))
        # avail procs at different times
        ap = []
        for l in self.xdata.log[4:]:
            if "+-" not in l:
                continue
            shtime = ' '.join(l.split(' ')[0:2])
            shtimex = time.mktime(time.strptime(shtime, "%y-%m-%d %H:%M:%S")) + time_align
            npr = float([k for k in l.split(' ') if len(k) != 0][7].split('/')[1])
            # time, avail procs, running procs
            ap.append([shtimex, npr, 0])
        self.job_time = ap[-1][0] - ap[0][0]
        if len(ap) >= 100:
            ap = ap[0:len(ap):len(ap) / 100]
        for x in ap:
            for tt in self.time_proc.values():
                if tt.has(x[0]): x[2] += 1.0
        self.avail_procs = np.array([a[1] for a in ap])
        self.run_procs = np.array([a[2] for a in ap])
        self.xtimes = { "walltime": self.job_time, 
            "nodetime": sum([p.running_time() for p in self.time_proc.values()]) }
        self.xtimes.update({ "cputime": self.cpu_time, 
            "costtime": self.factor * self.xtimes["nodetime"] })

    def times(self):
        self.cores = self.xdata.input["hosts"]["cores"]
        self.nodes = self.xdata.input["hosts"]["nodes"]
        # time item of the run
        self.time_item = TimeItem()
        # time item of all procs
        self.time_proc = {}
        # time and step of each structure
        self.time_step = {}
        self.cpu_time = 0.0
        if self.xdata.input["parallel"].get("proc_type", "para") == "para":
            self.factor = self.cores * self.nodes
        else:
            self.factor = 1.0
        ti = self.xdata.runtime["times"]
        traws = {}
        for j, tt in ti.items():
            self.time_proc[j] = TimeItem()
            for t in tt:
                if len(t) > 3:
                    self.time_proc[j].add(t[0], t[2])
                    self.cpu_time += t[3] * self.factor
            self.time_item.adds(self.time_proc[j].intervals)
            for t in tt:
                if len(t) < 3:
                    continue
                tinfo = t[-1]
                tminfo = None
                if isinstance(tinfo, list):
                    tminfo = tinfo[1]
                    tinfo = tinfo[0]
                if tinfo.startswith("nr.") or tinfo in ["no.consis.fail", "time.r.fail"]:
                    continue
                xname = str(t[1])
                if xname not in traws:
                    traws[xname] = []
                if tinfo.startswith("r.n.conv."):
                    xstep = 1
                elif tinfo in ["r.unk", "r.repeat", "r.fatal", "r.nostr"]:
                    xstep = 0
                elif tinfo.startswith("r."):
                    xstep = int(tinfo.split("~")[0].split(".")[2].split(":")[0])
                elif tinfo.startswith("neb.n.conv."):
                    xstep = 1
                elif tinfo in ["neb.unk", "neb.repeat", "neb.fatal"]:
                    xstep = 0
                elif tinfo.startswith("neb."):
                    xstep = int(tinfo.split(".")[2].split(":")[0])
                else:
                    xstep = 1
                mstep = -1
                mac = False
                if tminfo is not None and tminfo.startswith("mc."):
                    mstep = int(tminfo.split(">")[0].split(".")[-1]) - 1
                    mac = "ac" in tminfo or "init" in tminfo
                traws[xname].append(RawTimeRecord(tid=t[1], rstep=xstep, mstep=mstep,
                                                  rtime=t[3], stime=t[0], mac=mac))
        for _, tv in traws.items():
            ttmp = []
            tv.sort(key=lambda x: x.Stime)
            tx = 0
            for v in tv:
                ttmp.append(v)
                if v.Mstep != -1:
                    assert tx == v.Mstep
                    for m in ttmp:
                        m.Mstep = tx
                        m.Mac = v.Mac
                    tx += 1
                    ttmp = []
            for m in ttmp:
                m.Mstep = tx
            for v in tv:
                xname = "%d.%d" % (v.Tid, v.Mstep)
                if xname not in self.time_step:
                    self.time_step[xname] = StepAndTime(v.Tid, v.Mstep, v.Mac)
                self.time_step[xname].add(v.Rstep, v.Rtime)

    def params(self):
        self.args = self.xdata.input["parallel"]["arguments"][self.runt]
        self.task_type = self.args["task_type"]
        self.param = { \
            "program": self.args["program"],
            "basis": self.args.get("basis", ""),
            "functional": self.args["functional"],
            "others": self.args["others"],
            "charge": self.args["charge"],
            "step": self.args.get("step", 0),
            "max_step": self.args.get("max_step", 0),
            "freq_step": self.args.get("freq_step", 0.0),
            "name": self.xdata.input["creation"]["name"],
            "method": self.xdata.input["creation"]["method"],
            "number": self.xdata.input["creation"]["number"],
            "max_config": self.args.get("max_config", -1),
            "filtering": { \
                "create": self.xdata.input["filtering-create"],
                "run": self.xdata.input["filtering-parallel"],
                "final": self.xdata.input["filtering-report"]
            }
        }
        if self.xdata.input["creation"]["method"] == "blda":
            self.param["order"] = self.xdata.input["creation"]["order"]
            self.param["default_sigma"] = self.xdata.input["creation"]["default_sigma"]
            self.param["dist"] = self.xdata.input.get("dist", {"mu": {}, "sigma": {}})
        elif self.xdata.input["creation"]["method"] == "ck":
            self.param["drel"] = self.xdata.input["creation"]["drel"]
        self.param["surfname"] = None
        if "creation-surface" in self.xdata.input:
            xics = self.xdata.input["creation-surface"]
            if "surface" in self.args:
                xics = self.args["surface"]
            surfn = xics["surface"].split("/")[-1]
            self.param["surfname"] = surfn.split(".")[0]
            self.surf = copy.deepcopy(self.surf_locator[surfn])
            self.surf.unit_cell = xics["unit_cell"]
            self.surf.space_group = xics["space_group"]
            self.surf.space_group_ref = np.array(xics["space_group_ref"])
        if self.surf is not None:
            for key in ["sources", "local", "max", "dupl"]:
                if key in self.xdata:
                    self.xdata[key] = [ClusterAtSurface.from_cluster(l, self.surf) \
                        for l in self.xdata[key]]

        self.param["fname"] = formatted_name(self.param["name"], self.param["charge"], \
            surfn=self.param["surfname"])
        self.description = "[Run %s.%s.%s (M=%d)] %s - %s - %d; %s :: %s/%s" % \
            (self.parent.parent.parent.name, self.parent.parent.name, self.stage,
             self.m, self.param["fname"], self.param["method"], self.param["number"],
             self.param["program"], self.param["functional"], self.param["basis"])

# level 2
class StageData(object):
    def __init__(self, stage):
        self.stage = stage
        self.parent = None
        self.runs = []

    def add_run(self, run):
        self.runs.append(run)
        run.parent = self

    def focus(self):
        mts = list(set([r.m for r in self.runs]))
        for m in mts:
            xi = None
            for r in self.runs:
                if r.m != m:
                    continue
                if r.args["task_type"] in ["relax", "neb"] and (xi is None or \
                    (xi.runt < r.runt and (len(r.xdata.local) != 0 or len(xi.xdata.local) == 0))):
                    xi = r
            for r in self.runs:
                r.focus = xi
            if xi is not None:
                continue
            for r in self.runs:
                if r.m != m:
                    continue
                if r.args["task_type"] == "energy" and xi is None:
                    xi = r
            for r in self.runs:
                r.energy_focus = xi

    def analyze(self, order):
        for r in self.runs:
            r.analyze(order=order)
        if order == 0:
            self.focus()
            self.times()
            if self.runs[0].focus is not None:
                self.param = self.runs[0].focus.param
            else:
                self.param = self.runs[0].param
            self.fname = self.param["fname"]
            keys = ["name", "fname", "number", "dist", "order", "default_sigma", "drel", "method"]
            self.creation_param = {}
            for k in keys:
                if k in self.param:
                    self.creation_param[k] = self.param[k]
            felems, _, _ = understand_name(self.creation_param["name"])
            update_stat(felems, self.creation_param["dist"]["mu"], \
                self.creation_param["dist"]["sigma"], self.creation_param["default_sigma"])
            for k in self.creation_param["dist"]["mu"].keys():
                if AtomicWeight.x[k.split("-")[1]] > AtomicWeight.x[k.split("-")[0]]:
                    del self.creation_param["dist"]["mu"][k]
                    del self.creation_param["dist"]["sigma"][k]
            self.description = "[Stage %s.%s.%s] %s - %s - %d; %s :: %s/%s" % \
                (self.parent.parent.name, self.parent.name, self.stage, self.fname,
                 self.param["method"], self.param["number"],
                 self.param["program"], self.param["functional"], self.param["basis"])
        elif order == 3:
            self.conv()

    def conv(self):
        self.conv_data = {}
        self.conv_maxs = [0, 0]
        for r in self.runs:
            if not hasattr(r, "conv_data"):
                continue
            self.conv_data[str(r.m)] = r.conv_data
            ma = max([d[0] for d in r.conv_data])
            mb = max([d[1][1] for d in r.conv_data])
            if ma > self.conv_maxs[0]:
                self.conv_maxs[0] = ma
            if mb > self.conv_maxs[1]:
                self.conv_maxs[1] = mb

    def times(self):
        self.time_item = TimeItem()
        self.cpu_time = 0.0
        for r in self.runs:
            self.time_item.adds(r.time_item.intervals)
            self.cpu_time += r.cpu_time

    def get_group(self):
        runx = sorted(self.runs, key=lambda x: [x.runt, x.m])[0]
        return runx.xdata.input["parallel"]["secondary"]

# level 3
class GroupData(object):
    def __init__(self, group, parent):
        self.group = group
        self.stages = []
        self.global_min = None
        self.multi_locals = {}
        self.multi_finals = {}
        self.multi_neblocals = {}
        self.multi_nebmaxs = {}
        self.multi_stat = {}
        self.locals = []
        self.finals = []
        self.neblocals = []
        self.nebmaxs = []
        self.nebpaths = []
        self.nebmins = []
        self.stat = None
        self.nebgraph = None
        self.parent = parent

    def add_stage(self, stage):
        self.stages.append(stage)
        stage.parent = self

    def analyze(self):
        for i in range(0, 4):
            for s in self.stages:
                s.analyze(order=i)
        self.times()
        self.param = self.stages[0].param
        self.fname = self.param["fname"]
        self.surfname = self.param["surfname"]
        self.description = "[Group %s.%s] %s - %s - %d; %s :: %s/%s" % \
            (self.parent.name, self.name, self.fname, self.param["method"], self.param["number"],
             self.param["program"], self.param["functional"], self.param["basis"])
        self.filtering()
        self.similarity()
        self.initial_structs()

    def initial_structs(self):
        self.xsour_structi = {}
        for s in self.stages:
            for r in s.runs:
                if r.task_type == "relax":
                    k = r.sour_name.split(" ")[0]
                    if k not in self.xsour_structi:
                        self.xsour_structi[k] = []
                    for x in r.xdata.local:
                        self.xsour_structi[k].append(DictItem(energy=x.props["initial_energy"]))
        self.xsour_structi["final"] = self.finals

    def filtering(self):
        xcount = lambda k, g: len(filter(lambda x, k=k: x.minimum_type == k, g))
        for s in self.stages:
            for r in s.runs:
                if r.focus is not r and r.energy_focus is not r:
                    continue
                # if r.focus is not r and r.energy_focus is not r: continue
                if str(r.m) not in self.multi_locals:
                    self.multi_locals[str(r.m)] = []
                    self.multi_neblocals[str(r.m)] = []
                    self.multi_nebmaxs[str(r.m)] = []
                if r.focus is r and r.args["task_type"] == "relax":
                    self.multi_locals[str(r.m)] += r.xdata.local
                elif r.focus is r and r.args["task_type"] == "neb":
                    self.multi_neblocals[str(r.m)] += r.xdata.local
                    self.multi_nebmaxs[str(r.m)] += r.xdata.max
                elif r.energy_focus is r:
                    self.multi_locals[str(r.m)] += r.xdata.sources
        for k, l in self.multi_locals.items():
            # print k, len(l), np.array([c.energy for c in l]).min()
            l.sort(key=lambda x: x.energy)
            if len(l) != 0:
                filt = Filtering(**self.param["filtering"]["final"])
                filt.init_from_create(l, len(l), l[0].n)
                filt.flimit = False
                filt.filter(iprint=False)
                lf = filt.finals
            else:
                lf = []
            self.multi_finals[k] = lf
            self.finals += lf
            self.locals += l
            self.multi_stat[k] = [xcount(1, l), xcount(2, l), xcount(0, l), \
                xcount(1, lf), xcount(2, lf), xcount(0, lf)]
        for k, l in self.multi_neblocals.items():
            self.neblocals.extend(l)
        for k, l in self.multi_nebmaxs.items():
            self.nebmaxs.extend(l)
        self.stat = [xcount(i, self.locals) for i in [1, 2, 0]] + \
            [xcount(i, self.finals) for i in [1, 2, 0]]
        self.finals.sort(key=lambda x: x.energy)
        if len(self.finals) != 0 and isinstance(self.finals[0], ClusterAtSurface):
            from surface.surf_comp import surface_align
            for cc in self.finals:
                surface_align(cc)
        else:
            align = Align()
            align.clus = self.finals
            align.align()
        lmx = self.neblocals + self.nebmaxs
        self.nebpaths, self.nebmins = read_nebpaths(lmx, len(self.neblocals))
        if self.parent.parent.nebmins_ref is not None:
            print ("Used reference minina energies ...")
            cref = read_clusters("../" + self.parent.parent.nebmins_ref, iprint=True)
            for icr, cr in enumerate(cref):
                if str(icr) in self.nebmins:
                    self.nebmins[str(icr)].energy = cr.energy
            for np in self.nebpaths.values():
                np.structs[0].energy = self.nebmins[str(np.from_to[0])].energy
                np.structs[-1].energy = self.nebmins[str(np.from_to[1])].energy
                np.sort()
        rgraph = RGraph()
        for np in self.nebpaths.values():
            if np.focus and np.isdirect:
                rgraph.add_edge(REdge(*(np.from_to + list(np.barriers * httoev))))
        self.nebgraph = rgraph
        rgraph.mins = self.nebmins

    def similarity(self):
        mdff = self.param["filtering"]["final"]["max_diff_report"]
        if len(self.finals) != 0:
            c = self.finals[0]
            eles, ne = elem_num(c.elems)
        for i, c in enumerate(self.finals):
            mj = -1
            md = None
            for j, d in enumerate(self.finals[:i]):
                if isinstance(c, ClusterAtSurface):
                    v, _ = surface_compare(d, c, ne, eles, mdff)
                else:
                    v, _ = at_comp(d.atoms, c.atoms, ne, eles, mdff)
                if mj == -1 or md > v:
                    md = v
                    mj = j
            if mj != -1 and md != mdff:
                c.similar = [mj, md]
            else:
                c.similar = None

    def times(self):
        self.time_item = TimeItem()
        self.cpu_time = 0.0
        for s in self.stages:
            self.time_item.adds(s.time_item.intervals)
            self.cpu_time += s.cpu_time

# level 4
class InstanceData(object):
    def __init__(self, input_file, read_traj=True):
        self.input_file = input_file
        self.filename = os.path.basename(input_file)

        # read par_structs or not
        self.read_traj = read_traj
        self.xdata = DictItem()
        for key in ["log", "runtime", "sources", "local", "disp", "dupl", "max", \
            "props", "input", "structs", "filter"]:
            self.xdata[key] = DictItem()
        self.metas = []

        # output
        self.xglobal = {}
        self.xsystem = {}

    def read(self):
        zf = zipfile.ZipFile(self.input_file, 'r')
        print ("## %s" % self.filename)
        namel = zf.namelist()
        tree = FileItem("/")
        for l in namel:
            lf = FileItem(l)
            k = tree
            for i in lf.parents: k = k[i]
            k[lf.name] = lf
        tomasters = tree["tomaster"].path_names()
        masters = tree["master"].path_names()
        for m in masters:
            mfiles = tree["master"][m].file_names()
            for f in mfiles:
                fn = tree["master"][m][f].path
                idx = f.split('.')[2]
                xname = m + '.' + idx
                if f.startswith("par_structs") and not self.read_traj: continue
                if f.startswith("par_log"):
                    self.xdata.log[m + '.'] = zf.open(fn, 'rU').readlines()
                for key in ["props", "runtime"]:
                    if f.startswith("par_" + key):
                        self.xdata[key][xname] = json.loads(zf.open(fn, 'rU').read())
                for key in ["structs", "sources", "local", "disp", "dupl", "max", "filter"]:
                    if f.startswith("par_" + key):
                        self.xdata[key][xname] = read_clusters(fn, zf=zf, iprint=True)
                # fil_structs ignored

        for m in tomasters:
            mfiles = tree["tomaster"][m].file_names()
            for f in mfiles:
                fn = tree["tomaster"][m][f].path
                if f.startswith("para.json"):
                    xname = m.split('-')[1] + '.'
                    self.xdata.input[xname] = read_json(zf.open(fn, 'rU').read(),
                                                        cont=True, iprint=False)
                elif f.startswith("meta_info"):
                    self.metas.append(zf.open(fn, 'rU').readlines())

        self.surf_locator = {}
        self.cmd_hist = []
        for k in tree.subs.keys():
            if k.endswith(".xyz"):
                fn = tree[k].path
                self.surf_locator[k] = read_surface(fn, zf=zf, iprint=False)
            if k == "CMD-HISTORY":
                fn = tree[k].path
                self.cmd_hist = [m.strip() for m in zf.open(fn, 'rU').readlines()]
        print ("data loading finished")

    # after read
    def make_groups(self):
        self.locator = runl = {}
        for k, rs in self.xdata.runtime.items():
            runl[k] = RunData(k)
            runl[k].locator = self.locator
            runl[k].surf_locator = self.surf_locator
        for key in self.xdata.keys():
            # log file is shared by different idxs
            if key in ["log", "input"]:
                for k, ls in self.xdata[key].items():
                    for kk in runl.keys():
                        if kk.startswith(k): runl[kk].xdata[key] = ls
            else:
                for k, ls in self.xdata[key].items():
                    if k in runl: runl[k].xdata[key] = ls

        # runl is idx-level
        stagel = {}
        for k, r in runl.items():
            if "runtime" not in r.xdata: continue
            m = k.split(".")[0]
            if m not in stagel:
                stagel[m] = StageData(m)
            stagel[m].add_run(r)

        groupl = {}
        for k, s in sorted(stagel.items(), key=lambda x: tint(x[0])):
            d = s.get_group()
            if d == "":
                groupl[k] = GroupData(k, self)
                d = k
            groupl[d].add_stage(s)

        # sort results
        self.groups = sorted(groupl.values(), key=lambda x: tint(x.group))
        k = 'A'
        self.stages = []
        self.runs = []
        for g in self.groups:
            g.name = k
            k = chr(ord(k) + 1)
            g.stages.sort(key=lambda x: tint(x.stage))
            self.stages += g.stages
            for s in g.stages:
                s.group = k
                s.runs.sort(key=lambda x: [x.runt, x.m])
                self.runs += s.runs

    # after read, make groups
    def analyze(self):
        for g in self.groups:
            g.analyze()
        self.times()
        self.outputs()
        self.fnames = [g.fname for g in self.groups]
        self.lnames = [formatted_name(g.param["name"], g.param["charge"], rise=5, \
            surfn=g.param["surfname"]) for g in self.groups]

    def times(self):
        self.time_item = TimeItem()
        self.cpu_time = 0.0
        for g in self.groups:
            self.time_item.adds(g.time_item.intervals)
            self.cpu_time += g.cpu_time

    def outputs(self):
        # system info
        if len(self.metas) == 0:
            raise RuntimeError("No meta data found!")
        gmetas = []
        for meta in self.metas:
            metas = [[g.strip() for g in m.strip().split('=') if len(g) != 0] for m in meta if '=' in m]
            gmetas.append(metas)
        gmetasx = [[x[1] for x in xm] for xm in gmetas]
        self.xsystem = {m[0]: ", ".join(list(set(mg))) for m, mg in zip(gmetas[0], zip(*gmetasx))}

        self.xglobal["start_time"] = time.localtime(self.time_item.start())
        self.xglobal["finish_time"] = time.localtime(self.time_item.end())
        self.xglobal["running_time"] = self.time_item.running_time()
        self.xglobal["total_time"] = self.time_item.total_time()
        self.xglobal["idle_time"] = self.time_item.total_time() - self.time_item.running_time()
        self.xglobal["total_cpu_time"] = int(self.cpu_time)

class OptimizationData(object):
    def __init__(self, input_file, ipff, conv_cut=0.4, nebmins_ref=None):
        self.ipff = ipff
        self.conv_cut = conv_cut
        self.nebmins_ref = nebmins_ref
        self.input_file = input_file if isinstance(input_file, list) else [input_file]
        self.d = [InstanceData(i, read_traj=False) for i in self.input_file]
        for i, d in enumerate(self.d):
            d.parent = self
            d.read()
            d.make_groups()
            d.name = roman(i + 1)
            d.analyze()
        self.lnames = list(set(sum([d.lnames for d in self.d], [])))
        self.lnames.sort(key=lambda x: tint(x))
        self.lname = " / ".join(self.lnames)
        self.title = self.lname + " Global Optimization Report"
