
import os, sys, time, re, numpy as np, json, random
from formod.at_comp import at_comp
from cluster.base import Cluster, read_zip_clusters, read_clusters, elem_num
from surface.base import ClusterAtSurface, read_surface
from cluster.coval import AtomicWeight
from surface.surf_comp import surface_compare
from surface.surf_symm import to_cellmat, to_dcell, rev_num_list
from utils.io import read_json, write_json, new_file_name
from parallel.putils import RefStruct, RecArgs
from parallel.monte_carlo import make_movements, make_movements_simple, PreOptimization, ConnectivityCheck
from parallel.freqs import read_displacements, read_frequencies, make_displacements
from shutil import copyfile
from utils.base import ktoht
import copy, zipfile, re

class RecordWritter(object):

    def __init__(self, buffer_size=1024**2*50):
        self.buffer_size = buffer_size

    def write_json(self, json_data, fn, final=False):
        with open(fn, 'w', self.buffer_size) as f:
            if final: json.dump(json_data, f, indent=4)
            else: json.dump(json_data, f)

    def write_clusters(self, fn, clus, append=False):
        if len(clus) == 0: return
        elif isinstance(clus[0], Cluster):
            self.write_clusters_base(fn, clus, append)
        elif isinstance(clus[0], ClusterAtSurface):
            self.write_clusters_at_surf(fn, clus, append)

    def write_clusters_base(self, fn, clus, append=False):
        with open(fn, 'a' if append else 'w', self.buffer_size) as f:
            for c in clus:
                f.write('%d\n%s = %15.8f\n' % (c.n, c.label, c.energy))
                self.write_coords(f, c)

    def write_clusters_at_surf(self, fn, clus, append=False):
        with open(fn, 'a' if append else 'w', self.buffer_size) as f:
            for c in clus:
                f.write('%d\n%s = %15.8f SURF = %d\n' % (c.n + c.surf.n, c.label,
                                                         c.energy, c.surf.n))
                self.write_coords(f, c.surf)
                self.write_coords(f, c)
    
    def write_forces(self, fn, fstructs, append=False):
        with open(fn, 'a' if append else 'w', self.buffer_size) as f:
            for c in fstructs:
                f.write('%d\n%s\n' % (c.n, c.label))
                for y in c.data:
                    f.write('%18.10f%18.10f%18.10f\n' % (y[0], y[1], y[2]))
    
    def read_forces(self, fn):
        lg = []
        with open(fn, 'r', self.buffer_size) as f:
            fl = f.readlines()
            ii = 0
            while ii < len(fl):
                fil = fl[ii].strip()
                if fil == "":
                    break
                n = int(fil)
                g = ForcesStruct(n)
                g.label = fl[ii + 1].strip()
                for i in xrange(n):
                    a = [float(x) for x in fl[ii + i + 2].strip().split(" ") if x != ""]
                    g.data[i, :] = np.array(a)
                ii += n + 2
                lg.append(g)
        return lg

    def write_coords(self, f, c):
        if c.forces is not None:
            for x, y, z in zip(c.elems, c.atoms, c.forces):
                f.write('%5s%14.8f%14.8f%14.8f %12.5f %12.5f %12.5f\n' %
                        (x, y[0], y[1], y[2], z[0], z[1], z[2]))
        else:
            for x, y in zip(c.elems, c.atoms):
                f.write('%5s%14.8f%14.8f%14.8f\n' % (x, y[0], y[1], y[2]))

    def write_surfs(self, fnx, fnj, clus, append=False):
        with open(fnx, 'a' if append else 'w', self.buffer_size) as f:
            for c in clus:
                assert c.n == 0
                f.write(('%d\n CELL =' + ' %.8f' * len(c.surf.cell) + ' (%.8f)\n') %
                        ((c.surf.n, ) + tuple(c.surf.cell) + (c.surf.cellz, )))
                self.write_coords(f, c.surf)
        with open(fnj, 'w', self.buffer_size) as f:
            for c in clus:
                json.dump({"space_group": c.surf.space_group,
                           "space_group_ref": list(c.surf.space_group_ref),
                           "unit_cell": list(c.surf.unit_cell),
                           "fix": c.surf.fix, "citation": "relaxed"}, f, indent=4)

class MonteCarloParams(object):
    Temperature = 300.0
    StepLength = 0.4
    LightStepLength = 0.6
    CurrentStepLength = 0.0
    StepUpdateIter = 50
    StepUpdateFactor = 0.9
    AcceptTarget = 0.5
    ShortDistanceFactor = 0.7
    MaxIter = 0
    RefuseSeparate = True
    DetailedBalance = False
    LightShell = False
    KeepCH3 = False
    KeepCO = False
    SwapSite = False
    SwapSiteMakeSpace = []
    # only support when the first elem is unique in the solid molecule
    SolidMove = []
    SwapSiteRate = 1.0
    FinalFactor = 1.2
    # for debug only
    WriteSteps = False

    @staticmethod
    def read_from(ipmc):
        mc = MonteCarloParams()
        mc.MaxIter = ipmc.get("max-iter", 0)
        mc.StepLength = ipmc.get("step-length", 0.4)
        mc.LightStepLength = ipmc.get("light-step-length", 0.6)
        mc.Temperature = ipmc.get("temperature", 300.0)
        mc.AcceptTarget = ipmc.get("accept-target", 0.5)
        mc.StepUpdateFactor = ipmc.get("step-update-factor", 0.9)
        mc.StepUpdateIter = ipmc.get("step-update-iter", 50)
        mc.RefuseSeparate = ipmc.get("refuse-separate", True)
        mc.ShortDistanceFactor = ipmc.get("short-distance-factor", 0.7)
        mc.DetailedBalance = ipmc.get("detailed-balance", False)
        mc.LightShell = ipmc.get("light-shell", False)
        mc.FinalFactor = ipmc.get("final-factor", 1.2)
        mc.SwapSiteRate = ipmc.get("swap-site-rate", 1.0)
        mc.SwapSite = ipmc.get("swap-site", False)
        mc.SwapSiteMakeSapce = ipmc.get("swap-site-make-space", [])
        mc.KeepCH3 = ipmc.get("keep-ch3", False)
        mc.KeepCO = ipmc.get("keep-co", False)
        mc.SolidMove = ipmc.get("solid-move", [])
        if mc.MaxIter == 0:
            mc = None
        return mc

class FiniteDifferenceParams(object):
    StepLength = 0.01 # angstrom

    @staticmethod
    def read_from(ipfd):
        fd = FiniteDifferenceParams()
        fd.StepLength = ipfd.get("step-length", 0.01)
        return fd

class ForcesStruct(object):
    n = 0
    label = ""
    data = np.array([], dtype=float)
    def __init__(self, n):
        self.n = n
        self.data = np.zeros((n, 3), dtype=float)

class RecordKeeper(object):

    def __init__(self, parent, idx=0, fn_pat="./%s.%s"):

        # partial updates
        self.structs = [] # write-only
        self.nebstructs = [] # write-only
        self.max = []
        self.dupl = []
        self.disp = []

        # finite difference
        self.forces = []

        # monte carlo
        self.mcsources = []

        # all updates
        self.props = {}
        self.runtime = {}
        self.local = []
        self.sources = []
        self.surfs = []

        # no updates
        self.refs = []

        self.sources_changed = False
        self.props_changed = False
        self.runtime_changed = False

        self.p = parent
        self.surf = None
        self.idx = idx
        self.args = None
        self.next_rec = None
        self.ne, self.eles = None, None
        self.dmax, self.dmax_rep = None, None
        self.fn_pat = fn_pat + "." + str(idx)
        self.flags = {}
        self.update_flags()
        self.writter = RecordWritter()
    
    def update_flags(self):
        for nn in [ "structs", "max", "dupl", "disp", "local", "sources", 
            "mcsources", "surfs", "nebstructs", "forces" ]:
            self.flags[nn] = len(getattr(self, nn))
        self.sources_changed = False
        self.props_changed = False
        self.runtime_changed = False

    def save(self, final=False):
        for nn in [ "structs", "max", "dupl", "disp", "mcsources", "nebstructs" ]:
            if len(getattr(self, nn)) > self.flags[nn]:
                self.writter.write_clusters(self.fn_pat % (nn, 'xyz'), 
                    getattr(self, nn)[self.flags[nn]:], append=True)
        for nn in [ "local" ]:
            if len(getattr(self, nn)) > self.flags[nn]:
                self.writter.write_clusters(self.fn_pat % (nn, 'xyz'), clus=getattr(self, nn))
        for nn in [ "surfs" ]:
            if len(getattr(self, nn)) > self.flags[nn]:
                self.writter.write_surfs(self.fn_pat % (nn, 'xyz'), self.fn_pat % (nn, 'json'), 
                    clus=getattr(self, nn))
        for nn in [ "forces" ]:
            if len(getattr(self, nn)) > self.flags[nn]:
                self.writter.write_forces(self.fn_pat % (nn, 'xyz'), 
                    getattr(self, nn)[self.flags[nn]:], append=True)
        for nn in [ "sources" ]:
            if getattr(self, nn + "_changed"):
                self.writter.write_clusters(self.fn_pat % (nn, 'xyz'), clus=getattr(self, nn))
        for nn in [ "props", "runtime" ]:
            if getattr(self, nn + "_changed") or final:
                self.writter.write_json(getattr(self, nn), 
                    fn=self.fn_pat % (nn, 'json'), final=final)
        self.update_flags()
    
    def open(self):
        for nn in [ "local", "sources", "mcsources" ]:
            if os.path.isfile(self.fn_pat % (nn, 'xyz')):
                setattr(self, nn, read_zip_clusters(self.fn_pat % (nn, 'xyz'), \
                    bufsize=self.writter.buffer_size))
            if self.surf is not None:
                setattr(self, nn, [ClusterAtSurface.from_cluster(c, self.surf) \
                    for c in getattr(self, nn)])
            for c in getattr(self, nn):
                labels = c.label.split(":")
                if "." not in labels[1]:
                    c.tid = int(labels[1])
                    c.tidx = [int(labels[1]), 0]
                else:
                    c.tid = int(labels[1].split(".")[0])
                    c.tidx = [int(p) for p in labels[1].split(".")]
                c.stidx = ".".join([str(p) for p in c.tidx])
        for nn in [ "forces" ]:
            if os.path.isfile(self.fn_pat % (nn, 'xyz')):
                setattr(self, nn, self.writter.read_forces(self.fn_pat % (nn, 'xyz')))
        self.structs = [None] * len(self.local)
        for it, t in enumerate(self.local):
            self.structs[it] = t
        for nn in [ "props", "runtime" ]:
            if os.path.isfile(self.fn_pat % (nn, 'json')):
                setattr(self, nn, read_json(fn=self.fn_pat % (nn, 'json'), iprint=False))
        self.update_flags()

    def init(self, args):
        self.args = args
        if "surface" in self.args:
            ipcs = self.args["surface"]
            self.surf = read_surface(ipcs["surface"])
            self.surf.space_group = ipcs["space_group"]
            self.surf.space_group_ref = np.array(ipcs["space_group_ref"])
            self.surf.unit_cell = ipcs["unit_cell"]
            self.surf.fix = ipcs.get("fix", "")
        else:
            self.surf = self.p.surf
        if self.surf is not None:
            self.surf_dof = list(set(range(0, self.surf.n)) - set(rev_num_list(self.surf.fix)))
            self.surf_dof.sort()
        if self.args.get("do-monte-carlo", False):
            self.mc = self.p.mc
        else:
            self.mc = None
        if self.args.get("do-finite-diff", False):
            self.fd = self.p.fd
            self.args["rescue_bad"] = False
        else:
            self.fd = None
        self.neb = None
        if os.path.isfile(self.fn_pat % ("runtime", "json")):
            self.open()
            self.init_refs(args.get("refs", []), self.p.ref_dir, self.surf)
        else:
            self.init_refs(args.get("refs", []), self.p.ref_dir, self.surf)
            self.init_runtime()
            self.init_sources(args["sources"], args.get("max_config", -1), self.p.ref_dir, self.surf)
            self.save()

    def init_refs(self, xrefs, pre_dir, surf):
        assert len(self.refs) == 0
        if xrefs is not None:
            for x in xrefs:
                xx = RefStruct.solve(x, idx=self.idx, pre_dir=pre_dir)
                cs = read_zip_clusters(xx.xdir)
                if surf is not None:
                    cs = [ClusterAtSurface.from_cluster(c, surf) for c in cs]
                for c in cs:
                    labels = c.label.split(":")
                    if "." not in labels[1]:
                        ptid = int(labels[1])
                        ptidx = [int(labels[1]), 0]
                    else:
                        ptid = int(labels[1].split(".")[0])
                        ptidx = [int(p) for p in labels[1].split(".")]
                    c.tid = [ xx.xsta, xx.xmul, xx.xname, xx.xidx, ptid ]
                    c.tidx = [ xx.xsta, xx.xmul, xx.xname, xx.xidx, ptidx ]
                    self.refs.append(c)

    def init_sources(self, xsources, max_config, pre_dir, surf):
        assert len(self.sources) == 0
        for src in xsources:
            if src[0] == "inherit":
                sclus = []
            elif src[0] == "create":
                sclus = self.p.p.filter(create=True)
            elif src[0] == "read":
                sclus = read_zip_clusters(src[1])
                if len(src) >= 3:
                    sclus = sclus[:src[2]]
                if surf is not None:
                    sclus = [ClusterAtSurface.from_cluster(sc, surf) for sc in sclus]
            elif src[0] == "restart":
                restart_dir = os.path.dirname(self.p.restart_dir)
                xx = RefStruct.solve(src[1], idx=self.idx, pre_dir=pre_dir)
                sclus = read_zip_clusters(xx.xdir)
                if len(src) >= 3:
                    sclus = sclus[:src[2]]
                if surf is not None:
                    sclus = [ClusterAtSurface.from_cluster(sc, surf) for sc in sclus]
                for c in sclus:
                    labels = c.label.split(":")
                    c.resdirx = restart_dir + "/%s/restart.tar.gz.%d.%s" % (xx.xpre, \
                        xx.xidx, labels[1])
            else:
                raise RuntimeError('Unknown structure source {}!'.format(src[0]))
            self.sources += sclus
        # if inherit, then we cannot use new structures
        # otherwise the tid will be confused
        xxsour = [x[0] for x in xsources]
        assert "inherit" not in xxsour or len(self.sources) == 0
        assert "restart" not in xxsour or len(list(set(xxsour))) == 1
        if max_config != -1:
            self.sources.sort(key=lambda x: x.energy)
            self.sources = self.sources[:max_config]
        if self.args.get("sweep", None) is not None and self.args["sweep"]["copy"]:
            sourn = self.args["sweep"]["n"]
            sourx = []
            for c in self.sources:
                for _ in range(sourn):
                    sourx.append(copy.deepcopy(c))
            self.sources = sourx
        if len(self.sources) == 0 or not self.sources[0].label.startswith("CONN"):
            for ic, c in enumerate(self.sources):
                c.tid = ic
                c.tidx = [ic, 0]
                c.stidx = "%d.0" % ic
                if hasattr(c, 'mag') and c.mag is not None:
                    c.label = "MAG:%s:%.2f:0" % (c.stidx, c.mag)
                else:
                    c.label = "SRC:%s:0" % c.stidx
        else:
            # sources startswith CONN => NEB calculation
            pts = []
            ptn = [] # number of (all) images -1
            for c in self.sources:
                assert c.label.startswith("CONN")
                pid = c.label.split(":")[1].split(".")[0]
                if pid not in pts:
                    ic = len(pts)
                    pts.append(pid)
                    ptn.append(0)
                else:
                    ic = pts.index(pid)
                    ptn[ic] += 1
                c.tid = ic
                c.tidx = [ic, ptn[ic]]
                c.stidx = "%d.%d" % (ic, ptn[ic])
                c.label = "%s:%s:%.2f:0" % (c.label.split(":")[0], c.stidx,
                                            c.mag if c.mag is not None else 0.0)
            self.neb = True
        self.sources_changed = True
        for src in self.sources:
            if self.neb is not None and src.tidx[1] != 0:
                continue
            self.runtime["states"][str(src.tid)] = 0
            self.runtime["steps"][src.stidx] = 0
            if self.mc is not None:
                self.runtime["mc"][str(src.tid)] = [0, 0, 0, self.mc.StepLength, -1, -1] + [0] * 7
            if self.fd is not None:
                self.runtime["fd"][str(src.tid)] = self.find_degree_of_freedom(src)
            if hasattr(src, "resdirx"):
                self.runtime["restart"][src.stidx] = src.resdirx
            else:
                self.runtime["restart"][src.stidx] = None

    def find_source(self, tid):
        istart, iend = None, len(self.sources)
        for isrc, src in enumerate(self.sources):
            if src.tid == tid and istart is None:
                istart = isrc
            elif src.tid != tid and istart is not None:
                iend = isrc
                break
        assert istart is not None
        if iend == istart + 1:
            return istart
        else:
            return slice(istart, iend)
    
    def find_degree_of_freedom(self, stru):
        if self.surf is not None:
            x = self.surf_dof + range(self.surf.n, stru.n + self.surf.n)
        else:
            x = range(0, stru.n)
        y = range(0, len(x) * 6 + 1)
        return [y, x]

    def find_mcsource(self, tid, mcstep):
        for isrc, src in enumerate(self.mcsources):
            if src.tidx == [tid, mcstep]:
                return isrc
        return None

    # states: 0 not start or running; 1 converged; 2 bad struct;
    # 3 not converged; 4 max step; 5 converged duplicate; 6 non-converged dupl
    # mc: [ totalrej, currej, curiter, cursteplen, reftid, refmcstep ]
    # fd: [ curstep, totalstep ]
    def init_runtime(self):
        self.runtime = {"corr_list": {}, "times": {}, "states": {}, "steps": {}, "restart": {}}
        if self.mc is not None:
            self.runtime["mc"] = {}
        if self.fd is not None:
            self.runtime["fd"] = {}
        self.runtime_changed = True

    def inherit_source(self, struct, xtype):
        if self.next_rec is not None:
            for src in self.next_rec.args["sources"]:
                if src[0] == "inherit":
                    if xtype == src[1] or xtype in src[1]:
                        tids = [x.tid for x in self.next_rec.sources]
                        for tid in range(0, len(tids) + 1):
                            if tid not in tids:
                                break
                        stid = str(tid)
                        prev_tidx = [self.idx, struct.tidx]
                        prev_stidx = struct.stidx
                        struct = copy.deepcopy(struct)
                        struct.tid = tid
                        struct.tidx = [tid, 0]
                        struct.stidx = "%d.0" % tid
                        struct.label = "%s:%s:%s" % (struct.label.split(":")[0], struct.stidx, \
                            ":".join(struct.label.split(":")[2:]))
                        if self.next_rec.surf is not self.surf:
                            sc = ClusterAtSurface.from_cluster(struct.to_cluster(), \
                                self.next_rec.surf)
                            sc.tid = struct.tid
                            sc.tidx = struct.tidx
                            sc.stidx = struct.stidx
                            struct = sc
                        self.next_rec.sources.append(struct)
                        self.next_rec.sources.sort(key=lambda x: x.energy)
                        self.next_rec.runtime["states"][stid] = 0
                        self.next_rec.runtime["steps"][struct.stidx] = 0
                        if self.next_rec.mc is not None:
                            self.next_rec.runtime["mc"][stid] = [0, 0, 0, \
                                self.next_rec.mc.StepLength, -1, -1] + [0] * 7
                        if self.next_rec.fd is not None:
                            self.next_rec.runtime["fd"][stid] = self.next_rec.find_degree_of_freedom(struct)
                        self.next_rec.runtime["restart"][struct.stidx] = \
                            self.runtime["restart"][prev_stidx]
                        self.next_rec.props[struct.stidx] = {'prev_tidx': prev_tidx}
                        self.next_rec.props_changed = True
                        self.next_rec.runtime_changed = True
                        self.next_rec.sources_changed = True
                        break
            self.next_rec.check_sources()

    def check_sources(self):
        max_config = self.args.get("max_config", -1)
        if max_config != -1 and len(self.sources) > max_config:
            for isrc in range(len(self.sources) - 1, max_config, -1):
                src = self.sources[isrc]
                stid = str(src.tid)
                if self.runtime["states"][stid] == 0 and \
                    src.tid not in self.runtime["corr_list"].values():
                    del self.sources[isrc]
                    del self.runtime["states"][stid]
                    del self.runtime["steps"][src.stidx]
                    del self.runtime["restart"][src.stidx]
                    del self.props[src.stidx]
                    if self.mc is not None:
                        del self.runtime["mc"][stid]
                    if self.fd is not None:
                        del self.runtime["fd"][stid]

    def handle_forces(self, forces, sstr_id, xstr, dofl, forceslen):
        from utils.base import evtoht, angtobohr, amutoau
        spdc = 137.036
        auvtocm = 1.0 / 5.291772108E-9
        for_ref = None
        la = len(dofl)
        ll = la * 3
        for_dsp = [None] * len(dofl) * 6
        dofl = np.array(dofl)
        for ff in forces:
            a, b = ff.label.split(".")
            if a != sstr_id:
                continue
            b = int(b)
            if b == 0:
                for_ref = ff.data
            else:
                for_dsp[b - 1] = ff.data
        exx = xstr.elems
        if self.surf is not None:
            exx = np.array(list(self.surf.elems) + list(exx))
        masses = np.array([ np.sqrt(AtomicWeight.x[exx[x]]) for x in dofl ])
        msl = masses.reshape((la, 1))
        dmat = np.zeros((ll, ll))
        for i in xrange(0, ll):
            ddsp = (for_dsp[i * 2 + 1] - for_dsp[i * 2]) / msl
            dmat[i] = ddsp.flatten() / masses[i / 3]
        fac = 1.0 / (2 * forceslen) * evtoht / angtobohr / angtobohr / amutoau
        dmat *= fac
        for i in xrange(0, ll):
            for j in xrange(i + 1, ll):
                x = 0.5 * (dmat[i, j] + dmat[j, i])
                dmat[i, j] = x
                dmat[j, i] = x
        eigs, evec = np.linalg.eigh(dmat)
        imgs = eigs < 0
        eigs = np.sqrt(np.abs(eigs)) * (1.0 / (2.0 * np.pi * spdc) * auvtocm)
        eigs = np.array([-e if i else e for i, e in zip(imgs, eigs)])
        msll = (msl + np.zeros((la, 3))).flatten()
        facf = 1.0 / (angtobohr * np.sqrt(amutoau))
        for i in xrange(0, ll):
            x = evec[i] / msll
            evec[i] = x / np.linalg.norm(x) * facf
        return eigs.tolist(), evec.tolist()

    def handle_mc_newmove(self, acstr, str_id, mcst, accmp):
        sstr_id = str(str_id)
        mcst[2] += 1
        if mcst[2] == self.mc.MaxIter:
            self.runtime["states"][sstr_id] = 1
            yname = ["mc.end.%d>%d" % (mcst[2], mcst[0])]
        else:
            self.runtime["states"][sstr_id] = 0
            ac_stidx = "%d.%d" % (str_id, mcst[2])
            self.runtime["steps"][ac_stidx] = 0
            self.runtime["restart"][ac_stidx] = None
            self.mc.CurrentStepLength = mcst[3]
            mstr = make_movements(acstr, self.mc, self.p.p.get_create_opts(self.p.p.ip))
            mstr.tid = str_id
            mstr.tidx = [str_id, mcst[2]]
            mstr.stidx = "%d.%d" % (str_id, mcst[2])
            mlabel = acstr.label.split(":")
            mlabel[1] = mstr.stidx
            mstr.label = ":".join(mlabel)
            self.sources[self.find_source(str_id)] = mstr
            self.sources_changed = True
            yname = [ "mc.move.%s.%d>%d" % (accmp, mcst[2], mcst[0]) ]
        if self.mc.StepUpdateIter != 0 and mcst[2] % self.mc.StepUpdateIter == 0:
            rejrt = float(mcst[1]) / self.mc.StepUpdateIter
            rmcstp = mcst[3]
            if rejrt > self.mc.AcceptTarget / self.mc.StepUpdateFactor:
                mcst[3] = mcst[3] * self.mc.StepUpdateFactor
            elif rejrt < self.mc.AcceptTarget * self.mc.StepUpdateFactor:
                mcst[3] = mcst[3] / self.mc.StepUpdateFactor
            if rmcstp != mcst[3]:
                yname.append("mc.step.update.%.2f=>%.2f" % (rmcstp, mcst[3]))
            mcst[1] = 0
        return yname

    def handle_structs(self, structs, str_id, str_step, str_mcstep):
        for i, l in enumerate(structs):
            l.tid = str_id
            l.tidx = [str_id, str_mcstep]
            l.stidx = "%d.%d" % (str_id, str_mcstep)
            if hasattr(l, 'mag') and l.mag is not None:
                l.label = 'MAG:%s:%.2f:%d' % (l.stidx, l.mag, str_step + i)
            else:
                l.label = 'STR:%s:%d' % (l.stidx, str_step + i)
        # only two types; none and conv
        if self.args.get("filter_type", "none") == "none":
            return None
        if self.ne is None:
            self.eles, self.ne = elem_num(structs[-1].elems)
            self.dmax = self.p.ipfl["max_diff"]
            self.dmax_rep = self.p.ipfl["max_diff_report"]
        l = structs[-1]
        for cc in self.structs + self.refs:
            if cc.tidx == [str_id, str_mcstep]: continue
            if isinstance(cc.tid, int):
                if cc.stidx in self.props and "dupl" in self.props[cc.stidx]:
                    continue
                if self.runtime["states"][str(cc.tid)] == 0 and (self.mc is None or
                    self.runtime["mc"][str(cc.tid)][2] == cc.tidx[1]):
                    continue
            if l.energy + 1E-6 < cc.energy: continue
            if self.surf is not None:
                v, _ = surface_compare(l, cc, self.ne, self.eles, self.dmax_rep)
            else:
                v, _ = at_comp(l.atoms, cc.atoms, self.ne, self.eles, self.dmax_rep)
            if v < self.dmax:
                return cc.tidx
        return None

    def handle_neb_structs(self, lproc, pstep, cmtl, str_id):
        all_structs = []
        zf = zipfile.ZipFile(lproc + '/trajectory.zip', 'r')
        namel = zf.namelist()
        # sort filenames
        name_idx = [re.findall(r'[0-9]+', name) for name in namel]
        name_idx = [[int(i) for i in n] for n in name_idx]
        idx = sorted(range(len(name_idx)), key=name_idx.__getitem__)
        namel = [namel[i] for i in idx]
        for name in namel:
            if re.findall(r'[0-9]+', name)[0] < pstep:
                continue
            xclul = read_clusters(name, zf=zf, iprint=False)
            if self.surf is not None:
                xclul = [ClusterAtSurface.from_cluster(sc, self.surf) for sc in xclul]
            for ii, l in enumerate(xclul):
                l.tid = str_id
                l.tidx = [str_id, ii]
                l.stidx = "%d.%d" % (str_id, ii)
                lstep = int(l.label.split("=")[0].strip().split(":")[-1])
                l.label = "%s:%s:%.2f:%d" % (cmtl[0], l.stidx, l.mag if l.mag is not None
                                             else 0.0, lstep)
            all_structs.extend(xclul)
        zf.close()
        return all_structs

    def handle_neb(self, str_id, lproc, xtype):
        assert xtype in ["step", "conv"]
        sstr_id = str(str_id)
        str_stidx = "%d.0" % str_id
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        structs = read_zip_clusters(lproc + '/final.xyz', iprint=False)
        if self.surf is not None:
            structs = [ClusterAtSurface.from_cluster(sc, self.surf) for sc in structs]
        sours = self.sources[self.find_source(str_id)]
        lstep = int(structs[0].label.split("=")[0].strip().split(":")[-1])
        nstep = lstep + 1 - self.runtime["steps"][str_stidx]
        cmtl = sours[0].label.split("=")[0].strip().split(":")[0:2]
        for ii, l in enumerate(structs):
            l.tid = str_id
            l.tidx = [str_id, ii]
            l.stidx = "%d.%d" % (str_id, ii)
            l.label = "%s:%s:%.2f:%d" % (cmtl[0], l.stidx, l.mag if l.mag is not None
                                         else 0.0, lstep)
        all_structs = self.handle_neb_structs(lproc, self.runtime["steps"][str_stidx], \
            cmtl, str_id)
        self.nebstructs.extend(all_structs)
        self.runtime["steps"][str_stidx] += nstep
        xstr_step = "%d:%d" % (nstep, self.runtime["steps"][str_stidx])
        self.sources[self.find_source(str_id)] = structs
        self.sources_changed = True
        if xtype == "conv":
            self.runtime["states"][sstr_id] = 1
            self.local.extend(structs)
            self.local.sort(key=lambda x: x.tidx)
            # self.inherit_source(xstr, "local")
            xname = "neb.conv.%s" % xstr_step
        elif self.runtime["steps"][str_stidx] >= self.args["max_step"]:
            self.runtime["states"][sstr_id] = 4
            self.max.extend(structs)
            # self.inherit_source(xstr, "max")
            xname = "neb.max.%s" % xstr_step
        else:
            xname = "neb.step.%s" % xstr_step
        self.props_changed = True
        self.runtime_changed = True
        return [xname, max([x.energy for x in structs]) - structs[0].energy]

    def handle_bad_neb(self, str_id, xtype):
        assert xtype in ["not_conv", "bad_struct"]
        sstr_id = str(str_id)
        str_stidx = "%d.0" % str_id
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        if xtype == 'not_conv':
            self.runtime["states"][sstr_id] = 3
        elif xtype == 'bad_struct':
            self.runtime["states"][sstr_id] = 2
        if xtype == 'bad_struct':
            xname = "neb.bad.%d" % self.runtime["steps"][str_stidx]
        else:
            xname = "neb.n.conv.%d" % self.runtime["steps"][str_stidx]
        self.runtime_changed = True
        self.props_changed = True
        return xname
    
    def handle_relax(self, str_id, lproc, xtype, rfdstep):
        assert xtype in [ "energy", "step", "conv" ]
        sstr_id = str(str_id)
        str_mcstep = 0 if self.mc is None else self.runtime["mc"][sstr_id][2]
        str_fdstep = rfdstep
        if self.fd is not None:
            str_mcstep = str_fdstep
        str_stidx = "%d.%d" % (str_id, str_mcstep)
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        if os.path.isfile(lproc + '/trajectory.xyz'):
            structs = structsb = read_zip_clusters(lproc + '/trajectory.xyz', iprint=False)
            structs = structs[self.runtime["steps"][str_stidx]:]
        else:
            print ("no trajectory.xyz!")
            structs = []
        if len(structs) == 0:
            print ("no new structs! STR_ID = %s, LPROC = %s" % (str_stidx, lproc))
            if xtype == "conv" and len(structsb) != 0:
                structs = [ structsb[-1] ]
                if structs[0] is None: # probably filesystem error
                    self.runtime["steps"][str_stidx] = 0
                    return [ "r.nostr", None ]
            elif xtype == "energy":
                structs = [self.sources[self.find_source(str_id)].to_cluster()]
            else:
                return [ "r.nostr", None ]
        if self.surf is not None:
            structs = [ClusterAtSurface.from_cluster(sc, self.surf) for sc in structs]
        dupl = self.handle_structs(structs, str_id, self.runtime["steps"][str_stidx], str_mcstep)
        self.structs += structs
        xstr = structs[-1]
        self.props[str_stidx]["energy"] = xstr.energy
        if os.path.isfile(lproc + '/cell.txt'):
            xcf = open(lproc + '/cell.txt', "r").read().strip()
            xcell = [ float(x) for x in xcf.split(" ") if len(x) != 0 ]
            ori_cell = np.array(xstr.surf.cell)
            new_cell = np.array(xcell)
            xstr.surf.cell = new_cell
            xstr.surf.cellz *= new_cell[-1] / ori_cell[-1]
            orimat, newmat = to_cellmat(ori_cell), to_cellmat(new_cell)
            ucmat = to_cellmat(xstr.surf.unit_cell)
            ucmatx = ucmat * newmat / orimat
            ucmatx[np.isnan(ucmatx)] = ucmat[np.isnan(ucmatx)]
            xstr.surf.unit_cell = to_dcell(ucmatx)
            has_cell = True
        else:
            has_cell = False
        if os.path.isfile(lproc + '/bader.json'):
            bad = read_json(fn=lproc + '/bader.json', iprint=False)
            self.props[str_stidx].update(bad)
        if self.runtime["steps"][str_stidx] == 0:
            self.props[str_stidx]["initial_energy"] = structs[0].energy
        if xtype == "energy":
            self.runtime["states"][sstr_id] = 1
            self.inherit_source(xstr, "structs")
            xname = "energy.okay"
        else:
            # 0-empty 1-conv 2-bad.str 3-n.conv 4-max 5-dp.conv 6-dp.step
            self.runtime["steps"][str_stidx] += len(structs)
            xstr_step = "%d:%d" % (len(structs), self.runtime["steps"][str_stidx])
            self.sources[self.find_source(str_id)] = xstr
            self.sources_changed = True
            if dupl is not None:
                self.runtime["states"][sstr_id] = 5 if xtype == "conv" else 6
                self.props[str_stidx]["dupl"] = dupl
                self.dupl.append(xstr)
                if isinstance(dupl[0], int):
                    xname = "r.dupl.%s~%d.%d" % (xstr_step, dupl[0], dupl[1])
                else:
                    xname = "r.dupl.%s~ref.%s.%d.%d" % (xstr_step, dupl[0], dupl[-1][0], dupl[-1][1])
                self.inherit_source(xstr, "dupl")
            elif xtype == "conv":
                self.runtime["states"][sstr_id] = 1
                self.local.append(xstr)
                if has_cell and xstr.n == 0:
                    self.surfs.append(xstr)
                self.local.sort(key=lambda x: x.energy)
                self.inherit_source(xstr, "local")
                xname = "r.conv.%s" % xstr_step
            elif self.runtime["steps"][str_stidx] >= self.args["max_step"]:
                self.runtime["states"][sstr_id] = 4
                self.max.append(xstr)
                self.inherit_source(xstr, "max")
                xname = "r.max.%s" % xstr_step
            else:
                xname = "r.step.%s" % xstr_step
        # mc: [ totalrej, currej, curiter, cursteplen, reftid, refmcstep,  ]
        # mc: (5 +) 1 converged; 2 bad struct;
        # mc: 3 not converged; 4 max step; 5 converged duplicate; 6 non-converged dupl
        # mc: 7 separated
        if self.mc is not None and self.runtime["states"][sstr_id] != 0:
            mcst = self.runtime["mc"][sstr_id]
            mstu = self.runtime["states"][sstr_id]
            # first energy eval step
            if not self.mc.DetailedBalance:
                acstr = xstr
            else:
                acstr = self.mcsources[self.find_mcsource(str_id, str_mcstep)]
            acc = True
            accmp = "init"
            if mcst[2] != 0:
                rtid, rmcstep = mcst[4:6]
                rstidx = "%d.%d" % (rtid, rmcstep)
                rener = self.props[rstidx]["energy"]
                if rener is None:
                    rener = xstr.energy
                rd = np.random.random()
                rtemp = self.mc.Temperature * ktoht
                acc = np.exp(-(xstr.energy - rener) / rtemp) > rd
                accmp = "up" if xstr.energy > rener else "down"
                if self.mc.RefuseSeparate:
                    cck = ConnectivityCheck(xstr)
                    if not cck.check():
                        acc = False
                        accmp = "sprj"
                        mstu = 7
                if not acc:
                    acstr = self.mcsources[self.find_mcsource(rtid, rmcstep)]
                    mcst[0] += 1
                    mcst[1] += 1
                if accmp == "up":
                    accmp += "ac" if acc else "rj"
            if acc:
                mcst[4:6] = [str_id, mcst[2]]
                if not self.mc.DetailedBalance:
                    self.mcsources.append(xstr)
            mcst[5 + mstu] += 1
            yname = self.handle_mc_newmove(acstr, str_id, mcst, accmp)
            xname = [ xname ] + yname
        self.props_changed = True
        if self.fd is not None and self.runtime["states"][sstr_id] != 0:
            self.props_changed = False
            ntot = len(self.runtime["fd"][sstr_id][1]) * 6
            self.runtime["states"][sstr_id] = 0
            if "initial_energy" in self.props[str_stidx]:
                del self.props[str_stidx]["initial_energy"]
            if "energy" in self.props[str_stidx]:
                del self.props[str_stidx]["energy"]
            if str_fdstep == 0:
                yname = "fd.init.%d.%d" % (str_fdstep, ntot)
            elif len(self.runtime["fd"][sstr_id][0]) > 1:
                yname = "fd.step.%d.%d" % (str_fdstep, ntot)
            else:
                yname = "fd.end.%d.%d" % (str_fdstep, ntot)
                if self.runtime["fd"][sstr_id][0][0] == str_fdstep:
                    self.runtime["states"][sstr_id] = 1
            # store forces ...
            xforces = self.writter.read_forces(lproc + '/forces.xyz')[0]
            xforces.label = str_stidx
            self.forces.append(xforces)
            if str_fdstep in self.runtime["fd"][sstr_id][0]:
                self.runtime["fd"][sstr_id][0].remove(str_fdstep)
            zz_stidx = "%d.0" % (str_id)
            if self.runtime["states"][sstr_id] == 1:
                frs, disp = self.handle_forces(self.forces, sstr_id, xstr,
                    self.runtime["fd"][sstr_id][1], self.fd.StepLength)
                self.props[zz_stidx]["freqs"] = frs
                self.props[zz_stidx]["displacements"] = disp
                self.props_changed = True
            xname = [ xname, yname ]
        return [ xname, xstr.energy ]

    def handle_bad_relax(self, str_id, xtype):
        assert xtype in [ "not_conv", "bad_struct" ]
        sstr_id = str(str_id)
        str_mcstep = 0 if self.mc is None else self.runtime["mc"][sstr_id][2]
        str_mcstep = str_mcstep if self.fd is None else self.runtime["fd"][sstr_id][0]
        if self.fd is not None:
            if xtype == 'bad_struct':
                xname = "r.bad.0"
            else:
                xname = "r.n.conv.0"
            return xname
        str_stidx = "%d.%d" % (str_id, str_mcstep)
        if self.p.ippa.get("rescue_bad", False):
            self.runtime["states"][sstr_id] = 0
            xstr = self.sources[self.find_source(str_id)]
            mstr = make_movements_simple(xstr, self.p.ippa["rescue_length"])
            self.sources[self.find_source(str_id)] = mstr
            self.sources_changed = True
        else:
            if xtype == 'not_conv':
                self.runtime["states"][sstr_id] = 3
                if self.p.ippa.get("continue_scf", False):
                    self.runtime["states"][sstr_id] = 0
            elif xtype == 'bad_struct':
                self.runtime["states"][sstr_id] = 2
            if self.fd is not None:
                self.runtime["states"][sstr_id] = 0
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        self.props[str_stidx]["energy"] = None
        if xtype == 'bad_struct':
            xname = "r.bad.%d" % self.runtime["steps"][str_stidx]
        else:
            xname = "r.n.conv.%d" % self.runtime["steps"][str_stidx]
        if self.mc is not None:
            mcst = self.runtime["mc"][sstr_id]
            mstu = self.runtime["states"][sstr_id]
            mcst[5 + mstu] += 1
            if mcst[2] == 0:
                yname = [ "mc.init.failed" ]
            else:
                rtid, rmcstep = mcst[4:6]
                rstidx = "%d.%d" % (rtid, rmcstep)
                rtm = self.find_mcsource(rtid, rmcstep)
                while rtm is None:
                    rmcstep -= 1
                    rtm = self.find_mcsource(rtid, rmcstep)
                acstr = self.mcsources[rtm]
                mcst[0] += 1
                mcst[1] += 1
                accmp = "flrj"
                yname = self.handle_mc_newmove(acstr, str_id, mcst, accmp)
            xname = [ xname ] + yname
        self.runtime_changed = True
        self.props_changed = True
        return xname
    
    def handle_bad_freqs(self, str_id):
        sstr_id = str(str_id)
        str_mcstep = 0 if self.mc is None else self.runtime["mc"][sstr_id][2]
        str_stidx = "%d.%d" % (str_id, str_mcstep)
        self.runtime["states"][sstr_id] = 2
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        self.props[str_stidx]["freqs"] = None
        self.runtime_changed = True
        self.props_changed = True
        return "freq.bad"

    def handle_freqs(self, str_id, lproc):
        sstr_id = str(str_id)
        str_mcstep = 0 if self.mc is None else self.runtime["mc"][sstr_id][2]
        str_stidx = "%d.%d" % (str_id, str_mcstep)
        self.runtime["states"][sstr_id] = 1
        if str_stidx not in self.props:
            self.props[str_stidx] = {}
        frs = read_frequencies(lproc + '/vibspectrum')
        self.props[str_stidx]["freqs"] = frs
        disp = read_displacements(lproc + '/vib_normal_modes')
        self.props[str_stidx]["displacements"] = disp
        self.runtime_changed = True
        self.props_changed = True
        xstr = self.sources[self.find_source(str_id)]
        if self.args.get("freq_step", 0.0) != 0.0:
            md = make_displacements(xstr, frs, disp, self.args["freq_step"])
            self.disp += md
        self.inherit_source(xstr, "source")
        return "freq.okay"

class Event(object):

    def __init__(self, proc_id=-1, name="NTHG", str_id=-1, energy=None, idx=-1):
        self.proc_id = proc_id
        self.name = name
        self.str_id = str_id
        self.idx = idx
        self.energy = energy
        self.time = time.strftime("%y-%m-%d %H:%M:%S")
    
    def __repr__(self):
        if self.str_id == -1:
            return "%s | [%3d] %-15s" % (self.time, self.proc_id, self.name)
        else:
            if self.energy is not None:
                return "%s | [%3d] %-15s (%2d.%4d) E = %15.6f" % (self.time, self.proc_id, self.name, 
                    self.idx, self.str_id, self.energy)
            else:
                return "%s | [%3d] %-15s (%2d.%4d)" % (self.time, self.proc_id, self.name, 
                    self.idx, self.str_id)
    
class EventLogger(object):

    def __init__(self, fn):
        self.fn = fn
        self.events = []
    
    def add(self, event):
        if isinstance(event.name, list):
            for name in event.name:
                kevent = copy.deepcopy(event)
                kevent.name = name
                self.events.append(kevent)
        else:
            self.events.append(event)
    
    def add_header(self, rand_seed):
        if os.path.isfile(self.fn): return
        self.events.append("# random seed = %d" % rand_seed)
        lenx = 17 + 3 + 5 + 1 + 15 + 3 + 7 + 5 + 15 + 4
        self.events.append("-" * lenx)
        self.events.append("%17s   %5s %-15s (%7s) %-19s" % 
            ("time", "proc", "description", "struct", "energy"))
        self.events.append("-" * lenx)
    
    def add_summary(self, istr_rem, istr_tot, idx, idx_tot, iproc_act, iproc_tot, wait=False):
        self.events.append("%s +- STR %2d:%04d/%04d --- PROC %04d/%04d --- IDX %2d/%2d %s ---"
            % (time.strftime("%y-%m-%d %H:%M:%S"), idx, istr_rem, istr_tot, 
            iproc_act, iproc_tot, idx, idx_tot, "waiting" if wait else ""))

    def save(self):
        f = open(self.fn, "a", 1024**2/10)
        f.write("\n".join(map(str, self.events)) + "\n")
        f.close()
        self.events = []

class Shortlog(object):
    def __init__(self, fn):
        self.fn = fn
        self.args = ""
  
    def write(self, *args):
        sargs = " ".join(map(str, args))
        if sargs != self.args:
            self.args = sargs
            f = open(self.fn, "w")
            f.write(self.args)
            f.close()

class ParallelRun(object):

    def __init__(self, parent, ip):
        self.ip = ip
        self.p = parent
        self.ippa = self.ip["parallel"]
        self.ipfl = self.ip["filtering-parallel"]
        self.surf = None
        self.init_surface()
        self.mc = None
        self.init_mc()
        self.fd = None
        self.init_fd()
        self.idle_time = self.ippa.get("idle_time", 8.0)
        self.max_proc_time = self.ippa.get("max_no_repsonce_time", 300)
        self.max_master_time = self.ippa.get("max_time", None)
        self.proc_dir = self.ippa["proc_dir"] + "/" + self.ippa["proc_type"]
        self.restart_dir = self.ippa["restart_dir"]
        self.ref_dir = self.ippa["ref_dir"]
        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)
        self.args = self.ippa["arguments"]
        self.records = []
        for iarg, args in enumerate(self.args):
            self.records.append(RecordKeeper(parent=self, idx=iarg, fn_pat=self.p.para_x_name))
            self.records[iarg].init(args)
            if iarg != 0:
                self.records[iarg - 1].next_rec = self.records[iarg]
        self.log_fn = self.p.para_log_name + '.0'
        self.shortlog_fn = self.p.para_shlog_name + '.0'
        self.logger = EventLogger(fn=self.log_fn)
        self.shlogger = Shortlog(self.shortlog_fn)
    
    def init_surface(self):
        if "creation-surface" in self.ip:
            ipcs = self.ip["creation-surface"]
            self.surf = read_surface(ipcs["surface"])
            self.surf.space_group = ipcs["space_group"]
            self.surf.space_group_ref = np.array(ipcs["space_group_ref"])
            self.surf.unit_cell = ipcs["unit_cell"]
            self.surf.fix = ipcs.get("fix", "")
    
    def init_mc(self):
        if "monte-carlo" in self.ip:
            ipmc = self.ip["monte-carlo"]
            self.mc = MonteCarloParams.read_from(ipmc)
    
    def init_fd(self):
        if "finite-diff" in self.ip:
            ipfd = self.ip["finite-diff"]
            self.fd = FiniteDifferenceParams.read_from(ipfd)
    
    def has_finished(self):
        # check is it necessary to continue
        for rec in self.records:
            if 0 in rec.runtime["states"].values():
                return False
        return True
    
    def close_proc(self, name, proc_id, str_id=-1, resp_time=0, lproc=None, 
        energy=None, idx=-1):
        if lproc is not None and os.path.isfile(lproc + '/RESPONSE'):
            try:
                resp_name = new_file_name(lproc + '/RESPONSE.old.' + str(str_id))
                os.rename(lproc + '/RESPONSE', resp_name)
                f = open(lproc + '/FILELIST', 'a')
                f.write(resp_name[len(lproc) + 1:] + '\n')
                f.close()
            except IOError, OSError:
                pass
        sproc_id = str(proc_id)
        for rec in self.records:
            if sproc_id in rec.runtime["times"] and \
                len(rec.runtime["times"][sproc_id][-1]) == 2:
                rec.runtime["times"][sproc_id][-1] += [ time.time(), resp_time, name ]
                rec.runtime_changed = True
            if sproc_id in rec.runtime["corr_list"]:
                del rec.runtime["corr_list"][sproc_id]
                rec.runtime_changed = True
        self.logger.add(Event(proc_id=proc_id, name=name, str_id=str_id, 
            energy=energy, idx=idx))
    
    def assign_proc(self, xtype, lproc, xrec=None, xtid=-1, proc_id=-1, rstepid=-1):
        if xtype == "sleep":
            xtask = "SLEEP"
            self.logger.add(Event(proc_id=proc_id, name='sleep'))
        else:
            assert xrec is not None and xtid != -1 and proc_id != -1
            sproc_id = str(proc_id)
            xrec.runtime["corr_list"][sproc_id] = xtid
            if sproc_id not in xrec.runtime["times"]:
                xrec.runtime["times"][sproc_id] = []
            xrec.runtime["times"][sproc_id].append([ time.time(), xtid ])
            xrec.runtime_changed = True

            xxtid = str(xtid)
            xmcstep = 0 if xrec.mc is None else xrec.runtime["mc"][xxtid][2]
            if xrec.fd is not None:
                xmcstep = rstepid
                xrec.runtime["corr_list"][sproc_id] = (xtid, xmcstep)
            xfdstep = xmcstep
            xstidx = "%d.%d" % (xtid, xmcstep)

            if xrec.fd is not None and xfdstep != 0:
                xrec.runtime["restart"][xstidx] = xrec.runtime["restart"]["%d.0" % xtid]
                xrec.runtime["steps"][xstidx] = 0

            if xrec.runtime["restart"][xstidx] is not None:
                if xrec.runtime["steps"][xstidx] != 0 or xrec.args["read_restart"] or \
                    (xrec.fd is not None and xfdstep != 0):
                    xrdir = xrec.runtime["restart"][xstidx]
                    if xrdir[0] != '/': rdir = "../" + xrdir
                    else: rdir = xrdir
                    if os.path.islink(lproc + '/restart.tar.gz') or \
                        os.path.isfile(lproc + '/restart.tar.gz'):
                        os.remove(lproc + '/restart.tar.gz')
                    if os.path.isfile(xrdir):
                        # this is because the symlink itself is put in a deeper folder
                        os.symlink(rdir, lproc + '/restart.tar.gz')
                    else:
                        print ("no restart! RDIR = %s" % (rdir))
                        xrec.runtime["steps"][xstidx] = 0
                        if xrec.args["read_restart"]:
                            xrec.args["read_restart"] = False
            xstr = xrec.sources[xrec.find_source(xtid)]
            
            if xrec.mc is not None and xrec.mc.DetailedBalance and xrec.runtime["steps"][xstidx] == 0:
                if xrec.args["task_type"] == "relax" or xrec.args["task_type"] == "energy":
                    xrec.mcsources.append(xstr)

            start_type = "INIT"
            if xrec.args["read_restart"]: 
                start_type = "CONTINUE"
            elif xrec.fd is not None and xfdstep != 0:
                start_type = "CONTINUE"
            if xrec.args["task_type"] == "relax" or xrec.args["task_type"] == "energy":
                xtask = [ "RELAX", start_type if xrec.runtime["steps"][xstidx] == 0 else "CONTINUE" ]
            elif xrec.args["task_type"] == "neb":
                xtask = [ "NEB", start_type if xrec.runtime["steps"][xstidx] == 0 else "CONTINUE" ]
            elif xrec.args["task_type"] == "freq":
                xtask = [ "FREQ", start_type ]
            else:
                raise RuntimeError("Unknown task type: %s" % xrec.args["task_type"])
            
            xstrz = xstr[0] if isinstance(xstr, list) else xstr
            if isinstance(xstrz, ClusterAtSurface):
                dargs = RecArgs(xrec.args["others"])
                dargs.dict["cell"] = ":".join(map(str, xstrz.surf.cell))
                xrec.args["others"] = str(dargs)

            if xtask[1] == "INIT":
                if isinstance(xstr, list):
                    xrec.writter.write_clusters(lproc + '/coord.xyz', xstr, append=False)
                else:
                    xstr.write_xyz(fn=lproc + '/coord.xyz', append=False)
                xrargs = xrec.args["others"] + ";newstep"
            elif xrec.runtime["steps"][xstidx] == 0 and xtask[1] == "CONTINUE":
                xrargs = xrec.args["others"] + ";newstep"
            else:
                xrargs = xrec.args["others"]
            
            xxname = xrec.args["task_type"] + ".start"
            if xrec.fd is not None:
                if xfdstep == 0:
                    xrargs += ";forces=-1:-1"
                else:
                    xta, xtb = (xfdstep - 1) / 6, (xfdstep - 1) % 6
                    xrargs += ";forces=%d:%d;forceslen=%.5f" % (xrec.runtime["fd"][xxtid][1][xta], xtb,
                        xrec.fd.StepLength)
                xxname += "-%d" % xfdstep
            
            if xrec.args.get("sweep", None) is not None:
                xsn = xrec.args["sweep"]["n"]
                xsstart = xrec.args["sweep"]["start"]
                xsend = xrec.args["sweep"]["end"]
                xsname = xrec.args["sweep"]["name"]
                xsthis = 1.0 * (xsend - xsstart) / (xsn - 1) * (xtid % xsn) + xsstart
                xrargs = xrargs + (";%s=%g" % (xsname, xsthis))
            
            xtask += [
                xrec.args["step"] if xrec.args["task_type"] in ["relax", "neb"] else 0, 
                "%s/%s;%s" % (xrec.args["functional"], xrec.args.get("basis", ""), xrargs), 
            ]
            if isinstance(xrec.args["multiplicity"], int):
                xtask.append(xrec.args["multiplicity"])
            elif xrec.args["multiplicity"].startswith("read"):
                assert hasattr(xstr, 'mag') and xstr.mag is not None
                xamuls = "".join(xrec.args["multiplicity"].split("-")[1:])
                if xamuls.startswith("odd"):
                    xmulti = int(np.round(xstr.mag / 2)) * 2 + 1
                    xamulsx = xamuls[3:]
                elif xamuls.startswith("even"):
                    xmulti = int(np.round((xstr.mag - 1) / 2)) * 2 + 2
                    xamulsx = xamuls[4:]
                else:
                    raise RuntimeError("Must indicate odd/even multiplicity!")
                if xamulsx != "": xmulti += int(xamulsx)
                xtask.append(xmulti)
            xtask += [ xrec.args["charge"], xrec.args["program"].split(" ")[0] ]
            xtask = " ".join(map(str, xtask))

            self.logger.add(Event(proc_id=proc_id, name=xxname, str_id=xtid, idx=xrec.idx))
            
            f = open(lproc + '/REQUEST-INDEX', 'w')
            f.write("%d %d %d" % ( xrec.idx, xtid, xrec.runtime["steps"][xstidx] ))
            if xrec.fd is not None:
                f.write(" %d" % xfdstep)
            f.close()
        
        if os.path.isfile(lproc + '/REQUEST.nid'):
            try:
                os.rename(lproc + '/REQUEST.nid', new_file_name(lproc + '/REQUEST.nid.old'))
            except IOError:
                pass
        f = open(lproc + '/REQUEST', 'w')
        f.write(xtask)
        f.close()

    def move_restart(self, lproc, str_id, xrec, rfdstep):
        if os.path.isfile(lproc + '/restart.tar.gz'):
            xmcstep = 0 if xrec.mc is None else xrec.runtime["mc"][str(str_id)][2]
            if xrec.fd is not None:
                xmcstep = rfdstep
            xstidx = "%d.%d" % (str_id, xmcstep)
            new_tar = self.restart_dir + '/restart.tar.gz.' + str(xrec.idx) + '.' + xstidx
            os.rename(lproc + '/restart.tar.gz', new_tar)
            xrec.runtime["restart"][xstidx] = new_tar
            xrec.runtime_changed = True

    def run(self):
        self.logger.add_header(self.ip["random_seed"])
        stime = time.time()
        xblock = False
        while True:
            if self.has_finished():
                print ('all finished.')
                for rec in self.records: 
                    rec.save(final=True)
                break
            if self.max_master_time is not None and time.time() - stime > self.max_master_time:
                print ('max time reached.')
                break
            runnings = {}
            pprocs = {}
            res_found = False
            for proc in os.listdir(self.proc_dir):
                to_sleep = False
                lproc = os.path.join(self.proc_dir, proc)
                proc_id = int(re.findall(r'[0-9]+', proc)[0])
                sproc_id = str(proc_id)
                pprocs[sproc_id] = proc
                runnings[sproc_id] = False

                if os.path.isfile(lproc + '/BLOCK_RUNNING'):
                    xblock = True
                    runnings[sproc_id] = True
                elif os.path.isfile(lproc + '/RUNNING'):
                    last_proc_time = 0.0
                    for _ in range(4):
                        try:
                            last_proc_time = float(open(lproc + '/RUNNING').read().strip())
                        except (ValueError, IOError):
                            time.sleep(1.0)
                            continue
                        break
                    if time.time() - last_proc_time < self.max_proc_time:
                        runnings[sproc_id] = True
                if os.path.isfile(lproc + '/SLEEPING'):
                    runnings[sproc_id] = False
                if os.path.isfile(lproc + '/RESPONSE'):
                    try:
                        response = open(lproc + '/RESPONSE').read().strip().split(' ')
                        if response[0] != 'SUCCESS' and not runnings[sproc_id]:
                            self.close_proc("nr.resp.fail", proc_id, lproc=lproc)
                            continue
                        assert response[0] in [ 'SUCCESS', 'FAILED' ]
                    except (ValueError, IOError, AssertionError):
                        continue
                    try:
                        response_time = float(open(lproc + '/RESPONSE-TIME').read().strip())
                        rec_idx = map(int, open(lproc + '/REQUEST-INDEX-AC').read().strip().split(' '))
                        if len(rec_idx) == 3:
                            rec_idx.append(-1)
                    except (ValueError, IOError):
                        self.close_proc("time.r.fail", proc_id, lproc=lproc)
                        continue

                    str_id = None
                    for rec in self.records:
                        if sproc_id in rec.runtime["corr_list"]:
                            str_id = rec.runtime["corr_list"][sproc_id]
                            sub_step = -1
                            if isinstance(str_id, tuple) or isinstance(str_id, list):
                                sub_step = str_id[1]
                                str_id = str_id[0]
                            xrec = rec
                            xmcstep = 0 if xrec.mc is None else xrec.runtime["mc"][str(str_id)][2]
                            xstidx = "%d.%d" % (str_id, xmcstep)
                            str_step = rec.runtime["steps"][xstidx]
                            break
                    if str_id is None or [ xrec.idx, str_id, str_step, sub_step ] != rec_idx:
                        self.close_proc("no.consis.fail", proc_id, lproc=lproc, resp_time=response_time)
                        continue
                    rfdstep = rec_idx[3] if len(rec_idx) >= 4 else -1
                    
                    xener = None
                    if response[0] == 'SUCCESS':
                        # if move restart later, inherit will not get the right restart file
                        # this will influence freq calc
                        self.move_restart(lproc, str_id, xrec, rfdstep)
                        if response[1] == 'FREQ':
                            xname = xrec.handle_freqs(str_id, lproc)
                        elif response[1] == 'STEP_NEB':
                            xname, xener = xrec.handle_neb(str_id, lproc, "step")
                        elif response[1] == 'CONVERGED_NEB':
                            xname, xener = xrec.handle_neb(str_id, lproc, "conv")
                        elif response[1] == 'STEP':
                            xname, xener = xrec.handle_relax(str_id, lproc, "step", rfdstep)
                        elif response[1] == 'CONVERGED':
                            xname, xener = xrec.handle_relax(str_id, lproc, "conv", rfdstep)
                        elif response[1] == 'ENERGY':
                            xname, xener = xrec.handle_relax(str_id, lproc, "energy", rfdstep)
                        else:
                            raise RuntimeError("Unknown response type '%s'" % response[1])
                    else:
                        if response[1] == 'REPEAT_FREQ':
                            xname = 'freq.repeat'
                        elif response[1] == 'BAD_FREQ':
                            xname = xrec.handle_bad_freqs(str_id)
                        elif response[1] =='FATAL_FREQ':
                            xname = 'freq.fatal'
                            to_sleep = True
                        elif response[1] == 'UNKNOWN_FREQ':
                            xname = 'freq.unk'
                            to_sleep = True
                        elif response[1] == 'REPEAT':
                            xname = 'r.repeat'
                        elif response[1] == 'NOT_CONVERGED':
                            if self.ippa.get("continue_scf", False):
                                self.move_restart(lproc, str_id, xrec, rfdstep)
                            xname = xrec.handle_bad_relax(str_id, 'not_conv')
                        elif response[1] == 'BAD_STRUCTURE':
                            xname = xrec.handle_bad_relax(str_id, 'bad_struct')
                        elif response[1] =='FATAL':
                            xname = 'r.fatal'
                            to_sleep = True
                        elif response[1] == 'UNKNOWN':
                            xname = 'r.unk'
                            to_sleep = True
                        elif response[1] == 'REPEAT_NEB':
                            xname = 'neb.repeat'
                        elif response[1] == 'NOT_CONVERGED_NEB':
                            xname = xrec.handle_bad_neb(str_id, 'not_conv')
                        elif response[1] == 'BAD_STRUCTURE_NEB':
                            xname = xrec.handle_bad_neb(str_id, 'bad_struct')
                        elif response[1] =='FATAL_NEB':
                            xname = 'neb.fatal'
                            to_sleep = True
                        elif response[1] == 'UNKNOWN_NEB':
                            xname = 'neb.unk'
                            to_sleep = True
                    res_found = True
                    self.close_proc(xname, proc_id, lproc=lproc, str_id=str_id, 
                        resp_time=response_time, energy=xener, idx=xrec.idx)
                    
                if to_sleep:
                    if self.ippa.get("write_sleep", False):
                        self.assign_proc('sleep', lproc, proc_id=proc_id)
                    runnings[sproc_id] = False
                if not runnings[sproc_id] and not to_sleep and os.path.isfile(lproc + '/REQUEST'):
                    try:
                        os.rename(lproc + '/REQUEST', new_file_name(lproc + '/REQUEST.old'))
                    except IOError:
                        continue
                    self.close_proc("nr.req.del", proc_id, lproc=lproc)
                    continue
            
            for rec in self.records:
                for sproc_id in rec.runtime["corr_list"].keys():
                    if sproc_id not in pprocs:
                        continue
                    proc = pprocs[sproc_id]
                    lproc = os.path.join(self.proc_dir, proc)
                    proc_id = int(sproc_id)
                    if sproc_id not in runnings.keys() or not runnings[sproc_id]:
                        self.close_proc("nr.in.list", proc_id, lproc=lproc)
                    if sproc_id in runnings.keys() and runnings[sproc_id] and \
                        (os.path.isfile(lproc + '/REQUEST.nid') or \
                        (not os.path.isfile(lproc + '/REQUEST') and \
                        not os.path.isfile(lproc + '/RESPONSE'))):
                        try:
                            os.rename(lproc + '/REQUEST.nid', new_file_name(lproc + '/REQUEST.nid.old'))
                        except IOError:
                            pass
                        except OSError:
                            pass
                        self.close_proc("nr.in.list.nid", proc_id, lproc=lproc)
            
            waiting = False
            proc_act = 0
            for sproc_id, isrun in random.sample(runnings.items(), len(runnings)):
                if not isrun:
                    continue
                proc = pprocs[sproc_id]
                lproc = os.path.join(self.proc_dir, proc)
                proc_id = int(sproc_id)

                # check if proc has been used
                found = False
                for rec in self.records:
                    if sproc_id in rec.runtime["corr_list"]:
                        found = True
                        break
                if found:
                    proc_act += 1
                    continue

                if waiting: continue
                # find new structs
                found = None
                for rec in self.records:
                    rsour = rec.sources
                    if rec.mc is not None:
                        rsour.sort(key=lambda isrc: rec.runtime["mc"][str(isrc.tid)][2])
                    if rec.fd is not None:
                        # rsour.sort(key=lambda isrc: -len(rec.runtime["fd"][str(isrc.tid)][0]))
                        rsour.sort(key=lambda x: x.tid)
                        for isrc in rsour:
                            tid = isrc.tid
                            sta = rec.runtime["states"][str(tid)]
                            if sta == 0:
                                if 0 in rec.runtime["fd"][str(isrc.tid)][0]:
                                    if (tid, 0) not in rec.runtime["corr_list"].values():
                                        found = [ rec, tid, 0]
                                        break
                                else:
                                    for isp in rec.runtime["fd"][str(isrc.tid)][0]:
                                        if (tid, isp) not in rec.runtime["corr_list"].values():
                                            found = [ rec, tid, isp]
                                            break
                            if found is not None:
                                break
                    else:
                        for isrc in rsour:
                            tid = isrc.tid
                            sta = rec.runtime["states"][str(tid)]
                            if sta == 0 and tid not in rec.runtime["corr_list"].values():
                                found = [ rec, tid, -1 ]
                                break
                    if found is not None: break
                if found is None:
                    waiting = True
                else:
                    proc_act += 1
                    self.assign_proc(None, lproc, xrec=found[0], xtid=found[1], proc_id=proc_id, rstepid=found[2])
            
            xstates = None
            for rec in self.records:
                stv = rec.runtime["states"].values()
                iempt, iconv, ifail = stv.count(0), stv.count(1), stv.count(2) + stv.count(3)
                imax, idupl, itot = stv.count(4), stv.count(5) + stv.count(6), len(stv)
                iidx = rec.idx
                if iempt != 0: break
            
            for rec in self.records: rec.save()

            ipact, iptot = proc_act, runnings.values().count(True)
            iitot = len(self.records)
            self.logger.add_summary(iempt, itot, iidx, iitot, ipact, iptot, wait=waiting)
            self.logger.save()
            self.shlogger.write(*[ itot, iempt, iconv, ifail, imax, idupl, iidx, iitot, 
                ipact, iptot, "finished" if ipact == 0 and iempt == 0 and waiting else (
                    "waiting" if waiting else ( "paused" if ipact == 0 and itot != iempt else (
                    "init" if ipact == 0 else "running"))) ])
            
            if xblock and runnings.values().count(True) == 0:
                print ('all slept.')
                for rec in self.records: 
                    rec.save(final=True)
                break
            if not res_found:
                time.sleep(self.idle_time)
