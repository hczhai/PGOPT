
from __future__ import print_function
import numpy as np
from reportlab.platypus.flowables import Flowable
from reportlab.pdfbase.pdfmetrics import stringWidth
from formod.lbfgs import LBFGS

def to_hex(ints):
    return "#%02X%02X%02X" % tuple(ints)

def to_ints(hexx):
    return np.array([float(int(hexx[i:i+2], 16)) for i in range(1, 6, 2)])

def to_letters(intx):
    return chr(ord("A") + intx % 26) + "'" * (intx / 26)

def influ(ang):
    kk = np.abs(ang % (np.pi * 2))
    kk = min(kk, np.pi * 2 - kk)
    return kk

def draw_bh(canv, aa, ee, pos, cll):
    rot = aa / np.pi * 180.0
    cc = 1.5
    if rot > 90.0:
        rot -= 180.0
    elif rot < -90.0:
        rot += 180.0
    else:
        cc = 0.0
    canv.translate(pos[0], pos[1])
    canv.rotate(rot)
    x = "%.3f" % ee
    canv.setFillColor(cll)
    ww = stringWidth(x, fontName="Arial", fontSize=8)
    canv.setFont("Arial", 8.0)
    canv.drawString(-ww / 2, cc - 4.0, x)
    canv.rotate(-rot)
    canv.translate(-pos[0], -pos[1])

class DrawRGraph(Flowable, object):
    def __init__(self, rgraph, width=500):
        self.canv = None
        super(DrawRGraph, self).__init__()
        self.width = width
        self.rgraph = rgraph
        self.pos = rgraph.plot() * self.width
        print ("# of pathways = %d" % len(self.rgraph.es))
        self.n = len(self.pos) / 2
        self.height = self.pos[self.n:].max()
        from meta.report import Report
        from meta.projection import DrawCluster
        _, self.rto = Report.test_width(self.rgraph.mins.values()[0:5], clip_circle=True)
        self.rto *= 1.0
        self.dclus = {it: DrawCluster(t, simple=True, clip_circle=True, ratio=self.rto)
                 for it, t in self.rgraph.mins.items()}
        if len(self.dclus) != 0:
            hd = 12.0 + max([d.width / 2 for d in self.dclus.values()])
        else:
            hd = 12.0
        self.height += hd * 2
        self.width += hd * 2
        self.hd = hd

    def wrap(self, *opts):
        return (self.width, self.height)

    def draw(self):
        from meta.energetics import draw_flowable, httoev
        dclus = self.dclus
        self.canv.translate(self.hd, self.hd)
        ex = np.array([[e.wab, e.wba] for e in self.rgraph.es])
        if len(ex) == 0:
            return
        em, en = ex.min(), ex.max()
        egu = {k: np.pi / 3 for k in self.rgraph.mins.keys()}
        egd = {k: np.pi / 3 for k in self.rgraph.mins.keys()}
        srr = 20.0
        srl = []
        bhl = []
        for e in self.rgraph.es:
            ia = self.rgraph.vs.index(e.a)
            ib = self.rgraph.vs.index(e.b)
            wa = dclus[str(e.a)].width / 2
            wb = dclus[str(e.b)].width / 2
            atb = np.array([self.pos[ib] - self.pos[ia], self.pos[self.n+ib] - self.pos[self.n+ia]])
            atbn = atb / np.linalg.norm(atb)
            self.canv.setLineWidth(1.5)
            atbv = np.array([atbn[1], -atbn[0]])
            ra = np.array([self.pos[ia], self.pos[self.n+ia]]) + atbn * wa
            rb = np.array([self.pos[ib], self.pos[self.n+ib]]) - atbn * wb
            ta = ra + atbv * 1.0
            tb = rb + atbv * 1.0
            tc = tb - atbn * 8.0 + atbv * 1.8
            sa = ra - atbv * 1.0
            sb = rb - atbv * 1.0
            sc = sa + atbn * 8.0 - atbv * 1.8
            rd = np.array([self.pos[ia], self.pos[self.n+ia]]) + atb / 2
            td = rd + atbv * 5.0
            sd = rd - atbv * 5.0
            aa = np.arctan2(atbn[1], atbn[0])
            egu[str(e.a)] = min(egu[str(e.a)], influ(np.pi / 2 - aa))
            egd[str(e.a)] = min(egd[str(e.a)], influ(3 * np.pi / 2 - aa))
            ab = np.arctan2(-atbn[1], -atbn[0])
            egu[str(e.b)] = min(egu[str(e.b)], influ(np.pi / 2 - ab))
            egd[str(e.b)] = min(egd[str(e.b)], influ(3 * np.pi / 2 - ab))
            dobh = np.linalg.norm(ra - rb) > srr * 2
            for sx in srl:
                if np.linalg.norm(sx - rd) < srr:
                    dobh = False
                    break
            for ii, (i, j) in enumerate(zip(self.pos[:self.n], self.pos[self.n:])):
                ff = dclus[str(self.rgraph.vs[ii])]
                if np.linalg.norm(np.array([i, j]) - rd) < srr + ff.width / 2:
                    dobh = False
                    break
            if dobh:
                srl.append(rd)
            rto = (e.wab - em) / (en - em) if en != em else 0.5
            cll = to_hex((1.0 - rto * 0.8) * to_ints("#165766") + rto * 0.8 * to_ints("#FFFFFF"))
            self.canv.setStrokeColor(cll)
            self.canv.line(ta[0], ta[1], tb[0], tb[1])
            self.canv.line(tc[0], tc[1], tb[0], tb[1])
            if dobh:
                bhl.append([aa, e.wab, td, cll])
            rto = (e.wba - em) / (en - em) if en != em else 0.5
            cll = to_hex((1.0 - rto * 0.8) * to_ints("#165766") + rto * 0.8 * to_ints("#FFFFFF"))
            self.canv.setStrokeColor(cll)
            self.canv.line(sa[0], sa[1], sb[0], sb[1])
            self.canv.line(sa[0], sa[1], sc[0], sc[1])
            if dobh:
                bhl.append([ab, e.wba, sd, cll])
        for a, w, d, c in bhl:
            draw_bh(self.canv, a, w, d, c)
        self.canv.setStrokeColor("#165766")
        self.canv.setLineWidth(1.0)
        self.canv.setFillColor("#165766")
        erx = np.array([e.energy for e in self.rgraph.mins.values()])
        erm, ern = erx.min(), erx.max()
        for ii, (i, j) in enumerate(zip(self.pos[:self.n], self.pos[self.n:])):
            # self.canv.circle(i, j, 0.2, fill=0)
            ff = dclus[str(self.rgraph.vs[ii])]
            x = to_letters(self.rgraph.vs[ii])
            ww = stringWidth(x, fontName="Arial", fontSize=12)
            mxt, mxl = None, 0.0
            for _ in range(50):
                mxr = (np.random.random() * 0.5 + 0.5) * (ff.width - 10) / 2
                mxth = np.random.random() * 2 * np.pi
                mxx, mxy = mxr * np.cos(mxth), mxr * np.sin(mxth)
                mxlc = np.linalg.norm(ff.shs[:, :2] - np.array([mxx + ff.width / 2, \
                    mxy + ff.height / 2]), axis=1).min()
                if mxlc > mxl:
                    mxt, mxl = [mxx, mxy], mxlc
            draw_flowable(ff, self.canv, i - ff.width / 2, j + ff.height / 2)
            te = self.rgraph.mins[str(self.rgraph.vs[ii])].energy
            rto = (te - erm) / (ern - erm) if ern != erm else 0.5
            cll = to_hex((1.0 - rto * 0.7) * to_ints("#165766") + rto * 0.7 * to_ints("#FFFFFF"))
            self.canv.setFillColor(cll)
            self.canv.setStrokeColor(cll)
            self.canv.circle(i, j, ff.width / 2, fill=0)
            self.canv.setFont("Arial", 12.0)
            self.canv.drawString(i + mxt[0] - ww / 2, j + mxt[1] - 5, x)
            if egu[str(self.rgraph.vs[ii])] < egd[str(self.rgraph.vs[ii])]:
                epos = -ff.height / 2 - 5
            else:
                epos = ff.height / 2 + 5
            x = "[%.3f]" % ((te - erm) * httoev)
            ww = stringWidth(x, fontName="Arial", fontSize=8)
            self.canv.setFont("Arial", 8.0)
            self.canv.drawString(i - ww / 2, j + epos - 3, x)
        self.canv.translate(-self.hd, -self.hd)

class REdge(object):
    a = None
    b = None
    wab = None
    wba = None
    def __init__(self, a, b, wab, wba):
        self.a, self.b, self.wab, self.wba = a, b, wab, wba

class RGraph(object):
    def __init__(self):
        self.es = []
        self.vs = []
        self.matrix = None
        self.mins = None

    def add_edge(self, e):
        self.es.append(e)
        for x in [e.a, e.b]:
            if x not in self.vs:
                self.vs.append(x)

    def clear(self):
        self.es = []
        self.vs = []

    def sort(self):
        self.vs.sort()

    def __repr__(self):
        ss = ""
        for e in self.es:
            ss += "[%d.%d] %.3f / %.3f\n" % (e.a, e.b, e.wab, e.wba)
        return ss

    @staticmethod
    def ener(x, spath):
        n = len(x) / 2
        ee = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                ee += (np.sqrt((x[i] - x[j])**2 + (x[n+i] - x[n+j])**2) / spath[i, j] - 1)**2
        return ee

    @staticmethod
    def dener(x, spath):
        n = len(x) / 2
        dee = np.zeros(n * 2)
        for i in range(n):
            for j in range(i + 1, n):
                r = np.sqrt((x[i] - x[j])**2 + (x[n+i] - x[n+j])**2)
                dep = 2 * (r / spath[i, j] - 1) / spath[i, j]
                dee[i] += dep * (x[i] - x[j]) / r
                dee[j] += -dep * (x[i] - x[j]) / r
                dee[n+i] += dep * (x[n+i] - x[n+j]) / r
                dee[n+j] += -dep * (x[n+i] - x[n+j]) / r
        return dee

    @staticmethod
    def adjust(x, y):
        ine = np.zeros((2, 2))
        x = x - x.mean()
        y = y - y.mean()
        ine[1, 0] = ine[0, 1] = np.dot(x, y)
        ine[0, 0] = np.dot(x, x)
        ine[1, 1] = np.dot(y, y)
        _, v = np.linalg.eigh(ine)
        v = v[::-1]
        xc = np.tensordot(np.array(zip(x, y)), v, axes=[1, 1])
        if np.dot(xc[:, 0], abs(xc[:, 0])).sum() < 0:
            v = -v
        ca = np.tensordot(np.array(zip(x, y)), v, axes=[1, 1])
        return ca[:, 0], ca[:, 1]

    def plot(self):
        # shortest-path Floyd-Warshall
        nn = len(self.vs)
        mm = len(self.es)
        spath = np.zeros((nn, nn), dtype=float)
        evs = np.zeros(mm, dtype=float)
        inf, isrt, ilng = 5.0, 1.0, 3.0
        if len(self.es) == 0:
            return np.array([1.0, 1.0])
        evog = np.array([min(x.wab, x.wba) for x in self.es])
        evx, evm = evog.max(), evog.min()
        for k in range(0, mm):
            evs[k] = (evog[k] - evm) / (evx - evm) * (ilng - isrt) + isrt
        spath[:, :] = inf
        for i in range(0, nn):
            spath[i, i] = 0
        for k in range(0, mm):
            ia = self.vs.index(self.es[k].a)
            ib = self.vs.index(self.es[k].b)
            spath[ia, ib] = evs[k]
            spath[ib, ia] = evs[k]
        for k in range(0, nn):
            for i in range(0, nn):
                for j in range(0, nn):
                    d = spath[i, k] + spath[k, j]
                    if d < spath[i, j]:
                        spath[i, j] = d
        # init
        xx = np.zeros(nn * 2)
        ang = 2.0 * np.pi / nn / 2
        ar = 0.5 / np.sin(ang)
        ax = 0.0
        for i in range(0, nn):
            xx[i] = ar * np.cos(ax)
            xx[nn+i] = ar * np.sin(ax)
            ax += ang
        # minimization
        task = LBFGS(nn * 2)
        task.log_file = 0
        task.max_iter = 500
        task.p.eval = lambda x, spathx=spath: RGraph.ener(x, spathx)
        task.p.evald = lambda x, spathx=spath: RGraph.dener(x, spathx)
        task.start(xx)
        task.opt()
        x, y = RGraph.adjust(task.x[:nn], task.x[nn:])
        x -= min(x)
        y -= min(y)
        xy = np.array([x, y]).flatten()
        xy /= max(x)
        return xy

class RPath(object):
    def __init__(self, tid, ismax=False, eeps=1E-3):
        self.tid = tid
        self.structs = []
        self.energies = []
        self.props = {}
        self.nimages = 0
        self.surf = None
        self.from_to = None
        self.diffs = []
        self.ismax = ismax
        self.isdirect = True
        self.barriers = np.zeros(2)
        self.focus = False
        self.eeps = eeps

    def add(self, x):
        self.structs.append(x)

    def sort(self):
        self.structs.sort(key=lambda x: x.tidx[1])
        self.energies = [x.energy for x in self.structs]
        if hasattr(self.structs[0], "props"):
            self.props["time"] = self.structs[0].props.get("time", 0)
            self.props["step"] = self.structs[0].props.get("step", 0)
        self.nimages = len(self.structs) - 2
        if hasattr(self.structs[0], "surf"):
            self.surf = self.structs[0].surf
        ll = self.structs[0].label.split(":")[0].split("-")
        self.from_to = [int(x) for x in ll[1:]]
        self.diffs = []
        for i in range(0, len(self.structs)):
            v = np.linalg.norm(self.structs[i].atoms - self.structs[i - 1].atoms, axis=1).mean()
            self.diffs.append(v)
        for i in range(1, len(self.structs) - 1):
            if self.energies[i] < self.energies[i-1] - self.eeps and \
                self.energies[i] < self.energies[i+1] - self.eeps:
                self.isdirect = False
                break
        mh = max(self.energies)
        self.barriers = np.array([mh - self.energies[0], mh - self.energies[-1]])
