
from __future__ import print_function
import numpy as np
from reportlab.pdfbase.pdfmetrics import stringWidth

class Ptn(object):
    x = 0.0
    y = 0.0
    idx = 0
    charge = None
    def __init__(self, x, y, charge, idx):
        self.x, self.y, self.charge, self.idx = x, y, charge, idx

    def __repr__(self):
        return "[%.2f, %.2f]" % (self.x, self.y)

class ChargeItem(object):
    def __init__(self, charge, cpos, tpos, wid):
        self.charge, self.cpos, self.tpos, self.wid = charge, cpos, tpos, wid

def print_charge(c):
    return (u"%+.2f" % c).replace("-", u"\u2013")

# dcross = lambda a, b, c: np.cross(b - a, c - a)
# din = lambda a, b, c: np.dot(c - a, c - b)
# edges = lambda x: [[x[i], x[(i + 1) % 4]] for i in range(4)]
# def line_cross(a, b, c, d):
#     r = j, k, l, m = [[a, b, c], [a, b, d], [c, d, a], [c, d, b]]
#     return dcross(*j) * dcross(*k) < 0 and dcross(*l) * dcross(*m) < 0 \
#         or True in [dcross(*i) == 0 and din(*i) <= 0 for i in r]

# def rect_cross(apos, bpos, awid, bwid, ahei, bhei):
#     a = [apos, apos + np.array([awid, 0.0]), apos + np.array([awid, ahei])]
#     a += [apos + np.array([0.0, ahei])]
#     b = [bpos, bpos + np.array([bwid, 0.0]), bpos + np.array([bwid, bhei])]
#     b += [bpos + np.array([0.0, bhei])]
#     return True in [line_cross(*(i + j)) for i in edges(a) for j in edges(b)]
def cross_one(a, b, c, d):
    r = [[a, b, c], [a, b, d], [c, d, a], [c, d, b]]
    return True in [(n >= l and n <= m) or (n <= l and n >= m) for l, m, n in r]

def rect_cross(apos, bpos, awid, bwid, ahei, bhei):
    return cross_one(apos[0], apos[0] + awid, bpos[0], bpos[0] + bwid) and \
        cross_one(apos[1], apos[1] + ahei, bpos[1], bpos[1] + bhei)

def ctheta(x, y):
    m = np.linalg.norm([x, y])
    return [-x / m, m]

def idx_adj(i, j, n):
    return (j - i) % n

class ChargePlotter(object):
    atom_radius = 1.3 * 1.5
    def __init__(self, coords, charges, fontsize, pte):
        pts = []
        for idx, (co, ch) in enumerate(zip(coords, charges)):
            pts.append(Ptn(co[0], co[1], ch, idx))
        self.pts = pts
        self.hull = []
        self.find_hull()
        self.cpts = []
        self.fontsize = fontsize
        for ip, p in enumerate(self.pts):
            if p.charge is None:
                self.cpts.append(None)
            else:
                cpos = self.charge_point(ip)
                wid = stringWidth(print_charge(p.charge) + "|%s" % pte[ip], fontName="Arial",
                                  fontSize=self.fontsize)
                if p.x > cpos[0]:
                    wid = -wid
                self.cpts.append([p.charge, np.array(cpos), np.array([p.x, p.y]), wid])
        hei = self.fontsize
        for _ in range(10):
            adjed = False
            for ia, iag in enumerate(self.cpts):
                if iag is None:
                    continue
                _, cpos, apos, wid = iag
                ahei = hei if cpos[1] > apos[1] else -hei
                akked = False
                for ib, ibg in enumerate(self.cpts[:ia]):
                    if ibg is None:
                        continue
                    _, pcpos, papos, pwid = ibg
                    bhei = hei if pcpos[1] > papos[1] else -hei
                    l = np.linalg.norm(cpos - pcpos)
                    if l < ChargePlotter.atom_radius and l != 0:
                        k = (ChargePlotter.atom_radius - l) / 2
                        a = (cpos - pcpos) / l * k
                        self.cpts[ib][1] = pcpos - a
                        self.cpts[ia][1] = cpos + a
                        adjed = True
            if not adjed:
                break
        while True:
            adjed = False
            for ia, iag in enumerate(self.cpts):
                if iag is None:
                    continue
                _, cpos, apos, wid = iag
                ahei = hei if cpos[1] > apos[1] else -hei
                akked = False
                for ibg in self.cpts[:ia]:
                    if ibg is None:
                        continue
                    _, pcpos, papos, pwid = ibg
                    bhei = hei if pcpos[1] > papos[1] else -hei
                    if rect_cross(cpos, pcpos, wid, pwid, ahei, bhei):
                        akked = True
                        break
                if akked:
                    rl = cpos - apos
                    grl = np.linalg.norm(rl)
                    rl += rl / grl * ChargePlotter.atom_radius
                    self.cpts[ia][1] = apos + rl
                    adjed = True
            if not adjed:
                break

    def find_hull(self):
        # delete same points
        kxys = self.pts[:]
        xys = []
        for kxy in kxys:
            fd = False
            for xy in xys:
                if xy.x == kxy.x and xy.y == kxy.y:
                    fd = True
                    break
            if not fd:
                xys.append(kxy)

        # base point
        idm = -1
        for ixy, xy in enumerate(xys):
            if idm is -1 or xy.y < xys[idm].y or (xy.y == xys[idm].y and xy.x < xys[idm].x):
                idm = ixy

        pts2 = xys[:]
        del pts2[idm]

        pts2.sort(key=lambda x: ctheta(x.x - xys[idm].x, x.y - xys[idm].y))
        res = [xys[idm]]
        for ik, _ in enumerate(pts2):
            for ij in range(len(res) - 1, 0, -1):
                v1 = [res[ij].x - res[ij - 1].x, res[ij].y - res[ij - 1].y]
                v2 = [pts2[ik].x - res[ij].x, pts2[ik].y - res[ij].y]
                if np.cross(v1, v2) > 0.0:
                    break
                else:
                    del res[ij]
            res.append(pts2[ik])
        for ir, r in enumerate(res):
            r.hidx = ir
        # print ("pts = ", self.pts)
        self.hull = res

    def charge_point(self, idx):
        pt = self.pts[idx]
        if len(self.hull) == 1:
            return (pt.x, pt.y + ChargePlotter.atom_radius)
        elif pt in self.hull or len(self.hull) < 3:
            ct = np.array([[h.x, h.y] for h in self.hull]).mean(axis=0)
            cpt = np.array([pt.x, pt.y])
            cmt = cpt - ct
            return tuple(cpt + cmt / np.linalg.norm(cmt) * ChargePlotter.atom_radius)
        ptr = self.hull[:]
        ptr.sort(key=lambda x: np.linalg.norm([pt.x - x.x, pt.y - x.y]))
        pa = ptr[0]
        for ipb in ptr[1:]:
            l = (ipb.hidx - pa.hidx) % len(ptr)
            if l == 1 or l == len(ptr) - 1:
                pb = ipb
                break
        ltb = np.array([pb.x - pt.x, pb.y - pt.y])
        lta = np.array([pa.x - pt.x, pa.y - pt.y])
        lab = np.array([pb.x - pa.x, pb.y - pa.y])
        gtb, gta, gab = [np.linalg.norm(x) for x in [ltb, lta, lab]]
        ctab = (gab**2 + gta**2 - gtb**2) / (2 * gta * gab)
        lam = lab / gab * gta * ctab
        rl = lta + lam
        grl = np.linalg.norm(rl)
        rl += rl / grl * ChargePlotter.atom_radius
        return tuple(rl + np.array([pt.x, pt.y]))
