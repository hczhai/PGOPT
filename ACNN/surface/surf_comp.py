
if __name__ == "__main__":
    import sys
    sys.path.insert(0, sys.path[0] + '/..')

import numpy as np
from surface.base import ClusterAtSurface, Surface
from surface.spacegroupdata import SymOpsHall
from surface.surf_symm import apply_trans, to_direct, to_cartesian
from surface.surf_symm import to_cellmat
from formod.at_comp import at_comp
from formod.co_comp import co_comp
from formod.km_comp import kmc_comp, kmc_comp_best, weights

eps = np.finfo(dtype=float).eps * 100
TOLS = 0.50 - 0.02

class DisjointSet(object):
    def __init__(self, n):
        self.n = n
        self.parent = np.zeros(self.n, dtype=int)
        self.rank = np.zeros(self.n, dtype=int)
        self.clear()

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def clear(self):
        for i in range(0, len(self.rank)):
            self.rank[i] = 0
            self.parent[i] = i

    def union(self, x, y):
        x = self.find(x)
        y = self.find(y)
        if x == y:
            return
        elif self.rank[x] < self.rank[y]:
            self.parent[x] = y
        elif self.rank[x] > self.rank[y]:
            self.parent[y] = x
        else:
            self.parent[y] = x
            self.rank[x] += 1

def num_imagex(x, l, tol):
    k = np.round(np.abs(x) / l) * l
    return np.abs(np.abs(x) - k) < tol

def num_adjust(x, l):
    return np.round(x / l) * l

def check_imagex(ata, ela, atb, elb, cell, tol, cl):
    if ela != elb:
        return False
    for dx in range(0, 3):
        if cell[dx] is None:
            if np.abs(ata[dx] - atb[dx]) > tol / cl[dx]:
                return False
        else:
            if not num_imagex(ata[dx] - atb[dx], cell[dx], tol / cl[dx]):
                return False
    return True

# move periodic cluster atoms (at surface) for better interpolation
# (atoms should be matched)
def cas_adjust(a, b):
    assert isinstance(a, ClusterAtSurface)
    assert isinstance(b, ClusterAtSurface)
    adir = to_direct(a.atoms, a.surf.cell)
    bdir = to_direct(b.atoms, b.surf.cell)
    for i in range(a.n):
        for k in range(2):
            adir[i, k] += num_adjust(bdir[i, k] - adir[i, k], 1.0)
    a.atoms = to_cartesian(adir, a.surf.cell)

# move surface atoms for better interpolation
# (atoms should be matched)
def surface_adjust(a, b):
    assert isinstance(a, Surface)
    assert isinstance(b, Surface)
    adir = to_direct(a.atoms, a.cell)
    bdir = to_direct(b.atoms, b.cell)
    for i in range(a.n):
        for k in range(2):
            adir[i, k] += num_adjust(bdir[i, k] - adir[i, k], 1.0)
    a.atoms = to_cartesian(adir, a.cell)

# match surface atoms as much as possible


# match surface atoms
def surface_match(a, b, dmax, dmin=1E-2, nd=20):
    from cluster.base import elem_num
    assert isinstance(a, Surface)
    assert isinstance(b, Surface)
    adir = to_direct(a.atoms, a.cell)
    bdir = to_direct(b.atoms, a.cell)
    if len(a.elems) == len(b.elems):
        # efficient way
        assert (a.elems == b.elems).all()
        cellmat = to_cellmat(a.cell)
        eles, ne = elem_num(a.elems)
        ww = weights(a.atoms, b.atoms, ne, eles, cellmat)
        v, vx = kmc_comp(ww, ne, eles)
        vg = np.zeros(vx.shape, dtype=int)
        for i, ii in enumerate(vx):
            vg[ii - 1] = i + 1
        return v, vg
    else:
        # old way
        cl = np.linalg.norm(to_cellmat(a.cell), axis=1)
        dstep = (dmax/dmin)**(1.0/nd)
        d = dmin
        mt = [-1] * a.n
        mx = range(b.n)
        cell = [1.0, 1.0, None]
        for _ in range(nd):
            for i in range(a.n):
                if mt[i] == -1:
                    for ij, j in enumerate(mx):
                        if check_imagex(adir[i], a.elems[i], bdir[j], b.elems[j], cell, d, cl):
                            mt[i] = j
                            del mx[ij]
                            break
            if len(mx) == 0:
                break
            d *= dstep
        # diff, march indices
        return d, np.array(mt) + 1

# ne and eles only for cluster
# if best is not none, but a integer k, produce best k marchings
def surface_compare(a, b, ne, eles, _, best=None):
    assert isinstance(a, ClusterAtSurface)
    assert isinstance(b, ClusterAtSurface)
    ufactor = 0.3
    adir = to_direct(a.atoms, a.surf.cell)
    bdir = to_direct(b.atoms, b.surf.cell)
    bcent = bdir.mean(axis=0)
    bxcar = b.atoms
    surf = a.surf
    cellmat = to_cellmat(surf.cell)
    ucellmat = to_cellmat(surf.unit_cell)
    ucell = np.diag(to_direct(ucellmat, surf.cell))
    # space group transition reference atoms must be taken from space_group_ref
    atrref = surf.space_group_ref
    xvmins = []
    sg = SymOpsHall[surf.space_group]
    sg = [[x.strip() for x in s] for s in sg if s[2].strip() == 'z']
    for s in sg:
        ax = np.array(adir)
        ax[:, 0:2] -= atrref
        ax = apply_trans(ax, s, ucell)
        ax[:, 0:2] += atrref
        # for unit-cell periodicity, cluster must be moved as a whole
        # for supercell periodicity, atoms can be moved singlely and dim by dim
        acentpre = ax.mean(axis=0)
        dss = [[0.0], [0.0]]
        for k in range(2):
            acentpre[k] += num_adjust(bcent[k] - acentpre[k], ucell[k])
            if bcent[k] > acentpre[k] + ufactor * ucell[k]:
                dss[k].append(1.0)
            elif bcent[k] < acentpre[k] - ufactor * ucell[k]:
                dss[k].append(-1.0)
        for dx in dss[0]:
            for dy in dss[1]:
                dxy = np.array([dx * ucell[0], dy * ucell[1], 0.0])
                acent = acentpre + dxy
                axcar = to_cartesian(ax - ax.mean(axis=0) + acent, surf.cell)
                ww = weights(axcar, bxcar, ne, eles, cellmat)
                if best is None:
                    v, vx = kmc_comp(ww, ne, eles)
                    opmin = [s, vx, -ax.mean(axis=0) + acent, axcar]
                    xvmins.append([v, opmin])
                else:
                    v, vx = kmc_comp_best(best, ww, ne, eles)
                    for ij in range(best):
                        opmin = [s, vx[:, ij], -ax.mean(axis=0) + acent, axcar]
                        xvmins.append([v[ij], opmin])
    xvmins.sort(key=lambda x: x[0])
    if best is not None:
        xvmins = xvmins[:best]
    else:
        xvmins = xvmins[:1]
    for v, opmin in xvmins:
        if opmin is not None:
            asdir = to_direct(a.surf.atoms, surf.cell)
            asdir[:, 0:2] -= atrref
            asdir = apply_trans(asdir, opmin[0], ucell)
            asdir[:, 0:2] += atrref
            asdir += opmin[2]
            for ii, _ in enumerate(asdir):
                for k in range(2):
                    asdir[ii, k] += num_adjust(TOLS - asdir[ii, k], 1.0)
            ascar = to_cartesian(asdir, surf.cell)
            opmin.append(ascar)
    if best is not None:
        return xvmins
    else:
        return xvmins[0]

def surface_align(cc, deep=True):
    assert isinstance(cc, ClusterAtSurface)
    if cc.n == 0:
        return
    # part 1: put cluster fragments together (using supercell periodicity)
    # atoms having a distance less than 1/3 cell are considered connected
    # method: using disjointset finding all fragments
    # move all other atoms closer to the center of the biggest fragment
    if cc.n > 1:
        cdir = to_direct(cc.atoms, cc.surf.cell)
        disj = DisjointSet(len(cdir))
        for ii, i in enumerate(cdir):
            for jj, j in enumerate(cdir[:ii]):
                if np.linalg.norm(i - j) < 1.0 / 3.0:
                    disj.union(ii, jj)
        roots, rootn = [], []
        for ii, i in enumerate(cdir):
            ik = disj.find(ii)
            if ik not in roots:
                roots.append(ik)
                rootn.append(1)
            else:
                rootn[roots.index(disj.find(ik))] += 1
        rootm = max(zip(roots, rootn), key=lambda x: x[1])
        rootc = np.zeros(3)
        for ii, i in enumerate(cdir):
            if disj.find(ii) == rootm[0]:
                rootc += i
        rootc /= rootm[1]
        for ii, i in enumerate(cdir):
            if disj.find(ii) == rootm[0]:
                continue
            for k in range(2):
                cdir[ii, k] += num_adjust(rootc[k] - cdir[ii, k], 1.0)
        cc.atoms = to_cartesian(cdir, cc.surf.cell)
    # part 2: put cluster closer to center of surface (using unitcell periodicity)
    # can make the surface atoms indeices disordered
    surf = cc.surf
    if deep:
        ucellmat = to_cellmat(surf.unit_cell)
        ucell = np.diag(to_direct(ucellmat, surf.cell))
        cdir = to_direct(cc.atoms, surf.cell)
        # space group transition reference atoms must be taken from space_group_ref
        ctrref = surf.space_group_ref
        scent = np.array([0.50, 0.50])
        cxdmin = []
        sg = SymOpsHall[surf.space_group]
        sg = [[x.strip() for x in s] for s in sg if s[2].strip() == 'z']
        for s in sg:
            cx = np.array(cdir)
            cx[:, 0:2] -= ctrref
            cx = apply_trans(cx, s, ucell)
            cx[:, 0:2] += ctrref
            # for unit-cell periodicity, cluster must be moved as a whole
            # for supercell periodicity, atoms can be moved singlely and dim by dim
            ccent = cx.mean(axis=0)
            for k in range(2):
                ccent[k] += num_adjust(scent[k] - ccent[k], ucell[k])
            cxd = np.linalg.norm(scent[0:2] - ccent[0:2])
            cxcar = to_cartesian(cx - cx.mean(axis=0) + ccent, surf.cell)
            if len(cxdmin) == 0 or cxd < cxdmin[0]:
                cxdmin = [cxd, cxcar, s, -cx.mean(axis=0) + ccent]
        cc.atoms = cxdmin[1]
    # part 3: move surface atoms accodingly and also make surface boundary atoms better
    csdir = to_direct(surf.atoms, surf.cell)
    if deep:
        csdir[:, 0:2] -= ctrref
        csdir = apply_trans(csdir, cxdmin[2], ucell)
        csdir[:, 0:2] += ctrref
        csdir += cxdmin[3]
    for ii, _ in enumerate(csdir):
        for k in range(2):
            csdir[ii, k] += num_adjust(TOLS - csdir[ii, k], 1.0)
    surf.atoms = to_cartesian(csdir, surf.cell)

if __name__ == "__main__":
    from surface.base import read_surface
    surf = read_surface('../PGOPT/surfaces/mgo.xyz')
    surf.space_group = "-F 4 2 3"
    surf.unit_cell = np.array([4.257, 4.257, 4.257])
    from cluster.base import read_clusters, elem_num
    clus = read_clusters("/Users/fszix/Projects/test/master/1.1/par_local.xyz.0")
    cla = []
    nx = 72
    for c in clus:
        aa = ClusterAtSurface(c.n - nx, surf)
        surfmin = c.atoms[:nx].min(axis=0).reshape(1, 3)
        aa.atoms = c.atoms[nx:] - surfmin
        aa.elems = c.elems[nx:]
        cla.append(aa)
    eles, ne = elem_num(aa.elems)
    print surface_compare(cla[8], cla[9], ne, eles, 1.0)
    