
from __future__ import print_function
import copy
import numpy as np

class Surface(object):
    def __init__(self, n):
        self.n = n
        self.atoms = np.zeros((n, 3), dtype=float)
        self.elems = np.array([''] * n, dtype='|S20')
        self.cell = np.zeros(3, dtype=float)
        self.cellz = 0.0
        self.unit_cell = np.zeros(3)
        self.space_group = "P 1"
        self.space_group_ref = None # np.zeros(2) # in direct coordinates
        self.forces = None

# read surface file
def read_surface(fn, zf=None, iprint=True):
    if zf is None:
        fs = open(fn, 'r').readlines()
    else:
        fs = zf.open(fn, 'rU').readlines()
    fs = [[g for g in f.replace('\n', '').split(' ') if len(g) != 0] for f in fs]
    cc = 0
    if not fs[cc][0].isdigit():
        return None
    cn = int(fs[cc][0])
    clu = Surface(cn)
    i = 0
    if fs[cc + 1][0] != "CELL":
        return None
    if "(" in fs[cc + 1][-1]:
        clu.cellz = float(fs[cc + 1][-1].replace("(", "").replace(")", ""))
        fs[cc + 1] = fs[cc + 1][:-1]
    # can be of length 3 (cubic) or 5 (ax, ay, bx, by, cz)
    clu.cell = np.array([float(x) for x in fs[cc + 1][2:]])
    for f in fs[cc + 2:cc + 2 + cn]:
        clu.elems[i] = f[0]
        ar = np.asarray([float(g) for g in f[1:4]])
        clu.atoms[i] = ar
        i = i + 1
    if iprint:
        print ('atoms loaded: %d' % clu.n)
    # clu.atoms -= clu.atoms.min(axis=0).reshape(1, 3)
    return clu

def read_surface_data(sname):
    from utils.io import read_json
    surf = read_surface(sname + ".xyz")
    ipcs = read_json(sname + ".json")
    surf.space_group = ipcs["space_group"]
    surf.space_group_ref = ipcs["space_group_ref"]
    surf.unit_cell = ipcs["unit_cell"]
    surf.fix = ipcs.get("fix", "")
    return surf

def write_surface_data(sname, surf):
    from utils.main import RecordWritter
    RecordWritter().write_surfs(sname + ".xyz", sname + ".json", [ClusterAtSurface(0, surf=surf)])

# sample uniformly in the (half) unit sphere
def random_direction(half=False, quarter=False):
    d = np.random.random(size=2)
    phi = 2 * np.pi * d[0]
    if quarter:
        st = np.sqrt(2.0) / 2.0
        theta = np.arccos(d[1] * (1 - st) + st) # st~1
    elif half:
        theta = np.arccos(d[1]) # 0~1
    else:
        theta = np.arccos(2 * d[1] - 1) # -1~1
    return np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta),
                     np.cos(theta)])

# make sure cluster atoms are not below surface
def surface_clus_rel(cc):
    assert isinstance(cc, ClusterAtSurface)
    if cc.surf.n > 0:
        zmin = cc.surf.atoms[:, 2].min()
        for i in range(cc.n):
            if cc.atoms[i, 2] < zmin:
                cc.atoms[i, 2] += cc.surf.cell[-1]

class ClusterAtSurface(object):

    eps = np.finfo(dtype=float).eps * 100
    def __init__(self, n, surf):
        self.n = n
        self.atoms = np.zeros((n, 3))
        self.elems = np.array([''] * n, dtype='|S20')
        self.energy = 0.0
        self.label = ''
        self.mag = None
        self.forces = None
        self.surf = surf

    def to_cluster(self):
        from cluster.base import Cluster
        cc = Cluster(self.n + self.surf.n)
        cc.surfnum = self.surf.n
        cc.atoms = np.array(list(self.surf.atoms) + list(self.atoms))
        cc.elems = np.array(list(self.surf.elems) + list(self.elems))
        if self.surf.forces is not None and self.forces is not None:
            cc.forces = np.array(list(self.surf.forces) + list(self.forces))
        cc.energy = self.energy
        cc.mag = self.mag
        cc.label = self.label
        return cc

    def pointgroup(self):
        return 'C1'

    @staticmethod
    def from_cluster(clus, surf):
        c = ClusterAtSurface(clus.n - clus.surfnum, surf=surf)
        c.atoms = clus.atoms[clus.surfnum:]
        c.elems = clus.elems[clus.surfnum:]
        if clus.forces is not None:
            c.forces = clus.forces[clus.surfnum:]
        c.mag = clus.mag
        c.label = clus.label
        c.energy = clus.energy
        if surf.n != clus.surfnum:
            # surface replacement rule: the original two surfaces should be
            # A belongs to B. During all kinds of operations in the program,
            # the surface will not be shifted.
            # Therefore, we can first match the atoms here (directly),
            # and then add other atoms not included.
            # This does not require the all not included atoms to be after the atoms
            # of the small surface.
            print ("transforming %s ..." % c.label)
            from surface.surf_comp import surface_match
            slsurf = Surface(clus.surfnum)
            slsurf.atoms = clus.atoms[:clus.surfnum]
            slsurf.elems = clus.elems[:clus.surfnum]
            slsurf.cell = surf.cell
            dmax = 0.5
            for itx in range(0, 100):
                assert itx != 99
                _, mt = surface_match(surf, slsurf, dmax)
                mt = list(mt - 1)
                if mt.count(-1) == max(surf.n - slsurf.n, 0):
                    break
                dmax += 0.2
            c.surf = copy.deepcopy(surf)
            # the order will be the same as new (larger surface)
            # otherwise the fix will be problematic
            for ii, i in enumerate(mt):
                if i != -1:
                    c.surf.atoms[ii] = clus.atoms[i]
                    c.surf.elems[ii] = clus.elems[i]
        else:
            c.surf = copy.deepcopy(surf)
            c.surf.atoms = clus.atoms[:clus.surfnum]
            c.surf.elems = clus.elems[:clus.surfnum]
            if clus.forces is not None:
                c.surf.forces = clus.forces[:clus.surfnum]
        surface_clus_rel(c)
        return c

    # write coordinates to file
    def write_xyz(self, fn, append=True):
        f = open(fn, 'a' if append else 'w')
        if self.label == '':
            title = '%d\nE = %15.8f' % (self.n + self.surf.n, self.energy)
        else:
            title = '%d\n%s = %15.8f' % (self.n + self.surf.n, self.label, self.energy)
        title += " SURF = %d\n" % self.surf.n
        f.write(title)
        if self.surf.forces is not None:
            for x, y, z in zip(self.surf.elems, self.surf.atoms, self.surf.forces):
                f.write('%7s%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\n' %
                        (x, y[0], y[1], y[2], z[0], z[1], z[2]))
        else:
            for x, y in zip(self.surf.elems, self.surf.atoms):
                f.write('%7s%15.8f%15.8f%15.8f\n' % (x, y[0], y[1], y[2]))
        if self.forces is not None:
            for x, y, z in zip(self.elems, self.atoms, self.forces):
                f.write('%7s%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\n' %
                        (x, y[0], y[1], y[2], z[0], z[1], z[2]))
        else:
            for x, y in zip(self.elems, self.atoms):
                f.write('%7s%15.8f%15.8f%15.8f\n' % (x, y[0], y[1], y[2]))
        f.close()

if __name__ == "__main__":
    x = read_surface("tests/surface/newcoord.xyz")
    n = 5
    cl = ClusterAtSurface(n, x)
    xmean = {}
    xsigma = {}
    import sys
    sys.path.insert(0, sys.path[0] + '/..')
    from cluster.base import update_stat
    update_stat(["Pt", "Mg", "O"], xmean, xsigma, 0.04)
    cl.create(["Pt"] * n, xmean, xsigma)
    cl.write_xyz(fn='test.xyz', append=False)
    for i in range(20):
        print (i)
        cl.create(["Pt"] * n, xmean, xsigma)
        cl.write_xyz(fn='test.xyz', append=True)