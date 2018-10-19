
from __future__ import print_function

if __name__ == "__main__":
    import sys
    sys.path.insert(0, sys.path[0] + '/..')

import numpy as np
from surface.create_cluster import random_direction
from surface.surf_symm import to_direct, to_cartesian, to_cellmat
from surface.surf_comp import num_adjust

def enlarge(c, size, cell):
    satoms = []
    catoms = []
    selems = []
    celems = []
    sforces = []
    cforces = []
    diratoms = to_direct(c.atoms, cell)
    for x in range(size[0]):
        for y in range(size[1]):
            for z in range(size[2]):
                xatoms = diratoms + np.array([[x, y, z]], dtype=float)
                satoms += list(xatoms[:c.surfnum])
                catoms += list(xatoms[c.surfnum:])
                selems += list(c.elems[:c.surfnum])
                celems += list(c.elems[c.surfnum:])
                if c.forces is not None:
                    sforces += list(c.forces[:c.surfnum])
                    cforces += list(c.forces[c.surfnum:])
    c.atoms = to_cartesian(np.array(satoms + catoms), cell)
    c.elems = np.array(selems + celems)
    if c.forces is not None:
        c.forces = np.array(sforces + cforces)
    c.n = len(c.atoms)
    c.surfnum = c.surfnum * size[0] * size[1] * size[2]

class CreatePeriodic(object):
    eps = np.finfo(dtype=float).eps * 100
    def __init__(self, n=0, atoms=None, elems=None, surf=None):
        self.n = 0
        self.atoms = np.zeros((n, 3), dtype=float) if atoms is None else atoms
        self.elems = np.array([''] * n, dtype='|S20') if elems is None else elems
        self.surf = surf
        self.ext_range = None

    def to_cluster(self):
        from surface.base import ClusterAtSurface
        from copy import deepcopy
        c = ClusterAtSurface(self.n, surf=deepcopy(self.surf))
        c.elems = self.elems
        c.atoms = self.atoms
        return c

    @staticmethod
    def generator(surf, **opts):
        while True:
            c = CreatePeriodic(surf=surf)
            c.create(**opts)
            yield c.to_cluster()

    # first order atom
    def add_atom_f(self, telem, clelems, cura, mean, sigma, suelems, suatoms,
                   soelems=None, soatoms=None, atst=None, quarter=False, ext_range=None):
        elemr = np.array(list(clelems[:cura]) + list(suelems))
        atomr = np.array(list(self.atoms[:cura]) + list(suatoms))
        cellmat, zmax, start_height, hmax = atst
        while True:
            atst = np.array([0.0, 0.0, zmax + np.random.random() * start_height])
            atst += np.random.random() * cellmat[0]
            atst += np.random.random() * cellmat[1]
            x = 0.0
            kdir = random_direction(quarter=quarter, ext_range=ext_range)
            for j in range(0, len(elemr)):
                m = telem + '-' + elemr[j]
                l = np.random.normal(mean[m][0], sigma[m][0])
                va = atomr[j] - atst
                b = -2 * np.dot(va, kdir)
                c = np.square(va).sum() - l ** 2
                delta = b**2 - 4*c + CreatePeriodic.eps
                if delta > 0.0:
                    xx = (-b + np.sqrt(delta)) / 2
                    if xx > x:
                        x = xx
            xatom = atst + x * kdir
            if xatom[2] > hmax + zmax:
                continue
            okay = True
            if soelems is not None:
                for j in range(0, len(soelems)):
                    m = telem + '-' + soelems[j]
                    kl = np.random.normal(mean[m][0], sigma[m][0])
                    kla = np.linalg.norm(soatoms[j] - xatom)
                    if kla < kl:
                        okay = False
                        break
            if okay:
                break
        return xatom

    # make sure the new atom is in the primary cell range
    def image_ex(self, atomx, cell):
        xdir = to_direct(np.array([atomx]), cell)[0]
        for k in range(2):
            xdir[k] += num_adjust(0.5 - xdir[k], 1.0)
        return to_cartesian(np.array([xdir]), cell)[0]

    def image_cr(self, atomx, elemx, cell):
        adir = []
        xdir = to_direct(np.array([atomx]), cell)[0]
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                if x == 0 and y == 0:
                    continue
                aadd = np.array([[x, y, 0.0]])
                adir.append(xdir + aadd)
        atoms = to_cartesian(np.array(adir), cell)
        return atoms, np.array([elemx] * len(adir))

    def create(self, elems, mean, sigma, layer_height=0.1, start_height=1.0, max_height=3.5):
        self.elems = elems
        self.n = len(elems)
        self.atoms = np.zeros((self.n, 3))
        # surface part
        # enlarge surface to avoid drop off
        surfabak = np.array(self.surf.atoms)
        surfebak = np.array(self.surf.elems)
        atomsu = np.zeros((surfabak.shape[0] * 18, 3))
        elemsu = []
        cellmat = to_cellmat(self.surf.cell)
        cellmat[2, 2] = self.surf.cellz
        for ix, x in enumerate([-1, 0, 1]):
            for iy, y in enumerate([-1, 0, 1]):
                for iz, z in enumerate([-1, 0]):
                    ik = ix * 6 + iy * 2 + iz
                    surfadd = np.array([x * cellmat[0] + y * cellmat[1] + z * cellmat[2]])
                    atomsu[ik * surfabak.shape[0]:(ik + 1) * surfabak.shape[0]] = \
                        surfabak + surfadd
                    elemsu += list(surfebak)
        elemsu = np.array(elemsu)
        zmax = self.surf.atoms.max(axis=0)[2]
        xs_core = [ia for ia, a in enumerate(self.surf.atoms) if a[2] >= zmax - layer_height]
        xs_other = [ia for ia, a in enumerate(atomsu) if ia not in xs_core]
        soelems, soatoms = elemsu[xs_other], atomsu[xs_other]
        suelems, suatoms = self.surf.elems[xs_core], self.surf.atoms[xs_core]
        cura = 0 # current atom
        while cura != self.n:
            if cura == 0:
                qua, exr = True, None
            else:
                qua, exr = False, self.ext_range
            xatom = self.add_atom_f(self.elems[cura], self.elems, cura, mean, sigma, suelems, \
                suatoms, soelems=soelems, soatoms=soatoms, atst=[cellmat, zmax, start_height, \
                max_height], quarter=qua, ext_range=exr)
            xatom = self.image_ex(xatom, self.surf.cell)
            self.atoms[cura] = np.array(xatom)
            agatom, agelem = self.image_cr(xatom, self.elems[cura], self.surf.cell)
            soatoms = np.array(list(soatoms) + list(agatom))
            soelems = np.array(list(soelems) + list(agelem))
            cura += 1
        return self

if __name__ == "__main__":
    from cluster.base import elem_char, update_stat
    from surface.base import read_surface
    surf = read_surface('../PGOPT/surfaces/wcut.xyz')
    surf.unit_cell = np.array([3.172, 4.4859, 4.485877])
    surf.space_group = "-C 2b 2"
    xmean = {}
    xsigma = {}
    np.random.seed(100)
    k, _ = elem_char('B16')
    update_stat(np.array(list(k) + ["W"]), xmean, xsigma, 0.100)
    print (xmean, xsigma)
    app = False
    for i in range(100):
        print (i)
        c = CreatePeriodic()
        c.surf = surf
        c.create(k, xmean, xsigma).to_cluster().write_xyz('test.xyz', append=app)
        app = True
    