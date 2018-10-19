
from __future__ import print_function

if __name__ == "__main__":
  import sys
  sys.path.insert(0, sys.path[0] + '/..')

import numpy as np
from surface.spacegroupdata import SymOpsHall
from cluster.base import Cluster

TOLERANCE = 0.05
TOLX = 1E-6

def apply_trans(xx, tr, unit_cell):
    xx = xx / np.array(unit_cell).reshape(1, 3)
    x = np.array(xx[:, 0])
    y = np.array(xx[:, 1])
    z = np.array(xx[:, 2])
    for ic, c in enumerate(tr):
        c = c.replace("1/", "1.0/")
        xx[:, ic] = eval(c)
    return xx * np.array(unit_cell).reshape(1, 3)

def rev_num_list(l):
    if l == "-1":
        return []
    xl = []
    for il in l.split(","):
        if "-" in il:
            a, b = il.split("-")
            xl += range(int(a), int(b) + 1)
        else:
            xl += [int(il)]
    return xl

def num_list(l):
    kl = []
    for il in sorted(l):
        if len(kl) != 0:
            if il == kl[-1][1] + 1:
                kl[-1][1] = il
                continue
        kl.append([il, il])
    return ",".join([("%d" % i[0]) if i[0] == i[1] else ("%d-%d" % tuple(i)) for i in kl])

# return cluster, 1d-cell, fixed list
def read_contcar(fn):
    f = open(fn, 'r').readlines()
    ff = [l.strip() for l in f]
    cell = [[float(g) for g in f.split(' ') if len(g) != 0] for f in ff[2:5]]
    cell = np.array(cell)
    if np.linalg.norm(np.array([cell[0, 1], cell[1, 0]])) < TOLX:
        dcell = np.diag(cell)
        print (' CELL  = %15.7f%15.7f%15.7f' % tuple(dcell))
    else:
        dcell = np.array([cell[0, 0], cell[0, 1], cell[1, 0], cell[1, 1], cell[2, 2]])
        print (' CELL  = [%15.7f, %15.7f] [%15.7f, %15.7f] %15.7f' % tuple(dcell))
    elems = np.array([g for g in ff[5].split(' ') if len(g) != 0])
    elemsn = np.array([int(g) for g in ff[6].split(' ') if len(g) != 0])
    print (' ELEMS = ' + " ".join(["%s %d" % (i, j) for i, j in zip(elems, elemsn)]))
    print (' TOTAL = %d' % elemsn.sum())
    gelem = []
    for el, en in zip(elems, elemsn):
        gelem += [el] * en
    if ff[7].startswith("Selective"):
        lx = 9
    else:
        lx = 8
    assert ff[lx - 1] in ["Direct", "Cartesian"]
    gf = ff[lx - 1] == "Direct"
    gff = [[g for g in f.split(' ') if len(g) != 0] for f in ff[lx:lx+elemsn.sum()]]
    gff = [f[:6] for f in gff]
    coords = np.array([[float(g) if g not in ['T', 'F'] else 0.0 for g in f] for f in gff])
    rf = np.array([["relax" if g == 'T' else "fix" for g in f] for f in gff])
    fxl = [iif for iif, f in enumerate(rf) if len(f) >= 4 and f[5] == "fix"]
    lf = num_list(fxl)
    print (' FIX   = %s' % lf)
    if gf:
        coords = coords[:, 0:3].dot(cell)
    else:
        coords = coords[:, 0:3]
    assert len(gelem) == len(coords)
    c = Cluster(elemsn.sum(), np.array(gelem), coords)
    c.label = "MAG:0:0.0:1"
    return c, dcell, fxl

def to_cellmat(dcell):
    if len(dcell) == 3:
        return np.diag(dcell)
    else:
        return np.array([[dcell[0], dcell[1], 0.0], [dcell[2], dcell[3], 0.0],
                         [0.0, 0.0, dcell[4]]])

def to_dcell(cellmat):
    if np.linalg.norm(np.array([cellmat[0, 1], cellmat[1, 0]])) < TOLX:
        return np.diag(cellmat)
    else:
        return np.array([cellmat[0, 0], cellmat[0, 1], cellmat[1, 0],
                         cellmat[1, 1], cellmat[2, 2]])

def to_direct(coords, dcell):
    m = to_cellmat(dcell)
    return coords.dot(np.linalg.inv(m))

def to_cartesian(directs, dcell):
    m = to_cellmat(dcell)
    return directs.dot(m)

def num_image(x, l, cl):
    k = np.round(np.abs(x) / l) * l
    return np.abs(np.abs(x) - k) < TOLERANCE / cl

def num_move_in(x, l):
    if x < -TOLX or x >= l - TOLX:
        xi = np.floor(x / l) * l
        x -= xi
    if x >= l - TOLX:
        x -= l
    return x

def check_image(ata, ela, atb, elb, lcell, d, cl):
    if ela != elb:
        return False
    for dx in range(0, 3):
        if dx != d:
            if not num_image(ata[dx] - atb[dx], 1.0, cl[dx]):
                return False
        else:
            if not num_image(ata[dx] - atb[dx], lcell, cl[dx]):
                return False
    return True

def check_exist(cat, cel, at, el, lcell, cl):
    for i, j in zip(cat, cel):
        if j != el:
            continue
        rep = True
        for dx in range(0, 3):
            if lcell[dx] is None:
                if np.abs(i[dx] - at[dx]) > TOLERANCE / cl[dx]:
                    rep = False
                    break
            else:
                if not num_image(i[dx] - at[dx], lcell[dx], cl[dx]):
                    rep = False
                    break
        if rep:
            return True
    return False

def cell_move_in(cat, ucell):
    for i in cat:
        for k in range(0, 3):
            if ucell[k] is not None:
                i[k] = num_move_in(i[k], ucell[k])

def get_unit_cell(surf, ucellmat):
    sa, sb = np.array(surf.atoms), np.array(surf.elems)
    sa = to_direct(sa, surf.cell)
    ucell = np.diag(to_direct(ucellmat, surf.cell))
    cl = np.linalg.norm(to_cellmat(surf.cell), axis=1)
    gcell = [ucell[0], ucell[1], None]
    xa, xb = -TOLX + sa.min(axis=0)[0], gcell[0] + 2 * TOLERANCE + sa.min(axis=0)[0]
    ya, yb = -TOLX + sa.min(axis=0)[1], gcell[1] + 2 * TOLERANCE + sa.min(axis=0)[1]
    cell_at = []
    cell_el = []
    for i, j in zip(sa, sb):
        if i[0] > xa and i[0] < xb and i[1] > ya and i[1] < yb:
            if not check_exist(cell_at, cell_el, i, j, gcell, cl):
                cell_at.append(i)
                cell_el.append(j)
    cell_move_in(cell_at, gcell)
    if len(cell_at) * int(np.round(1.0 / ucell[0])) * int(np.round(1.0 / ucell[1])) != surf.n:
        raise RuntimeError("Surface x-y shifted! Cannot get unit cell!")
    return [to_cartesian(np.array(cell_at), surf.cell), cell_el]

# return matrix unit cell size
def determine_unit_cell(surf):
    sa, sb = np.array(surf.atoms), np.array(surf.elems)
    cmat = to_cellmat(surf.cell)
    cl = np.linalg.norm(cmat, axis=1)
    sa = to_direct(sa, surf.cell)
    sa = sa - sa.min(axis=0).reshape((1, 3))
    dx = np.array(cmat)
    for d in range(0, 2):
        for du in range(1, 11):
            gcell = [1.0, 1.0, None]
            gcell[d] = 1.0 / du
            if surf.n % du != 0:
                continue
            xa, xb = -TOLX, 1.0 / du + 2 * TOLERANCE / cl[d]
            cell_at = []
            cell_el = []
            other_at = []
            other_el = []
            for i, j in zip(sa, sb):
                lin = False
                if i[d] > xa and i[d] < xb:
                    if not check_exist(cell_at, cell_el, i, j, gcell, cl):
                        lin = True
                if lin:
                    cell_at.append(i)
                    cell_el.append(j)
                else:
                    other_at.append(i)
                    other_el.append(j)
            if len(cell_at) != surf.n / du:
                continue
            okay = True
            for i, j in zip(other_at, other_el):
                found = False
                for ip, jp, in zip(cell_at, cell_el):
                    if check_image(ip, jp, i, j, 1.0 / du, d, cl):
                        found = True
                        break
                if not found:
                    okay = False
                    break
            if okay:
                dx[d] = cmat[d] / du
    return dx

def num_distance(x, l):
    k = np.round(x / l) * l
    return x - k

def distance_in_cell(ata, ela, atb, elb, lcell):
    # if ela != elb:
    #     return np.linalg.norm(lcell)
    d = np.zeros((3, ))
    for dx in range(0, 3):
        d[dx] = num_distance(ata[dx] - atb[dx], 1.0) * lcell[dx]
    return np.linalg.norm(d)

# result will be exact lcell dist from dirs
# geometry similar to ref
def move_ref(dirs, ref, lcell):
    dirs = np.array(dirs)
    for dx in range(0, 3):
        dirs[dx] = ref[dx] + num_distance(dirs[dx] - ref[dx], lcell[dx])
    return dirs

class EqGroup(object):
    main = -1
    others = {}
    def __init__(self, main):
        self.main = main
        self.others = {}
    def __repr__(self):
        return "EgGroup(main = %d, others = %s)\n" % (self.main, self.others.__repr__())

def eqgroup_find(groups, x):
    for ig, g in enumerate(groups):
        if x == g.main or x in g.others:
            return ig
    return None

def inv_trans(x):
    if x == ['-y', 'x-y', 'z']:
        return ['-x+y', '-x', 'z']
    elif x == ['-x+y', '-x', 'z']:
        return ['-y', 'x-y', 'z']
    else:
        raise NotImplementedError

def inv_trans_chain(x):
    return map(inv_trans, x[::-1])

def eqgroup_merge(groups, im, io, op):
    imain = eqgroup_find(groups, im)
    iother = eqgroup_find(groups, io)
    g = groups[imain]
    go = groups[iother]
    chain = []
    if im != g.main:
        chain += g.others[im]
    chain += [op]
    if io != go.main:
        chain += inv_trans_chain(go.others[io])
    g.others[go.main] = chain
    for k, v in go.others.items():
        g.others[k] = chain + v
    del groups[iother]

def check_space_group_accuracy(surf):
    ucellmat = to_cellmat(surf.unit_cell)
    ucell = np.diag(to_direct(ucellmat, surf.cell))
    cellmat = to_cellmat(surf.cell)
    at, el = [np.array(i) for i in get_unit_cell(surf, ucellmat)]
    atbase = to_direct(at, surf.cell)
    direct = to_direct(surf.atoms, surf.cell)
    lcell = np.array([np.linalg.norm(c) for c in cellmat])
    ratiox = np.array([np.linalg.norm(ucellmat[i]) / np.linalg.norm(cellmat[i]) for i in range(0, 2)])
    nratio = map(int, map(np.round, 1 / ratiox))
    lrucell = np.array([1.0 / nratio[0], 1.0 / nratio[1], 1.0])
    lucell = np.array([lcell[0] / nratio[0], lcell[1] / nratio[1], 1.0])
    print('nratio = ', nratio)
    # print('atbase = ', atbase)
    # print('direct = ', direct)
    selected = [False] * surf.n
    atnew = np.zeros_like(atbase)
    maxd = 0.0
    for ix in range(0, nratio[0]):
        for iy in range(0, nratio[1]):
            for i in range(0, len(atbase)):
                x = atbase[i] + np.array([ix * ratiox[0], iy * ratiox[1], 0.0])
                d = [distance_in_cell(x, el[i], y, surf.elems[jy], lcell) for jy, y in enumerate(direct)]
                idx = np.argmin(d)
                # print("%3d %3d %3d %5d %10.5f" % (ix, iy, i, idx, np.min(d)))
                maxd = max(maxd, np.min(d))
                atnew[i] += move_ref(direct[idx], atbase[i], lrucell)
                assert selected[idx] == False
                selected[idx] = True
    print('maxd = ', maxd)
    maxd = 0.0
    atbase = atnew / (nratio[0] * nratio[1])
    atx = np.array([(x % lrucell + lrucell) % lrucell for x in atbase])
    ref = np.array(atx[np.argmin(np.linalg.norm(atx[:, 0:2], axis=1))])
    ref[2] = 0.0
    print('ref = ', ref)
    atx = np.array([x - ref for x in atx])
    atx = to_direct(to_cartesian(atx, surf.cell), surf.unit_cell)
    spg = [[x.strip() for x in s] for s in SymOpsHall[surf.space_group] if s[2].strip() == 'z']
    groups = [EqGroup(x) for x in range(0, len(atx))]
    spgfix_base = [False] * len(atx)
    for i in spg:
        si = ":".join(i)
        xat = apply_trans(atx, i, ucell)
        okay = True
        for ii, (xi, xj) in enumerate(zip(xat, el)):
            d = [distance_in_cell(xi, xj, aa, ee, lucell) for aa, ee in zip(atx, el)]
            # print(i, ii, np.argmin(d), min(d))
            idx = np.argmin(d)
            if ii == idx and (i == ['-x+y', '-x', 'z'] or i == ['-y', 'x-y', 'z']):
                spgfix_base[ii] = True
                pp = np.array(atx[idx])
                atx[idx] = move_ref(np.array([0.0, 0.0, atx[idx][2]]), atx[idx],
                    np.array([1.0 / 3, 1.0 / 3, 1.0]))
                print(pp, '->', atx[idx])
    for i in spg:
        si = ":".join(i)
        xat = apply_trans(atx, i, ucell)
        okay = True
        for ii, (xi, xj) in enumerate(zip(xat, el)):
            d = [distance_in_cell(xi, xj, aa, ee, lucell) for aa, ee in zip(atx, el)]
            print(i, ii, np.argmin(d), min(d))
            idx = np.argmin(d)
            if eqgroup_find(groups, idx) != eqgroup_find(groups, ii):
                gi = eqgroup_find(groups, ii)
                gj = eqgroup_find(groups, idx)
                # print ('merge : %d - %d, %s' % (groups[gi].main, groups[gj].main, str(i)))
                eqgroup_merge(groups, ii, idx, i)
    print(groups)
    for g in groups:
        for k, v in g.others.items():
            gat = np.array([atx[g.main]])
            for iv in v:
                gat = apply_trans(gat, iv, ucell)
            atx[k] = move_ref(gat[0], atx[k], lrucell)
    for i in spg:
        si = ":".join(i)
        xat = apply_trans(atx, i, ucell)
        okay = True
        for ii, (xi, xj) in enumerate(zip(xat, el)):
            d = [distance_in_cell(xi, xj, aa, ee, lucell) for aa, ee in zip(atx, el)]
            print(i, ii, np.argmin(d), min(d))

            # if ii <= idx:
            #     atx[idx] = (move_ref(xi, atx[idx], [1.0, 1.0, 1.0]) + atx[idx]) / 2
    atbase = to_direct(to_cartesian(atx, surf.unit_cell), surf.cell)
    selected = [False] * surf.n
    for ix in range(0, nratio[0]):
        for iy in range(0, nratio[1]):
            for i in range(0, len(atbase)):
                x = atbase[i] + np.array([ix * ratiox[0], iy * ratiox[1], 0.0])
                d = [distance_in_cell(x, el[i], y, surf.elems[jy], lcell) for jy, y in enumerate(direct)]
                idx = np.argmin(d)
                maxd = max(maxd, np.min(d))
                # print("%3d %3d %3d %5d %10.5f" % (ix, iy, i, idx, np.min(d)))
                # print (move_ref(x, direct[idx], lrucell) - x)
                direct[idx] = move_ref(x, direct[idx], lrucell)
                direct[idx] = x
                assert selected[idx] == False
                selected[idx] = True
    print('maxd = ', maxd)
    maxd = 0.0
    selected = [False] * surf.n
    spgfix = [False] * surf.n
    for ix in range(0, nratio[0]):
        for iy in range(0, nratio[1]):
            for i in range(0, len(atbase)):
                x = atbase[i] + np.array([ix * ratiox[0], iy * ratiox[1], 0.0])
                d = [distance_in_cell(x, el[i], y, surf.elems[jy], lcell) for jy, y in enumerate(direct)]
                idx = np.argmin(d)
                maxd = max(maxd, np.min(d))
                # print("%3d %3d %3d %5d %10.5f" % (ix, iy, i, idx, np.min(d)))
                assert selected[idx] == False
                if spgfix_base[i]:
                    spgfix[idx] = True
                selected[idx] = True
    print('maxd = ', maxd)
    spgfixnums = [i for i in range(0, surf.n) if spgfix[i]]
    print('spgfix = ' + num_list(spgfixnums))
    surf.atoms = to_cartesian(direct, surf.cell)

def determine_space_group(surf, ucellmat):
    at, el = [np.array(i) for i in get_unit_cell(surf, ucellmat)]
    atbase = to_direct(at, surf.cell)
    # when determine space group, must move one atom to origin
    # it should be determined. (which one (`for u`))
    ucell = np.diag(to_direct(ucellmat, surf.cell))
    cl = np.linalg.norm(to_cellmat(surf.cell), axis=1)
    gcell = [ucell[0], ucell[1], None]
    glist = []
    sg_trials = {}
    for k, v in SymOpsHall.items():
        sg = [[x.strip() for x in s] for s in v if s[2].strip() == 'z']
        # for not square unit cell, x and y cannot be mixed
        # if np.abs(cl[0] * ucell[0] - cl[1] * ucell[1]) > TOLX:
        #     sg = [[x.strip() for x in s] for s in sg \
        #         if "y" not in s[0] and "x" not in s[1]]
        sg_trials[k] = (sg, v)
    for u in range(-1, len(atbase)):
        print ('det space group ...', u, '/', len(atbase))
        at = np.array(atbase)
        if u == -1:
            atref = np.zeros_like(at[0, 0:2])
        else:
            atref = np.array(at[u, 0:2])
            at[:, 0:2] -= at[u, 0:2]
        tdict = {}
        for k, (sg, v) in sg_trials.items():
            ix = 0
            for i in sg:
                si = ":".join(i)
                if si not in tdict:
                    xat = apply_trans(at, i, ucell)
                    okay = True
                    for xi, xj in zip(xat, el):
                        if not check_exist(at, el, xi, xj, gcell, cl):
                            okay = False
                            break
                    tdict[si] = okay
                if tdict[si]:
                    ix += 1
            if ix == len(sg):
                glist.append([k, ix, len(v), atref, len(atbase)])
    glist.sort(key=lambda x: [-x[1], x[2]])
    return [g for g in glist if g[1] == glist[0][1]]

if __name__ == "__main__":
    from surface.base import read_surface
    surf = read_surface('../PGOPT/surfaces/alumina.xyz')
    ucell = determine_unit_cell(surf)
    sp = determine_space_group(surf, ucell)
    print (sp)
    sp = sp[0]
    print (SymOpsHall[sp[0]])
