#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import math
import time
import shlex
import json
import thread
import psutil
import socket
from psutil import Popen
from subprocess import PIPE, STDOUT
from sys import argv
import numpy as np
from shutil import copyfile
import zipfile


class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


sys.stdout = Unbuffered(sys.stdout)


class QMInput(object):
    Method = 'pbe'
    Coords = ''
    Title = ''
    Charge = 0
    Multiplicity = 0
    Para = 1
    Check = 'OUT'
    Opt = 0  # opt step
    NebOpt = 0  # neb step
    Scf = []
    Opts = {}
    Line = ''
    RunType = None
    mirablocks = ['', '', '']  # block, corner, shape

    def __init__(self):
        fin = sys.stdin
        while True:
            line = fin.readline().strip()
            if len(line) == 0:
                break
            if line[0] == '!':
                pass
            if line[0] == '%':
                xs = line[1:].split('=')
                xsn = xs[0].strip().lower()
                if xsn == 'workers' or xsn == 'lindaworkers':
                    nodes = xs[1].strip().split(',')
                    self.Para = 0
                    for n in nodes:
                        ns = n.strip().split(':')
                        if len(ns) == 1:
                            self.Para += 1
                        else:
                            num = int(ns[1])
                            self.Para += num
                elif xsn == 'mirablocks':
                    self.mirablocks = xs[1].strip().split(',')
                elif xsn == 'check' or xsn == 'chk':
                    self.Check = xs[1].strip()
                elif xsn == 'nproc' or xsn == 'nprocshared':
                    xs = int(xs[1].strip())
                    self.Para = xs
            if line[0] == '#':
                self.Line = line.strip()
                xs = line[1:].strip()
                xsr = self._opts(xs)
                for x in xsr:
                    if x[1] == '' and '/' in x[0]:
                        xl = x[0].split('/')
                        self.Method = xl[0].strip().lower()
                    elif x[1] != '' and x[0].lower() == 'coords':
                        self.Coords = x[1].strip()
                    elif x[1] != '' and x[0].lower() == 'scf':
                        self.Scf = self._opts(x[1])
                    elif x[0].lower() in ['opt', 'neb']:
                        xgr = self._opts(x[1])
                        fp = {}
                        for g in xgr:
                            if g[0] == 'iter':
                                if x[0].lower() == "opt":
                                    self.Opt = int(g[1])
                                else:
                                    self.NebOpt = int(g[1])
                            else:
                                fp[g[0]] = g[1]
                        if len(fp) != 0:
                            self.Opts[x[0].lower()] = fp
                    else:
                        self.Opts[x[0].lower()] = x[1]

        while True:
            line = fin.readline().strip()
            if len(line) == 0:
                break
            else:
                self.Title = line

        while True:
            line = fin.readline().strip()
            if len(line) == 0:
                break
            else:
                self.Charge = int(line.split(' ')[0])
                self.Multiplicity = int(line.split(' ')[1])

    @staticmethod
    def _opts(xs):
        xsr = []
        xcur = ['', '']
        pcnt = 0
        pp = 0
        for x in xs:
            if x == '(' and ('/' not in xcur[0]):
                pcnt += 1
                if pcnt > 1:
                    xcur[pp] += x
                else:
                    pp = 1
            elif x == ')' and ('/' not in xcur[0]):
                pcnt -= 1
                if pcnt == 0:
                    xsr.append(xcur)
                    xcur = ['', '']
                    pp = 0
                else:
                    xcur[pp] += x
            elif x == '=':
                if pp == 0:
                    pp = 1
                else:
                    xcur[1] += x
            elif x == ' ' or x == ',':
                if pcnt == 0 and len(xcur[0]) != 0:
                    xsr.append(xcur)
                    xcur = ['', '']
                    pp = 0
                elif pp == 1:
                    xcur[1] += x
            else:
                xcur[pp] += x
        if pcnt == 0 and len(xcur[0]) != 0:
            xsr.append(xcur)
        return xsr

# PBE/ continue encut=500 scf(iter=200) opt(iter=300) cell=15 fix(Mg,O) nocenter
# norun
# use nocenter to avoid auto centering for clusters


httoev = 27.21138505


def solve_list(k):
    p = []
    kk = k.split(",")
    for i in kk:
        if '-' in i and i[0] != '-':
            a, b = map(int, i.split('-'))
            for x in range(a, b + 1):
                p.append(x)
        elif i.isdigit() or (i[0] == '-' and i[1:].isdigit()):
            p.append(int(i))
        else:
            p.append(i)
    return p


class Struct(object):
    Coords = []
    Forces = []
    Energy = 0.0
    Mag = 0.0
    Id = -1


class StructTraj(object):
    Structs = []
    ScfNum = []
    RelaxConv = False
    ScfConv = False

    def write_structs(self, fn, elems, pstep, surfnum, cont, last, zf=None):
        if last:
            strs = self.Structs[-1:]
        elif cont:
            strs = self.Structs[1:]
        else:
            strs = self.Structs
        assert (not last) or (not cont)
        gl = []
        exx = [ep + [ipx] for ipx, ep in enumerate(elems)]
        exx.sort(key=lambda x: x[0])
        for t in strs:
            if t.Forces == []:
                t.Forces = np.zeros(np.array(t.Coords).shape, dtype=float)
            gl.append("%d\n" % len(elems))
            cmt = "MAG::%.2f:%d = %.8f" % (t.Mag, t.Id + pstep, t.Energy)
            if surfnum != 0:
                cmt += " SURF = %d" % surfnum
            gl.append(cmt + "\n")
            for _, j, k in exx:
                gl.append("%7s%18.10f%18.10f%18.10f%18.10f%18.10f%18.10f\n" %
                          ((j, ) + tuple(t.Coords[k]) + tuple(t.Forces[k])))
        if zf is None:
            with open(fn, "a" if cont else "w") as g:
                for l in gl:
                    g.write(l)
        else:
            if fn not in zf.namelist():
                zf.writestr(fn, "".join(gl))

    def write_scfnum(self, fn, cont, zf=None):
        if cont:
            scfn = self.ScfNum[1:]
        else:
            scfn = self.ScfNum
        if zf is None:
            with open(fn, "a" if cont else "w") as g:
                for i in scfn:
                    g.write("%5d " % i)
        else:
            if fn not in zf.namelist():
                zf.writestr(fn, "".join(["%5d " % i for i in scfn]))

    @staticmethod
    def read_outcar(f, nocenter, dcell):
        pmag = 0.0
        state = 0
        kk = None  # coords
        kf = None  # forces
        st = None
        scfn = 0
        traj = StructTraj()
        traj.Structs = []
        traj.ScfNum = []
        traj.RelaxConv = False
        traj.ScfConv = False
        for line in f:
            line_ori = line
            line = line.strip()
            if line.startswith("number of electron"):
                ll = [x for x in line.split(" ") if len(x) != 0]
                pmag = float(ll[5]) if len(ll) > 5 else 0.0
            elif "- Iteration" in line:
                scfn = int(line.split("(")[1].split(")")[0].strip())
            elif line.startswith("POSITION"):
                state = 1
            elif line.startswith("------------") and state == 1:
                state = 2
                kk = []
                kf = []
            elif line.startswith("------------") and state == 2:
                state = 3
                st = Struct()
                if not nocenter:  # "nocenter" not in ip.Opts:
                    kk = np.array(kk) - np.array([dcell / 2])
                st.Coords = np.array(kk)
                st.Forces = np.array(kf)
                st.Mag = pmag
                traj.Structs.append(st)
                traj.ScfNum.append(scfn)
            elif state == 2:
                ll = [x for x in line.split(" ") if len(x) != 0][:6]
                if len(ll) != 6:
                    st = [0, 13, 26, 39, 56, 70, 84]
                    ll = [line_ori[a:b] for a, b in zip(st, st[1:])]
                ll = [float(l) if '*' not in l else 0.0 for l in ll]
                kk.append(ll[:3])
                kf.append(ll[3:6])
            elif "energy(sigma->0)" in line and state == 3:
                ll = [x.strip() for x in line.split("=") if len(x) != 0]
                state = 0
                st.Energy = float(ll[2]) / httoev
            if "reached required accuracy" in line:
                traj.RelaxConv = True
            elif "EDIFF  " in line:
                ediff = float(line.split()[2])
            elif 'total energy-change' in line:
                if 'MIXING' in line:
                    continue
                sx = line.split(':')
                if "**********" in line:
                    traj.ScfConv = False
                    continue
                a = float(sx[1].split('(')[0])
                b = sx[1].split('(')[1].replace(')', '')
                if 'e' not in b.lower():
                    # replace last occurrence of - (assumed exponent) with -e
                    bsplit = b.split('-')
                    bsplit[-1] = 'e' + bsplit[-1]
                    b = '-'.join(bsplit).replace('-e', 'e-')
                b = float(b)
                if [abs(a), abs(b)] < [ediff, ediff]:
                    traj.ScfConv = True
                else:
                    traj.ScfConv = False
        for it, t in enumerate(traj.Structs):
            t.Id = it
        return traj


class FreqI(object):
    freq = 0.0
    disps = []
    disp2 = []

# avoid overwritting


def new_file_name(x):
    i = 0
    y = x + '.' + str(i)
    while os.path.isfile(y):
        i += 1
        y = x + '.' + str(i)
    return y


def num_list(l):
    kl = []
    for il in sorted(l):
        if len(kl) != 0:
            if il == kl[-1][1] + 1:
                kl[-1][1] = il
                continue
        kl.append([il, il])
    return ",".join([("%d" % i[0]) if i[0] == i[1] else ("%d-%d" % tuple(i)) for i in kl])


def read_contcar(fn):
    TOLX = 1E-6
    f = open(fn, 'r').readlines()
    ff = [l.strip() for l in f]
    cell = [[float(g) for g in f.split(' ') if len(g) != 0] for f in ff[2:5]]
    cell = np.array(cell)
    if np.linalg.norm(np.array([cell[0, 1], cell[1, 0]])) < TOLX:
        dcell = np.diag(cell)
        print (' CELL  = %15.7f%15.7f%15.7f' % tuple(dcell))
    else:
        dcell = np.array(
            [cell[0, 0], cell[0, 1], cell[1, 0], cell[1, 1], cell[2, 2]])
        print (
            ' CELL  = [%15.7f, %15.7f] [%15.7f, %15.7f] %15.7f' % tuple(dcell))
    elems = np.array([g for g in ff[5].split(' ') if len(g) != 0])
    elemsn = np.array([int(g) for g in ff[6].split(' ') if len(g) != 0])
    print (' ELEMS = ' + " ".join(["%s %d" % (i, j)
                                   for i, j in zip(elems, elemsn)]))
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
    coords = np.array([[float(g) if g not in ['T', 'F'] else 0.0 for g in f.split(' ')
                        if len(g) != 0] for f in ff[lx:lx + elemsn.sum()]])
    rf = np.array([["relax" if g == 'T' else "fix" for g in f.split(' ') if len(g) != 0]
                   for f in ff[lx:lx + elemsn.sum()]])
    fxl = [iif for iif, f in enumerate(rf) if len(f) >= 4 and f[5] == "fix"]
    lf = num_list(fxl)
    print (' FIX   = %s' % lf)
    if gf:
        coords = coords[:, 0:3].dot(cell)
    else:
        coords = coords[:, 0:3]
    assert len(gelem) == len(coords)
    return [elemsn.sum(), np.array(gelem), coords], dcell, fxl

def proc_stat(pid, isneb):
    pp = psutil.Process(pid)
    with open("proc-stat.txt", "a") as g:
        g.write("# pid = %d, host = %s\n" % (pid, socket.gethostname()))
        g.flush()
        cnt = 0
        psmem = 0
        mcnt = 0
        while pp.is_running():
            ps = [pp] + pp.children(True)
            rpcnt, cpu, cps, cpp, smem = 0, 0.0, 0.0, 0.0, 0
            for p in ps:
                if p.status() == "running":
                    rpcnt += 1
                with p.oneshot():
                    cpu += p.cpu_times().user
                    cps += p.cpu_times().system
                    cpp += p.cpu_percent()
                    smem += p.memory_info().rss / 1024
            g.write("%03d %s %15.0f %15.0f %10.2f %10d KB\n" % (rpcnt, time.strftime("%y-%m-%d %H:%M:%S"), cpu, 
                cps, cpp, smem))
            g.flush()
            time.sleep(30)
            cnt += 1
            if psmem == smem:
                mcnt += 1
            else:
                psmem = smem
                mcnt = 0
            if cnt > 10 and (rpcnt == 0 or (mcnt > 2 and not isneb)):
                g.write("PROC ERROR !!!\n")
                g.flush()
                print ("PROC ERROR!!!")
                # for p in ps:
                #     p.kill()
                # time.sleep(2)
                # psutil.Process().kill()
                return
            if smem > 30 * 1024 ** 2:
                g.write("MEMORY ERROR !!!\n")
                g.flush()
                print("MEMORY ERROR !!!\n")
                for p in ps:
                    p.kill()
                time.sleep(2)
                psutil.Processes().kill()

class RunVASP(object):
    VaspDIR = ''
    VaspPPDIR = ''
    Funcs = ['pw91', 'lda', 'pbe', 'tpss', 'rpbe', 'pbesol', 'rtpss',
             'm06l', 'mbj', 'pbe0', 'tpssh', 'hse03', 'hse06']
    InCar = {
        "SYSTEM": "UNK", "ISTART": 0,
        "ISMEAR": 0, "SIGMA": 0.1, "NELM": 200,
        "ISPIN": 2, "LWAVE": ".FALSE.", "LCHARG": ".FALSE.",
        "EDIFF": 1e-6, "ALGO": "Fast"
    }

    ExInCar = ["ENCUT", "PREC", "EDIFFG", "POTIM", "IDIPOL", "LDIPOL",
               "IBRION", "NSW", "ICHARG", "NUPDOWN", "LSOL", "EB_K",
               "LVDW", "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "IVDW",
               "ISIF", "ISYM", "LREAL", "IOPT", "SPRING", "LASPH", "LMAXMIX",
               "LCLIMB", "EFIELD", "LORBIT", "AMIX", "BMIX", "NBANDS", "NPAR",
               "ICHAIN", "TEBEG", "TEEND", "SMASS", "KBLOCK", "NBLOCK",
               "IWAVPR", "DIPOL"]

    def prepare(self, ip):
        self.pstep = 0
        if ip.Multiplicity != 0 and "NUPDOWN" not in self.InCar:
            self.InCar["NUPDOWN"] = ip.Multiplicity - 1

        if "smass" in ip.Opts:
            self.md = True
        else:
            self.md = False
        
        if "continue" in ip.Opts:
            if "istart" not in ip.Opts:
                self.InCar["ISTART"] = 1
            if os.path.isfile("./final.xyz"):
                copyfile("./final.xyz", "./coord.xyz")
                cmtl = open("./final.xyz", "r").readlines()[1]
                cmtl = cmtl.split("=")[0].strip().split(":")[-1]
                self.pstep = int(cmtl)
                if ip.Opt != 0:
                    ip.Opt += 1
                if ip.NebOpt != 0:
                    ip.NebOpt += 1
            else:
                raise RuntimeError("Cannot continue without final.xyz!")
        else:
            if os.path.isfile('../' + ip.Coords):
                copyfile('../' + ip.Coords, "./coord.xyz")
            else:
                raise RuntimeError("Cannot start without coord file!")

        for s in ip.Scf:
            if s[0].lower() == 'iter':
                self.InCar["NELM"] = s[1]

        if 'freq' in ip.Opts:
            if "ibrion" not in ip.Opts:
                self.InCar["IBRION"] = 5
            if "nsw" not in ip.Opts:
                self.InCar["NSW"] = 1
            # self.InCar["ALGO"] = "Very_Fast"
        elif ip.Opt <= 1 and ip.NebOpt <= 0:
            if "ibrion" not in ip.Opts:
                self.InCar["IBRION"] = -1
            if "nsw" not in ip.Opts:
                self.InCar["NSW"] = 0
            if 'bader' in ip.Opts:
                self.keep_chg = False
                if "lcharg" in ip.Opts:
                    self.keep_chg = ip.Opts["lcharg"] in [ '.TRUE.', 'T' ]
                self.InCar["LCHARG"] = '.TRUE.'
                self.InCar["LAECHG"] = '.TRUE.'
                # elf.InCar["PREC"] = 'Accurate'
        elif 'dimer' in ip.Opts:
            if "ibrion" not in ip.Opts:
                self.InCar["IBRION"] = 3
            if "iopt" not in ip.Opts:
                self.InCar["IOPT"] = 2
            self.InCar["POTIM"] = 0
            self.InCar["ICHAIN"] = 2
            if "nsw" not in ip.Opts:
                self.InCar["NSW"] = ip.Opt
        elif ip.NebOpt > 0:
            if "ibrion" not in ip.Opts:
                self.InCar["IBRION"] = 1
            if "nsw" not in ip.Opts:
                self.InCar["NSW"] = ip.NebOpt
        else:
            if "ibrion" not in ip.Opts:
                self.InCar["IBRION"] = 2  # 1 quasi-newton, 2 CG
            if "nsw" not in ip.Opts:
                self.InCar["NSW"] = ip.Opt

        for k in self.InCar.keys() + self.ExInCar:
            if k.lower() in ip.Opts:
                self.InCar[k] = ip.Opts[k.lower()]

        if ":" in ip.Opts["cell"]:
            self.Cell = np.array([float(x)
                                  for x in ip.Opts["cell"].split(":")])
        else:
            self.Cell = np.array([float(ip.Opts["cell"])] * 3)

        self.ExeVersion = "gam"
        if "std" in ip.Opts:
            self.ExeVersion = "std"
        if "kgrid" in ip.Opts:
            self.KGrid = np.array([int(x)
                                   for x in ip.Opts["kgrid"].split(":")])
            self.ExeVersion = "std"
        else:
            self.KGrid = np.array([1, 1, 1])
        if "gpu" in ip.Opts:
            self.ExeVersion = "gpu"
            self.InCar["LREAL"] = "A"

        self.KMethod = "Monkhorst Park"
        if "kmethod" in ip.Opts:
            if ip.Opts["kmethod"].startswith("M"):
                self.KMethod = "Monkhorst Park"
            elif ip.Opts["kmethod"].startswith("G"):
                self.KMethod = "Gamma"
            else:
                raise RuntimeError("Unknown KPoint Method!")

        fu = open("./coord.xyz", 'r').readlines()
        gxf = [g for g in [x.strip() for x in fu][1].split(' ') if len(g) != 0]
        if "SURF" in gxf:
            self.surfnum = int(gxf[gxf.index("SURF") + 2])
        else:
            self.surfnum = 0
        nf = int(fu[0].strip())
        f = [x.strip() for x in fu][2:2 + nf]
        fx = [[i for i in x.split(' ') if len(i) != 0] for x in f]
        elems = [[ix, x[0]] for ix, x in enumerate(fx)]
        elems.sort(key=lambda x: x[1])
        self.Elems = elems
        self.UElems = list(set([x[1] for x in self.Elems]))
        self.UElems.sort()
        f = open('sort.dat', 'w')
        f.write(' '.join([str(x[0]) for x in elems]) + '\n')
        f.write(' '.join([x for x in self.UElems]) + '\n')
        f.close()

        if 'dipol' in ip.Opts:
            self.InCar["DIPOL"] = ip.Opts.replace(':', ' ')
        elif 'LDIPOL' in self.InCar:
            assert self.surfnum != 0
            surf_coords = np.array([[float(i) for i in l[1:4]] for l in fx[:self.surfnum]])
            surf_center = surf_coords.mean(axis=0)
            if len(self.Cell) == 3:
                cellmat = np.diag(self.Cell)
            else:
                cellmat = np.array([[self.Cell[0], self.Cell[1], 0.0],
                    [self.Cell[2], self.Cell[3], 0], [0, 0, self.Cell[4]]])
            surf_direct = surf_center.dot(np.linalg.inv(cellmat))
            self.InCar["DIPOL"] = "%10.6f %10.6f %10.6f" % tuple(surf_direct)

        if 'forces' in ip.Opts:
            fca, fcb = map(int, ip.Opts["forces"].split(":"))
            if fca != -1:
                fci, fcd = fcb / 2, fcb % 2
                fdd = float(ip.Opts["forceslen"])
                fdd = fdd if fcd == 0 else -fdd
                if "fix" in ip.Opts and ip.Opts["fix"] != "-1":
                    pp = solve_list(ip.Opts["fix"])
                else:
                    pp = []
                assert fca not in pp
                fpp = np.array(map(float, fx[fca][1:4]))
                fpp[fci] += fdd
                fx[fca][1:4] = map(str, fpp)
                self.InCar["LWAVE"] = "F"
                self.InCar["LCHARG"] = "F"
                self.InCar["ICHARG"] = "1"
                self.save_wave = False
            else:
                self.InCar["LWAVE"] = "T"
                self.InCar["LCHARG"] = "T"
                self.save_wave = True

        if ip.NebOpt > 0:
            self.OrigStr = []
            self.Coords = []
            fur = [x.strip() for x in fu]
            ifur = 0
            while ifur < len(fur):
                fl = [g for g in fur[1 + ifur].split(" ") if len(g) != 0]
                f = fur[2 + ifur:2 + nf + ifur]
                fx = [[i for i in x.split(' ') if len(i) != 0] for x in f]
                coords = [map(float, fx[ex[0]][1:4]) for ex in elems]
                self.Coords.append(np.array(coords))
                xstr = Struct()
                xstr.Coords = self.Coords[-1]
                if len(fl) >= 3:
                    xstr.Energy = float(fl[2])
                if fl[0].startswith("MAG"):
                    xstr.Mag = float(fl[0].split("~")[0].split(":")[2])
                self.OrigStr.append(xstr)
                ifur += 2 + nf
            self.images = len(self.Coords)
            assert self.images > 2
            self.InCar["IMAGES"] = self.images - 2
            self.DModes = None
        else:
            coords = [map(float, fx[ex[0]][1:4]) for ex in elems]
            self.Coords = np.array(coords)
            if 'dimer' in ip.Opts or 'dimerpre' in ip.Opts:
                dmodes = [map(float, fx[ex[0]][4:7]) for ex in elems]
                self.DModes = np.array(dmodes)
            else:
                self.DModes = None

        cpucnt = psutil.cpu_count()
        if 'freq' not in ip.Opts and "gpu" not in ip.Opts and "NPAR" not in self.InCar:
            ipp = ip.Para if ip.NebOpt == 0 else ip.Para / self.InCar["IMAGES"]
            ca, cb = 1, ipp
            for ca in range(2, ipp):
                if ipp % ca == 0 and np.abs(ca - cb) > np.abs(ipp / ca - ca):
                    ca, cb = ca, ipp / ca
            if ca > cb:
                ca, cb = cb, ca
            cc = ca
            self.InCar["NPAR"] = ipp / cc if ipp % cc == 0 else ipp
            # if self.InCar["NPAR"] > cpucnt and ipp % cpucnt == 0:
            #     self.InCar["NPAR"] = cpucnt
        
        if 'randmix' in ip.Opts:
            self.InCar["AMIX"] = "%.3f" % (0.2 + np.random.random() * 0.2 - 0.1)
            self.InCar["BMIX"] = "%.3f" % (2.0 + np.random.random() * 0.3 - 0.15)
        
        # if ip.Para >= 256:
        #     self.InCar["LPLANE"] = 'F'

        # xcm = np.array([np.array(coords).mean(axis=0)])
        # self.Coords = np.array(coords) - xcm

        exx = [ep + [ipx] for ipx, ep in enumerate(self.Elems)]
        exx.sort(key=lambda x: x[0])
        self.Fixed = [False] * len(elems)
        hasfix = False
        if "fixz" in ip.Opts:
            fixz = float(ip.Opts["fixz"])
            for p in range(0, len(elems)):
                if self.Coords[p][2] < fixz:
                    self.Fixed[p] = True
                    hasfix = True
        elif "fix" in ip.Opts and ip.Opts["fix"] != "-1":
            pp = solve_list(ip.Opts["fix"])
            for p in pp:
                if isinstance(p, int):
                    self.Fixed[exx[p][2]] = True
                elif p in self.UElems:
                    for g in exx:
                        if g[1] == p:
                            self.Fixed[g[2]] = True
                else:
                    raise RuntimeError("Unknown Element: " % p)
                hasfix = True
        self.SpFixed = [False] * len(elems)
        if "spfixxy" in ip.Opts:
            pp = solve_list(ip.Opts["spfixxy"])
            for p in pp:
                self.SpFixed[exx[p][2]] = True
        if hasfix:
            fixnum = len([x for x in self.Fixed if x])

        if ip.Method not in self.Funcs:
            raise RuntimeError('Functional %s is not valid!' % ip.Method)
        elif ip.Method == 'pw91':
            pp = self.VaspPPDIR + "/potpaw_GGA/"
        elif ip.Method == 'lda':
            pp = self.VaspPPDIR + "/potpaw/"
        else:
            pp = self.VaspPPDIR + "/potpaw_PBE/"
            if ip.Method == 'pbe':
                pass
            elif ip.Method == 'rpbe':
                self.InCar["GGA"] = 'RP'
            elif ip.Method == 'pbesol':
                self.InCar["GGA"] = 'PS'
            elif ip.Method == "pbe0":
                self.InCar["LHFCALC"] = '.TRUE.'
            elif ip.Method == "tpssh":
                self.InCar["LHFCALC"] = '.TRUE.'
                self.InCar["AEXX"] = '0.10'
                self.InCar["METAGGA"] = 'TPSS'
            elif ip.Method == "hse03":
                self.InCar["LHFCALC"] = '.TRUE.'
                self.InCar["HFSCREEN"] = '0.3'
            elif ip.Method == "hse06":
                self.InCar["LHFCALC"] = '.TRUE.'
                self.InCar["HFSCREEN"] = '0.2'
            elif ip.Method in ['tpss', 'rtpss', 'm06l', 'mbj']:
                self.InCar["METAGGA"] = ip.Method.upper()
            else:
                raise RuntimeError('Unknown functional %s!' % ip.Method)

        self.PotCar = []
        self.ElemNum = []
        self.VElec = {}
        self.Elec = 0
        for u in self.UElems:
            if os.path.exists(pp + u):
                self.PotCar.append(open(pp + u + "/POTCAR", 'r').read())
                self.ElemNum.append(len([x for x in self.Elems if x[1] == u]))
                nelec = int(float(self.PotCar[-1].split('\n')[1]) + 1E-10)
                self.VElec[u] = nelec
                self.Elec += self.ElemNum[-1] * nelec

        self.InCar["NELECT"] = self.Elec - ip.Charge

        self.InCar = {ix: str(x).replace(",", " ")
                      for ix, x in self.InCar.items()}

        if "nocenter" not in ip.Opts:
            self.Coords = np.array(self.Coords) + np.array([self.Cell / 2])

        self.PosCar = [ip.Title, "%18.10f" % 1.0]
        if len(self.Cell) == 3:
            self.PosCar += ["%18.10f%18.10f%18.10f" % (self.Cell[0], 0.0, 0.0),
                            "%18.10f%18.10f%18.10f" % (0.0, self.Cell[1], 0.0),
                            "%18.10f%18.10f%18.10f" % (0.0, 0.0, self.Cell[2])]
        else:
            self.PosCar += ["%18.10f%18.10f%18.10f" % (self.Cell[0], self.Cell[1], 0.0),
                            "%18.10f%18.10f%18.10f" % (
                                self.Cell[2], self.Cell[3], 0.0),
                            "%18.10f%18.10f%18.10f" % (0.0, 0.0, self.Cell[4])]
        self.PosCar += [("%7s" * len(self.UElems)) % tuple(self.UElems),
                        ("%7d" * len(self.ElemNum)) % tuple(self.ElemNum)]

        if hasfix:
            self.PosCar.append("Selective dynamics")
        self.PosCar.append("Cartesian")

        if ip.NebOpt > 0:
            pos_head = list(self.PosCar)
            self.PosCar = []
            for coo in self.Coords:
                pos_car = list(pos_head)
                for c, k in zip(coo, self.Fixed):
                    cck = "%18.10f%18.10f%18.10f" % tuple(c)
                    if hasfix:
                        if k:
                            cck += " F F F"
                        elif "xyfix" in ip.Opts:
                            cck += " F F T"
                        else:
                            cck += " T T T"
                    pos_car.append(cck)
                self.PosCar.append(pos_car)
        else:
            # spfixxy: fix xy because space group requirements
            # if fix is true, spfixxy is false, set to T T F
            # if fix is true, spfixxy is true, set to F F F
            # if spfixxy is true set to F F T
            # if spfixxy is false set to T T T
            if "spfixxy" in ip.Opts:
                for c, k, spk in zip(self.Coords, self.Fixed, self.SpFixed):
                    cck = "%18.10f%18.10f%18.10f" % tuple(c)
                    if k and not spk:
                        cck += " T T F"
                    elif k and spk:
                        cck += " F F F"
                    elif not k and spk:
                        cck += " F F T"
                    elif "xyfix" in ip.Opts:
                        cck += " F F T"
                    else:
                        cck += " T T T"
                    self.PosCar.append(cck)
            else:
                for c, k in zip(self.Coords, self.Fixed):
                    cck = "%18.10f%18.10f%18.10f" % tuple(c)
                    if hasfix:
                        if k and "xyrelax" in ip.Opts:
                            cck += " T T F"
                        elif k:
                            cck += " F F F"
                        elif "xyfix" in ip.Opts:
                            cck += " F F T"
                        else:
                            cck += " T T T"
                    self.PosCar.append(cck)

        self.KPoints = [
            "K-Points", "0", self.KMethod, "  ".join(
                map(str, self.KGrid)), "0  0  0"
        ]

        if "kpoints" in ip.Opts and os.path.isfile(ip.Opts["kpoints"]):
            copyfile(ip.Opts["kpoints"], "KPOINTS")
        else:
            open("KPOINTS", "w").write("\n".join(self.KPoints))

        if "potcar" in ip.Opts and os.path.isfile(ip.Opts["potcar"]):
            copyfile(ip.Opts["potcar"], "POTCAR")
        else:
            open("POTCAR", "w").write("".join(self.PotCar))

        if "continue" in ip.Opts and self.md:
            copyfile("CONTCAR.prev", "POSCAR")
        elif ip.NebOpt > 0:
            for ipos, pos in enumerate(self.PosCar):
                ip_dir = "%02d" % ipos
                if not os.path.exists(ip_dir):
                    os.mkdir(ip_dir)
                open("%s/POSCAR" % ip_dir, "w").write("\n".join(pos))
        else:
            open("POSCAR", "w").write("\n".join(self.PosCar))

        open("INCAR", "w").write("\n".join([
            "%s = %s" % (ix, x) for ix, x in self.InCar.items()]))
        
        if self.DModes is not None:
            self.ModeCar = ["%18.10f%18.10f%18.10f" % tuple(c) for c in self.DModes]
            open("MODECAR", "w").write("\n".join(self.ModeCar))

    def __init__(self, ip):
        if "VASPHOME" in os.environ and len(os.environ["VASPHOME"]) != 0:
            self.VaspDIR = os.environ["VASPHOME"]
        else:
            raise RuntimeError("VASPHOME env variable must be set!")

        if "VASP_PP_PATH" in os.environ and len(os.environ["VASP_PP_PATH"]) != 0:
            self.VaspPPDIR = os.environ["VASP_PP_PATH"]
        else:
            raise RuntimeError("VASP_PP_PATH env variable must be set!")

        if len(ip.Check) != 0:
            self.sub_dir = ip.Check.replace(" ", "_")

        if not os.path.exists(self.sub_dir):
            os.mkdir(self.sub_dir)
        os.chdir(self.sub_dir)

        self.InCar["SYSTEM"] = ip.Title
        self.InCar["ISTART"] = 0

        self.prepare(ip)
        surfnum = self.surfnum
        pstep = self.pstep

        if "norun" not in ip.Opts:
            while True:
                if "post" not in ip.Opts:
                    if os.path.isfile("OUTCAR"):
                        os.rename("OUTCAR", new_file_name("OUTCAR"))
                    if os.path.isfile("OSZICAR"):
                        os.rename("OSZICAR", new_file_name("OSZICAR"))
                    if ip.NebOpt > 0:
                        for ipos, pos in enumerate(self.PosCar):
                            ip_dir = "%02d" % ipos
                            outf, oszf = ip_dir + "/OUTCAR", ip_dir + "/OSZICAR"
                            if os.path.exists(ip_dir):
                                if os.path.isfile(outf):
                                    os.rename(outf, new_file_name(outf))
                                if os.path.isfile(oszf):
                                    os.rename(oszf, new_file_name(oszf))
                
                prog_name = "vasp_" + self.ExeVersion
                if ip.Para == 1 and ip.RunType != "mpirunst":
                    cmd = self.VaspDIR + "/" + prog_name
                else:
                    if ip.RunType == "aprun":
                        cmd = ("aprun -n %d " % ip.Para) + \
                            self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "ibrun":
                        cmd = ("ibrun -np %d " % ip.Para) + \
                            self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "hydra":
                        cmd = ("mpiexec.hydra -np %d " % ip.Para) + \
                            self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "srun":
                        cmd = ("srun --mpi=pmi2 --wait=60 --kill-on-bad-exit -n %d -N %d "
                               % (ip.Para, ip.Para / 16)) + self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "sruncori":
                        cmd = ("srun -n %d -c %d " % (ip.Para, 1)) + self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "mpirunst":
                        cmd = ("mpirun -n %d -ppn %s -hosts %s "
                            % (ip.Para, ip.mirablocks[0], ip.mirablocks[1].replace(":", ","))) + self.VaspDIR + "/" + prog_name
                    elif ip.RunType == "mira":
                        cmd = ("runjob --block %s -n %d -p 16 --corner %s --shape %s : "
                               % (ip.mirablocks[0], ip.Para, ip.mirablocks[1],
                                  ip.mirablocks[2])) + self.VaspDIR + "/" + prog_name
                    else:
                        cmd = ("mpirun -n %d " % ip.Para) + \
                            self.VaspDIR + "/" + prog_name
                
                if "pre" in ip.Opts:
                    open("../RUNDIR", "w").write(self.sub_dir)
                    open("../PROGRAM", "w").write(self.VaspDIR + "/" + prog_name)
                    return
                if "post" in ip.Opts:
                    cmd = "date"
                self.shell = Popen(shlex.split(cmd), shell=False, universal_newlines=True, bufsize=1,
                                   stdout=PIPE, stderr=STDOUT)
                if "post" not in ip.Opts:
                    thread.start_new_thread(proc_stat, (self.shell.pid, ip.NebOpt > 0))
                need_fix = False
                while self.shell.poll() is None:
                    line = self.shell.stdout.readline().rstrip()
                    if "had an illegal value" in line or "segmentation fault occurred" in line or \
                        "*** glibc detected ***" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        if self.InCar["ALGO"] == "Fast":
                            need_fix = True
                            self.InCar["ALGO"] = "Normal"
                        if ip.RunType == "mira":
                            time.sleep(4.0)
                        break
                    elif "Error EDDDAV: Call to ZHEGV failed" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        if os.path.isfile("CHGCAR") and ip.NebOpt <= 0 and "freq" not in ip.Opts:
                            fb = open("OUTCAR", "r").readlines()
                            traj = StructTraj.read_outcar(
                                fb, "nocenter" in ip.Opts, self.Cell)
                            if len(traj.Structs) != 0 and os.path.isfile("CONTCAR") and self.InCar["ALGO"] != "All":
                                traj.write_structs("trajectory.xyz", self.Elems, self.pstep, surfnum,
                                                   "continue" in ip.Opts, False)
                                traj.write_structs(
                                    "final.xyz", self.Elems, self.pstep, surfnum, False, True)
                                traj.write_scfnum(
                                    "scf_num.txt", "continue" in ip.Opts)
                                ip.Opts["continue"] = ""
                                os.remove("CHGCAR")
                                copyfile("CONTCAR", "POSCAR")
                                self.InCar["ALGO"] = "All"
                                need_fix = True
                            break
                    elif "Fatal error in MPI_Allreduce" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        print ("ired error")
                        break
                    elif "internal error, the gradient is not orthogonal" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        print ("ired error")
                        break
                    elif "very serious problems" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        print ("ired error")
                        break
                    elif "call to ZHEGV failed" in line:
                        self.shell.kill()
                        print ("ERROR DETECTED!! >>%s<<" % line)
                        print ("timeout error")
                        break
                    if line != '':
                        print (line)
                if need_fix:
                    open("INCAR", "w").write("\n".join([
                        "%s = %s" % (ix, x) for ix, x in self.InCar.items()]))
                    continue
                else:
                    break

            if self.shell.returncode == 0 and "bader" in ip.Opts:
                if "VTSTHOME" in os.environ and len(os.environ["VTSTHOME"]) != 0:
                    chgsum_prog = os.environ["VTSTHOME"] + "/chgsum.pl"
                else:
                    chgsum_prog = self.VaspDIR + "/../../vtstscripts-927/chgsum.pl"
                if "BADERHOME" in os.environ and len(os.environ["BADERHOME"]) != 0:
                    bader_prog = os.environ["BADERHOME"] + "/bader"
                else:
                    bader_prog = self.VaspDIR + "/../../bader/bader"

                if not os.path.isfile(bader_prog):
                    raise RuntimeError(
                        'No bader program found! %s' % bader_prog)
                if not os.path.isfile(chgsum_prog):
                    raise RuntimeError(
                        'No vtst program found! %s' % chgsum_prog)

                if "post" in ip.Opts:
                    prefix = ""
                elif ip.RunType == "aprun":
                    prefix = "aprun -n 1 "
                elif ip.RunType == "ibrun":
                    prefix = "ibrun -np 1 "
                elif ip.RunType == "hydra":
                    prefix = "mpiexec.hydra -np 1 "
                elif ip.RunType == "srun":
                    prefix = "srun --mpi=pmi2 --wait=60 --kill-on-bad-exit -n 1 -N 1 "
                elif ip.RunType == "sruncori":
                    prefix = "srun -n 1 -c 1 "
                elif ip.RunType == "mpirunst":
                    # prefix = "mpirun -n 1 -ppn 1 -hosts %s " % ip.mirablocks[1].replace(":", ",")
                    prefix = ""
                elif ip.RunType == "mira":
                    # prefix = "runjob --block %s -n 1 -p 1 --corner %s --shape %s : " \
                    #     % tuple(ip.mirablocks)
                    prefix = ""
                else:
                    prefix = ""

                print ("Run Bader charge analysis ...")
                cmd = prefix + chgsum_prog + " AECCAR0 AECCAR2"
                self.shell = Popen(shlex.split(cmd), shell=False, universal_newlines=True,
                                   bufsize=1, stdout=PIPE, stderr=STDOUT)

                while self.shell.poll() is None:
                    line = self.shell.stdout.readline().rstrip()
                    print (line)

                cmd = prefix + bader_prog + " CHGCAR -ref CHGCAR_sum"
                self.shell = Popen(shlex.split(cmd), shell=False, universal_newlines=True,
                                   bufsize=1, stdout=PIPE, stderr=STDOUT)

                while self.shell.poll() is None:
                    line = self.shell.stdout.readline().rstrip()
                    print (line)

                arr = []
                for i in ['CHG', 'CHGCAR', 'CHGCAR_sum', 'AECCAR0', 'AECCAR1', 'AECCAR2']:
                    if self.keep_chg and i == 'CHGCAR':
                        continue
                    if os.path.isfile(i):
                        arr.append(i)
                        try:
                            os.remove(i)
                        except OSError:
                            pass
                print ("removed: " + " ".join(arr))
            
            if "post" not in ip.Opts:
                if self.shell.returncode == 0:
                    print ("Normal termination")
                else:
                    print ("Error termination")

        if (os.path.isfile("OUTCAR") or os.path.isfile("01/OUTCAR")) and self.shell.returncode == 0:
            if "newstep" in ip.Opts:
                for fl in ["trajectory.xyz", "scf_num.txt", "trajectory.zip", "scf_num.zip", "bader.json"]:
                    if os.path.isfile(fl):
                        os.rename(fl, new_file_name(fl))
            if "freq" in ip.Opts:
                f = open("OUTCAR", "r").readlines()
                self.Freqs = []
                state = 0
                for line in f:
                    line = line.strip()
                    if line.startswith("Eigenvectors and eigenvalues of"):
                        state = 1
                    elif (state == 1 or state == 4) and line.startswith("Finite differences"):
                        break
                    elif "=" in line and (state == 1 or state == 4):
                        ll = [x for x in line.split(" ") if len(x) != 0]
                        ff = FreqI()
                        self.Freqs.append(ff)
                        if "f/i" in line:
                            ff.freq = -float(ll[-4])
                        else:
                            ff.freq = float(ll[-4])
                        ff.disps = []
                        ff.disp2 = []
                        state = 2
                    elif "dx" in line and state == 2:
                        state = 3
                    elif state == 3 and len(line) == 0:
                        state = 4
                    elif state == 4 and len(line) == 0:
                        break
                    elif state == 3:
                        ll = [x for x in line.split(" ") if len(x) != 0]
                        ff.disps.append(map(float, ll[3:]))
                self.Freqs.sort(key=lambda x: x.freq)
                exx = [ep + [ipx] for ipx, ep in enumerate(self.Elems)]
                exx.sort(key=lambda x: x[0])
                for i, j, k in exx:
                    for f in self.Freqs:
                        f.disp2.append(f.disps[k])
                for f in self.Freqs:
                    f.disp2 = np.array(f.disp2)
                f = open("vibspectrum", "w")
                f.write("$vibrational spectrum\n")
                f.write("#  mode symmetry wave number\n")
                f.write("#                  cm**(-1) \n")
                for ix, fr in enumerate(self.Freqs):
                    f.write("%6d%10s%12.3f\n" % (ix + 1, 'a', fr.freq))
                f.write("$end\n")
                f.close()
                f = open("vib_normal_modes", "w")
                f.write("$vibrational normal modes\n")
                x = np.array([fr.disp2.reshape(fr.disp2.shape[0] *
                                               fr.disp2.shape[1]) for fr in self.Freqs]).T
                for ix, xx in enumerate(x):
                    for iy, yy in enumerate(range(0, len(xx), 5)):
                        dd = xx[yy:yy + 5]
                        f.write(("%4d%4d" + ("%15.10f" * len(dd)) + "\n") %
                                ((ix + 1, iy + 1, ) + tuple(dd)))
                f.write("$end\n")
                f.close()
                f = open("vibration.xyz", "w")
                for ix, fr in enumerate(self.Freqs):
                    f.write("%d\n" % len(self.Elems))
                    cmt = "FREQ:%d:%10.5f\n" % (ix + 1, fr.freq)
                    if surfnum != 0:
                        cmt += " SURF = %d" % surfnum
                    f.write(cmt + "\n")
                    for i, j, k in exx:
                        f.write("%7s%18.10f%18.10f%18.10f%18.10f%18.10f%18.10f\n"
                                % ((j, ) + tuple(self.Coords[k]) + tuple(fr.disps[k])))
                f.close()
            elif ip.NebOpt > 0:
                trajs = []
                for ii in range(1, self.images - 1):
                    ii_dir = "%02d" % ii
                    f = open("%s/OUTCAR" % ii_dir, "r").readlines()
                    traj = StructTraj.read_outcar(
                        f, "nocenter" in ip.Opts, self.Cell)
                    trajs.append(traj)
                zf = zipfile.ZipFile("trajectory.zip", "a" if "continue" in ip.Opts else "w",
                                     compression=zipfile.ZIP_DEFLATED)
                zff = zipfile.ZipFile("scf_num.zip", "a" if "continue" in ip.Opts else "w",
                                      compression=zipfile.ZIP_DEFLATED)
                for ist in range(0, len(trajs[0].Structs)):
                    xtrajfs = [t.Structs[ist] for t in trajs]
                    xscfs = [t.ScfNum[ist] for t in trajs]
                    self.OrigStr[0].Id = self.OrigStr[-1].Id = xtrajfs[0].Id
                    nebt = StructTraj()
                    nebt.Structs = [self.OrigStr[0]] + \
                        xtrajfs + [self.OrigStr[-1]]
                    nebt.ScfNum = [0] + xscfs + [0]
                    nebt.write_structs("%d.xyz" % (xtrajfs[0].Id + pstep), self.Elems, pstep,
                                       surfnum, False, False, zf=zf)
                    nebt.write_scfnum("%d.txt" %
                                      (xtrajfs[0].Id + pstep), False, zf=zff)
                    if ist == len(trajs[0].Structs) - 1:
                        nebt.write_structs(
                            "final.xyz", self.Elems, pstep, surfnum, False, False)
                zf.close()
                zff.close()

                if trajs[0].RelaxConv:
                    print ("NEB_CONVERGED")
                else:
                    print ("NEB DID NOT CONVERGE WITHIN")
                if not trajs[0].ScfConv:
                    print ("ERROR: your energy calculation did not converge")
            else:
                f = open("OUTCAR", "r").readlines()
                traj = StructTraj.read_outcar(
                    f, "nocenter" in ip.Opts, self.Cell)
                if 'dimer' in ip.Opts:
                    df = open("NEWMODECAR", "r").readlines()
                    mds = [[float(g) for g in f.strip().split(" ") if len(g) != 0] for f in df]
                    traj.Structs[-1].Forces = np.array(mds)
                elif 'dimerpre' in ip.Opts:
                    traj.Structs[-1].Forces = np.array(self.DModes)
                cont, _, _ = read_contcar('CONTCAR')
                traj.Structs[-1].Coords = cont[2]
                traj.write_structs("trajectory.xyz", self.Elems, pstep, surfnum,
                                   "continue" in ip.Opts and "forces" not in ip.Opts and ip.Opt > 1, False)
                traj.write_structs("final.xyz", self.Elems,
                                   pstep, surfnum, False, True)
                traj.write_scfnum("scf_num.txt", "continue" in ip.Opts)

                if traj.RelaxConv:
                    print ("GEO_OPT_CONVERGED")
                else:
                    print ("OPTIMIZATION DID NOT CONVERGE WITHIN")
                if not traj.ScfConv:
                    print ("ERROR: your energy calculation did not converge")
            if "forces" in ip.Opts:
                print ("Extract forces ...")
                if not os.path.isfile("vasprun.xml"):
                    print ("NO XML FILE!!")
                else:
                    with open("vasprun.xml", "r", 1024**2) as v:
                        ff = []
                        edt = -1
                        for pl in v:
                            if "</varray>" in pl:
                                edt = -1
                            elif 'varray name="forces"' in pl:
                                ff = []
                                edt = 0
                            elif edt != -1:
                                if not self.Fixed[edt]:
                                    ff.append(pl.split(">")[1].split("<")[0])
                                else:
                                    ff.append(None)
                                edt += 1
                        ff2 = []
                        exx = [ep + [ipx] for ipx, ep in enumerate(self.Elems)]
                        exx.sort(key=lambda x: x[0])
                        for _, _, k in exx:
                            if ff[k] is not None:
                                ff2.append(ff[k])
                        assert len(ff2) != 0
                        ff2.append("")
                        with open("forces.xyz", "w", 1024**2) as f:
                            f.write(("%d\n\n" % (len(ff2) - 1)) + "\n".join(ff2))
                if not self.save_wave:
                    arr = []
                    for i in ['WAVECAR', 'CHG', 'CHGCAR']:
                        if os.path.isfile(i):
                            arr.append(i)
                            try:
                                os.remove(i)
                            except OSError:
                                pass
                    print ("removed: " + " ".join(arr))

        else:
            if "norun" not in ip.Opts:
                print ("NO OUTCAR!!")

        if "ISIF" in self.InCar and self.InCar["ISIF"] == "3" and \
                os.path.isfile("CONTCAR") and self.shell.returncode == 0:
            f = open("CONTCAR", "r").readlines()
            fl = [[float(i) for i in x.strip().split(" ") if len(i) != 0]
                  for x in f[2:5]]
            fn = open("cell.txt", "w")
            if np.abs(fl[0][1]) < 1E-3 and np.abs(fl[1][0]) < 1E-3:
                fn.write("%15.10f %15.10f %15.10f\n" %
                         (fl[0][0], fl[1][1], fl[2][2]))
            else:
                fn.write("%15.10f %15.10f %15.10f %15.10f %15.10f\n" %
                         (fl[0][0], fl[0][1], fl[1][0], fl[1][1], fl[2][2]))
            fn.close()
        
        if self.md:
            copyfile("CONTCAR", "CONTCAR.prev")

        if "bader" in ip.Opts and os.path.isfile("ACF.dat") and self.shell.returncode == 0:
            f = open("ACF.dat", "r").readlines()
            fl = [[float(i) for i in x.strip().split(" ") if len(i) != 0]
                  for x in f[2:-4]]
            fl = np.array(fl)
            elecs_orig = fl[:, 4]
            exx = [ep + [ipx] for ipx, ep in enumerate(self.Elems)]
            exx.sort(key=lambda x: x[0])
            elecs = []
            charges = []  # minus # of electron
            for i, j, k in exx:
                elecs.append(elecs_orig[k])
                charges.append(self.VElec[j] - elecs_orig[k])
            nelecg = float(f[-1].strip().split(" ")[-1])
            vacv = float(f[-2].strip().split(" ")[-1])
            vacc = float(f[-3].strip().split(" ")[-1])
            baderd = {"bader_electrons": elecs, "bader_charges": charges,
                      "bader_vacuum_charge": vacc,
                      "bader_vacuum_volume": vacv, "bader_nelectrons": nelecg,
                      "bader_nelec_diff": nelecg - (self.Elec - ip.Charge)}
            with open('bader.json', 'w') as f:
                json.dump(baderd, f, indent=4)


def read_opts(argv):
    opts = {}
    for i in argv:
        if not i.startswith('--'):
            raise RuntimeError('Unknown argument: %s' % i)
        if '=' in i:
            a, b = i.split('=')
        else:
            a, b = i, ''
        opts[a[2:]] = b
    return opts


if __name__ == "__main__":
    if len(argv) > 1 and 'clean' in read_opts(argv[1:]):
        arr = []
        for i in ['CHG', 'CHGCAR', 'CHGCAR_sum', 'CONTCAR', 'DOSCAR', 'EIGENVAL',
                  'IBZKPT', 'OSZICAR', 'OUTCAR', 'PCDAT', 'vasprun.xml',
                  'WAVECAR', 'XDATCAR', 'PROCAR', 'LOCPOT', 'AECCAR0',
                  'AECCAR1', 'AECCAR2', 'REPORT']:
            if os.path.isfile(i):
                arr.append(i)
                try:
                    os.remove(i)
                except OSError:
                    pass
        print ("removed: " + " ".join(arr))
    elif len(argv) > 1 and 'example' in read_opts(argv[1:]):
        with open('relax.in', 'w') as f:
            f.write("% workers=unknown:32\n")
            f.write("% chk=NAME.chk\n")
            f.write("# PBE/ fix=0-0 cell=1:0:0:1:1 coords=coords.xyz scf(iter=300) " +
                "encut=400 algo=Fast ediff=1E-4 ediffg=1E-3 lwave=F nocenter sigma=0.1 " +
                "lreal=A kmethod=MP kgrid=3:3:1 opt(iter=300)\n\nNAME\n\n0 0\n\n")
    elif len(argv) > 1 and 'traj' in read_opts(argv[1:]):
        assert os.path.isfile('CONTCAR') and os.path.isfile('OUTCAR')
        c, dcell, fixl = read_contcar('CONTCAR')
        if os.path.isfile('sort.dat'):
            ielems = [int(x) for x in open('sort.dat', 'r').readline().split()]
            xelems = [[ix, x] for ix, x in zip(ielems, c[1])]
        else:
            xelems = [[ix, x] for ix, x in enumerate(c[1])]
        nct = np.abs(dcell[0] - dcell[1]) > 1E-6 or np.abs(dcell[0] - dcell[2]) > 1E-6
        traj = StructTraj.read_outcar(open('OUTCAR', 'r').readlines(), nct, dcell)
        traj.Structs[-1].Coords = c[2]
        traj.write_structs('trajectory.xyz', xelems, 0, 0, False, False)
    elif len(argv) > 1 and 'npar' in read_opts(argv[1:]):
        ipp = int(read_opts(argv[1:])['npar'])
        ca, cb = 1, ipp
        for ca in range(2, ipp):
            if ipp % ca == 0 and np.abs(ca - cb) > np.abs(ipp / ca - ca):
                ca, cb = ca, ipp / ca
        if ca > cb:
            ca, cb = cb, ca
        cc = ca
        print (ipp / cc if ipp % cc == 0 else ipp)
    else:
        ipx = QMInput()
        for rt in ['aprun', 'srun', 'ibrun', 'hydra', 'mira', 'mpirunst', 'sruncori']:
            if len(argv) > 1 and rt in read_opts(argv[1:]):
                ipx.RunType = rt
        if 'norun' in read_opts(argv[1:]):
            ipx.Opts['norun'] = ''
        if 'pre' in read_opts(argv[1:]):
            ipx.Opts['pre'] = ''
        if 'post' in read_opts(argv[1:]):
            ipx.Opts['post'] = ''
        vp = RunVASP(ipx)
