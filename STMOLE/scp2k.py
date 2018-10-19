#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import json
import copy
import numpy as np
import itertools
from shutil import copyfile
from sys import argv


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
    Basis_Set = []
    Potential = []
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
                    elif x[1] != '' and x[0].lower() == 'basis_set':
                        self.Basis_Set = self._opts(x[1])
                    elif x[1] != '' and x[0].lower() == 'potential':
                        self.Potential = self._opts(x[1])
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


def trans_cp2k_input(js, name=None):
    r = []
    pre = ""
    if name is not None:
        if "_" in js:
            r.append("&" + name.split(' ')[0] + " " + str(js["_"]))
        else:
            r.append("&" + name.split(' ')[0])
        pre = "  "
    for k, v in sorted(js.items(), key=lambda (k, v):
                       int(k.split()[-1]) if k.split()[-1].isdigit() else k):
        k = k.split(' ')[0]
        if k == '_':
            continue
        if isinstance(v, dict):
            r.extend([pre + x for x in trans_cp2k_input(v, k)])
        else:
            r.append(pre + k + " " + str(v))
    if name is not None:
        r.append("&END " + name.split(' ')[0])
    return r


class RunCP2K(object):
    CP2kDIR = ''
    CP2kDataDIR = ''
    STMoleDIR = ''
    BasicPart = ''

    Input = None

    def __init__(self, ip):
        if "CP2KHOME" in os.environ and len(os.environ["CP2KHOME"]) != 0:
            self.CP2kDIR = os.environ["CP2KHOME"]
        else:
            raise RuntimeError("CP2KHOME env variable must be set!")

        if "STMOLE_HOME" in os.environ and len(os.environ["STMOLE_HOME"]) != 0:
            self.STMoleDIR = os.environ["STMOLE_HOME"]
        else:
            raise RuntimeError("STMOLE_HOME env variable must be set!")

        if "CP2KDATA" in os.environ and len(os.environ["CP2KDATA"]) != 0:
            self.CP2kDataDIR = os.environ["CP2KDATA"]
        else:
            raise RuntimeError("CP2KDATA env variable must be set!")

        with open(os.environ['STMOLE_HOME'] + "/cp2k-template.json", "r") as f:
            self.BasicPart = json.load(f)

        if len(ip.Check) != 0:
            self.sub_dir = ip.Check.replace(" ", "_")

        if not os.path.exists(self.sub_dir):
            os.mkdir(self.sub_dir)
        os.chdir(self.sub_dir)

        self.Input = copy.deepcopy(self.BasicPart)
        self.Input["GLOBAL"]["PROJECT"] = ip.Title

        self.Input["FORCE_EVAL"]["DFT"]["BASIS_SET_FILE_NAME"] = self.CP2kDataDIR + "/BASIS_MOLOPT"
        self.Input["FORCE_EVAL"]["DFT"]["POTENTIAL_FILE_NAME"] = self.CP2kDataDIR + \
            "/GTH_POTENTIALS"

        self.prepare(ip)

    def prepare(self, ip):
        x = ["%14.8f" % float(x) for x in ip.Opts["cell"].split(":")]
        z = "%14.8f" % 0.0
        self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"] = {}
        if len(x) == 1:
            self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"]["ABC"] = " ".join(x * 3)
        elif len(x) == 3:
            self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"]["ABC"] = " ".join(x)
        elif len(x) == 5:
            self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"]["A"] = " ".join([
                                                                       x[0], x[1], z])
            self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"]["B"] = " ".join([
                                                                       x[2], x[3], z])
            self.Input["FORCE_EVAL"]["SUBSYS"]["CELL"]["C"] = " ".join([
                                                                       z, z, x[4]])
        else:
            assert False

        if "continue" in ip.Opts:
            assert False
        else:
            assert os.path.isfile('../' + ip.Coords)

        with open('../' + ip.Coords, 'r') as f:
            natoms = int(f.readline())
            f.readline()
            atoms = np.zeros((natoms, 3))
            elems = []
            for i, l in enumerate(f):
                l = l.split()
                elems.append(l[0])
                atoms[i] = np.array(l[1:4], dtype=float)

        if "fix" in ip.Opts and ip.Opts["fix"] != "-1":
            fix_list = solve_list(ip.Opts["fix"])

        self.Input["FORCE_EVAL"]["SUBSYS"]["COORD"] = {}
        elem_fix = {}
        for i in range(0, natoms):
            if i in fix_list:
                kk = elems[i] + "FIX %d" % i
                vv = ("%14.8f" * 3) % tuple(atoms[i])
                self.Input["FORCE_EVAL"]["SUBSYS"]["COORD"][kk] = vv
                elem_fix[elems[i] + "FIX"] = elems[i]
            else:
                kk = elems[i] + " %d" % i
                vv = ("%14.8f" * 3) % tuple(atoms[i])
                self.Input["FORCE_EVAL"]["SUBSYS"]["COORD"][kk] = vv

        uelems = [k for k, g in itertools.groupby(elems)]
        basis = "DZVP-MOLOPT-SR-GTH"
        potential = "GTH-PBE"
        for k, v in ip.Basis_Set:
            if v == '':
                basis = k
        for k, v in ip.Potential:
            if v == '':
                potential = k
        kind = {"BASIS_SET": basis, "POTENTIAL": potential}
        for k in self.Input["FORCE_EVAL"]["SUBSYS"]:
            if k.startswith("KIND"):
                del self.Input["FORCE_EVAL"]["SUBSYS"][k]
        for e in uelems:
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " +
                                               e] = copy.deepcopy(kind)
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " + e]["_"] = e
        for fe, e in elem_fix.items():
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " +
                                               fe] = copy.deepcopy(kind)
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " + fe]["_"] = fe
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " + fe]["ELEMENT"] = e

        for k, v in ip.Basis_Set:
            if v == '':
                continue
            bs = self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " +
                                                    k]["BASIS_SET"].split('-')
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " +
                                               k]["BASIS_SET"] = '-'.join([v] + bs[1:])
        for k, v in ip.Potential:
            if v == '':
                continue
            self.Input["FORCE_EVAL"]["SUBSYS"]["KIND " + k]["POTENTIAL"] = v

        self.Input["MOTION"]["CONSTRAINT"]["FIXED_ATOMS"] = {}
        if "fix" in ip.Opts and ip.Opts["fix"] != "-1":
            pp = fix_list
            g = []
            for p in pp:
                p += 1
                assert isinstance(p, int)
                if len(g) != 0 and g[-1][-1] == p - 1:
                    g[-1][-1] = p
                else:
                    g.append([p, p])
            for i, j in g:
                self.Input["MOTION"]["CONSTRAINT"]["FIXED_ATOMS"]["LIST " +
                                                                  str(i)] = str(i) + ".." + str(j)

        if "cutoff" in ip.Opts:
            self.Input["FORCE_EVAL"]["DFT"]["MGRID"]["CUTOFF"] = "[Ry] " + \
                ip.Opts["cutoff"]

        if "add_mos" in ip.Opts:
            self.Input["FORCE_EVAL"]["DFT"]["SCF"]["ADDED_MOS"] = int(
                ip.Opts["add_mos"])

        if "max_scf" in ip.Opts:
            self.Input["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"] = int(
                ip.Opts["max_scf"])
        if "eps_scf" in ip.Opts:
            self.Input["FORCE_EVAL"]["DFT"]["SCF"]["EPS_SCF"] = ip.Opts["eps_scf"]

        if "sp" in ip.Opts:
            self.Input["GLOBAL"]["RUN_TYPE"] = "ENERGY"
            del self.Input["MOTION"]
        else:
            self.Input["GLOBAL"]["RUN_TYPE"] = "GEO_OPT"
        
        if "kmethod" in ip.Opts and "kgrid" in ip.Opts:
            if ip.Opts["kmethod"].startswith("M"):
                self.Input["FORCE_EVAL"]["DFT"]["KPOINTS"] = {}
                self.Input["FORCE_EVAL"]["DFT"]["KPOINTS"]["SCHEME"] = "MONKHORST-PACK " + \
                    ip.Opts["kgrid"].replace(':', ' ')
        
        if "fftw_plan" in ip.Opts:
            self.Input["GLOBAL"]["FFTW_PLAN_TYPE"] = ip.Opts["fftw_plan"]
        
        if "elpa" in ip.Opts:
            self.Input["GLOBAL"]["PREFERRED_DIAG_LIBRARY"] = "ELPA"

        open("input.in", "w").write("\n".join(trans_cp2k_input(self.Input)))


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
        for i in ['????']:
            if os.path.isfile(i):
                arr.append(i)
                try:
                    os.remove(i)
                except OSError:
                    pass
        print ("removed: " + " ".join(arr))
    else:
        ipx = QMInput()
        for rt in ['aprun', 'srun', 'ibrun', 'hydra', 'mira', 'mpirunst']:
            if len(argv) > 1 and rt in read_opts(argv[1:]):
                ipx.RunType = rt
        vp = RunCP2K(ipx)
