#!/usr/bin/env python

import os
import sys
from subprocess import Popen, PIPE
import signal
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

class Alarm(Exception):
    pass


def alarm_handler(signum, frame):
    raise Alarm


class QMInput:
    Method = 'HF'
    Basis = 'def2-SV(P)'
    Ecp = ''
    Coords = ''
    Title = ''
    Charge = 0
    Multiplicity = 1
    Para = None
    Nodes = []
    Check = 'OUT'
    Opt = 0
    Skr = ''
    Scf = []
    Opts = {}
    Line = ''

    def __init__(self, doreduce=False):
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
                    self.Nodes = []
                    nodes = xs[1].strip().split(',')
                    if len(nodes) == 1 and len(nodes[0]) == 0:
                        self.Para = "SMP"
                    else:
                        self.Para = "MPI"
                    for n in nodes:
                        ns = n.strip().split(':')
                        xhost = os.popen('hostname').readline().strip() 
                        if ns[0] == 'unknown': ns[0] = xhost
                        if len(ns) == 1:
                            self.Nodes.append(ns[0])
                        else:
                            num = int(ns[1])
                            self.Nodes += [ns[0]] * num
                    if doreduce and len(self.Nodes) >= 2 and self.Nodes[0] == self.Nodes[1]:
                      self.Nodes = self.Nodes[1:]
                    if len(self.Nodes) == 1:
                        self.Para = None
                elif xsn == 'check' or xsn == 'chk':
                    self.Check = xs[1].strip()
                elif xsn == 'nproc' or xsn == 'nprocshared':
                    xs = int(xs[1].strip())
                    self.Nodes = [''] * xs
                elif xsn == 'skr':
                    xs = xs[1].strip()
                    self.Skr = xs
            if line[0] == '#':
                self.Line = line.strip()
                xs = line[1:].strip()
                xsr = self._opts(xs)
                for x in xsr:
                    if x[1] == '' and '/' in x[0]:
                        xl = x[0].split('/')
                        self.Method = xl[0].strip()
                        self.Basis = xl[1].strip()
                        if len(xl) == 3:
                            self.Ecp = xl[2].strip()
                    elif x[1] != '' and x[0].lower() == 'coords':
                        self.Coords = x[1].strip()
                    elif x[1] != '' and x[0].lower() == 'scf':
                        self.Scf = self._opts(x[1])
                    elif x[1] != '' and x[0].lower() == 'frag':
                        xgr = self._opts(x[1])
                        fp = []
                        for g in xgr:
                            fp.append([g[0].strip(), [h.strip() for h in g[1].split(',')]])
                        self.Opts[x[0].lower()] = fp
                    elif x[0].lower() == 'ri':
                        xgr = self._opts(x[1])
                        fp = {}
                        for g in xgr:
                            fp[g[0]] = g[1]
                        self.Opts[x[0].lower()] = fp
                    elif x[0].lower() == 'opt':
                        xgr = self._opts(x[1])
                        fp = {}
                        for g in xgr:
                            if g[0] == 'iter':
                                self.Opt = int(g[1])
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
            elif ' ' in line and len(line.split(' ')) == 2:
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


class RunTM:
    # TurboDIR = '/home/hcpc/Desktop/UCLA-2015/ProgramsS/TURBOMOLE'
    # TurboDIR = '/u/home/a/ana/TURBOMOLE'
    # TurboDIR = '/u/home/h/hczhai/project-ana/TURBOMOLE'
    # TurboDIR = '/p/home/ana/hczhai/program/TURBOMOLE'
    TurboDIR = ''
    old_version = False
    Funcs = ['slater-dirac-exchange', 's-vwn', 'vwn', 's-vwn_Gaussian',
             'pwlda', 'becke-exchange', 'b-lyp', 'b-vwn', 'lyp',
             'b-p', 'pbe', 'tpss', 'bh-lyp', 'b3-lyp', 'b3-lyp_Gaussian',
             'pbe0', 'tpssh', 'lhf', 'b97-d', 'b2-plyp']
    blank_quit = ['scf', 'ff']
    timeout = 10
    timeout_err = 'timeout error'
    shell = None
    state = ''
    errstate = ''
    stdline = []
    errline = []
    sub_dir = ''
    ctrl_file = "control"
    nodes_file = "nodes.txt"
    def_file = "define.txt"
    def_out = sys.stdout
    post_def = []
    is_ri = False

    tasks = {}

    def __init__(self, ip):
        if "TURBOMOLE_HOME" in os.environ and len(os.environ["TURBOMOLE_HOME"]) != 0:
            self.TurboDIR = os.environ["TURBOMOLE_HOME"]
        elif "TURBODIR" in os.environ and len(os.environ["TURBODIR"]) != 0:
            self.TurboDIR = os.environ["TURBODIR"]
        else:
            raise RuntimeError("Turbomole root directory must be set in TURBODIR or TURBOMOLE_HOME!")
        if len(ip.Check) != 0:
            self.sub_dir = ip.Check.replace(" ", "_")

        if not os.path.exists(self.sub_dir):
            os.mkdir(self.sub_dir)
        os.chdir(self.sub_dir)

        os.environ["TURBODIR"] = self.TurboDIR
        if ip.Para == "SMP" or ip.Para == "MPI":
            os.environ["PARA_ARCH"] = ip.Para
            os.environ["PARNODES"] = str(len(ip.Nodes))
        else:
            os.environ["PARA_ARCH"] = ""
        if ip.Para == "MPI":
            fnodes = open(self.nodes_file, mode='w')
            fnodes.write('\n'.join(ip.Nodes) + '\n')
            fnodes.flush()
            fnodes.close()
            os.environ["HOSTS_FILE"] = os.curdir + '/' + self.nodes_file
            if 'PBS_NODEFILE' in os.environ and not os.path.isfile(os.environ['PBS_NODEFILE']):
                os.environ['PBS_NODEFILE'] = os.environ["HOSTS_FILE"]
        os.environ['PATH'] = os.environ['TURBODIR'] + '/scripts:' + os.environ['PATH']
        xsysname = os.popen('sysname').readline().strip().replace("sgi", "unknown")
        os.environ['PATH'] = os.environ['TURBODIR'] + '/bin/' + xsysname + \
            ':' + os.environ['PATH']
        if 'LD_LIBRARY_PATH' not in os.environ: os.environ['LD_LIBRARY_PATH'] = ''
        os.environ['LD_LIBRARY_PATH'] = os.environ['TURBODIR'] + '/libso/' + \
            os.popen('sysname').readline().strip() + ':' + os.environ['LD_LIBRARY_PATH']
        if len(ip.Skr) != 0:
            os.environ['TURBOTMPDIR'] = ip.Skr
        if 'continue' not in ip.Opts.keys() and os.path.isfile(self.ctrl_file):
            os.rename(self.ctrl_file, self.ctrl_file + '.bak')
        self.def_out = open(self.def_file, mode='a')
        self.def_out.write(ip.Line + '\n')
        self.shell = Popen("define", shell=True, universal_newlines=True,
                           stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1)

        self.stdline.append(("IF YOU WANT TO READ DEFAULT-DATA FROM ANOTHER control-TYPE FILE,", 2, "control"))
        self.stdline.append(("INPUT TITLE OR", 1, "title"))
        self.stdline.append(("HIT >return< TO ACCEPT DEFAULT TITLE OR", 1, "title"))
        self.stdline.append(("TERMINATE MOLECULAR GEOMETRY SPECIFICATION", 5, "geometry"))
        self.stdline.append(("GOBACK=& (TO GEOMETRY MENU !)", 1, "atomic"))
        self.stdline.append(("THE COMMANDS  use  OR  eht  OR  *  OR q(uit) TERMINATE THIS MENU", 2, "occupation"))
        self.stdline.append(("LIST OF MO-SHELL INDICES (LIKE  1-5,7-8,11)", 2, "occ_num"))
        self.stdline.append(("GO BACK TO OCCUPATION/ORBITAL ASSIGNMENT MENU", 1, "general"))
        if self.old_version:
            self.stdline.append(("pdiag: PREDIAGONALIZATION", 1, "scf"))
        else:
            self.stdline.append(("soghf: SPIN ORBIT GENERALIZED SCF", 1, "scf"))
        self.stdline.append(("ENTER NEW VALUE FOR MAXIMUM NUMBER OF SCF-ITERATIONS", 2, "scf_iter"))
        self.stdline.append(("off: TO SWITCH off COUNTERPOISE CORRECTION", 1, "cp_on"))
        self.stdline.append(("on:   TO SWITCH on COUNTERPOISE CORRECTION", 1, "cp_off"))
        self.stdline.append(("on:   TO SWITCH ON  DFT", 1, "dft_off"))
        self.stdline.append(("off:  TO SWITCH OFF DFT", 1, "dft_on"))
        self.stdline.append(("on: TO SWITCH ON  RI", 1, "ri_off"))
        self.stdline.append(("off: TO SWITCH OFF RI", 1, "ri_on"))
        self.stdline.append(("old :  to switch DFT-D2 correction on", 1, "dsp"))
        self.stdline.append(("Please input fragment # of atom", 0, "fragment"))
        self.stdline.append(("Blank input starts calculation", 0, "ff"))
        self.stdline.append(("MULTIPOLE APPROXIMATION ALREADY SWITCHED ON", 1, "marij_on"))
        self.stdline.append(("truncated RI ALREADY SWITCHED ON", 1, "trunc_on"))
        self.stdline.append(("YOUR INPUT HAS TO BE AN INTEGER IN THE RANGE 4-9", 2, "scf_conv"))
        self.stdline.append(("3) none -- no DIIS is performed", 2, "scf_diis"))

        self.stdline.append(("DO YOU WANT TO CHANGE THE GEOMETRY DATA ?", 1, "ask:geometry"))
        self.stdline.append(("ATOMIC ATTRIBUTE DATA (BASES,CHARGES,MASSES,ECPS) HAVE BEEN", 3, "ask:atomic"))
        self.stdline.append(("MOLECULAR ORBITAL DATA (OCCUPATION NUMBERS,MOS) HAVE BEEN", 3, "ask:occupation"))

        self.stdline.append(("CONFIRM REMOVAL OF THESE ATOMS BY TYPING  y", 1, "y"))
        self.stdline.append(("DO YOU WANT THE DEFAULT PARAMETERS FOR THE EXTENDED HUECKEL CALCULATION ?", 2, "y"))
        self.stdline.append(("IF YOU WANT, I SKIP TO THE LIBRARY", 1, "y"))
        self.stdline.append(("ENTER THE MOLECULAR CHARGE", 1, str(ip.Charge)))
        self.stdline.append(("ENTER THE ATOMIC CHARGE", 1, str(ip.Charge)))
        self.stdline.append(("DO YOU ACCEPT THIS OCCUPATION ?", 1, "acc_occupation"))
        self.stdline.append(("DO YOU WANT THE DEFAULT OCCUPATION ASSIGNMENT FOR ATOMS ?", 1, "n"))
        # delete previous datagroups
        self.stdline.append(("LEFT OVER FROM PREVIOUS CALCULATIONS ?", 0, 'y' if 'clear' in ip.Opts.keys() else 'n'))
        self.stdline.append(("DO YOU REALLY WANT TO USE $forceapprox", 0, 'n' if 'clear' in ip.Opts.keys() else 'y'))
        self.stdline.append(("delete data group $", 0, "y"))
        self.stdline.append(("DO YOU REALLY WANT TO WRITE OUT NATURAL ORBITALS ?", 1, "n"))  # natural orbitals
        self.stdline.append(("TO CONTINUE, ENTER <return>", 0, ''))
        self.stdline.append(("Enter the number to change a value or <return> to accept all.", 0, ''))

        self.errline.append(("nick    - REPEAT NICKNAME INPUT", 8, "basis set name error"))
        self.errline.append(("CHOSEN OCCUPATION IMPOSSIBLE", 0, "multiplicity error"))
        self.errline.append(("NO MORE DATA AVAILABLE FOR EHT !!", 3, "eht error"))
        self.errline.append(("Please define more than one fragment!", 3, "frag error"))
        self.errline.append(("NO ATOMS, NO MOLECULE, NOTHING !", 3, "no atoms error"))
        self.errline.append(("WARNING! Lowest Eigenvalue", 0, "ired error"))
        self.errline.append(("ired failed, retry with cartesian coords", 0, "ired error"))
        self.errline.append(("REDCOR: STRANGE LOCAL GEOMETRY?", 0, "ired error"))
        # self.errline.append(("WARNING! Lowest Eigenvalue for BEND at Center", 0, "ired error"))
        self.errline.append((self.timeout_err, 1, self.timeout_err))

        self.tasks = { 'geometry': [],
                      'atomic': [],
                      'occupation': [],
                      'occ_num': [],
                      'general': [],
                      'scf': [],
                      'dft_on': [], 'dft_off': [],
                      'ri_on': [], 'ri_off': [],
                      'cp_on': [], 'cp_off': [], 
                      'dsp': [] }

        if ip.Coords != '' and 'continue' not in ip.Opts.keys():
            self.tasks['geometry'].append('del all')
            if ip.Coords[0] != '!':
                os.popen("x2t ../" + ip.Coords + " > coord").close()
                self.tasks['geometry'].append('a coord')
            else:
                self.tasks['geometry'].append('a ! ' + ip.Coords[1:])

            if 'desy' in ip.Opts.keys():
                if ip.Opts['desy'] == '':
                    self.tasks['geometry'].append('desy')
                else:
                    self.tasks['geometry'].append('desy ' + ip.Opts['desy'])

            self.tasks['geometry'].append('ired')

            if 'frag' in ip.Opts.keys():
                fl = {}
                for i in ip.Opts['frag']:
                    for j in i[1]:
                        fl[j] = i[0]
                fll = ['0'] * len(fl)
                for i in range(0, len(fll)):
                    fll[i] = fl[str(i + 1)]
                self.tasks['geometry'].append('frag')
                self.tasks['cp_off'] = ['on']
                self.tasks['fragment'] = fll
            else:
                self.tasks['geometry'].append('frag')
                self.tasks['cp_on'] = ['off']

            if 'ff' in ip.Opts.keys():
                self.tasks['geometry'].append('ff')
                self.tasks['ff'] = ['c ' + str(ip.Charge)]
                if ip.Opts['ff'] != '':
                    self.tasks['ff'].append('m ' + ip.Opts['ff'])
                self.tasks['geometry'].append('w ffopt')
                self.post_def.append('t2x ffopt > ffopt.xyz')

        if ip.Basis != '' and 'continue' not in ip.Opts.keys():
            self.tasks['atomic'].append('b all ' + ip.Basis)

        if ip.Ecp != '' and 'continue' not in ip.Opts.keys():
            self.tasks['atomic'].append('ecp all ' + ip.Ecp)

        if ip.Method in self.Funcs:
            if 'ri' in ip.Opts.keys():
                self.is_ri = True

        if ip.Method in self.Funcs:
            self.tasks['general'].append('dft')
            self.tasks['dft_off'].append('on')
            self.tasks['dft_on'].append('func ' + ip.Method)
            if 'dftgrid' in ip.Opts.keys():
                self.tasks['dft_on'].append('grid ' + ip.Opts['dftgrid'])

            if 'ri' in ip.Opts.keys():
                self.tasks['general'].append('ri')
                self.tasks['ri_off'].append('on')
                if 'mem' in ip.Opts['ri'].keys():
                    self.tasks['ri_on'].append('m ' + ip.Opts['ri']['mem'])
                if 'marij' in ip.Opts['ri'].keys():
                    self.tasks['general'].append('marij')
                    self.tasks['marij_on'] = ['no']
                else:
                    self.tasks['marij_on'] = ['yes']
                if 'trunc' in ip.Opts['ri'].keys():
                    self.tasks['general'].append('trunc')
                    self.tasks['trunc_on'] = ['no']
                else:
                    self.tasks['trunc_on'] = ['yes']
            else:
                self.tasks['general'].append('ri')
                self.tasks['ri_on'].append('off')
        else:
            self.tasks['general'].append('dft')
            self.tasks['dft_on'].append('off')
        
        if 'dsp' in ip.Opts.keys():
            self.tasks['general'].append('dsp')
            self.tasks['dsp'].append('on')
        else:
            self.tasks['general'].append('dsp')
            self.tasks['dsp'].append('off')

        if len(ip.Scf) != 0:
            self.tasks['general'].append('scf')
            for s in ip.Scf:
                if s[0].lower() == 'iter':
                    self.tasks['scf'].append('iter')
                    self.tasks['scf_iter'] = [s[1]]
                if s[0].lower() == 'conv':
                    self.tasks['scf'].append('conv')
                    self.tasks['scf_conv'] = [s[1]]
                if s[0].lower() == 'diis':
                    self.tasks['scf'].append('diis')
                    if s[1].lower() == 'fds': self.tasks['scf_diis'] = ['1']
                    elif s[1].lower() == 'sfds': self.tasks['scf_diis'] = ['2']
                    elif s[1].lower() == 'none': self.tasks['scf_diis'] = ['3']

        if 'continue' not in ip.Opts.keys():
            self.tasks['occupation'].append("eht")
            if ip.Multiplicity == 0:
                self.tasks['acc_occupation'] = ['y']
            else:
                self.tasks['acc_occupation'] = ['n']
                self.tasks['occ_num'].append("u " + str(ip.Multiplicity - 1))

        while self._read():
            if self.state == "control":
                self._write("")
            elif self.state == "title":
                self._write(ip.Title)
            elif self.state in self.tasks.keys():
                if len(self.tasks[self.state]) != 0:
                    self._write(self.tasks[self.state].pop(0))
                else:
                    self._write('*' if self.state not in self.blank_quit else '')
            elif len(self.state) > 4 and self.state[0:4] == 'ask:' and self.state[4:] in self.tasks.keys():
                if len(self.tasks[self.state[4:]]) != 0:
                    self._write('y')
                else:
                    self._write('n')
            else:
                self._write(self.state)

        if len(self.errstate) != 0:
            self.def_out.write("!!! ERROR | " + self.errstate + '\n')
            print self.errstate

        self.def_out.flush()
        self.def_out.close()

        if len(self.errstate) == 0:

            for p in self.post_def:
                os.popen(p).close()
            if len(ip.Skr) != 0:
                os.popen("sed -i 's/file=.*twoint/file=twoint/g' control").close()
                fl = open('control', 'r').read().replace('file=twoint', 'file=' + ip.Skr + '/twoint')
                flx = open('control', 'w')
                flx.write(fl)
                flx.flush()
                flx.close()

            print ''
            if 'freq' in ip.Opts.keys():
                exitcode = self._run('aoforce')
            elif ip.Opt == 0:
                exitcode = self._run('ridft' if self.is_ri else 'dscf')
            else:
                exitcode = self._run('jobex' + (' -ri' if self.is_ri else '') + 
                  ' -relax -c ' + str(ip.Opt))
                print os.popen('echo GEO_*').readline().strip()
                if exitcode == 0:
                    os.popen('t2x coord > final.xyz').close()
                    os.popen('t2x > trajectory.xyz').close()

            print 'Normal termination' if exitcode == 0 else 'Error termination'

        else:
            print 'Define Error termination'

    def _run(self, cmd):
        self.shell = Popen(cmd, shell=True, universal_newlines=True, bufsize=1)
        return self.shell.wait()

    def __del__(self):
        if not (self.shell is None) and self.shell.poll() is None:
            self.shell.kill()

    def _write(self, x):
        assert self.shell.poll() is None
        self.def_out.write(">>> " + x + '\n')
        self.shell.stdin.write(x + '\n')

    def _read(self):
        while self.shell.poll() is None:
            line = self._readline()
            self.def_out.write(line + '\n')
            if line == '':
                continue
            for i in self.stdline:
                if i[0] in line:
                    for j in range(0, i[1]):
                        line = self._readline()
                        self.def_out.write(line + '\n')
                    self.state = i[2]
                    return True
            for i in self.errline:
                if i[0] in line:
                    self.errstate = i[2]
                    self._write('qq')
                    while self.shell.poll() is None:
                        line = self._readline()
                        self.def_out.write(line + '\n')
                    return False
        return False

    def _readline(self):
        if not self.timeout:
            return self.shell.stdout.readline().strip()
        signal.alarm(self.timeout)
        try:
            line = self.shell.stdout.readline().strip()
            return line
        except Alarm:
            return self.timeout_err
        finally:
            signal.alarm(0)

def read_opts(argv):
    opts = {}
    for i in argv:
        if not i.startswith('--'):
            raise RuntimeError('Unknown argument: %s' % i)
        if '=' in i: a, b = i.split('=')
        else: a, b = i, ''
        opts[a[2:]] = b
    return opts

if __name__ == "__main__":
    signal.signal(signal.SIGALRM, alarm_handler)
    if len(argv) > 1 and 'reduce' in read_opts(argv[1:]):
        ipx = QMInput(doreduce=True)
    else:
        ipx = QMInput()
    tm = RunTM(ipx)
