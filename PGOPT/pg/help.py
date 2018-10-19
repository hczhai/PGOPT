
class PGHelp(object):
  def __init__(self):
    pass
  
  def init(self):
    """
    pgopt init name [--opt[=value]]
    pgopt init name [number [method]] [--opt[=value]]
      --no-scratch:  store all files in current directory
      [--number=]200:  default number of structures to generate
      --charge=0:    default charge of the system
      [--method=]blda: default creation method
      --func=tpssh:  default DFT functional (for vasp, pbe)
      --basis=def2-SV(P):  default basis set (for vasp, '')
      --program=tm:  default DFT program (turbomole) (or vasp)
      --model=<auto>:  override auto determined model
      --cores=<auto>:  cores per node
      --nodes=1: number of node used by each processor
      --surface=: [mgo] surface name
      --encut=: will be added to (vasp) incar
      --periodic: periodic generation mode
      this cmd will also update host information (name, model, cores)
    """
  
  def torun(self):
    """
    pgopt torun [--runt=]relax/freq [from to] [--opt[=value]]
      this will create both torun and procs directories
      [--from=]0: start index of processors
      [--to=]50:  end index of processors (not included)
      --time=23:00:00: time limit of each processor
      --clean:   clean corresponding torun submission scripts
      --x:   use job array in hoffman, no effect for non-hoffman
    """
  
  def relax(self):
    """
    pgopt relax [--stage=]1 [--multi=]1-9 [--opt[=value]]
      --sec=:   the stage will be treated as secondary stage
      --time=23:00:00: time limit of master process
      --source=create: structure source (comma separated)
        format: create or [stage[.multi]-]name.idx[:max_number]
      --ref=:   reference structures to avoid duplicates (comma separated)
        format: [stage[.multi]-]name.idx, name can be local or structs
      --rseed=: random seed (comma separated for each multi)
      --nprev=30: structures from previous multiplicity
      --max-config=-1: max number of inital configs, -1 means no limit
      --idx=0:  index in master directory
      --sp:     single point calculation
      --sp=n:   n step relax, step will be n, max_step will be n / 2
      --batch:  only generate batch script
      --copy:   if source/ref files in local but not in remote, will copy 
                from local to remote
      [examples]
        append new structures: 
          pgopt relax 2 --sec=1 --source=create --ref=1-structs.0
        relax transition states:
          pgopt relax 4 --sec=1 --source=1-structs.1 --ref=1-structs.0 --nprev=0
        relax in new basis: 
          pgopt relax 3 --source=1-local.0,4-local.0 --nprev=0 --max-config=..
    """
  
  def freq(self):
    """
    pgopt freq [--stage=]1 [[--multi=]1-9] [--opt[=value]]
      implies: 
        pgopt relax --idx=1 --source=local.0 --nprev=0 --runt=freq --sp
      if no multi, will be determined by searching relax dirs in tomaster
    """
  
  def submit(self):
    """
    pgopt submit relax/freq [--stage=]1 [[--multi=]1-9]
      if multi is not a comma expression, will submit run-master.sh
      otherwise, will submit run-batch.sh
      if no multi, will submit most recent run-batch.sh
    pgopt submit torun [--runt=]relax/freq [from to]
    """
  
  def depend(self):
    """
    pgopt depend relax/freq [--stage=]1 [[--multi=]1-9] --jobid=1
      if multi is not a comma expression, will change run-master.sh
      otherwise, will change run-batch.sh
      if no multi, will change most recent run-batch.sh
      change the batch/master scripts so that it will excute
        after (the jobid of a batch/master) ends
    pgopt depend torun [--runt=]relax/freq [from to] --jobid=1
      --jobid=: the job id. if empty (default), will remove dependence
      change the torun scripts so that it will excute
        after (the jobid of a batch/master) begins
    """
  
  def sync(self):
    """
    pgopt sync [--restarts]
      --restarts: will also sync restarts files
      within run-batch.sh or run-master.sh, pgopt sync is executed every 2 min
      if init with --no-scratch, this will generate error
    """
  
  def show(self):
    """
    pgopt show [[--dir=]local/remote]
      show local/remote dir path
      if no local/remote given, will show the path other than current one
    """
  
  def report(self):
    """
    pgopt report
      this will first pgopt sync and then generate the pdf report
    """
  
  def final(self):
    """
    pgopt final
      this will first pgopt sync and then zip master and tomaster
    """
  
  def log(self):
    """
    pgopt log
      show short log
    """
  
  def create(self):
    """
    pgopt create [--stage=]1 [--multi=]1-9 [--number=]5000
    """
  
  def check(self):
    """
    pgopt check [--runt=]relax/freq [[--opt=]value]
      --runt=relax: default checking relax
      --maxt=60: max time of no writting output
      --delete:  if presented, will delete job and RUNNING file
        for noresponse procs and sleeping procs
      check if proc is in good state
      autocheck (per 3600) will be enabled in scratch mode
    """
  
  def surface(self):
    """
    [examples-surface]
    pgopt init Pt5 100 --program=vasp --surface=mgo --encut=320
    pgopt relax 1 1 --time=23:00:00 --nprev=0
    pgopt torun relax 0 50 --time=23:00:00 --x
    (--x will be efficient for cascade/hoffman)
    pgopt submit relax 1 1
    pgopt submit torun relax

    <acnnmain filter.json>
    pgopt init Pt5 100 --program=vasp --surface=mgo4 --encut=400 --cores=32
    <step = 10>
    pgopt relax 2 1 --time=24:00:00 --source=F1-filter.0 --nprev=0 --max-config=24
    pgopt submit relax 2 1
    """
  
  def periodic(self):
    """
    [examples-periodic]
    pgopt init B16 200 --program=vasp --surface=wcut --periodic --encut=280 --nodes=4
    pgopt relax 1 1 --time=24:00:00 --nprev=0
    pgopt torun relax 0 24 --time=24:00:00
    pgopt submit relax 1 1
    pgopt submit torun relax
    """

  def example(self):
    """
    [examples-mac]
    pgopt init Pt9 500 blda --func=tpssh --no-scratch --model=mac --cores=4

    [examples-small basis] 
    pgopt init Pt6 50 blda --func=tpssh
    pgopt init --basis=def2-SV\\(P\\)
    pgopt torun relax 0 30 --time=24:00:00 --x
    pgopt relax 1 1-9 --time=24:00:00
    pgopt submit relax 1
    pgopt depend torun relax --jobid=1
    pgopt submit torun relax
    pgopt log

    [examples-big basis]
    pgopt init --basis=def2-TZVP
    pgopt relax 2 1-9 --time=24:00:00 --source=1-local.0 --nprev=0 --max-config=..
    pgopt submit relax 2

    [examples-big basis-freq]
    pgopt freq 2 --time=10:00:00
    pgopt torun freq 0 48 --time=10:00:00
    pgopt submit freq 2
    pgopt depend torun freq --jobid=1
    pgopt submit torun freq

    [examples-big basis-displacement re-relax]
    pgopt relax 3 1-11 --sec=2 --source=2-structs.1 --ref=2-structs.0 \\
      --nprev=0 --time=24:00:00
    pgopt submit relax 3

    [examples-pre-fit 10 step single points]
    pgopt init Pt9 1000 blda --func=tpssh --basis=def2-TZVP
    pgopt relax 1 1-9 --time=24:00:00 --sp=10 --nprev=0

    [examples-pre-fit vasp]
    pgopt init Pt13 50 blda --program=vasp
    pgopt relax 1 1 --time=23:00:00 --sp=10 --nprev=0
    pgopt torun relax 0 10 --time=10:00:00 --x
    pgopt submit relax 1 1
    pgopt submit torun relax

    [examples-after-fit relax]
    pgopt relax 2 1 --time=24:00:00 --source=G1-filter.0:500 --nprev=0 --copy
    """

  def help(self, args):
    """
    pgopt
    pgopt help [@@]
      show this help message
    """
    keys = [ "help", "init", "torun", "relax", "freq", "submit", "depend", 
      "sync", "show", "report", "log", "example", "check", "final", "create", 
      "surface", "periodic" ]
    getattr(PGHelp.help, "__func__").__doc__ = \
      self.help.__doc__.replace("@@", "/".join(keys))
    for k in keys:
      if len(args) == 0 or args[0] == k:
        print (getattr(self, k).__doc__)
