
from pg.driver import BaseDriver, solve_multi
import os, time
from pg.utils import read_json, write_json, read_opts, optcopy
from pg.hosts import get_host, get_model
from gpu.hosts import gpu_models, gpu_host_cores

class GPUDriver(BaseDriver):
    def __init__(self, argv):
        super(GPUDriver, self).__init__(argv)
        self.tasks = [ "init", "fit", "opt", "recv", "send", "pretest" ]
        if "TMP_HOST" not in os.environ:
            raise RuntimeError("must set TMP_HOST environ variable!")
        self.tmpdir = os.environ["TMP_HOST"]
    
    def send(self, args):
        opts = { "id": "0", "stage": "1", "dir": self.tmpdir }
        optl = [ ] + opts.keys()
        def_pos = { "0": "id", "1": "stage" }
        opts.update(read_opts(args, def_pos, optl))
        if "." in opts["stage"]:
            cmd = "zip %s.zip %s" % (opts["id"], opts["stage"])
        elif opts["stage"].startswith("G"):
            cmd = "zip %s.zip master/%s.*/*" % (opts["id"], opts["stage"])
        else:
            cmd = "zip %s.zip master/%s.*/par_structs.xyz.0" % (opts["id"], opts["stage"])
        print (os.popen(cmd).read().strip())
        cmd = "scp %s.zip %s" % (opts["id"], opts["dir"])
        for l in os.popen(cmd):
            print (l.strip())
        cmd = "rm %s.zip" % (opts["id"])
        print (os.popen(cmd).read().strip())
    
    def recv(self, args):
        opts = { "id": "0", "dir": self.tmpdir }
        optl = [ ] + opts.keys()
        def_pos = { "0": "id" }
        opts.update(read_opts(args, def_pos, optl))
        cmd = "scp %s/%s.zip ." % (opts["dir"], opts["id"])
        for l in os.popen(cmd):
            print (l.strip())
        cmd = "unzip %s.zip" % (opts["id"])
        print (os.popen(cmd).read().strip())
        cmd = "rm %s.zip" % (opts["id"])
        print (os.popen(cmd).read().strip())

    def pre_info(self):
        if os.path.isfile("./net-template.json"):
            net_temp = read_json("./net-template.json")
        else:
            raise RuntimeError("please run 'gpuopt init ...' first!")
        xmodel = net_temp["gpu_hosts"]["model"]
        xhost = net_temp["gpu_hosts"]["name"]
        self.scripts_render = ScriptsRender(self.scripts_dir, self.scripts_templates_dir,
                                            self.scripts_spec_dir, xhost, xmodel)
        # self.model_dir = self.scripts_dir + "." + xmodel
        return net_temp
    
    def fit(self, args):
        return self.fit_or_opt(args, task="fit")
    
    def opt(self, args):
        return self.fit_or_opt(args, task="opt")
    
    def pretest(self, args):
        pre = self.pre_info()
        opts = { "ecut": "0.20" }
        optl = [ ] + opts.keys()
        def_pos = { "0": "stage", "1": "multi" }
        pre["tasks"] = [ "fit" ]
        pre["output_dir"] = "./ptest"
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "stage", "multi" ]:
            if k not in opts:
                raise RuntimeError("no %s argument found!" % k)
        pre["fitting"]["input_file"] = "./master/%s.%s/par_structs.xyz.0" % (opts["stage"], 
            opts["multi"])
        pre["outer_loop"] = 0
        pre["write_summary"] = False
        pre["fitting"]["energy_cut"] = "+%.2f" % float(opts["ecut"])
        write_json(pre, "./ptest.json")
        cmd = "$ACNNHOME/acnnmain ./ptest.json"
        for l in os.popen(cmd, 'r', 1):
            print (l.strip())
        os.remove("./ptest.json")
        if os.path.exists("./ptest"):
            os.rmdir("./ptest")

    def fit_or_opt(self, args, task):
        pre = self.pre_info()
        if task == "fit":
            opts = { "time": "7:00:00" }
            pre["tasks"] = [ "fit" ]
        else:
            opts = { "time": "2:00:00" }
            pre["tasks"] = [ "opt", "filter" ]
        optl = [ "ecut", "step", "epochs", "ref" ] + opts.keys()
        def_pos = { "0": "stage", "1": "multi" }
        opts.update(read_opts(args, def_pos, optl))
        for k in [ "stage", "multi" ]:
            if k not in opts:
                raise RuntimeError("no %s argument found!" % k)
        if "ecut" in opts:
            pre["fitting"]["energy_cut"] = "+%.2f" % float(opts["ecut"])
        if "step" in opts:
            pre["fitting_net"]["step"] = float(opts["step"])
        if "epochs" in opts:
            pre["fitting_net"]["epochs"] = int(opts["epochs"])
        
        # multi and random seed list solve
        mk = solve_multi(opts["multi"])
        rs = None
        if "rseed" in opts:
            rs = [int(x) for x in opts["rseed"].split(",")]
            while len(rs) < len(mk): rs += [ rs[-1] ]
        
        if "ref" in opts:
            pre["output_dir"] = "../../master/G@RUNS.REF"
            pre["fitting_net"]["load_network"] = -1
            if len(mk) > 1:
                raise RuntimeError("Reference cannot be multiple multiplicities!")
        else:
            pre["output_dir"] = "../../master/G@RUNS.@MULT"
            pre["fitting_net"]["load_network"] = \
                "../../master/G%s.REF/fit_network.dill.0" % opts["stage"]
        
        # for each multiplicity
        mpr = [ pre["output_dir"], pre["fitting"]["input_file"], 
            pre["optimization"]["input_file"] ]
        for im, m in enumerate(mk):
            if rs is None:
                pre["random_seed"] = int(time.time() * 99991) % (2**32 - 1)
            else:
                pre["random_seed"] = rs[im]
            sopts = { "@RUNS": opts["stage"], "@MULT": str(m) }
            mpx = []
            for x in mpr:
                for k, v in sopts.items():
                    x = x.replace(k, v)
                mpx.append(x)
            [ pre["output_dir"], pre["fitting"]["input_file"], 
                pre["optimization"]["input_file"] ] = mpx
            
            if im != 0:
                pre["fitting_net"]["load_network"] = \
                    "../../master/G%s.%d/fit_network.dill.0" % (opts["stage"], mk[im - 1])

            # gtomaster dir
            if "ref" in opts:
                gmdir = "./gtomaster/GR-%s" % (opts["stage"])
            else:
                gmdir = "./gtomaster/G-%s.%s" % (opts["stage"], str(m))
            print ("create %s" % gmdir)
            if not os.path.exists(gmdir): os.makedirs(gmdir)
            os.chdir(gmdir)
            write_json(pre, "./net.json")
            # run-gpu.sh
            ropts = {
                "@NPROCS": str(pre["gpu_hosts"]["cores"]), 
                "@TIME": opts["time"], 
                "@PROJ": self.project_name
            }
            optcopy(self.scripts_render.get("run-gpu.sh"), "./run-gpu.sh", ropts)
            # gmeta_info.txt
            f = open("./gmeta_info.txt", "w")
            f.write("hostname = %s\n" % pre["gpu_hosts"]["name"])
            f.write("cores = %d\n" % pre["gpu_hosts"]["cores"])
            xosname = os.popen("uname -s").read().strip()
            if xosname == "Darwin":
                cmd = "sysctl -n hw.cpufrequency | xargs -I % expr % / 1000000"
                f.write("MHz = %s\n" % os.popen(cmd).read().strip())
                cmd = "sysctl -n machdep.cpu.brand_string"
                f.write("CPU = %s\n" % os.popen(cmd).read().strip())
                cmd = "echo '-productName' '-productVersion' | xargs -n 1 sw_vers | xargs"
                f.write("os = %s\n" % os.popen(cmd).read().strip())
            else:
                cmd = "lscpu | grep 'MHz' | awk '{print $3}'"
                f.write("MHz = %s\n" % os.popen(cmd).read().strip())
                cmd = "lscpu | grep 'Vendor' | awk '{print $3}'"
                f.write("CPU = %s\n" % os.popen(cmd).read().strip())
                cmd = "lsb_release -d | awk -F : '{print $2}'"
                f.write("os = %s\n" % os.popen(cmd).read().strip())
            f.close()
            os.chdir("../..")

    def init(self, args):
        opts = { }
        optl = [ "model", "cores", "gmodel" ] + opts.keys()
        if os.path.isfile("./net-template.json"):
            net_temp = read_json("./net-template.json")
        else:
            net_temp = read_json(self.scripts_dir + "/net-template.json")
        def_pos = { }
        opts.update(read_opts(args, def_pos, optl))
        xhost = get_host()
        if xhost != "hoffman":
            if "PROJECT_NAME" not in os.environ:
                raise RuntimeError("must set PROJECT_NAME environ variable!")
            self.project_name = os.environ["PROJECT_NAME"]
        if "gmodel" in opts:
            xgmodel = opts["gmodel"]
        elif "gpu_hosts" in net_temp:
            xgmodel = net_temp["gpu_hosts"]["gmodel"]
        else:
            xgmodel = gpu_models[xhost]
        if "model" in opts:
            xmodel = opts["model"]
        elif "gpu_hosts" in net_temp:
            xmodel = net_temp["gpu_hosts"]["model"]
        else:
            xmodel = get_model(xhost)
        if "cores" in opts:
            xgcores = int(opts["cores"])
        elif "gpu_hosts" in net_temp:
            xgcores = net_temp["gpu_hosts"]["cores"]
        else:
            xgcores = gpu_host_cores[xhost]
        net_temp["gpu_hosts"] = { "name": xhost, "cores": xgcores, "gmodel": xgmodel, 
            "model": xmodel }
        write_json(net_temp, "./net-template.json")
    
    def help(self, args):
        print """
        gpuopt init [--opt[=value]]
            --gmodel=<auto>: override auto determined gpu acc model
            --model=<auto>:    override auto determined host model
            --cores=<auto>:    cores per gpu node
        
        gpuopt fit [--stage=]1 [--multi=]1
            --rseed=: random seed
            --time=7:00:00
            --step=0.1
            --epochs=1400
            --ecut=0.20
            --ref
        
        gpuopt fit 1 1 --time=7:00:00 --step=0.1 --epochs=1400 --ecut=0.20 --ref
        gpuopt fit 1 3-9,1 --time=3:00:00 --step=0.01 --epochs=500 --ecut=0.20
        
        gpuopt pretest [--stage=]1 [--multi=]1 --ecut=0.20
        
        gpuopt send [[--id=]0] [[--stage=]1] [--dir=..]
        gpuopt recv [[--id=]0] [--dir=..]

        [examples-mac]
        gpuopt init --gmodel=mac --cores=4 --model=mac
        """
    
