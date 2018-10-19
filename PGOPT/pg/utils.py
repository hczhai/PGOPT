
import re, os, json

# copy files, and replace some option keys with values
def optcopy(src_cnt, dst, opts):
    kk = src_cnt
    for k, v in opts.items():
        if '\n' in v:
            c = re.findall(r'\n(\s*)%s' % k, kk)
            if len(c) != 0:
                v = v.replace('\n', '\n' + c[0])
        kk = kk.replace(k, v)
    f = open(dst, 'w')
    f.write(kk)
    f.close()
    if dst.endswith('.sh'):
        os.popen('chmod +x \"%s\"' % dst).read()

# read options of the form --key=value
def read_opts(argv, def_pos=None, optl=None):
    opts = {}
    if def_pos is None:
        def_pos = {}
    for ii, i in enumerate(argv):
        if not i.startswith('--'):
            if str(ii) not in def_pos:
                raise RuntimeError('Unknown argument: %s' % i)
            else:
                i = "--%s=%s" % (def_pos[str(ii)], i)
        c = i.split('=')
        if len(c) == 1:
            a, b = c[0], ""
        else: a, b = i.split('=')
        if a[2:] in def_pos.values() or optl is None or a[2:] in optl:
            opts[a[2:]] = b
        else:
            raise RuntimeError("Unknown argument: %s" % a[2:])
    return opts

# reverse of read_opts
def write_opts(opts, def_pos=None):
    argv = []
    if def_pos == None:
        def_pos = {}
    ropts = {}
    ropts.update(opts)
    for k, v in sorted(def_pos.items(), key=lambda x: int(x[0])):
        if int(k) != len(argv) or v not in opts:
            break
        else:
            argv.append(opts[v])
            del ropts[v]
    for k, v in ropts.items():
        argv.append("--%s=%s" % (k, v))
    return argv

# read json input
def read_json(fn):
    json_data = open(fn).read()
    json_data = re.sub(r'//.*\n', '\n', json_data)
    return json.loads(json_data)

# write json summary
def write_json(json_data, fn):
    f = open(fn, 'w')
    json.dump(json_data, f, indent=4)
    f.close()
