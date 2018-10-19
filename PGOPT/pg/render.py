
# for rendering the scripts for different hosts
import copy
import os

class ScriptPart(object):
    def __init__(self, name, text, aint=0):
        self.text = text
        self.name = name
        self.aint = aint

def search(parts, name):
    for p in parts:
        if p.name == name:
            return p.text
    return None

def render(base, target):
    base = copy.deepcopy(base)
    for t in target:
        if t.name is None:
            continue
        for b in base:
            if b.name == t.name:
                b.text = t.text
                ci = b.aint - t.aint
                if ci != 0 and b.text != "":
                    b.text = (" " * ci) + b.text.replace("\n", "\n" + " " * ci)
    return base

def tostring(parts):
    return "\n".join([p.text for p in parts])

def decompose(fn):
    parts = []
    with open(fn, "r") as f:
        aname = None
        atext = ""
        aint = 0
        for line in f:
            if line.lstrip().startswith("### start:"):
                assert aname is None
                aname = line.split(":")[1].strip()
                aint = len(line) - len(line.lstrip())
            elif line.lstrip().startswith("### end:"):
                assert aname == line.split(":")[1].strip()
                parts.append(ScriptPart(aname, atext, aint))
                aname, atext = None, ""
            elif line.lstrip().startswith("### ::"):
                assert aname is None
                aint = len(line) - len(line.lstrip())
                if "=" not in line:
                    parts.append(ScriptPart(line.split("::")[1].strip(), "", aint))
                else:
                    c = line.split("::")[1].strip()
                    a = c.split("=")[0].strip()
                    b = c.split("=")[1].strip()
                    btext = search(parts, b)
                    assert btext is not None
                    parts.append(ScriptPart(a, btext, aint))
            else:
                if aname is None:
                    if len(parts) != 0 and parts[-1].name is None:
                        parts[-1].text += "\n" + line.rstrip()
                    else:
                        parts.append(ScriptPart(None, line.rstrip()))
                else:
                    if atext == "":
                        atext += line.rstrip()
                    else:
                        atext += "\n" + line.rstrip()
    return parts

class ScriptsRender(object):
    def __init__(self, scri_dir, temp_dir, spec_dir, hostname, modelname):
        self.temp_dir = temp_dir
        self.spec_dir = spec_dir
        self.scri_dir = scri_dir
        self.scripts = {}
        self.templates = {}
        self.specs = {}
        self.hostname = hostname
        self.modelname = modelname
        self.read_dirs()

    def read_dirs(self):
        for k in os.listdir(self.scri_dir):
            if k.endswith(".sh") or k.endswith(".in"):
                self.scripts[k] = open(self.scri_dir + "/" + k, "r").read()
        for k in os.listdir(self.temp_dir):
            if k.endswith(".sh"):
                self.templates[k] = decompose(self.temp_dir + "/" + k)
        for k in os.listdir(self.spec_dir):
            if k.endswith(".spec.sh"):
                l = k[:-8]
                self.specs[l] = decompose(self.spec_dir + "/" + k)

    def get(self, scriptname):
        hostname = self.hostname
        modelname = self.modelname
        if scriptname in self.templates:
            base = self.templates[scriptname]
            if modelname in self.specs:
                base = render(base, self.specs[modelname])
            if hostname != modelname and hostname in self.specs:
                base = render(base, self.specs[hostname])
            return tostring(base)
        else:
            return self.scripts[scriptname]
