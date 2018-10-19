
import os, json, dill, re

# avoid overwritting
def new_file_name(x):
  i = 0
  y = x + '.' + str(i)
  while os.path.isfile(y):
    i += 1
    y = x + '.' + str(i)
  return y

# avoid overwritting
def new_path_name(x):
  i = 0
  y = x + '.' + str(i)
  while os.path.exists(y):
    i += 1
    y = x + '.' + str(i)
  return y

# read json input
def read_json(fn, i=None, iprint=True, cont=False):
  if i is not None:
    fn += '.' + str(i)
  if iprint: print ('read json: ' + fn)
  json_data = open(fn).read() if not cont else fn
  json_data = re.sub(r'//.*\n', '\n', json_data)
  return json.loads(json_data)

# write json summary
def write_json(json_data, fn, new=True, i=None, iprint=True):
  if new:
    fn = new_file_name(fn)
  if isinstance(i, int): 
    fn += '.' + str(i)
  if iprint: print ('write json: ' + fn)
  f = open(fn, 'w')
  json.dump(json_data, f, indent=4)
  f.close()

def dump_data(name, obj, i=None, iprint=True):
  if i is None: 
    name = new_file_name(name)
  elif isinstance(i, int): 
    name += '.' + str(i)
  if iprint: print ('dump data: ' + name)
  with open(name, 'wb') as f:
    dill.dump(obj, f)

def load_data(name, i=None, iprint=True):
  if i is not None:
    name += '.' + str(i)
  if iprint: print ('load data: ' + name)
  with open(name, 'rb') as f:
    return dill.load(f)

def load_name(name, i=None):
  if not i is None: name += '.' + str(i)
  return name