
from acnn.layers import BasicLayerSize as B, CoLayerSize as C, \
  FeatureLayerSize as F, MixedLayerSize as M
from acnn.layers import Reshape, Average, Minimum, Maximum
from acnn.layers import Sigmoid, Tanh, TanhSig, Softmax, Softplus, Linear, Relu

ldic = { 'average': Average, 'minimum': Minimum, 'maximum': Maximum, 
  'reshape': Reshape, 'sigmoid': Sigmoid, 'tanh': Tanh, 'tanhsig': TanhSig, 
  'softmax': Softmax, 'softplus': Softplus, 'linear': Linear, 'relu': Relu }
sdic = { 'B': B, 'C': C, 'F': F, 'M': M }

# translate layers short notation for json file
# to layers object

class LayerList(object):
  def __init__(self, layers):
    self.layers = layers
  
  def __repr__(self):
    rep = '\n'
    for layer in self.layers:
      rep += ('{:3d} -- {}'.format(layer.layerid, layer)) + '\n'
    rep += ('{:3s} -- {}'.format('', 'END')) + '\n'
    return rep

def trans_layers(list, ipt, diff=False):
  layers = []
  prev_x = [C(i) for i in ipt]
  prev_x[-1] = B(prev_x[-1].n)
  prev = [i for i in ipt]
  for l in list:
    lc = [m for m in l.split(':') if len(m) != 0]
    if len(lc) == 1:
      layers.append(ldic[lc[0]]())
      prev = layers[-1].size_output(prev)
    else:
      lr = [m for m in lc[1].split('>') if len(m) != 0]
      if len(lr) == 1:
        ls = [m for m in lr[0].split(' ') if len(m) != 0]
        if len(ls) == 1 and ls[0] == '..':
          ls = prev
        else:
          for xi in range(len(ls)):
            if ls[xi] == '.': ls[xi] = prev[xi]
            elif ls[xi][0] == '.':
              ls[xi] = prev[int(ls[xi][1:])]
            else:
              ls[xi] = int(ls[xi])
        if len(ls) != 1: ls = tuple(ls)
        else: ls = ls[0]
        layers.append(ldic[lc[0]](ls))
        prev = layers[-1].size_output(prev)
      else:
        par = []
        for g in lr:
          ls = [m for m in g.split(' ') if len(m) != 0]
          if len(ls) == 1 and ls[0] == '..':
            ls = prev_x[:]
          else:
            for xi in range(len(ls)):
              if ls[xi] == '.': ls[xi] = prev_x[xi]
              elif ls[xi][0] == '.':
                ls[xi] = prev_x[int(ls[xi][1:])]
              elif ls[xi][0] in ['B', 'C', 'F'] and ls[xi][1] == '.':
                if len(ls[xi][2:]) == 0: ls[xi] = sdic[ls[xi][0]](prev_x[xi].n)
                else: ls[xi] = sdic[ls[xi][0]](prev_x[int(ls[xi][2:])].n)
              elif ls[xi][0] in ['B', 'C', 'F']:
                ls[xi] = sdic[ls[xi][0]](int(ls[xi][1:]))
              else:
                ls[xi] = B(int(ls[xi]))
          par.append(ls[:])
          prev_x = ls[:]
        layers.append(ldic[lc[0]](*par))
        prev = [i.n for i in prev_x]
  if diff:
    for layer in layers:
      if isinstance(layer, ParameterLayer):
        isize = [C(2)] + layer.input_size
        osize = [C(2)] + layer.output_size
        layer.init_size(isize, osize)
      elif isinstance(layer, Reshape):
        shape = (2, ) + layer.shape if layer.shape != None else (2, -1)
        layer.shape = shape
      elif isinstance(layer, AxisPoolLayer):
        layer.axis += 1
    layers.append(Diff(0))
  print ('layers resolved: {}'.format(LayerList(layers)))
  return layers
