
amutoau = 1822.88839
angtobohr = 1.0 / 0.52917721092
ktoht = 1.0 / 315774.646
cmitoht = 1.0 / 219474.6313705
evtoht = 1.0 / 27.21138505

def cast_float(value):
    import theano
    import numpy as np
    from theano.tensor.var import TensorVariable
    from theano.tensor.sharedvar import TensorSharedVariable
    import theano.tensor as T
    float_type = theano.config.floatX
    if isinstance(value, (np.matrix, np.ndarray)):
        if value.dtype != np.dtype(float_type):
            return value.astype(float_type)
        else:
            return value
    elif isinstance(value, (TensorVariable, TensorSharedVariable)):
        return T.cast(value, float_type)
    float_x_type = np.cast[float_type]
    return float_x_type(value)

def reproducible(seed=0):
    import random, numpy as np
    np.random.seed(seed)
    random.seed(seed)

def print_dict(d, group=4, pre='## '):
    dc = []
    dt = d.items()
    dt.sort()
    for i in range(len(dt) / group + 1):
        x = dt[i * group:(i + 1) * group]
        if len(x) != 0: dc.append(x)
    for x in dc: print (pre + ', '.join(["{} = {}".format(i, j) for i, j in x]))