
import os

def get_host():
    xhost = os.popen("hostname -f | grep '\\.' || hostname -a").read().strip()
    keys = ['hoffman', 'gordon', 'conrad', 'garnet', 'armstrong', 'haise',
            'kilrain', 'excalibur', 'shepard', 'topaz', 'thunder', 'mac',
            'lightning', 'mira', 'bridges', 'kalk', 'onyx', 'centennial']
    for k in keys:
        if k in xhost:
            return k
    keys = [('emsl', 'cascade'), ('nd200', 'kalk')]
    for k, v in keys:
        if k in xhost:
            return v
    xhost = os.popen("hostname -f").read().strip()
    for k in keys:
        if k in xhost:
            return k
    return 'unknown'

host_cores = {'hoffman': 16, 'gordon': 32, 'conrad': 32, 'garnet': 32,
              'armstrong': 24, 'haise': 16, 'kilrain': 16, 'excalibur': 32,
              'topaz': 36, 'thunder': 36, 'cascade': 16, 'mac': 4,
              'lightning': 24, 'shepard': 24, 'mira': 16, 'bridges': 28,
              'kalk': 8, 'onyx': 44, 'centennial': 40}

# do not use this in computational nodes
# not useful for sync
def get_model(thost):
    if thost in ['gordon', 'conrad', 'garnet', 'armstrong', 'excalibur',
                 'shepard', 'lightning', 'onyx']:
        return 'cray'
    elif thost in ['haise', 'kilrain', 'topaz', 'thunder', 'centennial']:
        return 'idataplex'
    elif thost in ['hoffman']:
        return 'hoffman'
    elif thost in ['mac']:
        return 'mac'
    elif thost in ['cascade']:
        return 'cascade'
    elif thost in ['mira']:
        return 'mira'
    elif thost in ['bridges']:
        return 'bridges'
    elif thost in ['kalk']:
        return 'kalk'
    else:
        raise RuntimeError('Unknown host! %s' % thost)

def get_job_cmd(xmodel):
    job_cmd = {'cascade': ('msub', 'canceljob'), 'bridges': ('sbatch', 'scancel')}
    if xmodel in job_cmd:
        return job_cmd[xmodel]
    else:
        return ('qsub', 'qdel')
