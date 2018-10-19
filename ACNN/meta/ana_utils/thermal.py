
import numpy as np
from cluster.coval import AtomicWeight
from cluster.base import Cluster
from surface.base import ClusterAtSurface
from parallel.freqs import freqs_net, imag
from utils.base import amutoau, angtobohr, ktoht, cmitoht

def rot_order(pgroup):
    x = { "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6, "C7": 7, "C8": 8, 
        "D2": 4, "D3": 6, "D4": 8, "D5": 10, "D6": 12, 
        "C2v": 2, "C3v": 3, "C4v": 4, "C5v": 5, "C6v": 6, 
        "C2h": 2, "C3h": 3, "C4h": 4, "C5h": 5, "C6h": 6, 
        "D2h": 4, "D3h": 6, "D4h": 8, "D5h": 10, "D6h": 12, 
        "D2d": 4, "D3d": 6, "D4d": 8, "D5d": 10, "D6d": 12, 
        "S4": 2, "S6": 3, "S8": 4, "Coov": 9999, "Dooh": 9999, 
        "T": 12, "Td": 12, "Th": 12, "O": 24, "Oh": 24, "I": 60, "Ih": 60 }
    if pgroup in x: return x[pgroup]
    else: return 1

def moments_of_inertia(c):
    i = np.zeros((3, 3))
    for el, at in zip(c.elems, c.atoms):
        mi = AtomicWeight.x[el] * amutoau * angtobohr ** 2
        i[0, 0] += mi * (at[1] ** 2 + at[2] ** 2)
        i[1, 1] += mi * (at[0] ** 2 + at[2] ** 2)
        i[2, 2] += mi * (at[0] ** 2 + at[1] ** 2)
        i[0, 1] += -mi * at[0] * at[1]
        i[0, 2] += -mi * at[0] * at[2]
        i[1, 0] += -mi * at[1] * at[0]
        i[1, 2] += -mi * at[1] * at[2]
        i[2, 0] += -mi * at[2] * at[0]
        i[2, 1] += -mi * at[2] * at[1]
    return list(np.linalg.eigvalsh(i))

def gas_phase_zpe(c):
    assert isinstance(c, Cluster)
    if "freqs" not in c.props: return 0.0
    ff = np.array(freqs_net(c.props["freqs"])) * cmitoht
    zpe = 0.0
    for f in ff:
        if f > 0.0: zpe += f / 2.0
    return zpe

def surface_zpe(c):
    assert isinstance(c, ClusterAtSurface)
    if "freqs" not in c.props: return 0.0
    ff = np.array(c.props["freqs"]) * cmitoht
    zpe = 0.0
    for f in ff:
        if f > 0.0: zpe += f / 2.0
    return zpe

def surface_partition(c, t, rel_energy=0.0, vib=True, iprint=False):
    assert isinstance(c, ClusterAtSurface)
    tt = t * ktoht
    vib_part = np.float128(1.0)
    if vib:
        if "freqs" not in c.props:
            if iprint: print ("no freqs found for %s!" % c.aname)
            return 0.0
        ff = np.array(c.props["freqs"]) * cmitoht
        for f in ff:
            f = np.float128(f)
            if f > 0.0:
                assert (1 - np.exp(-f / tt)) != np.float128(0.0)
                vib_part *= np.exp(-f / (2 * tt)) / (1 - np.exp(-f / tt))
    elec_part = c.multiplicity * np.exp(-(c.energy - rel_energy) / tt)
    if not (vib_part * elec_part != 0.0 or c.multiplicity == 0):
        print("surface partition wrong!!")
        return 1.0
    return vib_part * elec_part

def gas_phase_partition(c, t, rel_energy=0.0, vib=True, rot=True, iprint=False):
    assert isinstance(c, Cluster)
    tt = t * ktoht
    rot_part = (2.0 * tt) **1.5 * \
      np.sqrt(np.prod(c.props["moments.of.inertia"]) * np.pi) / c.props["rot.order"]
    vib_part = 1.0
    if vib:
        if "freqs" not in c.props:
            if iprint: print ("no freqs found for %s!" % c.aname)
            return 0.0
        ff = np.array(freqs_net(c.props["freqs"])) * cmitoht
        if imag(c.props["freqs"]) > 0:
            if iprint: print ("ignore TS for %s!" % c.aname)
            return 0.0
        for f in ff:
            if f > 0.0:
                vib_part *= np.exp(-f / (2 * tt)) / (1 - np.exp(-f / tt))
    elec_part = c.multiplicity * c.props["spatial.degeneracy"] * \
        np.exp(-(c.energy - rel_energy) / tt)
    if vib and rot: return rot_part * vib_part * elec_part
    elif vib: return vib_part * elec_part
    elif rot: return rot_part * elec_part
    else:
        return elec_part

def gas_phase_partition_d(c, t, rel_energy=0.0, vib=True, rot=True):
    assert isinstance(c, Cluster)
    beta = 1 / (t * ktoht)
    part = - (c.energy - rel_energy)
    if vib:
        ff = np.array(freqs_net(c.props["freqs"])) * cmitoht
        for f in ff:
            if f > 0.0:
                part -= 0.5 * f + f / (np.exp(beta * f) - 1)
    if rot: part -= 3.0 / (2 * beta)
    return part

def gas_phase_partition_dd(c, t, vib=True, rot=True):
    assert isinstance(c, Cluster)
    beta = 1 / (t * ktoht)
    part = 0
    if vib:
        ff = np.array(freqs_net(c.props["freqs"])) * cmitoht
        for f in ff: 
            if f > 0.0:
                part += f**2 * np.exp(beta * f) / (np.exp(beta * f) - 1) ** 2
    if rot: part += 3.0 / (2 * beta**2)
    return part