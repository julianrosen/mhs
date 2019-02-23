def set_vars(*var):
    global _max_weight, _num, _mhs_basis, _weights, _mzv_alg
    global _mzv_mono, _mzv_to_mono, _mono_to_basis, _mhs_to_mono
    if len(var) == 9:
        (_max_weight, _num, _mhs_basis, _weights, _mzv_alg,
             _mzv_mono, _mzv_to_mono, _mono_to_basis, _mhs_to_mono) = var

def max_weight():
    return _max_weight

def num():
    return _num

def mhs_basis():
    return _mhs_basis

def weights():
    return _weights

def mzv_alg():
    return _mzv_alg

def mzv_mono():
    return _mzv_mono

def mzv_to_mono():
    return _mzv_to_mono

def mono_to_basis():
    return _mono_to_basis

def mhs_to_mono():
    return _mhs_to_mono