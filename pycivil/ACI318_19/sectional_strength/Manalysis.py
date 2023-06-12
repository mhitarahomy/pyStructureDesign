from functools import partial
from typing import Tuple

import numpy as np
from numpy.typing import NDArray
from shapely import Point

from scipy.optimize import root_scalar

from pycivil.ACI318_19.designProps import DesignData
import pycivil.ACI318_19.sectional_strength.PMManalysis as PMM
import pycivil.ACI318_19.assumptions as assump


def get_neg_design_data(data: DesignData) -> DesignData:
    _data = data
    _data.section = rotate_section(_data.section, np.float32(180))
    _data.coords = rotate_rebar_coords(_data.coords, np.float32(180))
    return _data


def calc_d(data: DesignData) -> float:
    _, miny, _, maxy = data.section.bounds
    return maxy-miny-data.clear_cover-data.max_conf_rebar_size-(data.max_rebar_size/2)


def calc_dp(data: DesignData) -> float:
    return data.clear_cover+data.max_conf_rebar_size+(data.max_rebar_size/2)


def set_pos_As(data: DesignData, As: float, _d: float|None = None) -> DesignData:
    d = calc_d(data) if _d==None else _d
    minx, _, maxx, maxy= data.section.bounds
    _data = data
    _data.coords=np.array([Point((maxx+minx)/2, maxy-d)])
    _data.As = np.array([As])
    return _data


def set_neg_As(data: DesignData, As: float, _dp: float|None = None) -> DesignData:
    dp = calc_dp(data) if _dp==None else _dp
    minx, _, maxx, maxy= data.section.bounds
    _data = data
    _data.coords=np.array([Point((maxx+minx)/2, maxy-dp)])
    _data.As = np.array([As])
    return _data


# * Only for 0, 90, 180, 270 degrees
rotate_section = PMM.rotate_section

# * Only for 0, 90, 180, 270 degrees
rotate_rebar_coords = PMM.rotate_rebar_coords

calc_neutral_axis = partial(PMM.calc_neutral_axis, angle=np.float32(0))

calc_neutral_region = partial(PMM.calc_neutral_region, angle=np.float32(0))

calc_max_pressure_point = partial(PMM.calc_max_pressure_point, angle=np.float32(0))

calc_pressure_axis = partial(PMM.calc_pressure_axis, angle=np.float32(0))

calc_pressure_region = partial(PMM.calc_pressure_region, angle=np.float32(0))

calc_es = partial(PMM.calc_es, angle=np.float32(0))

calc_ec = partial(PMM.calc_ec, angle=np.float32(0))

calc_fs = partial(PMM.calc_fs, angle=np.float32(0))

calc_Fs = partial(PMM.calc_Fs, angle=np.float32(0))

calc_Cc = partial(PMM.calc_Cc, angle=np.float32(0))


def calc_Fsz(data: DesignData, c: np.float32) -> np.float32:
    Fs = calc_Fs(data, c)
    yCoords = np.array([point.y for point in data.coords])
    return np.sum(Fs * yCoords)


def calc_pos_M(data: DesignData, c: np.float32, IsPhi: bool = True) -> np.float32:
    Fszy = calc_Fsz(data, c)
    Cc = calc_Cc(data, c)
    pressure_region = calc_pressure_region(data.section, data.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data.section, data.coords, c)
    phi = assump.PHI_MOMENT_AXIAL(data.fy, data.Es, min(es))
    M = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    return M


def calc_neg_M(data: DesignData, c: np.float32, IsPhi: bool = True) -> np.float32:
    _data = get_neg_design_data(data)
    return calc_pos_M(_data, c, IsPhi)


calc_P = partial(PMM.calc_P, angle=np.float32(0))


def calc_pos_c(data: DesignData) -> np.float32:
    def _optim_c(x):
        return calc_P(data, x)
    _, miny, _, maxy = data.section.bounds
    return root_scalar(_optim_c, bracket=[0.0001, 2*(maxy-miny)]).root


def calc_neg_c(data: DesignData) -> np.float32:
    neg_data = get_neg_design_data(data)
    return calc_pos_c(neg_data)


def calc_As_pos_max(data: DesignData, _d: float|None=None, As_interval:float = 50):
    As_list = np.array([], dtype=np.float32)
    _As = 1 
    while True:
        _data = set_pos_As(data, _As, _d)
        _c = calc_pos_c(_data)
        if abs(np.min(calc_es(_data.section, _data.coords, _c))) < 0.005: break
        As_list = np.append(As_list, _As)
        _As += As_interval
    return np.max(As_list)


def calc_As_neg_max(data: DesignData, _dp: float|None=None, As_interval:float = 50):
    As_list = np.array([], dtype=np.float32)
    _As = 1
    while True:
        _data = set_neg_As(data, _As, _dp)
        _c = calc_neg_c(_data)
        if abs(np.min(calc_es(_data.section, _data.coords, _c))) < 0.005: break
        As_list = np.append(As_list, _As)
        _As += As_interval
    return np.max(As_list)


def calc_Mx_max(data: DesignData, _d: float|None=None, _dp: float|None=None) -> Tuple[np.float32, np.float32]:
    As_max_pos = calc_As_pos_max(data, _d)
    c_p_max = calc_pos_c(set_pos_As(data, As_max_pos, _d))
    As_max_neg = calc_As_neg_max(data, _dp)
    c_n_max = calc_neg_c(set_neg_As(data, As_max_neg, _dp))
    return calc_pos_M(data, c_p_max), calc_neg_M(data, c_n_max) 


def calc_Mn(data: DesignData) -> Tuple[np.float32, np.float32]:
    c_pos = calc_pos_c(data)
    c_neg = calc_neg_c(data)
    return calc_pos_M(data, c_pos), calc_neg_M(data, c_neg)


def calc_M_ratio(data: DesignData, Mux: np.float32) -> np.float32:
    Mn_pos, Mn_neg = calc_Mn(data)
    Mn = Mn_pos if Mux>=0 else Mn_neg
    return Mux/Mn


def calc_percent(data: DesignData, Mux: np.float32, _d: float|None=None) -> np.float32:
    def _optim_As(x):
        _data = set_pos_As(data, x, _d) if Mux>=0 else set_neg_As(data, x, _d)
        return calc_M_ratio(_data, Mux) - 1

    As_max = calc_As_pos_max(data, _d) if Mux>=0 else calc_As_neg_max(data, _d)
    M_pos_max, M_neg_max = calc_Mx_max(data, _d, _d)
    M_max = M_pos_max if Mux>=0 else M_neg_max
    if Mux > M_max: raise Exception(f"Mu is greater than Mn max ({M_max}N.mm)")
    output = root_scalar(_optim_As, bracket=[1, As_max])
    return output.root
