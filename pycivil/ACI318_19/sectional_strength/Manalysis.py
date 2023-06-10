from functools import partial
from typing import Tuple

import numpy as np
from numpy.typing import NDArray
from shapely import Point

from scipy.optimize import root_scalar

from pycivil.ACI318_19.designProps import DesignData
import pycivil.ACI318_19.sectional_strength.PMManalysis as PMM
import pycivil.ACI318_19.assumptions as assump

def set_As(data: DesignData, As: NDArray[np.float32]) -> DesignData:
    """create design data with new array of rebar area

    Args:
        data (DesignData): design data
        As (NDArray[np.float32]): array of rebar area, [mm2]

    Returns:
        DesignData: design data
    """
    return DesignData(section=data.section, bw=data.bw, fy= data.fy, 
                      fyt=data.fyt, fc=data.fc, coords=data.coords, As=As,
                      Es=data.Es, Av=data.Av, conf_dist=data.conf_dist, 
                      clear_cover=data.clear_cover, conf_type=data.conf_type, 
                      max_conf_rebar_size=data.max_conf_rebar_size,
                      max_rebar_size=data.max_rebar_size)


def get_As_percent(data: DesignData) -> np.float32:
    """get percent of rebars area to concrete area

    Args:
        data (DesignData): design data

    Returns:
        np.float32: percent, %
    """
    return (np.sum(data.As)/data.section.area)*100


def set_As_percent(data: DesignData, percent: float) -> DesignData:
    """create new design data with percent of rebars area to concrete area

    Args:
        data (DesignData): design data
        percent (float): percent, %

    Returns:
        DesignData: new design data
    """
    totalAs = data.section.area * (percent/100)
    return set_As(data, np.array([totalAs/ len(data.As) for i in range(len(data.As))]))


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


def calc_M(data: DesignData, c: np.float32, IsPhi: bool = True) -> Tuple[np.float32, np.float32]:
    Fszy = calc_Fsz(data, c)
    Cc = calc_Cc(data, c)
    pressure_region = calc_pressure_region(data.section, data.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data.section, data.coords, c)
    phi = assump.PHI_MOMENT_AXIAL(data.fy, data.Es, min(es))
    M_pos = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    
    data_180 = data
    data_180.section = rotate_section(data_180.section, np.float32(180))
    data_180.coords = rotate_rebar_coords(data_180.coords, np.float32(180))
    Fszy = calc_Fsz(data_180, c)
    Cc = calc_Cc(data_180, c)
    pressure_region = calc_pressure_region(data_180.section, data_180.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data_180.section, data_180.coords, c)
    phi = assump.PHI_MOMENT_AXIAL(data_180.fy, data_180.Es, min(es))
    M_neg = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    return M_pos, M_neg


calc_P = partial(PMM.calc_P, angle=np.float32(0))


def calc_d(data: DesignData) -> float:
    _, miny, _, maxy = data.section.bounds
    return maxy-miny-data.clear_cover-data.max_conf_rebar_size-(data.max_rebar_size/2)


def calc_As_max(data: DesignData, _d: float|None=None, As_interval:float = 50):
    d = calc_d(data) if _d==None else _d
    minx, _, maxx, maxy= data.section.bounds
    _data = data
    _data.coords=np.array([Point((maxx+minx)/2, maxy-d)])
    _data.As = np.array([0])
    As_list = np.array([], dtype=np.float32)
    # c_list = np.array([], dtype=np.float32)
    _As = 1 
    while True:
        _data = set_As(_data, np.array([_As]))
        _c = calc_c(_data)
        if abs(np.min(calc_es(_data.section, _data.coords, _c))) < 0.005: break
        As_list = np.append(As_list, _As)
        # c_list = np.append(c_list, _c)
        _As += As_interval
    return np.max(As_list)


def calc_Mx_max(data: DesignData, _d: float|None=None):
    As_max = calc_As_max(data, _d)
    


# def calc_c_max(data: DesignData) -> np.float32:
#     ety = data.fy / data.Es
#     _, _, _, maxy = data.section.bounds
#     miny_rebar = np.min([point.y for point in data.coords])
#     dt = maxy - miny_rebar
#     return dt / (1-(ety/assump.ECU))


def calc_c(data: DesignData) -> np.float32:
    def _optim_c(x):
        return calc_P(data, x)
    _, miny, _, maxy = data.section.bounds
    return root_scalar(_optim_c, bracket=[0.0001, 2*(maxy-miny)]).root


def calc_Mn(data: DesignData) -> Tuple[np.float32, np.float32]:
    c = calc_c(data)
    return calc_M(data, c)


def calc_M_ratio(data: DesignData, Mux: np.float32) -> np.float32:
    Mn_pos, Mn_neg = calc_Mn(data)
    Mn = Mn_pos if Mux>=0 else Mn_neg
    return Mux/Mn



# ! min and max reinforcement must add
def calc_percent(data: DesignData, Mux: np.float32, _d: float|None=None) -> np.float32:
    d = calc_d(data) if _d==None else _d
    minx, _, maxx, maxy= data.section.bounds
    _data = data
    _data.coords=np.array([Point((maxx+minx)/2, maxy-d)])

    def _optim_As(x):
        __data = set_As(_data, np.array([x]))
        return calc_M_ratio(__data, Mux) - 1

    As_max, c_max = calc_As_c_max(data, _d)
     



    pass
    # def _optim_percent(x):
    #     data_percent = set_As_percent(data, x)
    #     return calc_M_ratio(data_percent, Mux) - 1
    # data_one_percent = set_As_percent(data, 1)
    # data_eight_percent = set_As_percent(data, 8)
    # if calc_M_ratio(data_one_percent, Mux) <= 1: 
    #     output_percent = 1
    # elif calc_M_ratio(data_eight_percent, Mux) > 1:
    #     raise ValueError("section is weak")
    # else:
    #     output_percent = root_scalar(_optim_percent, bracket=[1, 8]).root
    # return np.float32(output_percent)
