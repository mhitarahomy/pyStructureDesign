from functools import partial
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from scipy.optimize import root_scalar

from pycivil.ACI318_19.designProps import DesignData
import pycivil.ACI318_19.PMManalysis as PMM
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
                      fyt=data.fyt, fc=data.fc, Coords=data.Coords, As=As,
                      Es=data.Es, Av=data.Av, conf_dist=data.conf_dist, 
                      cover=data.cover, conf_type=data.conf_type)


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

calc_neutral_axis = partial(PMM.calc_neutral_axis, angle=0)

calc_neutral_region = partial(PMM.calc_neutral_region, angle=0)

calc_max_pressure_point = partial(PMM.calc_max_pressure_point, angle=0)

calc_pressure_axis = partial(PMM.calc_pressure_axis, angle=0)

calc_pressure_region = partial(PMM.calc_pressure_region, angle=0)

calc_es = partial(PMM.calc_es, angle=0)

calc_ec = partial(PMM.calc_ec, angle=0)

calc_fs = partial(PMM.calc_fs, angle=0)

calc_Fs = partial(PMM.calc_Fs, angle=0)

calc_Cc = partial(PMM.calc_Cc, angle=0)

def calc_Fsz(data: DesignData, c: float) -> np.float32:
    Fs = calc_Fs(data, c)
    yCoords = np.array([point.y for point in data.Coords])
    return np.sum(Fs * yCoords)


def calc_M(data: DesignData, c: float, IsPhi: bool = True) -> Tuple[np.float32, np.float32]:
    Fszy = calc_Fsz(data, c)
    Cc = calc_Cc(data, c)
    pressure_region = calc_pressure_region(data.section, data.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data.section, data.Coords, c)
    phi = assump.phif(data.fy, data.Es, min(es))
    M_pos = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    
    data_180 = data
    data_180.section = rotate_section(data_180.section, 180)
    data_180.Coords = rotate_rebar_coords(data_180.Coords, 180)
    Fszy = calc_Fsz(data_180, c)
    Cc = calc_Cc(data_180, c)
    pressure_region = calc_pressure_region(data_180.section, data_180.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data_180.section, data_180.Coords, c)
    phi = assump.phif(data_180.fy, data_180.Es, min(es))
    M_neg = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    return M_pos, M_neg


calc_P = partial(PMM.calc_P, angle=0)


def calc_c_max(data: DesignData) -> np.float32:
    ety = data.fy / data.Es
    _, _, _, maxy = data.section.bounds
    miny_rebar = np.min([point.y for point in data.Coords])
    dt = maxy - miny_rebar
    return dt / (1-(ety/assump.ecu))


def calc_c(data: DesignData) -> float:
    def _optim_c(x):
        return calc_P(data, x)
    c_max = calc_c_max(data)
    return root_scalar(_optim_c, bracket=[0.0001, c_max]).root


def calc_Mn(data: DesignData) -> Tuple[np.float32, np.float32]:
    c = calc_c(data)
    return calc_M(data, c)


def calc_M_ratio(data: DesignData, Mux: float) -> np.float32:
    Mn_pos, Mn_neg = calc_Mn(data)
    Mn = Mn_pos if Mux>=0 else Mn_neg
    return Mux/Mn

# ! min and max reinforcement must add
def calc_percent(data: DesignData, Mux: float) -> np.float32: 
    def _optim_percent(x):
        data_percent = set_As_percent(data, x)
        return calc_M_ratio(data_percent, Mux) - 1
    data_one_percent = set_As_percent(data, 1)
    data_eight_percent = set_As_percent(data, 8)
    if calc_M_ratio(data_one_percent, Mux) <= 1: 
        output_percent = 1
    elif calc_M_ratio(data_eight_percent, Mux) > 1:
        raise ValueError("section is weak")
    else:
        output_percent = root_scalar(_optim_percent, bracket=[1, 8]).root
    return np.float32(output_percent)