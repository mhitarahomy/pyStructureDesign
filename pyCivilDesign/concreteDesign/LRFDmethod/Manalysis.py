from functools import partial

import numpy as np

from scipy.optimize import root #,least_squares, minimize

from pyCivilDesign.concreteDesign.designProps import DesignData
import pyCivilDesign.concreteDesign.LRFDmethod.PMManalysis as PMM
import pyCivilDesign.concreteDesign.LRFDmethod.assumptions as assump

# def setAs(data: DesignData, As: NDArray[np.float32]) -> DesignData:
#     return DesignData(section=data.section, bw=data.bw, d=data.d, fy= data.fy, 
#                       fyt=data.fyt, fc=data.fc, Coords=data.Coords, As=As, Es=data.Es)


# def setAsPercent(data: DesignData, percent: float) -> DesignData:
#     totalAs = data.section.area * (percent/100)
#     return setAs(data, np.array([totalAs/ len(data.As) for i in range(len(data.As))]))


# def AsPercent(data: DesignData) -> np.float32:
#     return (np.sum(data.As)/data.section.area)*100


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


def calc_M(data: DesignData, c: float, IsPhi: bool = True) -> np.float32:
    Fszy = calc_Fsz(data, c)
    Cc = calc_Cc(data, c)
    pressure_region = calc_pressure_region(data.section, data.fc, c)
    zcy = abs(pressure_region.centroid.y)
    es = calc_es(data.section, data.Coords, c)
    phi = assump.phif(data.fy, data.Es, min(es))
    Mx = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    return Mx


calc_P = partial(PMM.calc_P, angle=0)


def _optim_F(x, *args):
    c = x[0]
    data = args[0]
    return calc_P(data, c)


def calc_c(data: DesignData) -> np.float32:
    _, miny, _, maxy = data.section.bounds
    return root(_optim_F, ((maxy-miny)/2,), args=(data,)).x[0]
