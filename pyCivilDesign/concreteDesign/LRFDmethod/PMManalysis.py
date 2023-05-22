from typing import Tuple

from numpy.typing import NDArray
import numpy as np

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

from scipy.optimize import least_squares, minimize, root

from pyCivilDesign.concreteDesign.designProps import DesignData
import pyCivilDesign.concreteDesign.LRFDmethod.assumptions as assump


def set_As(data: DesignData, As: NDArray[np.float32]) -> DesignData:
    return DesignData(section=data.section, bw=data.bw, d=data.d, fy= data.fy, 
                      fyt=data.fyt, fc=data.fc, Coords=data.Coords, As=As, Es=data.Es)


def set_As_percent(data: DesignData, percent: float) -> DesignData:
    totalAs = data.section.area * (percent/100)
    return set_As(data, np.array([totalAs/ len(data.As) for i in range(len(data.As))]))


def As_percent(data: DesignData) -> np.float32:
    return (np.sum(data.As)/data.section.area)*100


def calc_P0(data: DesignData) -> np.float32:
    return 0.65 * (0.85 * data.fc * (data.section.area - sum(data.As))\
                   + sum(data.As)*data.fy)


def calc_Pn_max(data: DesignData) -> np.float32:
    return 0.8 * calc_P0(data)


def calc_Pt_max(data: DesignData) -> np.float32:
    return -0.9 * (sum(data.As) * data.fy)


def calc_alpha(Mx: float, My: float) -> np.float32:
    alpha = np.degrees(np.arctan(abs(My/Mx))) if Mx != 0 else 90
    alpha = alpha if (Mx>=0 and My>=0) else 180-alpha if (Mx<0 and My>0) else \
        180+ alpha if (Mx<0 and My<0) else 360-alpha if (Mx>0 and My<0) else alpha
    return np.float32(alpha)


def _optim_angle(x, *args):
    _angle = x[0]
    _P = args[0]
    _alpha = args[1]
    _data = args[2]
    return abs(calc_Mn(_data, _P, _angle)[3]-_alpha)


def calc_angle_from_forces(data: DesignData, P: float, Mx: float, My: float) -> np.float32:
    alpha = calc_alpha(Mx, My)
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<=180 else\
         270 if 180<alpha<=270 else 360
    output = least_squares(_optim_angle, (alpha), bounds=((lAngle), (uAngle)),
                           args=(P, alpha, data))
    return output.x[0]


def calc_angle_from_alpha(data:DesignData, P:float, alpha:float) -> np.float32:
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<= 180 else\
         270 if 180<alpha<=270 else 360
    output = least_squares(_optim_angle, (alpha), bounds=((lAngle), (uAngle)),
                           args=(P, alpha, data))
    return output.x[0]


def rotate_section(section: Polygon, angle: float) -> Polygon:
    return  rotate(section, angle, origin=Point([0, 0])) if angle!=0 else section


def rotate_rebar_coords(coords: NDArray[Point], angle: float) -> NDArray[Point]:
    return np.array([rotate(coord, angle, origin=Point([0, 0])) \
                           for coord in coords]) if angle !=0 else coords


def calc_neutral_axis(section: Polygon, c: float, angle: float) -> LineString:
    rot_section = rotate_section(section, angle) if angle!=0 else section 
    minx, _, maxx, maxy = rot_section.bounds
    return LineString([(maxx+10, maxy-c), (minx-10, maxy-c)])


def calc_neutral_region(section: Polygon, c: float, angle: float) -> Polygon:
    rot_section = rotate_section(section, angle) if angle!=0 else section
    minx, _, maxx, maxy = rot_section.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-c),
                      (minx-10, maxy-c), (minx-10, maxy)])
    NL = calc_neutral_axis(section, c, angle)
    unioned = rot_section.boundary.union(NL)
    NR = [poly for poly in polygonize(unioned) if \
        poly.representative_point().within(topArea)]
    return NR[0]


def calc_max_pressure_point(section: Polygon, c: float, angle: float) -> np.float32:
    NL = calc_neutral_axis(section, c, angle)
    NR = calc_neutral_region(section, c, angle)
    return np.max([NL.distance(Point(p)) for p in list(NR.exterior.coords)])


def calc_pressure_axis(section: Polygon, fc: np.float32, c: float, angle: float) -> LineString:
    MaxPrPoint = calc_max_pressure_point(section, c, angle)
    NL = calc_neutral_axis(section, c, angle)
    return NL.parallel_offset(distance=MaxPrPoint*(1-assump.beta1(fc)), side="right")


def calc_pressure_region(section: Polygon, fc: np.float32, c: float, angle: float) -> Polygon:
    rot_section = rotate_section(section, angle) if angle!=0 else section
    minx, _, maxx, maxy = rot_section.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-(0.85*c)),
                      (minx-10, maxy-(0.85*c)), (minx-10, maxy)])
    PL = calc_pressure_axis(section, fc, c, angle)
    unioned = rot_section.boundary.union(PL)
    PR = [poly for poly in polygonize(unioned) if \
          poly.representative_point().within(topArea)]
    return PR[0]


def calc_es(section: Polygon, coords:NDArray[Point], c: float, angle: float) -> NDArray[np.float32]:
    rot_coords = rotate_rebar_coords(coords, angle) if angle!=0 else coords
    NL = calc_neutral_axis(section, c, angle)
    MaxPrPoint = calc_max_pressure_point(section, c, angle)
    NR = calc_neutral_region(section, c, angle)
    esSign = np.array([1 if NR.contains(point) else -1 for point in rot_coords])
    return np.array([((esSign[i]*NL.distance(rot_coords[i]))/MaxPrPoint)*assump.ecu \
        for i in range(len(rot_coords))])


def calc_ec(section: Polygon, c: float, angle: float, point: Point) -> np.float32:
    NL = calc_neutral_axis(section, c, angle)
    MaxPrPoint = calc_max_pressure_point(section, c, angle)
    NR = calc_neutral_region(section, c, angle)
    ecSign = 1 if NR.contains(point) else -1
    return ((ecSign * NL.distance(point))/MaxPrPoint)*assump.ecu


def calc_fs(data: DesignData, c: float, angle: float) -> NDArray[np.float32]:
    es = calc_es(data.section, data.Coords, c, angle)
    fs = np.minimum(np.abs(es) * data.Es, data.fy)
    return np.copysign(fs, es)


def calc_Fs(data: DesignData, c: float, angle: float) -> NDArray[np.float32]:
    fs = calc_fs(data, c, angle)
    return np.where(data.As*fs <= 0, data.As*fs, data.As*(fs-0.85*data.fc))


def calc_Cc(data: DesignData, c: float, angle: float) -> np.float32:
    return 0.85 * data.fc * calc_pressure_region(data.section, data.fc, c, angle).area


def calc_Fsz(data: DesignData, c: float, angle: float) -> Tuple[np.float32, np.float32]:
    Fs = calc_Fs(data, c, angle)
    xCoords = np.array([point.x for point in data.Coords])
    yCoords = np.array([point.y for point in data.Coords])
    return np.sum(Fs * xCoords), np.sum(Fs * yCoords)


def calc_M(data: DesignData, c: float, angle: float, IsPhi: bool = True) \
        -> Tuple[np.float32, np.float32, np.float32]:
    signX = 1 if (0<=angle<=90) or (270<=angle<=360) else -1
    signY = 1 if (0<=angle<=180) else -1
    Fszx, Fszy = calc_Fsz(data, c, angle)
    Cc = calc_Cc(data, c, angle)
    Pr = calc_pressure_region(data.section, data.fc, c, angle)
    rot_pressure_region = rotate(Pr, -angle, Point([0, 0])) if angle!=0 else Pr
    zcy = abs(rot_pressure_region.centroid.y)
    zcx = abs(rot_pressure_region.centroid.x)
    es = calc_es(data.section, data.Coords, c, angle)
    phi = assump.phif(data.fy, data.Es, min(es))
    Mx = phi*(Cc*zcy + abs(Fszy)) if IsPhi else Cc*zcy + abs(Fszy)
    My = phi*(Cc*zcx + abs(Fszx)) if IsPhi else Cc*zcx + abs(Fszx)
    M = pow(Mx**2+My**2, 0.5)
    return M, signX*Mx, signY*My


def calc_P(data: DesignData, c: float, angle: float, IsPhi: bool = True) -> np.float32:
    Fs = calc_Fs(data, c, angle)
    Cc = calc_Cc(data, c, angle) 
    es = calc_es(data.section, data.Coords, c, angle)
    phi = assump.phif(data.fy, data.Es, min(es))
    return phi*(Cc+sum(Fs)) if IsPhi else Cc+sum(Fs)


def _optim_F(x, *args):
    c = x[0]
    P = args[0]
    angle = args[1]
    data = args[2]
    return P - calc_P(data, c, angle)


def calc_c(data: DesignData, P: float, angle: float) -> np.float32:
    rot_section = rotate_section(data.section, angle) if angle!=0 else data.section
    _, miny, _, maxy = rot_section.bounds
    return least_squares(_optim_F, ((maxy-miny)*2.5), 
                         bounds=((0.00001), (5*(maxy-miny))), 
                         args=(P, angle, data)).x[0]


def _optim_max_M(x, *args):
    _P = x[0]
    _angle = args[0]
    _data = args[1]
    rot_section = rotate_section(_data.section, _angle) if _angle!=0 else _data.section
    _, miny, _, maxy = rot_section.bounds
    _c = least_squares(_optim_F, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
                       args=(_P, _angle, _data)).x[0]
    return -calc_M(_data, _c, _angle)[0]


def calc_PM_max(data: DesignData, angle: float) -> Tuple[np.float32, np.float32]:
    output = minimize(_optim_max_M, ((calc_P0(data)/2)), method="L-BFGS-B", args=(angle, data))
    return output.x[0], -output.fun


def _optim_M(x, *args):
    _P = x[0]
    _angle = args[0]
    _e0 = args[1]
    _data = args[2]
    rot_section = rotate_section(_data.section, _angle) if _angle!=0 else _data.section
    _, miny, _, maxy = rot_section.bounds
    #[ ] TODO: must be faster
    # _c = least_squares(OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
                    #    args=(_P, _angle, _data, _assump)).x[0]
    _c = root(_optim_F, ((maxy-miny)/2), args=(_P, _angle, _data)).x[0]
    return abs(calc_M(_data, _c, _angle)[0]/_P - _e0)


def calc_PM_ratio(data: DesignData, P: float, Mx: float, My: float) -> np.float32:
    angle = calc_angle_from_forces(data, P, Mx, My)
    M = pow(Mx**2 + My**2, 0.5)
    #[ ] TODO: must be faster
    # _P = least_squares(OptimM, ((P0(data)+PtMax(data))/2), 
    #                    bounds=((PtMax(data)), (P0(data))), 
    #                    args=(angle, M/P, data, assump)).x[0] if P != 0 else 1
    _P = root(_optim_M,(calc_P0(data)+calc_Pt_max(data))/2, args=(angle, M/P, data)).x[0] if P!=0 else 1
    _M = calc_Mn(data, _P, angle)[0] # type: ignore
    return np.float32(P/_P if (P/_P) != 0 else M/_M)


def calc_Mn(data: DesignData, P: float, angle: float) -> Tuple[np.float32, np.float32,
                                                          np.float32, np.float32]:
    c = calc_c(data, P, angle)
    _M, Mx, My = calc_M(data, float(c), angle)
    alpha = calc_alpha(Mx, My) # type: ignore
    return _M, Mx, My, alpha


def _optim_percent(x, *args):
    percent = x[0]
    P = args[0]
    Mx = args[1]
    My = args[2]
    data = args[3]
    Data = set_As_percent(data, percent)
    return calc_PM_ratio(Data, P, Mx, My) - 1


def calc_As_percent(data:DesignData, P: float, Mx: float, My: float) -> np.float32:
    maxPercent = np.float32(10)
    minPercent = np.float32(0.1)
    e = 1
    percent = np.float32(0)
    while abs(e) > 0.01:
        percent = np.float32((maxPercent + minPercent)/2)
        data = set_As_percent(data, percent) # type: ignore
        e = calc_PM_ratio(data, P, Mx, My) - 1
        maxPercent = (maxPercent + minPercent)/2 if e<0 else maxPercent
        minPercent = minPercent if e<0 else (maxPercent + minPercent)/2
    return percent
    # return root(OptimPercent, (4,), args=(P, Mx, My, data, assump)).x[0]
