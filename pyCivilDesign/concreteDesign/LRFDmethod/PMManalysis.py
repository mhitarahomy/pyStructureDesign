from typing import Tuple

from numpy.typing import NDArray
import numpy as np

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

import scipy.optimize as opt

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
        180+alpha if (Mx<0 and My<0) else 360-alpha if (Mx>0 and My<0) else alpha
    return np.float32(alpha)


def _optim_angle(x, *args):
    _angle = x[0]
    _P = args[0]
    _alpha = args[1]
    _data = args[2]
    _, _Mx, _My = calc_Mn(_data, _P, _angle)
    return abs(calc_alpha(_Mx, _My)-_alpha) # type: ignore


def calc_angle_from_forces(data: DesignData, P: float, Mx: float, My: float) -> np.float32:
    alpha = calc_alpha(Mx, My)
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<=180 else\
         270 if 180<alpha<=270 else 360
    output = opt.root(fun=_optim_angle, method="lm", x0=(alpha, ), args=(P, alpha, data))
    # output = opt.least_squares(_optim_angle, (alpha), bounds=((lAngle), (uAngle)),
    #                        args=(P, alpha, data))
    return output.x[0]


def calc_angle_from_alpha(data:DesignData, P:float, alpha:float) -> np.float32:
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<= 180 else\
         270 if 180<alpha<=270 else 360
    output = opt.least_squares(_optim_angle, (alpha), bounds=((lAngle), (uAngle)),
                           args=(P, alpha, data))
    return output.x[0]


def rotate_section(section: Polygon, angle: float) -> Polygon:
    """rotate cross section

    Args:
        section (Polygon): cross section shape
        angle (float): angle of rotate, degree

    Returns:
        Polygon: rotated cross section
    """
    return  rotate(section, angle, origin=Point([0, 0])) if angle!=0 else section


def rotate_rebar_coords(coords: NDArray[Point], angle: float) -> NDArray[Point]:
    """rotate rebar point coordinates

    Args:
        coords (NDArray[Point]): list of point coordinates
        angle (float): angle of rotate, degree

    Returns:
        NDArray[Point]: list of rotated point coordinates
    """
    return np.array([rotate(coord, angle, origin=Point([0, 0])) \
                           for coord in coords]) if angle !=0 else coords


def calc_neutral_axis(section: Polygon, c: float, angle: float) -> LineString:
    """calculate neutral axis for section that rotated

    Args:
        section (Polygon): cross section shape
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate section, degree

    Returns:
        LineString: netural line
    """
    c = c if c>=0.0001 else 0.0001
    rot_section = rotate_section(section, angle)
    minx, _, maxx, maxy = rot_section.bounds
    return LineString([(maxx+10, maxy-c), (minx-10, maxy-c)])


def calc_neutral_region(section: Polygon, c: float, angle: float) -> Polygon:
    """calculate neutral region for section that rotated

    Args:
        section (Polygon): cross section shape
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        Polygon: neutral region shape
    """
    c = c if c>=0.0001 else 0.0001
    rot_section = rotate_section(section, angle)
    minx, _, maxx, maxy = rot_section.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-c),
                      (minx-10, maxy-c), (minx-10, maxy)])
    neutral_line = calc_neutral_axis(section, c, angle)
    unioned = rot_section.boundary.union(neutral_line)
    neutral_region = [poly for poly in polygonize(unioned) if \
        poly.representative_point().within(topArea)]
    return neutral_region[0]


def calc_max_pressure_point(section: Polygon, c: float, angle: float) -> np.float32:
    """calculate maximum distance of neutral region for rotated section

    Args:
        section (Polygon): cross section shape
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        np.float32: maximum distance of pressure point
    """
    c = c if c>=0.0001 else 0.0001
    neutral_line = calc_neutral_axis(section, c, angle)
    neutral_region = calc_neutral_region(section, c, angle)
    return np.max([neutral_line.distance(Point(p)) for p in list(neutral_region.exterior.coords)])


def calc_pressure_axis(section: Polygon, fc: np.float32, c: float, angle: float) -> LineString:
    """calculate pressure line based on ACI code

    Args:
        section (Polygon): cross section shape
        fc (np.float32): specified compressive strength of concrete, MPa
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        LineString: pressure line
    """
    c = c if c>=0.0001 else 0.0001
    max_pressure_point = calc_max_pressure_point(section, c, angle)
    neutral_line = calc_neutral_axis(section, c, angle)
    return neutral_line.parallel_offset(distance=max_pressure_point*
                                        (1-assump.beta1(fc)), side="right")


def calc_pressure_region(section: Polygon, fc: np.float32, c: float, angle: float) -> Polygon:
    """calculate pressure region based on ACI code

    Args:
        section (Polygon): cross section shape
        fc (np.float32): specified compressive strength of concrete, MPa
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        Polygon: pressure region
    """
    c = c if c>=0.0001 else 0.0001
    rot_section = rotate_section(section, angle) if angle!=0 else section
    minx, _, maxx, maxy = rot_section.bounds
    top_area = Polygon([(maxx+10, maxy), (maxx+10, maxy-(0.85*c)),
                      (minx-10, maxy-(0.85*c)), (minx-10, maxy)])
    pressure_line = calc_pressure_axis(section, fc, c, angle)
    unioned = rot_section.boundary.union(pressure_line)
    pressure_region = [poly for poly in polygonize(unioned) if \
          poly.representative_point().within(top_area)]
    return pressure_region[0]


def calc_es(section: Polygon, coords:NDArray[Point], c: float, angle: float) -> NDArray[np.float32]:
    """calculate strain of rebars for rotated section

    Args:
        section (Polygon): cross section shape
        coords (NDArray[Point]): coordinates of rebars, mm
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        NDArray[np.float32]: _description_
    """
    rot_coords = rotate_rebar_coords(coords, angle) if angle!=0 else coords
    neutral_line = calc_neutral_axis(section, c, angle)
    max_pressure_point = calc_max_pressure_point(section, c, angle)
    neutral_region = calc_neutral_region(section, c, angle)
    es_sign = np.array([1 if neutral_region.contains(point) else -1 for point in rot_coords])
    return np.array([((es_sign[i]*neutral_line.distance(rot_coords[i]))/max_pressure_point)*assump.ecu \
        for i in range(len(rot_coords))])


def calc_ec(section: Polygon, c: float, angle: float, point: Point) -> np.float32:
    """calculate strain of concrete for each point on section

    Args:
        section (Polygon): cross section shape
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree
        point (Point): one point on section

    Returns:
        np.float32: strain of concrete
    """
    neutral_line = calc_neutral_axis(section, c, angle)
    max_pressure_point = calc_max_pressure_point(section, c, angle)
    neutral_region = calc_neutral_region(section, c, angle)
    ec_sign = 1 if neutral_region.contains(point) else -1
    return ((ec_sign * neutral_line.distance(point))/max_pressure_point)*assump.ecu


def calc_fs(data: DesignData, c: float, angle: float) -> NDArray[np.float32]:
    """calculate stress for each rebar for rotated section

    Args:
        data (DesignData): design data
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        NDArray[np.float32]: stress of rebars
    """
    es = calc_es(data.section, data.Coords, c, angle)
    fs = np.minimum(np.abs(es) * data.Es, data.fy)
    return np.copysign(fs, es)


def calc_Fs(data: DesignData, c: float, angle: float) -> NDArray[np.float32]:
    """calculate force for each rebar for rotated section

    Args:
        data (DesignData): design data
        c (float): distance from extreme compression fiber to neutral axis, mm
        angle (float): angle of rotate, degree

    Returns:
        NDArray[np.float32]: force of rebars
    """
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
    # c = c if c>=0.0001 else 0.0001
    Fs = calc_Fs(data, c, angle)
    Cc = calc_Cc(data, c, angle) 
    es = calc_es(data.section, data.Coords, c, angle)
    phi = assump.phif(data.fy, data.Es, min(es))
    _P = phi*(Cc+sum(Fs)) if IsPhi else Cc+sum(Fs)
    return _P if c>=0.0001 else 0


# def _optim_F(x, *args):
#     c = x[0]
#     P = args[0]
#     angle = args[1]
#     data = args[2]
#     return P - calc_P(data, c, angle)


def calc_c(data: DesignData, P: float, angle: float):
    rot_section = rotate_section(data.section, angle)
    _, miny, _, maxy = rot_section.bounds

    def _optim_F(x):
        c = x[0]
        return P - calc_P(data, c, angle)

    result = opt.root(fun=_optim_F, x0=((maxy-miny)/2, ))
    return result.x[0]
    # return result.x[0]
    # return opt.least_squares(_optim_F, ((maxy-miny),), 
    #                      bounds=((0.00001), (5*(maxy-miny))), 
    #                      args=(P, angle, data)).x[0]


def _optim_max_M(x, *args):
    _P = x[0]
    _angle = args[0]
    _data = args[1]
    _c = calc_c(_data, _P, _angle)
    # _c = least_squares(_optim_F, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
    #                    args=(_P, _angle, _data)).x[0]
    return -calc_M(_data, _c, _angle)[0] # type: ignore


def calc_PM_max(data: DesignData, angle: float) -> Tuple[np.float32, np.float32]:
    output = opt.minimize(_optim_max_M, ((calc_P0(data)/2)), method="L-BFGS-B", args=(angle, data))
    return output.x[0], -output.fun


def _optim_M(x, *args):
    _P = x[0]
    _angle = args[0]
    _e0 = args[1]
    _data = args[2]
    _c = calc_c(_data, _P, _angle)
    # _c = opt.least_squares(_optim_F, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
    #                    args=(_P, _angle, _data)).x[0]
    return abs(calc_M(_data, _c, _angle)[0]/_P - _e0) # type: ignore


def calc_PM_ratio(data: DesignData, P: float, Mx: float, My: float) -> np.float32:
    angle = calc_angle_from_forces(data, P, Mx, My)
    M = pow(Mx**2 + My**2, 0.5)
    # result = opt.root(fun=_optim_M, method="lm", x0=((calc_P0(data)+calc_Pt_max(data))/2,), args=(angle, M/P, data))
    # _P = result.x[0]
    _P = opt.least_squares(_optim_M, ((calc_P0(data)+calc_Pt_max(data))/2), 
                       bounds=((calc_Pt_max(data)), (calc_P0(data))), 
                       args=(angle, M/P, data)).x[0] if P != 0 else 1
    _M = calc_Mn(data, _P, angle)[0] # type: ignore
    return np.float32(P/_P if (P/_P) != 0 else M/_M)


def calc_Mn(data: DesignData, P: float, angle: float) -> Tuple[np.float32, np.float32, np.float32]:
    c = calc_c(data, P, angle)
    M, Mx, My = calc_M(data, float(c), angle)
    return M, Mx, My


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
    # return least_squares(_optim_percent, (4,), args=(P, Mx, My, data)).x[0]
