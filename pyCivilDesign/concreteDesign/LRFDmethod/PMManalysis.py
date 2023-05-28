from typing import Tuple

from numpy.typing import NDArray
import numpy as np

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

from scipy.optimize import root_scalar

from pyCivilDesign.concreteDesign.designProps import DesignData
import pyCivilDesign.concreteDesign.LRFDmethod.assumptions as assump


def set_As(data: DesignData, As: NDArray[np.float32]) -> DesignData:
    return DesignData(section=data.section, bw=data.bw, d=data.d, fy= data.fy, 
                      fyt=data.fyt, fc=data.fc, Coords=data.Coords, As=As, Es=data.Es)


def get_As_percent(data: DesignData) -> np.float32:
    return (np.sum(data.As)/data.section.area)*100


def set_As_percent(data: DesignData, percent: float) -> DesignData:
    totalAs = data.section.area * (percent/100)
    return set_As(data, np.array([totalAs/ len(data.As) for i in range(len(data.As))]))


def calc_P0(data: DesignData) -> np.float32:
    return (0.85 * data.fc * (data.section.area - sum(data.As)))\
                   + (sum(data.As)*data.fy)


def calc_phi_P0(data: DesignData) -> np.float32:
    return 0.65 * calc_P0(data)


def calc_Pn_max(data: DesignData) -> np.float32:
    return 0.8 * calc_P0(data)


def calc_Pt_max(data: DesignData) -> np.float32:
    return -0.9 * (sum(data.As) * data.fy)


def calc_alpha(Mx: float, My: float) -> np.float32:
    """angle between Mx & My on M-M chart

    Args:
        Mx (float): x direction moment
        My (float): y direction moment

    Returns:
        np.float32: alpha
    """
    alpha = np.degrees(np.arctan(abs(My/Mx))) if Mx != 0 else 90
    alpha = alpha if (Mx>=0 and My>=0) else 180-alpha if (Mx<0 and My>0) else \
        180+alpha if (Mx<0 and My<0) else 360-alpha if (Mx>0 and My<0) else alpha
    return np.float32(alpha)


def calc_angle(data: DesignData, P: float, Mx: float, My: float) -> np.float32:
    """calculate angle of rotation for section on PMM force

    Args:
        data (DesignData): design data
        P (float): axial force
        Mx (float): x direction moment
        My (float): y direction moment

    Returns:
        np.float32: angle of rotation
    """
    alpha = calc_alpha(Mx, My)
    def _optim_angle(x):
        c = calc_c(data, P, x)
        _, _Mx, _My = calc_M(data, c, x)
        return calc_alpha(_Mx, _My) - alpha # type: ignore
    
    lbound = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    ubound = 90 if 0<=alpha<=90 else 180 if 90<alpha<=180 else\
         270 if 180<alpha<=270 else 360
    return root_scalar(_optim_angle, bracket=[lbound, ubound]).root


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
    return LineString([(maxx+10, maxy-c), (minx-10, maxy-c)]) # type: ignore


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
    rot_pressure_region = rotate(Pr, angle, Point([0, 0])) if angle!=0 else Pr # [] TODO: check sign of angle rotation
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
    return _P


def calc_c_max(data, angle):
    ety = data.fy / data.Es
    rot_section = rotate_section(data.section, angle)
    _, _, _, maxy = rot_section.bounds
    rot_Coords = rotate_rebar_coords(data.Coords, angle)
    miny_rebar = np.min([point.y for point in rot_Coords])
    dt = maxy - miny_rebar
    return dt / (1-(ety/assump.ecu))


def calc_c(data, P, angle):
    def _optim_c(x):
        return calc_P(data, x, angle) - P
    c_max = calc_c_max(data, angle)
    return root_scalar(_optim_c, bracket=[0.0001, c_max]).root


def calc_Pc_list(data, angle, num=20, is_phi=True):
    c_max = calc_c_max(data, angle)
    c_list = np.linspace(1e-6, c_max, num=num, dtype=np.float32)
    P_list = np.array([calc_P(data, _c, angle, is_phi) for _c in c_list], dtype=np.float32)
    return P_list, c_list


def calc_PM_list(data, angle, num=20, is_phi=True):
    P_list, c_list = calc_Pc_list(data, angle, num, is_phi)
    M_list = np.array([calc_M(data, _c, angle, is_phi) for _c in c_list], dtype=np.float32)
    M_list = M_list.reshape(-1, 3)
    return P_list, M_list


def calc_Mn(data, P, angle):
    c = calc_c(data, P, angle)
    return calc_M(data, c, angle)


def calc_PM_max(data: DesignData, angle: float) -> np.float32:
    _, M_list = calc_PM_list(data, angle, num=50)
    return np.max(M_list[:,0])


def calc_PM_ratio(data, P, Mx, My, angle:float|None=None):
    if angle == None: angle = calc_angle(data, P, Mx, My)
    M = pow(Mx**2 + My**2, 0.5)
    e = M/P
    P_list, M_list = calc_PM_list(data, angle)
    M_custom = np.max(M_list[:,0]) * 1.1
    P_custom = M_custom/e

    PM_line = LineString(list(zip(M_list[:,0], P_list)))
    PMx_line = LineString([(0, 0), (M_custom, P_custom)])
    
    inter_point = PM_line.intersection(PMx_line)
    _M = inter_point.x
    _P = inter_point.y
    return P/_P if P!=0 else M/_M


def calc_percent(data, P, Mx, My, ratio=1, num=8):
    angle = calc_angle(data, P, Mx, My)
    def _optim_percent(x):
        data_percent = set_As_percent(data, x)
        c = calc_c(data_percent, P, angle)
        return calc_M(data_percent, c, angle)[0] - pow(Mx**2+My**2, 0.5)
    
    data_one_percent = set_As_percent(data, 1)
    data_eight_percent = set_As_percent(data, 8)
    if calc_PM_ratio(data_one_percent, P, Mx, My, angle) <= 1: 
        output_percent = 1
    elif calc_PM_ratio(data_eight_percent, P, Mx, My, angle) > 1:
        raise ValueError("section is weak")
    else:
        output_percent = root_scalar(_optim_percent, bracket=[1, 8])
    return output_percent

