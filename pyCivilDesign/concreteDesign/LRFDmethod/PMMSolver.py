from dataclasses import dataclass, field
from math import sqrt
from typing import Tuple, List
from matplotlib import pyplot as plt

from numpy.typing import NDArray
import numpy as np

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

from scipy.optimize import least_squares, minimize, root

from pyCivilDesign.concreteDesign.designProps import Assumptions, defaultAssumption, DesignData


def setAs(data: DesignData, As: NDArray[np.float32]) -> DesignData:
    return DesignData(section=data.section, bw=data.bw, d=data.d, fy= data.fy, 
                      fyt=data.fyt, fc=data.fc, Coords=data.Coords, As=As, Es=data.Es)


def setAsPercent(data: DesignData, percent: float) -> DesignData:
    totalAs = data.section.area * (percent/100)
    return setAs(data, np.array([totalAs/ len(data.As) for i in range(len(data.As))]))


def AsPercent(data: DesignData) -> np.float32:
    return (np.sum(data.As)/data.section.area)*100


def onePercentData(data: DesignData) -> DesignData:
    return setAsPercent(data, 1)


def eightPercentData(data: DesignData) -> DesignData:
    return setAsPercent(data, 8)


def P0(data: DesignData) -> np.float32:
    return 0.65 * (0.85 * data.fc * (data.section.area - sum(data.As))\
                   + sum(data.As)*data.fy)


def PnMax(data: DesignData) -> np.float32:
    return 0.8 * P0(data)


def PtMax(data: DesignData) -> np.float32:
    return -0.9 * (sum(data.As) * data.fy)


def Alpha(data: DesignData, Mx: float, My: float,
          assump:Assumptions=defaultAssumption) -> np.float32:
    alpha = np.degrees(np.arctan(abs(My/Mx))) if Mx != 0 else 90
    alpha = alpha if (Mx>=0 and My>=0) else 180-alpha if (Mx<0 and My>0) else \
        180+ alpha if (Mx<0 and My<0) else 360-alpha if (Mx>0 and My<0) else alpha
    return np.float32(alpha)


def OptimAngle(x, *args):
    _angle = x[0]
    _P = args[0]
    _alpha = args[1]
    _data = args[2]
    _assump = args[3]
    return abs(CalcMn(_data, _P, _angle, _assump)[3]-_alpha)


def AngleFromForces(data: DesignData, P: float, Mx: float, My: float,
                    assump: Assumptions=defaultAssumption) -> np.float32:
    alpha = Alpha(data, Mx, My, assump)
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<=180 else\
         270 if 180<alpha<=270 else 360
    output = least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)),
                           args=(P, alpha, data, assump))
    return output.x[0]


def AngleFromAlpha(data:DesignData, P:float, alpha:float,
                   assump:Assumptions=defaultAssumption) -> np.float32:
    lAngle = 0 if 0<=alpha<=90 else 90 if 90<alpha<=180 else\
         180 if 180<alpha<=270 else 270
    uAngle = 90 if 0<=alpha<=90 else 180 if 90<alpha<= 180 else\
         270 if 180<alpha<=270 else 360
    output = least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)),
                           args=(P, alpha, data, assump))
    return output.x[0]


def rotateSection(data: DesignData, angle: float) -> Polygon:
    return rotate(data.section, angle, origin=Point([0, 0]))


def rotateRebarCoords(data: DesignData, angle: float) -> NDArray[np.float32]:
    points: List[Point] = [rotate(Point(coord), angle, origin=Point([0, 0])) \
                           for coord in data.Coords]
    return np.array([(point.x, point.y) for point in points])


def NeutralAxis(data: DesignData, c: float, angle: float) -> NDArray[np.float32]:
    rSection: Polygon = rotateSection(data, angle)
    minx, _, maxx, maxy = rSection.bounds
    line: LineString = LineString([(maxx+10, maxy-c), (minx-10, maxy-c)])
    return np.array(line.coords)


def NeutralRegion(data: DesignData, c: float, angle: float) -> Polygon:
    rSection = rotateSection(data, angle)
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-c),
                      (minx-10, maxy-c), (minx-10, maxy)])
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    unioned = rSection.boundary.union(NL)
    NR: List[Polygon] = [poly for poly in polygonize(unioned)\
                         if poly.representative_point().within(topArea)]
    return NR[0]


def MaxPressurePoint(data: DesignData, c: float, angle: float) -> np.float32:
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    NR: Polygon = NeutralRegion(data, c, angle)
    return np.max([NL.distance(Point(p)) for p in list(NR.exterior.coords)])


def PressureAxis(data: DesignData, c: float, angle: float,
                 assump:Assumptions=defaultAssumption) -> NDArray[np.float32]:
    MaxPrPoint: np.float32 = MaxPressurePoint(data, c, angle)
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    PL: LineString = NL.parallel_offset(distance=MaxPrPoint*(1-assump.beta1(data)),
                                        side="right")
    return np.array(PL.coords)


def PressureRegion(data: DesignData, c: float, angle: float,
                   assump:Assumptions=defaultAssumption) -> Polygon:
    rSection: Polygon = rotateSection(data, angle)
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-(0.85*c)),
                      (minx-10, maxy-(0.85*c)), (minx-10, maxy)])
    PL: LineString = LineString(PressureAxis(data, c, angle, assump))
    unioned = rSection.boundary.union(PL)
    PR: List[Polygon] = [poly for poly in polygonize(unioned) if \
          poly.representative_point().within(topArea)]
    return PR[0]


def es(data: DesignData, c: float, angle: float, 
       assump:Assumptions=defaultAssumption) -> NDArray[np.float32]:
    rCoords: NDArray[np.float32] = rotateRebarCoords(data, angle)
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    MaxPrPoint: np.float32 = MaxPressurePoint(data, c, angle)
    NR: Polygon = NeutralRegion(data, c, angle)
    esSign = [1 if NR.contains(Point(point)) else -1 for point in rCoords]
    return np.array([((esSign[i]*NL.distance(Point(rCoords[i])))/MaxPrPoint)*assump.ecu \
        for i in range(len(rCoords))])


def ec(data: DesignData, c: float, angle: float, point: Tuple[float, float],
       assump:Assumptions=defaultAssumption) -> np.float32:
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    MaxPrPoint: np.float32 = MaxPressurePoint(data, c, angle)
    NR: Polygon = NeutralRegion(data, c, angle)
    ecSign = 1 if NR.contains(Point(point)) else -1
    return ((ecSign * NL.distance(Point(point)))/MaxPrPoint)*assump.ecu


def fs(data: DesignData, c: float, angle: float,
       assump:Assumptions) -> NDArray[np.float32]:
    _es = es(data, c, angle, assump)
    _fs = np.minimum(np.abs(_es) * data.Es, data.fy)
    return np.copysign(_fs, _es)


def Fs(data: DesignData, c: float, angle: float, 
       assump:Assumptions=defaultAssumption) -> NDArray[np.float32]:
    _fs = fs(data, c, angle, assump)
    return np.where(data.As*_fs <= 0, data.As*_fs, data.As*(_fs-0.85*data.fc))


def Cc(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> np.float32:
    _Cc = 0.85 * data.fc * PressureRegion(data, c, angle, assump).area
    return _Cc


def Fsz(data: DesignData, c: float, angle: float, 
        assump:Assumptions=defaultAssumption) -> Tuple[np.float32, np.float32]:
    _Fs = Fs(data, c, angle, assump)
    xCoords = [point.x for point in data.Coords] #data.Coords[:,0]
    yCoords = [point.y for point in data.Coords] #data.Coords[:,1]
    return np.sum(_Fs * xCoords), np.sum(_Fs * yCoords)


def M(data: DesignData, c: float, angle: float, 
      assump:Assumptions=defaultAssumption, 
      IsPhi: bool = True) -> Tuple[np.float32, np.float32, np.float32]:
    signX = 1 if (0<=angle<=90) or (270<=angle<=360) else -1
    signY = 1 if (0<=angle<=180) else -1
    _Fszx, _Fszy = Fsz(data, c, angle, assump)
    _Cc = Cc(data, c, angle, assump)
    rPr = rotate(PressureRegion(data, c, angle, assump), -angle, Point([0, 0]))
    _zcy = abs(rPr.centroid.y)
    _zcx = abs(rPr.centroid.x)
    _es = es(data, c, angle, assump)
    _phi = assump.phif(data, min(_es))
    _Mx = _phi*(_Cc*_zcy + abs(_Fszy)) if IsPhi else _Cc*_zcy + abs(_Fszy)
    _My = _phi*(_Cc*_zcx + abs(_Fszx)) if IsPhi else _Cc*_zcx + abs(_Fszx)
    _M = pow((_Mx)**2+(_My)**2, 0.5)
    return _M, signX*_Mx, signY*_My


def P(data: DesignData, c: float, angle: float, 
      assump:Assumptions=defaultAssumption, IsPhi: bool = True) -> np.float32:
    _Fs = Fs(data, c, angle, assump)
    _Cc = Cc(data, c, angle, assump) 
    _es = es(data, c, angle, assump)
    _phi = assump.phif(data, min(_es))
    return _phi*(_Cc+sum(_Fs)) if IsPhi else _Cc+sum(_Fs)


def OptimF(x, *args):
    _c = x[0]
    _P = args[0]
    _angle = args[1]
    _data = args[2]
    _assump = args[3]
    return abs(_P - P(_data, _c, _angle, _assump))


def C(data: DesignData, P: float, angle: float, 
      assump:Assumptions=defaultAssumption) -> np.float32:
    _, miny, _, maxy = rotateSection(data, angle).bounds
    return least_squares(OptimF, ((maxy-miny)*2.5), 
                         bounds=((0.00001), (5*(maxy-miny))), 
                         args=(P, angle, data, assump)).x[0]


def OptimMaxM(x, *args):
    _P = x[0]
    _angle = args[0]
    _data = args[1]
    _assump = args[2]
    rSection = rotateSection(_data, _angle)
    _, miny, _, maxy = rSection.bounds
    _c = least_squares(OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
                       args=(_P, _angle, _data, _assump)).x[0]
    return -M(_data, _c, _angle, _assump)[0]


def CalcPMmax(data: DesignData, angle: float, 
              assump:Assumptions=defaultAssumption) -> Tuple[np.float32, np.float32]:
    output = minimize(OptimMaxM, ((P0(data)/2)), method="L-BFGS-B", args=(angle, data, assump))
    return output.x[0], -output.fun


def OptimM(x, *args):
    _P = x[0]
    _angle = args[0]
    _e0 = args[1]
    _data = args[2]
    _assump = args[3]
    rSection = rotateSection(_data, _angle)
    _, miny, _, maxy = rSection.bounds
    #[ ] TODO: must be faster
    # _c = least_squares(OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), 
                    #    args=(_P, _angle, _data, _assump)).x[0]
    _c = root(OptimF, ((maxy-miny)/2), args=(_P, _angle, _data, _assump)).x[0]
    return abs(M(_data, _c, _angle, _assump)[0]/_P - _e0)


#* c, angle, alpha, es, fs, Fs, Cc, P, M, Mx, My, ratio
def CalcPMRatio(data: DesignData, P: float, Mx: float, My: float, 
                assump:Assumptions=defaultAssumption) -> np.float32:
    angle = AngleFromForces(data, P, Mx, My, assump)
    M = pow(Mx**2 + My**2, 0.5)
    #[ ] TODO: must be faster
    # _P = least_squares(OptimM, ((P0(data)+PtMax(data))/2), 
    #                    bounds=((PtMax(data)), (P0(data))), 
    #                    args=(angle, M/P, data, assump)).x[0] if P != 0 else 1
    _P = root(OptimM,(P0(data)+PtMax(data))/2, args=(angle, M/P, data, assump)).x[0] if P!=0 else 1
    _M = CalcMn(data, _P, angle, assump)[0] # type: ignore
    return np.float32(P/_P if (P/_P) != 0 else M/_M)


#* c, angle, alpha, es, fs, Fs, Cc, P, M, Mx, My 
def CalcMn(data: DesignData, P: float, angle: float, 
           assump:Assumptions=defaultAssumption) -> Tuple[np.float32, np.float32,
                                                          np.float32, np.float32]:
    c = C(data, P, angle, assump)
    _M, Mx, My = M(data, float(c), angle, assump)
    alpha = Alpha(data, Mx, My, assump) # type: ignore
    return _M, Mx, My, alpha


def OptimPercent(x, *args):
    percent = x[0]
    P = args[0]
    Mx = args[1]
    My = args[2]
    data = args[3]
    assump = args[4]
    Data = setAsPercent(data, percent)
    return CalcPMRatio(Data, P, Mx, My, assump) - 1


#* c, angle, alpha, es, fs, Fs, Cc, P, M, Mx, My, percent
def CalcAsPercent(data:DesignData, P: float, Mx: float, My: float, # [x] BUG: not work!
                assump: Assumptions=defaultAssumption) -> np.float32:
    maxPercent = np.float32(10)
    minPercent = np.float32(0.1)
    e = 1
    percent = np.float32(0)
    while abs(e) > 0.01:
        percent = np.float32((maxPercent + minPercent)/2)
        data = setAsPercent(data, percent) # type: ignore
        e = CalcPMRatio(data, P, Mx, My, assump) - 1
        maxPercent = (maxPercent + minPercent)/2 if e<0 else maxPercent
        minPercent = minPercent if e<0 else (maxPercent + minPercent)/2
    return percent
    # return root(OptimPercent, (4,), args=(P, Mx, My, data, assump)).x[0]
