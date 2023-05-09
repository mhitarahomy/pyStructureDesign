from dataclasses import dataclass, field
from typing import Protocol, Tuple, List, Callable
from math import copysign, degrees, atan

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

from scipy.optimize import least_squares, minimize

from pyCivilDesign.sections.section import ListOfPoints
from pyCivilDesign.concreteDesign.designAssumptions import Assumptions, defaultAssumption, DesignData

   
def setAs(data: DesignData, As: List[float]) -> DesignData:
    return DesignData(data.section, data.fy, data.fc, data.Coords, As, data.Es)

def setAsPercent(data: DesignData, percent: float) -> DesignData:
    totalAs = Polygon(data.section).area * (percent/100)
    return setAs(data, [totalAs/ len(data.Coords) for i in range(len(data.As))])

def AsPercent(data: DesignData) -> float:
    return (sum(data.As)/Polygon(data.section).area)*100

def onePercentData(data: DesignData) -> DesignData:
    return setAsPercent(data, 1)

def eightPercentAnalysis(data: DesignData) -> DesignData:
    return setAsPercent(data, 8)

def P0(data: DesignData):
    return 0.65 * (0.85 * data.fc * (Polygon(data.section).area-sum(data.As))+sum(data.As)*data.fy)

def PnMax(data: DesignData):
    return 0.8 * P0(data)

def PtMax(data: DesignData):
    return -0.9 * (sum(data.As) * data.fy)

def Alpha(data: DesignData, Mx: float, My: float, assump:Assumptions=defaultAssumption) -> float:
    alpha = degrees(atan(abs(My/Mx))) if Mx != 0 else 90
    alpha = alpha if (Mx >= 0 and My >= 0) else 180-alpha if (Mx < 0 and My > 0) else 180 + \
        alpha if (Mx < 0 and My < 0) else 360 - \
        alpha if (Mx > 0 and My < 0) else alpha
    return alpha

def OptimAngle(x, *args):
    _angle = x[0]
    _P = args[0]
    _alpha = args[1]
    _data = args[2]
    _assump = args[3]
    return abs(CalcMn(_data, _P, _angle, _assump)[3]-_alpha)

def AngleFromForces(data: DesignData, P: float, Mx: float, My: float, assump: Assumptions=defaultAssumption) -> float:
    alpha = Alpha(data, Mx, My, assump)
    lAngle = 0 if 0 <= alpha <= 90 else 90 if 90 < alpha <= 180 else 180 if 180 < alpha <= 270 else 270
    uAngle = 90 if 0 <= alpha <= 90 else 180 if 90 < alpha <= 180 else 270 if 180 < alpha <= 270 else 360
    return least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)), args=(P, alpha, data, assump)).x[0]

def AngleFromAlpha(data:DesignData, P:float, alpha:float, assump:Assumptions=defaultAssumption) -> float:
    lAngle = 0 if 0 <= alpha <= 90 else 90 if 90 < alpha <= 180 else 180 if 180 < alpha <= 270 else 270
    uAngle = 90 if 0 <= alpha <= 90 else 180 if 90 < alpha <= 180 else 270 if 180 < alpha <= 270 else 360
    return least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)), args=(P, alpha, data, assump)).x[0]

def rotateSection(data: DesignData, angle: float) -> ListOfPoints:
    rsection: Polygon = rotate(Polygon(data.section), angle, origin=Point([0, 0]))
    return ListOfPoints(rsection.exterior.coords)

def rotateRebarCoords(data: DesignData, angle: float) -> ListOfPoints:
    points: List[Point] = [rotate(Point(coord), angle, origin=Point([0, 0])) for coord in data.Coords]
    return ListOfPoints([(point.x, point.y) for point in points])

def NeutralAxis(data: DesignData, c: float, angle: float) -> ListOfPoints:
    rSection: Polygon = Polygon(rotateSection(data, angle))
    minx, _, maxx, maxy = rSection.bounds
    line: LineString = LineString([(maxx+10, maxy-c), (minx-10, maxy-c)])
    return ListOfPoints(line.coords)

def NeutralRegion(data: DesignData, c: float, angle: float) -> ListOfPoints:
    rSection = Polygon(rotateSection(data, angle))
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-c),
                      (minx-10, maxy-c), (minx-10, maxy)])
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    unioned = rSection.boundary.union(NL)
    NR: List[Polygon] = [poly for poly in polygonize(unioned) if poly.representative_point().within(topArea)]
    return ListOfPoints(NR[0].exterior.coords)

def MaxPressurePoint(data: DesignData, c: float, angle: float) -> float:
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    NR: Polygon = Polygon(NeutralRegion(data, c, angle))
    return max([NL.distance(Point(p)) for p in list(NR.exterior.coords)])

def PressureAxis(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> ListOfPoints:
    MaxPrPoint: float = MaxPressurePoint(data, c, angle)
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    PL: LineString = NL.parallel_offset(distance=MaxPrPoint*(1-assump.beta1(data)), side="right")
    return ListOfPoints(PL.coords)

def PressureRegion(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> ListOfPoints:
    rSection: Polygon = Polygon(rotateSection(data, angle))
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-(0.85*c)),
                      (minx-10, maxy-(0.85*c)), (minx-10, maxy)])
    PL: LineString = LineString(PressureAxis(data, c, angle, assump))
    unioned = rSection.boundary.union(PL)
    PR = [poly for poly in polygonize(unioned) if poly.representative_point().within(topArea)]
    return ListOfPoints(PR[0].exterior.coords)

def es(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> List[float]:
    rCoords: ListOfPoints = rotateRebarCoords(data, angle)
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    MaxPrPoint: float = MaxPressurePoint(data, c, angle)
    NR: Polygon = Polygon(NeutralRegion(data, c, angle))
    esSign = [1 if NR.contains(Point(point)) else -1 for point in rCoords]
    return [((esSign[i]*NL.distance(rCoords[i]))/MaxPrPoint)*assump.ecu for i in range(len(rCoords))]

def ec(data: DesignData, c: float, angle: float, point: Tuple[float, float], assump:Assumptions=defaultAssumption) -> float:
    NL: LineString = LineString(NeutralAxis(data, c, angle))
    MaxPrPoint: float = MaxPressurePoint(data, c, angle)
    NR: Polygon = Polygon(NeutralRegion(data, c, angle))
    ecSign = 1 if NR.contains(Point(point)) else -1
    return ((ecSign * NL.distance(Point(point)))/MaxPrPoint)*assump.ecu

def fs(data: DesignData, c: float, angle: float, assump:Assumptions) -> List[float]:
    _es = es(data, c, angle, assump)
    _fs = [min(abs(e)*data.Es, data.fy) for e in _es]
    return [copysign(_fs[i], _es[i]) for i in range(len(_fs))]

def Fs(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> List[float]:
    _fs = fs(data, c, angle, assump)
    _Fs =[data.As[i]*_fs[i] if data.As[i]*_fs[i] <= 0 else data.As[i]*(_fs[i]-0.85*data.fc) for i in range(len(_fs))]
    return _Fs

def Cc(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> float:
    _Cc = 0.85 * data.fc * Polygon(PressureRegion(data, c, angle, assump)).area
    return _Cc

def Fsz(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption) -> Tuple[float, float]:
    _Fs = Fs(data, c, angle, assump)
    rebarXcoords = [p[0] for p in data.Coords]
    rebarYcoords = [p[1] for p in data.Coords]
    return sum([_Fs[i]*rebarXcoords[i] for i in range(len(_Fs))]), sum([_Fs[i]*rebarYcoords[i] for i in range(len(_Fs))])

def M(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption, IsPhi: bool = True) -> Tuple[float, float, float]:
    signX = 1 if (0 <= angle <= 90) or (270 <= angle <= 360) else -1
    signY = 1 if (0 <= angle <= 180) else -1
    _Fszx, _Fszy = Fsz(data, c, angle, assump)
    _Cc = Cc(data, c, angle, assump)
    rPr = rotate(Polygon(PressureRegion(data, c, angle, assump)), -angle, Point([0, 0]))
    _zcy = abs(rPr.centroid.y)
    _zcx = abs(rPr.centroid.x)
    _es = es(data, c, angle, assump)
    _phi = assump.phif(data, min(_es))
    _Mx = _phi*(_Cc*_zcy + abs(_Fszy)) if IsPhi else _Cc*_zcy + abs(_Fszy)
    _My = _phi*(_Cc*_zcx + abs(_Fszx)) if IsPhi else _Cc*_zcx + abs(_Fszx)
    _M = pow((_Mx)**2+(_My)**2, 0.5)
    return _M, signX*_Mx, signY*_My

def P(data: DesignData, c: float, angle: float, assump:Assumptions=defaultAssumption, IsPhi: bool = True) -> float:
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

def C(data: DesignData, P: float, angle: float, assump:Assumptions=defaultAssumption) -> float:
    _, miny, _, maxy = Polygon(rotateSection(data, angle)).bounds
    return least_squares(OptimF, ((maxy-miny)*2.5), bounds=((0.00001), (5*(maxy-miny))), args=(P, angle, data, assump)).x[0]

def OptimMaxM(x, *args):
    _P = x[0]
    _angle = args[0]
    _data = args[1]
    _assump = args[2]
    rSection = Polygon(rotateSection(_data, _angle))
    _, miny, _, maxy = rSection.bounds
    _c = least_squares(OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), args=(_P, _angle, _data, _assump)).x[0]
    return -M(_data, _c, _angle, _assump)[0]

def CalcPMmax(data: DesignData, angle: float, assump:Assumptions=defaultAssumption) -> Tuple[float, float]:
    output = minimize(OptimMaxM, ((P0(data)/2)), method="L-BFGS-B", args=(angle, data, assump))
    return output.x[0], -output.fun

def OptimM(x, *args):
    _P = x[0]
    _angle = args[0]
    _e0 = args[1]
    _data = args[2]
    _assump = args[3]
    rSection = Polygon(rotateSection(_data, _angle))
    _, miny, _, maxy = rSection.bounds
    _c = least_squares(OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), args=(_P, _angle, _data, _assump)).x[0]
    return abs(M(_data, _c, _angle, _assump)[0]/_P - _e0)

def CalcPMRatio(data: DesignData, P: float, Mx: float, My: float, assump:Assumptions=defaultAssumption) -> float:
    angle = AngleFromForces(data, P, Mx, My, assump)
    M = pow(Mx**2 + My**2, 0.5)
    _P = least_squares(OptimM, ((P0(data)+PtMax(data))/2), bounds=((PtMax), (P0)), args=(angle, M/P, data, assump)).x[0] if P != 0 else 1
    _M = CalcMn(data, _P, angle, assump)[0]
    return P/_P if (P/_P) != 0 else M/_M

def CalcMn(data: DesignData, P: float, angle: float, assump:Assumptions=defaultAssumption) -> Tuple[float, float, float, float]:
    c = C(data, P, angle, assump)
    _M, Mx, My = M(data, float(c), angle, assump)
    alpha = Alpha(data, Mx, My, assump)
    return _M, Mx, My, alpha

def OptimPercent(x, *args):
    percent = x[0]
    P = args[0]
    Mx = args[1]
    My = args[2]
    data = args[3]
    assump = args[4]
    Data = setAsPercent(data, percent)
    return abs(CalcPMRatio(Data, P, Mx, My, assump) - 1)

def CalcPercent(data:DesignData, P: float, Mx: float, My: float, assump: Assumptions=defaultAssumption) -> float:
    return least_squares(OptimPercent, (1,), bounds=((1,), (8,)), args=(P, Mx, My, data, assump)).x[0]
