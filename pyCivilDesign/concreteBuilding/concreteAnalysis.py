from dataclasses import dataclass, field
from typing import Tuple, List

import numpy as np
import numpy.typing as npt

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate
from shapely.ops import polygonize

from scipy.optimize import least_squares, minimize

from pyCivilDesign.sections.section import ListOfPoints


@dataclass
class PMMDesignData():
    section: Polygon
    fy: float
    fc: float   
    Es: float = field(default=2e5) 
    ecu: float = field(default=0.003)


beta1 = lambda fc: 0.85 if fc<=28 else max(0.85-(0.05*(fc-28)/7), 0.65)

# def setAs(self, _As: npt.NDArray[np.float32]) -> None:
#     self.As = _As

# def setAsPercent(self, percent: float) -> None:
#     self.setAs(np.ones_like(self.rebarCoords) *
#                (self.section.area * (percent/100) / len(self.rebarCoords)))

def rotateSection(section: ListOfPoints, angle: float) -> ListOfPoints:
    rsection: Polygon = rotate(Polygon(section), angle, origin=Point([0, 0]))
    return ListOfPoints(rsection.exterior.coords)

def rotateRebarCoords(coords: ListOfPoints, angle: float) -> ListOfPoints:
    points: List[Point] = [rotate(Point(coord), angle, origin=Point([0, 0])) for coord in coords]
    return ListOfPoints([(point.x, point.y) for point in points])

def NeutralAxis(section: ListOfPoints, c: float, angle: float) -> ListOfPoints:
    rSection: Polygon = Polygon(rotateSection(section, angle))
    minx, _, maxx, maxy = rSection.bounds
    line: LineString = LineString([(maxx+10, maxy-c), (minx-10, maxy-c)])
    return ListOfPoints(line.coords)

def NeutralRegion(section: ListOfPoints, c: float, angle: float) -> ListOfPoints:
    rSection = Polygon(rotateSection(section, angle))
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-c),
                      (minx-10, maxy-c), (minx-10, maxy)])
    NL: LineString = LineString(NeutralAxis(section, c, angle))
    unioned = rSection.boundary.union(NL)
    NR: List[Polygon] = [poly for poly in polygonize(unioned) if poly.representative_point().within(topArea)]
    return ListOfPoints(NR[0].exterior.coords)

def MaxPressurePoint(section: ListOfPoints, c: float, angle: float) -> float:
    NL: LineString = LineString(NeutralAxis(section, c, angle))
    NR: Polygon = Polygon(NeutralRegion(section, c, angle))
    return max([NL.distance(Point(p)) for p in list(NR.exterior.coords)])

def PressureAxis(section: ListOfPoints, fc:float, c: float, angle: float) -> ListOfPoints:
    MaxPrPoint: float = MaxPressurePoint(section, c, angle)
    NL: LineString = LineString(NeutralAxis(section, c, angle))
    PL: LineString = NL.parallel_offset(distance=MaxPrPoint*(1-beta1(fc)), side="right")
    return ListOfPoints(PL.coords)

def PressureRegion(section: ListOfPoints, fc:float, c: float, angle: float) -> ListOfPoints:
    rSection: Polygon = Polygon(rotateSection(section, angle))
    minx, _, maxx, maxy = rSection.bounds
    topArea = Polygon([(maxx+10, maxy), (maxx+10, maxy-(0.85*c)),
                      (minx-10, maxy-(0.85*c)), (minx-10, maxy)])
    PL: LineString = LineString(PressureAxis(section, fc, c, angle))
    unioned = rSection.boundary.union(PL)
    PR = [poly for poly in polygonize(unioned) if poly.representative_point().within(topArea)]
    return ListOfPoints(PR[0].exterior.coords)

def es(section: ListOfPoints, coords:ListOfPoints, ecu: float, c: float, angle: float) -> List[float]:
    rCoords: ListOfPoints = rotateRebarCoords(coords, angle)
    NL: LineString = LineString(NeutralAxis(section, c, angle))
    MaxPrPoint: float = MaxPressurePoint(section, c, angle)
    NR: Polygon = Polygon(NeutralRegion(section, c, angle))
    esSign = [1 if NR.contains(Point(point)) else -1 for point in rCoords]
    return [((esSign[i]*NL.distance(rCoords[i]))/MaxPrPoint)*ecu for i in range(len(rCoords))]

def ec(section: ListOfPoints, ecu: float, c: float, angle: float, point: Tuple[float, float]) -> float:
    NL: LineString = LineString(NeutralAxis(section, c, angle))
    MaxPrPoint: float = MaxPressurePoint(section, c, angle)
    NR: Polygon = Polygon(NeutralRegion(section, c, angle))
    ecSign = 1 if NR.contains(Point(point)) else -1
    return ((ecSign * NL.distance(Point(point)))/MaxPrPoint)*ecu

# def fs(self, c: float, angle: float) -> npt.NDArray[np.float32]:
#     _es = self.es(c, angle)
#     _fs = np.array([min(abs(e)*self.Es, self.fy) for e in _es], dtype=np.float32)
#     # return np.where(abs(_es)*self.Es < self.fy, abs(_es)*self.Es, self.fy)*np.sign(_es)
#     return _fs * np.sign(_es)
# def Fs(self, c: float, angle: float) -> npt.NDArray[np.float32]:
#     _fs = self.fs(c, angle)
#     _Fs = np.array([self.As*f if self.As*f <= 0 else self.As*(f-0.85*self.fc) for f in _fs], dtype=np.float32)
#     # return np.where(self.As*_fs <= 0, self.As*_fs, self.As*(_fs-0.85*self.fc))
#     return _Fs
# def Cc(self, c: float, angle: float) -> np.float32:
#     _Cc = 0.85 * self.fc * self.PressureRegion(c, angle).area
#     return _Cc
# def Fsz(self, c: float, angle: float) -> Tuple[np.float32, np.float32]:
#     return np.sum(self.Fs(c, angle) * self.rebarXcoords), np.sum(self.Fs(c, angle) * self.rebarYcoords)
# def phi(self, c: float, angle: float) -> np.float32:
#     ety = self.fy / self.Es
#     et = abs(np.min(self.es(c, angle)))
#     return np.float32(0.65) if et <= ety else 0.65+0.25*((et-ety)/0.003) if ety < et < ety+0.003 else np.float32(0.9)
# def M(self, c: float, angle: float, IsPhi: bool = True) -> Tuple[np.float32, np.float32, np.float32]:
#     signX = 1 if (0 <= angle <= 90) or (270 <= angle <= 360) else -1
#     signY = 1 if (0 <= angle <= 180) else -1
#     _Fszx, _Fszy = self.Fsz(c, angle)
#     _Cc = self.Cc(c, angle)
#     _zcy = abs(rotate(self.PressureRegion(c, angle), -
#                angle, Point([0, 0])).centroid.y)
#     _zcx = abs(rotate(self.PressureRegion(c, angle), -
#                angle, Point([0, 0])).centroid.x)
#     _phi = self.phi(c, angle)
#     _Mx = _phi*(_Cc*_zcy + abs(_Fszy)) if IsPhi else _Cc*_zcy + abs(_Fszy)
#     _My = _phi*(_Cc*_zcx + abs(_Fszx)) if IsPhi else _Cc*_zcx + abs(_Fszx)
#     _M = pow((_Mx)**2+(_My)**2, 0.5)
#     return _M, signX*_Mx, signY*_My
# def P(self, c: float, angle: float, IsPhi: bool = True) -> np.float32:
#     _Fs = self.Fs(c, angle)
#     _Cc = self.Cc(c, angle)
#     _phi = self.phi(c, angle)
#     return _phi*(_Cc+np.sum(_Fs)) if IsPhi else _Cc+np.sum(_Fs)
# def C(self, P: float, angle: float) -> np.float32:
#     def OptimF(x, *args):
#         c = x[0]
#         P = args[0]
#         angle = args[1]
#         return abs(P - self.P(c, angle))
#     _, miny, _, maxy = self.rotateSection(angle).bounds
#     return least_squares(OptimF, ((maxy-miny)*2.5), bounds=((0.00001), (5*(maxy-miny))), args=(P, angle)).x[0]
# def CalcPMmax(self, angle: float) -> Tuple[np.float32, np.float32]:
#     def OptimF(x, *args):
#         c = x[0]
#         P = args[0]
#         angle = args[1]
#         return abs(P - self.P(c, angle))
#     def OptimMaxM(x, *args):
#         P = x[0]
#         angle = args[0]
#         rSection = self.rotateSection(angle)
#         _, miny, _, maxy = rSection.bounds
#         c = least_squares(
#             OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), args=(P, angle)).x[0]
#         return -self.M(c, angle)[0]
#     output = minimize(OptimMaxM, ((self.P0/2)),
#                       method="L-BFGS-B", args=(angle,))
#     return output.x[0], -output.fun
# def CalcPMRatio(self, P: float, Mx: float, My: float) -> np.float32:
#     def OptimF(x, *args):
#         c = x[0]
#         P = args[0]
#         angle = args[1]
#         return abs(P - self.P(c, angle))
#     def OptimM(x, *args):
#         P = x[0]
#         angle = args[0]
#         e0 = args[1]
#         rSection = self.rotateSection(angle)
#         _, miny, _, maxy = rSection.bounds
#         c = least_squares(
#             OptimF, ((maxy-miny)/2), bounds=((0.00001), (5*(maxy-miny))), args=(P, angle)).x[0]
#         return abs(self.M(c, angle)[0]/P - e0)
#     angle = self.AngleFromForces(P, Mx, My)
#     M = pow(Mx**2 + My**2, 0.5)
#     _P = least_squares(OptimM, ((self.P0+self.PtMax)/2), bounds=(
#         (self.PtMax), (self.P0)), args=(angle, M/P)).x[0] if P != 0 else 1
#     _M = self.CalcMn(_P, angle)[0]
#     return np.float32(P/_P) if (P/_P) != 0 else np.float32(M/_M)
# def CalcMn(self, P: float, angle: float) -> Tuple[float, float, float, float]:
#     c = self.C(P, angle)
#     M, Mx, My = self.M(float(c), angle)
#     alpha = self.Alpha(float(Mx), float(My))
#     return M, Mx, My, alpha
# def AngleFromForces(self, P, Mx, My) -> float:
#     def OptimAngle(x, *args):
#         angle = x[0]
#         P = args[0]
#         alpha = args[1]
#         return abs(self.CalcMn(P, angle)[3]-alpha)
#     alpha = self.Alpha(Mx, My)
#     lAngle = 0 if 0 <= alpha <= 90 else 90 if 90 < alpha <= 180 else 180 if 180 < alpha <= 270 else 270
#     uAngle = 90 if 0 <= alpha <= 90 else 180 if 90 < alpha <= 180 else 270 if 180 < alpha <= 270 else 360
#     return least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)), args=(P, alpha)).x[0]
# def AngleFromAlpha(self, P, alpha) -> float:
#     def OptimAngle(x, *args):
#         angle = x[0]
#         P = args[0]
#         alpha = args[1]
#         return abs(self.CalcMn(P, angle)[3]-alpha)
#     lAngle = 0 if 0 <= alpha <= 90 else 90 if 90 < alpha <= 180 else 180 if 180 < alpha <= 270 else 270
#     uAngle = 90 if 0 <= alpha <= 90 else 180 if 90 < alpha <= 180 else 270 if 180 < alpha <= 270 else 360
#     return least_squares(OptimAngle, (alpha), bounds=((lAngle), (uAngle)), args=(P, alpha)).x[0]
# def Alpha(self, Mx: float, My: float) -> float:
#     alpha = np.degrees(np.arctan(abs(My/Mx))
#                        ) if Mx != 0 else np.degrees(np.arctan(np.inf))
#     alpha = alpha if (Mx >= 0 and My >= 0) else 180-alpha if (Mx < 0 and My > 0) else 180 + \
#         alpha if (Mx < 0 and My < 0) else 360 - \
#         alpha if (Mx > 0 and My < 0) else alpha
#     return alpha
# def CalcPercent(self, P: float, Mx: float, My: float) -> float:
#     def OptimPercent(x, *args):
#         percent = x[0]
#         sctAnalysis = args[0]
#         P = args[1]
#         Mx = args[2]
#         My = args[3]
#         sctAnalysis.setAsPercent(percent)
#         return abs(sctAnalysis.CalcPMRatio(P, Mx, My) - 1)
#     sctAnalysis = PMMAnalysis(
#         self.section, float(self.fy), float(self.fc), list(self.rebarCoords), self.As)
#     return least_squares(OptimPercent, (1,), bounds=((1,), (8,)), args=(sctAnalysis, P, Mx, My)).x[0]
# @property
# def P0(self):
#     return 0.65 * (0.85 * self.fc * (self.section.area-np.sum(self.As))+np.sum(self.As)*self.fy)
# @property
# def PnMax(self):
#     return 0.8 * self.P0
# @property
# def PtMax(self):
#     return -0.9 * (np.sum(self.As) * self.fy)
# @property
# def AsPercent(self) -> np.float32:
#     return (np.sum(self.As)/self.section.area)*100
# @property
# def onePercentAnalysis(self):
#     sctAnalysis = PMMAnalysis(
#         self.section, float(self.fy), float(self.fc), list(self.rebarCoords), self.As)
#     sctAnalysis.setAs(np.ones_like(self.rebarCoords) *
#                       (self.section.area * (1/100) / len(self.rebarCoords)))
#     return sctAnalysis
# @property
# def eightPercentAnalysis(self):
#     sctAnalysis = PMMAnalysis(
#         self.section, float(self.fy), float(self.fc), list(self.rebarCoords), self.As)
#     sctAnalysis.setAs(np.ones_like(self.rebarCoords) *
#                       (self.section.area * (8/100) / len(self.rebarCoords)))
#     return sctAnalysis