from dataclasses import dataclass, field
from typing import Any, TypedDict, Optional, List, Tuple
from enum import StrEnum, auto
from math import sin, cos, pi
from matplotlib import pyplot as plt

from shapely import Point, Polygon, LineString
from shapely.affinity import rotate, translate
from shapely.plotting import plot_line, plot_points, plot_polygon

from pyCivil.materials import ConcreteMat, RebarMat, AIII, AII, C25
from pyCivil.sections.rebarSections import Rebar, GRebars


class ListOfPoints(List[Tuple[float, float]]):
    ...

class SectionType(StrEnum):
    Triangle = auto()
    Rectangle = auto()
    Trapezoid = auto()
    TShape = auto()
    LShape = auto()
    Box = auto()
    Pipe = auto()
    RegularPolygon = auto()
    Circle = auto()
    Ellipse = auto()
    Custom = auto()


class ConfType(StrEnum):
    Tie = auto()
    Spiral = auto()


@dataclass
class Cover():
    Top: float
    Bottom: float
    Right: float
    Left: float
    

@dataclass
class RebarCoords():
    point: tuple[float, float]
    rebar: Rebar|GRebars

C25def = lambda: C25
AIIdef = lambda: AII
AIIIdef = lambda: AIII

@dataclass()
class ConcreteSct():
    section: ListOfPoints
    sectionType: SectionType
    concMat: ConcreteMat
    lBarMat: RebarMat 
    cBarMat: RebarMat
    confType: ConfType
    rebarCoords: list[RebarCoords] = field(init=False, default_factory=list)
    
    @property
    def sectionCoords(self) -> ListOfPoints:
        return ListOfPoints(list(Polygon(self.section).exterior.coords))
    @property
    def rebarNums(self) -> int: return len(self.rebarCoords)

    @property
    def coords(self) -> list[tuple[float, float]]: return [self.rebarCoords[i].point for i in range(self.rebarNums)]
    
    @property
    def rebarList(self) -> list[Rebar|GRebars]: return [self.rebarCoords[i].rebar for i in range(self.rebarNums)]

    @property
    def As(self) -> list[float]: return [self.rebarList[i].Area for i in range(self.rebarNums)]


@dataclass()
class RectConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default = 600)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: ListOfPoints = field(init=False)
    sectionType: SectionType = field(init=False, default=SectionType.Rectangle)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    confType: ConfType = field(default=ConfType.Tie)
    rebarCoords: list[RebarCoords] = field(init=False, default_factory=list)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h":
            self.section =  RectangleSct(self.b, self.h)
    
    @property
    def TopEdge(self): return ListOfPoints([self.sectionCoords[0], self.sectionCoords[1]])
    
    @property
    def BottomEdge(self): return ListOfPoints([self.sectionCoords[2], self.sectionCoords[3]])
    
    @property
    def LeftEdge(self): return ListOfPoints([self.sectionCoords[1], self.sectionCoords[2]])
    
    @property
    def RightEdge(self): return ListOfPoints([self.sectionCoords[3], self.sectionCoords[0]])


# region Create Sections
def moveCentroidToOrigin(section: Polygon) -> Polygon:
    return translate(section, -section.centroid.x, -section.centroid.y)


def TriangleSct(b: float, h: float):
    sct = Polygon([(0, h), (-b/2, 0), (b/2, 0)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def RectangleSct(b: float, h: float) -> ListOfPoints:
    sct = Polygon([(b/2, h/2), (-b/2, h/2), (-b/2, -h/2), (b/2, -h/2)])
    return ListOfPoints(list(sct.exterior.coords))


def TrapzoidSct(b1: float, b2: float, h: float) -> list[tuple[float, float]]:
    sct = Polygon([(b1/2, h/2), (-b1/2, h/2), (-b2/2, -h/2), (b2/2, -h/2)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def TShapeSct(b: float, h: float, th1: float, tb1: float, th2: float|None=None, tb2: float|None=None) -> list[tuple[float, float]]:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(b/2, h), (-b/2, h), (-b/2, h-th1), (-tb2/2, h-th2),
                  (-tb1/2, 0), (tb1/2, 0), (tb2/2, h-th2), (b/2, h-th1), (b/2, h)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def LShapeSct(b: float, h: float, th1: float, tb1: float, th2: float|None=None, tb2: float|None=None) -> list[tuple[float, float]]:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(0, 0), (0, h), (b, h), (b, h-th1),
                  (tb2, h-th2), (tb1, 0), (0, 0)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def BoxSct(b: float, h: float, th: float) -> list[tuple[float, float]]:
    outerRect = Polygon(RectangleSct(b, h))
    innerRect = Polygon(RectangleSct(b-2*th, h-2*th))
    sct = outerRect.difference(innerRect)
    return list(sct.exterior.coords)


def CircleSct(r: float) -> list[tuple[float, float]]:
    sct = Point(0, 0).buffer(r)
    return list(sct.exterior.coords)


def PipeSct(r: float, th: float) -> list[tuple[float, float]]:
    outerCircle = Polygon(CircleSct(r))
    innerCircle = Polygon(CircleSct(r-2*th))
    sct = outerCircle.difference(innerCircle)
    return list(sct.exterior.coords)


def CreateEllipseSct(a: float, b: float) -> list[tuple[float, float]]:
    n = 100
    theta = [2 * pi * i / n for i in range(n)]
    sct = Polygon([(a * cos(t), b * sin(t)) for t in theta])
    return list(sct.exterior.coords)
# endregion

# region Create Rebar Shapes
def QuadrilateralRebarShape(xNum: int, yNum: int, rebarSct: Rebar|GRebars) -> list[list[Rebar|GRebars]]:
    return [[rebarSct for i in range(xNum)] if (j == 0 or j == yNum-1) else [rebarSct, rebarSct] for j in range(yNum)]


def LinearRebarShape(num: int, rebarSct:Rebar|GRebars) -> list[Rebar|GRebars]:
    return [rebarSct for i in range(num)]


def CircularRebarShape(num: int, rebarSct:Rebar|GRebars) -> list[Rebar|GRebars]:
    return [rebarSct for i in range(num)]
# endregion

# region Create & Edit Rebars
def setCover(cover:Cover|float|int):
    return cover if type(cover)==Cover else Cover(cover, cover, cover, cover) if (
        type(cover)==float or type(cover)==int) else Cover(0, 0, 0, 0)

def PointRebar(rebar: Rebar|GRebars, point: tuple[float, float]) -> list[RebarCoords]:
    return [RebarCoords(point=point, rebar=rebar)]

def LinearRebars(barShape: list[Rebar|GRebars], startCover:float, endCover:float , line: ListOfPoints, 
                 offsetDist: float=0, offsetSide: str = "left", 
                 barDistances: list[float]|None = None) -> list[RebarCoords]:
    barsNum = len(barShape)
    if (barDistances != None) and (len(barDistances) != len(barShape)):
        raise ValueError("barShape and barDistances must have same length.")
    _line = LineString(line)
    pStart = _line.parallel_offset(offsetDist, offsetSide).interpolate(startCover)
    pEnd = _line.parallel_offset(offsetDist, offsetSide).interpolate(_line.length-endCover)
    d = ([i / (barsNum-1) for i in range(barsNum)] if barsNum!=1 else [0.5]) if barDistances == None else barDistances
    return [RebarCoords(point=list(LineString([pStart, pEnd]).interpolate(d[i], normalized=True).coords)[0], rebar=barShape[i]) for i in range(barsNum)]

def QuadrilateralRebars(barShape: list[list[Rebar|GRebars]], cover: Cover|float|int, 
                              points: ListOfPoints, barDistances: list[list[float]]|None=None,
                              layerDistances: list[float]|None=None) -> list[RebarCoords]:
    _cover = setCover(cover)
    barsNumList: list[int] = [len(s) for s in barShape]
    if 5<len(points)<4: raise ValueError("Number of points must be 4")
    p00, p11, p22, p33 = Point(points[0]), Point(points[1]), Point(points[2]), Point(points[3])
    line0 = LineString([p00, p33])
    line1 = LineString([p11, p22])
    p0 = line0.interpolate(_cover.Top)
    p3 = line0.interpolate(line0.length-_cover.Bottom)
    p1 = line1.interpolate(_cover.Top)
    p2 = line1.interpolate(line1.length-_cover.Bottom)
    d = ([i / (len(barsNumList)-1) for i in range(len(barsNumList))]
         if len(barsNumList) != 1 else [0.5]) if layerDistances == None else layerDistances
    rebarCoords = []
    for i in range(len(barsNumList)):
        pStart: Point = LineString([p0, p3]).interpolate(d[i], normalized=True)
        pEnd: Point = LineString([p1, p2]).interpolate(d[i], normalized=True)
        if  barDistances == None:
            _rebarCoords = LinearRebars(barShape[i], _cover.Right, _cover.Left, [pStart.coords[0], pEnd.coords[0]])  # type: ignore
        else: 
            _rebarCoords = LinearRebars(barShape[i], _cover.Right, _cover.Left, [pStart.coords[0], pEnd.coords[0]], barDistances=barDistances[i]) # type: ignore
        rebarCoords.extend(_rebarCoords)
    return rebarCoords

def RectSectRebars(section: ListOfPoints, xNum: int, yNum: int, rebar: Rebar|GRebars, cover: Cover|float|int):
    rshape = QuadrilateralRebarShape(xNum, yNum, rebar)
    return QuadrilateralRebars(rshape, cover, section)

def CircularBars(barShape: list[Rebar|GRebars], cover: float, D: float) -> list[RebarCoords]:
    barsNum = len(barShape)
    r = (D/2) - cover
    alphas = [i * (360/barsNum) for i in range(barsNum)]
    return [RebarCoords(rebar=barShape[i], point=rotate(Point(0, r), alphas[i], Point(0, 0))) for i in range(barsNum)]

def DeleteRebars(rebarCoords: list[RebarCoords], *args: int) -> list[RebarCoords]:
    output = rebarCoords
    for index in sorted(args, reverse=True):
            del output[index]
    return output

def ChangeRebars(rebarCoords: list[RebarCoords], changeList: list[tuple[int, Rebar|GRebars]]):
    output = rebarCoords
    for item in changeList:
            output[item[0]].rebar = item[1]
    return output

def showSection(concSct: ConcreteSct) -> None:
    fig = plt.figure(dpi=90)
    ax = fig.add_subplot()  # type: ignore
    plot_polygon(Polygon(concSct.section), add_points=False, linewidth=1)
    if concSct.rebarCoords != None:
        plot_points([Point(rcoord.point) for rcoord in concSct.rebarCoords], color='black')
        for i in range(len(concSct.rebarCoords)):
            rcoord = concSct.rebarCoords[i]
            ax.annotate(f"{rcoord.rebar}\n{i}", xy=(rcoord.point[0], rcoord.point[1]), xycoords='data',
                    xytext=(1.5, -10), textcoords='offset points')
    plt.show()
# endregion
