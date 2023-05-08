from dataclasses import dataclass, field
from typing import Any, List, Tuple
from enum import StrEnum, auto
from math import sin, cos, pi

from shapely import Point, Polygon
from shapely.affinity import translate

from ..materials import ConcreteMat, RebarMat, AIII, AII, C25


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
    
    @property
    def sectionCoords(self) -> ListOfPoints:
        return ListOfPoints(list(Polygon(self.section).exterior.coords))
    

@dataclass()
class RectConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default = 600)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: ListOfPoints = field(init=False)
    sectionType: SectionType = field(init=False, default=SectionType.Rectangle)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h":
            self.section =  RectangleSct(self.b, self.h)
    

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


# this code must go to concreteAnalysis
# def showSection(concSct: ConcreteSct) -> None:
#     fig = plt.figure(dpi=90)
#     ax = fig.add_subplot()  # type: ignore
#     plot_polygon(Polygon(concSct.section), add_points=False, linewidth=1)
#     if concSct.rebarCoords != None:
#         plot_points([Point(rcoord.point) for rcoord in concSct.rebarCoords], color='black')
#         for i in range(len(concSct.rebarCoords)):
#             rcoord = concSct.rebarCoords[i]
#             ax.annotate(f"{rcoord.rebar}\n{i}", xy=(rcoord.point[0], rcoord.point[1]), xycoords='data',
#                     xytext=(1.5, -10), textcoords='offset points')
#     plt.show()
