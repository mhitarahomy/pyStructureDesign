from enum import StrEnum, auto
from typing import List, Tuple
from math import sin, cos, pi

from shapely import Point, LineString, Polygon
from shapely.affinity import translate, rotate


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


def moveCentroidToOrigin(section: Polygon) -> Polygon:
    return translate(section, -section.centroid.x, -section.centroid.y)


def TriangleSct(b: float, h: float):
    sct = Polygon([(0, h), (-b/2, 0), (b/2, 0)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def RectangleSct(b: float, h: float) -> ListOfPoints:
    sct = Polygon([(b/2, h/2), (-b/2, h/2), (-b/2, -h/2), (b/2, -h/2)])
    return ListOfPoints(list(sct.exterior.coords))


def TrapzoidSct(b1: float, b2: float, h: float) -> ListOfPoints:
    sct = Polygon([(b1/2, h/2), (-b1/2, h/2), (-b2/2, -h/2), (b2/2, -h/2)])
    msct = moveCentroidToOrigin(sct)
    return ListOfPoints(msct.exterior.coords)


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