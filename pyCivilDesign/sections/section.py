from enum import StrEnum, auto
from typing import List, Tuple
from math import sin, cos, pi

from shapely import Point, LineString, Polygon
from shapely.affinity import translate, rotate


ListOfPoints = List[Tuple[float, float]]

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
    return list(sct.exterior.coords)


def TrapzoidSct(b1: float, b2: float, h: float) -> ListOfPoints:
    sct = Polygon([(b1/2, h/2), (-b1/2, h/2), (-b2/2, -h/2), (b2/2, -h/2)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def TShapeSct(b: float, h: float, th1: float, tb1: float, th2: float|None=None, tb2: float|None=None) -> ListOfPoints:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(b/2, h), (-b/2, h), (-b/2, h-th1), (-tb2/2, h-th2),
                  (-tb1/2, 0), (tb1/2, 0), (tb2/2, h-th2), (b/2, h-th1), (b/2, h)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def LShapeSct(b: float, h: float, th1: float, tb1: float, th2: float|None=None, tb2: float|None=None) -> ListOfPoints:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(0, 0), (0, h), (b, h), (b, h-th1),
                  (tb2, h-th2), (tb1, 0), (0, 0)])
    msct = moveCentroidToOrigin(sct)
    return list(msct.exterior.coords)


def BoxSct(b: float, h: float, th: float) -> ListOfPoints:
    outerRect = Polygon(RectangleSct(b, h))
    innerRect = Polygon(RectangleSct(b-2*th, h-2*th))
    sct = outerRect.difference(innerRect)
    return list(sct.exterior.coords)


def CircleSct(r: float) -> ListOfPoints:
    sct = Point(0, 0).buffer(r)
    return list(sct.exterior.coords)


def PipeSct(r: float, th: float) -> ListOfPoints:
    outerCircle = Polygon(CircleSct(r))
    innerCircle = Polygon(CircleSct(r-2*th))
    sct = outerCircle.difference(innerCircle)
    return list(sct.exterior.coords)


def CreateEllipseSct(a: float, b: float) -> ListOfPoints:
    n = 100
    theta = [2 * pi * i / n for i in range(n)]
    sct = Polygon([(a * cos(t), b * sin(t)) for t in theta])
    return list(sct.exterior.coords)


def DistanceFrom(section: ListOfPoints, dist: float, position: str="top") -> ListOfPoints:
    sct = Polygon(section)
    minx, miny, maxx, maxy = sct.bounds
    if (position == "top" or position=="bottom") and (dist > maxy-miny): 
        raise ValueError("distance is greater than height of shape.")
    if (position == "right" or position=="left") and (dist > maxx-minx): 
        raise ValueError("distance is greater than width of shape.")
    line = LineString([(maxx+10, maxy-dist), (minx-10, maxy-dist)]) if position=="top" else\
            LineString([(maxx+10, miny+dist), (minx-10, miny+dist)]) if position=="bottom" else\
            LineString([(maxx-dist, maxy+10), (maxx-dist, miny-10)]) if position=="right" else\
            LineString([(minx+dist, maxy+10), (minx+dist, miny-10)]) if position=="left" else None
    return list(sct.intersection(line).coords)


def Edge(section: list, position: str = "top"):
    sct = Polygon(section)
    minx, miny, maxx, maxy = sct.bounds
    Points = [point for point in section if point[1]==maxy] if position=="top" else\
                [point for point in section if point[1]==miny] if position=="bottom" else\
                [point for point in section if point[0]==maxx] if position=="right" else\
                [point for point in section if point[0]==minx] if position=="left" else None
    if Points==None: raise ValueError("position is wrong")
    if len(Points) > 1:
        output = LineString(Points)
    else:
        firstPoint = Points[0]
        if section.index(firstPoint)==0:
            secondPoint = section[section.index(firstPoint)+1]
        elif section.index(firstPoint)==len(section)-1:
            secondPoint = section[section.index(firstPoint)-1]
        else:
            spoint1 = section[section.index(firstPoint)-1]
            spoint2 = section[section.index(firstPoint)+1]
            if position == "top":
                secondPoint = spoint1 if firstPoint[1]-spoint1[1] < firstPoint[1]-spoint2[1] else spoint2
            elif position == "bottom":
                secondPoint = spoint1 if spoint1[1]-firstPoint[1] < spoint2[1]-firstPoint[1] else spoint2
            elif position == "right":
                secondPoint = spoint1 if firstPoint[0]-spoint1[0] < firstPoint[0]-spoint2[0] else spoint2
            elif position == "left":
                secondPoint = spoint1 if spoint1[0]-firstPoint[0] < spoint2[0]-firstPoint[0] else spoint2
        output = LineString([firstPoint, secondPoint])
    return list(output.coords)