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


def TriangleSct(b: float, h: float) -> Polygon:
    sct = Polygon([(0, h), (-b/2, 0), (b/2, 0)])
    return moveCentroidToOrigin(sct)


def RectangleSct(b: float, h: float) -> Polygon:
    return Polygon([(b/2, h/2), (-b/2, h/2), (-b/2, -h/2), (b/2, -h/2)])


def TrapzoidSct(b1: float, b2: float, h: float) -> Polygon:
    sct = Polygon([(b1/2, h/2), (-b1/2, h/2), (-b2/2, -h/2), (b2/2, -h/2)])
    return moveCentroidToOrigin(sct)


def TShapeSct(b: float, h: float, th1: float, tb1: float,
              th2: float|None=None, tb2: float|None=None) -> Polygon:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(b/2, h), (-b/2, h), (-b/2, h-th1), (-tb2/2, h-th2),
                  (-tb1/2, 0), (tb1/2, 0), (tb2/2, h-th2), (b/2, h-th1), (b/2, h)])
    return moveCentroidToOrigin(sct)


def LShapeSct(b: float, h: float, th1: float, tb1: float,
              th2: float|None=None, tb2: float|None=None) -> Polygon:
    tb2 = tb1 if tb2 == None else tb2
    th2 = th1 if th2 == None else th2
    sct = Polygon([(0, 0), (0, h), (b, h), (b, h-th1),
                  (tb2, h-th2), (tb1, 0), (0, 0)])
    return moveCentroidToOrigin(sct)


def BoxSct(b: float, h: float, th: float) -> Polygon:
    outerRect = Polygon(RectangleSct(b, h))
    innerRect = Polygon(RectangleSct(b-2*th, h-2*th))
    return outerRect.difference(innerRect)


def CircleSct(d: float) -> Polygon:
    return Point(0, 0).buffer(d/2)


def PipeSct(d: float, th: float) -> Polygon:
    outerCircle = Polygon(CircleSct(d))
    innerCircle = Polygon(CircleSct(d-(2*th)))
    return outerCircle.difference(innerCircle)


def CreateEllipseSct(a: float, b: float) -> Polygon:
    n = 100
    theta = [2 * pi * i / n for i in range(n)]
    return Polygon([(a*cos(t), b*sin(t)) for t in theta])


def DistanceFrom(section: Polygon, dist: float, position: str="top") -> ListOfPoints:
    minx, miny, maxx, maxy = section.bounds
    if (position == "top" or position=="bottom") and (dist > maxy-miny): 
        raise ValueError("distance is greater than height of shape.")
    if (position == "right" or position=="left") and (dist > maxx-minx): 
        raise ValueError("distance is greater than width of shape.")
    line = LineString([(maxx+10, maxy-dist), (minx-10, maxy-dist)]) if position=="top" else\
            LineString([(maxx+10, miny+dist), (minx-10, miny+dist)]) if position=="bottom" else\
            LineString([(maxx-dist, maxy+10), (maxx-dist, miny-10)]) if position=="right" else\
            LineString([(minx+dist, maxy+10), (minx+dist, miny-10)]) if position=="left" else None
    return list(section.intersection(line).coords)


def Edge(section: Polygon, position: str = "top") -> LineString:
    coords = list(section.exterior.coords)
    del coords[-1]
    minx, miny, maxx, maxy = section.bounds
    Points = [point for point in coords if point[1]==maxy] if position=="top" else\
                [point for point in coords if point[1]==miny] if position=="bottom" else\
                [point for point in coords if point[0]==maxx] if position=="right" else\
                [point for point in coords if point[0]==minx] if position=="left" else None
    spoint1 = (0, 0)
    spoint2 = (0, 0)
    IsThreePoint: bool = False
    if Points==None: raise ValueError("position is wrong")
    if len(Points) > 1:
        output = LineString(Points)
    else:
        firstPoint = Points[0]
        secondPoint = (0, 0)
        if coords.index(firstPoint)==0:
            secondPoint = coords[coords.index(firstPoint)+1]
        elif coords.index(firstPoint)==len(coords)-1:
            secondPoint = coords[coords.index(firstPoint)-1]
        else:
            spoint1 = coords[coords.index(firstPoint)-1]
            spoint2 = coords[coords.index(firstPoint)+1]
            if position == "top":
                IsThreePoint = firstPoint[1]-spoint1[1] == firstPoint[1]-spoint2[1]
                secondPoint = spoint1 if firstPoint[1]-spoint1[1] < firstPoint[1]-spoint2[1] else spoint2
            elif position == "bottom":
                IsThreePoint = spoint1[1]-firstPoint[1] == spoint2[1]-firstPoint[1]
                secondPoint = spoint1 if spoint1[1]-firstPoint[1] < spoint2[1]-firstPoint[1] else spoint2
            elif position == "right":
                IsThreePoint = firstPoint[0]-spoint1[0] == firstPoint[0]-spoint2[0]
                secondPoint = spoint1 if firstPoint[0]-spoint1[0] < firstPoint[0]-spoint2[0] else spoint2
            elif position == "left":
                IsThreePoint = spoint1[0]-firstPoint[0] == spoint2[0]-firstPoint[0]
                secondPoint = spoint1 if spoint1[0]-firstPoint[0] < spoint2[0]-firstPoint[0] else spoint2
        
        output = LineString([spoint1, firstPoint, spoint2]) if IsThreePoint else LineString([firstPoint, secondPoint])
    return output
