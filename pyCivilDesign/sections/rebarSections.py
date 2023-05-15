from enum import StrEnum, auto
from dataclasses import dataclass
from typing import Callable, List, Tuple
from math import pi, ceil

from shapely import Polygon, LineString, Point
from shapely.affinity import rotate

from pyCivilDesign.sections.section import ListOfPoints

# region Rebar Sections
@dataclass
class Rebar():
    d: float
    Area = property(lambda self: CalcArea(self.d, 1))

    def __repr__(self) -> str:
        return f"\u03C6{self.d}"

    def __mul__(self, num: int):
        return GRebars({f"{self.d}": num})

    def __rmul__(self, num: int):
        return GRebars({f"{self.d}": num})

    def __add__(self, rebar):
        return GRebars({f"{self.d}": 2}) if rebar.d == self.d else GRebars({f"{self.d}": 1, f"{rebar.d}": 1})


@dataclass
class GRebars():
    listOfRebars: dict[str, float]
    Area = property(lambda self: CalcRebarsArea(**self.listOfRebars))

    def __repr__(self) -> str:
        rebarStr = ""
        for k, v in list(self.listOfRebars.items()):
            rebarStr += f"{v}\u03C6{k}+"
        return rebarStr[:-1]

    def __add__(self, rebar):
        output = GRebars(self.listOfRebars)
        if type(rebar) == Rebar:
            if str(rebar.d) in self.listOfRebars:
                output.listOfRebars[f"{rebar.d}"] += 1
            else:
                output.listOfRebars.update({f"{rebar.d}": 1})
        elif type(rebar) == GRebars:
            for k, v in list(rebar.listOfRebars.items()):
                if k in output.listOfRebars:
                    output.listOfRebars[k] += v
                else:
                    output.listOfRebars.update({k: v})
        output.listOfRebars = {
            k: output.listOfRebars[k] for k in sorted(output.listOfRebars)}
        return output

d8 = Rebar(8)
d10 = Rebar(10)
d12 = Rebar(12)
d14 = Rebar(14)
d16 = Rebar(16)
d18 = Rebar(18)
d20 = Rebar(20)
d22 = Rebar(22)
d25 = Rebar(25)
d28 = Rebar(28)
d32 = Rebar(32)
d36 = Rebar(36)
# endregion


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
    point: Point
    rebar: Rebar|GRebars


CalcArea: Callable[[float, int], float] = lambda d, num=1: num * (pi*d**2)/4

def CalcRebarsArea(**listOfRebars) -> float:
    area: float = 0
    for k, v in list(listOfRebars.items()):
        area += CalcArea(float(k), v)
    return round(area, 2)

def CalcNumOfBars(area: float, bar: float) -> int:
    return ceil(area/bar)

# region Create Rebar Shapes
def QuadrilateralRebarShape(xNum: int, yNum: int, rebarSct: Rebar|GRebars) -> List[List[Rebar|GRebars]]:
    return [[rebarSct for i in range(xNum)] if (j == 0 or j == yNum-1) else [rebarSct, rebarSct] for j in range(yNum)]


def LinearRebarShape(num: int, rebarSct:Rebar|GRebars) -> List[Rebar|GRebars]:
    return [rebarSct for i in range(num)]


def LinearRebarShapeDist(length: float, distance: int, rebarSct:Rebar|GRebars) -> List[Rebar|GRebars]:
    num = ceil(length / distance) + 1
    return [rebarSct for i in range(num)]


def CircularRebarShape(num: int, rebarSct:Rebar|GRebars) -> List[Rebar|GRebars]:
    return [rebarSct for i in range(num)]
# endregion

# region Create & Edit Rebars
def setCover(cover:Tuple[float, float, float, float]|float|int) -> Tuple[float, float, float, float]:
    return cover if type(cover)==Tuple[float, float, float, float] else (cover, cover, cover, cover) if (
        type(cover)==float or type(cover)==int) else (0, 0, 0, 0)

def CreateLineWithCover(line: LineString, startCover: float, endCover: float) -> LineString:
    _line = list(line.coords)
    startPoint = line.interpolate(startCover)
    endPoint = line.interpolate(line.length - endCover)
    _points =[_line[i] for i in range(1,len(_line)-1,1) if LineString(_line[:i+1]).length < (line.length - endCover)]
    _points.insert(0, _line[0])
    _points.extend([(endPoint.x, endPoint.y)])
    points = [_points[i+1] for i in range(len(_points)-1) if LineString(_points[:i+2]).length > startCover]
    points.insert(0, (startPoint.x, startPoint.y))
    return LineString(points)

def PointRebar(rebar: Rebar|GRebars, point: Point) -> RebarCoords:
    return RebarCoords(point=point, rebar=rebar)

def LinearRebars(barShape: List[Rebar|GRebars], cover: Tuple[float, float], 
                 line: LineString, *, offsetDist: float=0, 
                 offsetSide: str = "left", IsFirst: bool=True, IsEnd: bool=True,
                 barDistances: List[float]|None = None) -> List[RebarCoords]:
    barsNum = len(barShape)
    if (barDistances != None) and (len(barDistances) != len(barShape)):
        raise ValueError("barShape and barDistances must have same length.")
    lineStr = line.parallel_offset(offsetDist, offsetSide) \
        if offsetDist != 0 else line
    _line = CreateLineWithCover(lineStr, cover[0], cover[1])
    d = ([i / (barsNum-1) for i in range(barsNum)] if barsNum!=1 else [0.5]) \
        if barDistances == None else barDistances
    pointList = [PointRebar(barShape[i], 
                            _line.interpolate(d[i], normalized=True))\
                                for i in range(barsNum)]
    if not IsFirst: del pointList[0]
    if not IsEnd: del pointList[-1]
    return pointList


def QuadrilateralRebars(barShape: List[List[Rebar|GRebars]], cover: Tuple[float, float, float, float]|float|int, 
                              line0: LineString, line1: LineString, *, barDistances: List[List[float]]|None=None,
                              layerDistances: List[float]|None=None) -> List[RebarCoords]:
    _cover = setCover(cover)
    barsNumList: list[int] = [len(s) for s in barShape]
    coords0 = list(line0.coords)
    coords1 = list(line1.coords)
    line00 = LineString(coords0.sort(key=lambda x: x[1], reverse=True))
    line11 = LineString(coords1.sort(key=lambda x: x[1], reverse=True))
    sline = CreateLineWithCover(line00, _cover[0], _cover[1])
    eline = CreateLineWithCover(line11, _cover[0], _cover[1])
    d = ([i / (len(barsNumList)-1) for i in range(len(barsNumList))]
         if len(barsNumList) != 1 else [0.5]) if layerDistances == None else layerDistances
    rebarCoords = []
    for i in range(len(barsNumList)):
        pStart: Point = sline.interpolate(d[i], normalized=True)
        pEnd: Point = eline.interpolate(d[i], normalized=True)
        if  barDistances == None:
            _rebarCoords = LinearRebars(barShape[i], (_cover[2], _cover[3]), LineString([pStart.coords[0], pEnd.coords[0]])) 
        else: 
            _rebarCoords = LinearRebars(barShape[i], (_cover[2], _cover[3]), LineString([pStart.coords[0], pEnd.coords[0]]), barDistances=barDistances[i])
        rebarCoords.extend(_rebarCoords)
    return rebarCoords


# def RectSectRebars(section: ListOfPoints, xNum: int, yNum: int, rebar: Rebar|GRebars, cover: Cover|float|int):
#     rshape = QuadrilateralRebarShape(xNum, yNum, rebar)
#     return QuadrilateralRebars(rshape, cover, section)


def CircularBars(barShape: List[Rebar|GRebars], cover: float, D: float) -> List[RebarCoords]:
    barsNum = len(barShape)
    r = (D/2) - cover
    alphas = [i * (360/barsNum) for i in range(barsNum)]
    return [RebarCoords(rebar=barShape[i], point=rotate(Point(0, r), alphas[i], Point(0, 0))) for i in range(barsNum)]


def DeleteRebars(rebarCoords: List[RebarCoords], *args: int) -> List[RebarCoords]:
    output = rebarCoords
    for index in sorted(args, reverse=True):
            del output[index]
    return output


def ChangeRebars(rebarCoords: List[RebarCoords], changeList: List[tuple[int, Rebar|GRebars]]):
    output = rebarCoords
    for item in changeList:
            output[item[0]].rebar = item[1]
    return output
# endregion

