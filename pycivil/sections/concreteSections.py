from dataclasses import dataclass, field
from typing import Any, List, Tuple
from matplotlib import pyplot as plt

from shapely import Polygon, Point
from shapely.plotting import plot_points, plot_polygon
from pycivil.errors import RebarCoordsError, SectionError

from pycivil.materials import ConcreteMat, RebarMat, AIII, AII, C25
from pycivil.sections.rebarSections import GRebars, Rebar, RebarCoords, ConfType, ConfRebars
import pycivil.sections.section as Sct 


# * This lambdas is for use in field default_factory 
ListOfPoints = List[Tuple[float, float]]
C25def = lambda: C25
AIIdef = lambda: AII
AIIIdef = lambda: AIII
ConfRebarsdef = lambda: ConfRebars(0, [])


@dataclass()
class ConcreteSct():
    section: Polygon
    section_type: Sct.SectionType
    bw: float
    concMat: ConcreteMat = field(default_factory=C25def)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    cover: float = field(default=70)
    rebarCoords: List[RebarCoords] = field(default_factory=list)
    conf_rebars: ConfRebars = field(default_factory=ConfRebarsdef)
    conf_type: ConfType = field(default=ConfType.Tie)
    
    @property
    def sectionCoords(self) -> ListOfPoints:
        return list(self.section.exterior.coords)
    
    @property
    def Coords(self) -> List[Point]:
        return [rcoord.point for rcoord in self.rebarCoords]
    
    @property
    def As(self) -> List[float]:
        return [rcoord.rebar.Area for rcoord in self.rebarCoords]
    
    @property
    def Av(self) -> List[float]:
        return [rebar.Area for rebar in self.conf_rebars.rebars]
    
    @property
    def conf_dist(self) -> float:
        return self.conf_rebars.deistance


@dataclass()
class RectConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    bw: float = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.Rectangle)
    conf_type: ConfType = field(init=False, default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type" or __name=="conf_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h":
            self.section =  Sct.RectangleSct(self.b, self.h)
            self.bw = self.b


@dataclass()
class TrapzoidConcreteSct(ConcreteSct):
    b1: float = field(default=400)
    b2: float = field(default=300)
    bw: float|None = field(init=False)
    h: float = field(default=600)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.Trapezoid)
    conf_type: ConfType = field(init=False, default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type" or __name=="conf_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="b1" or __name=="b2" or __name=="h":
            self.section =  Sct.TrapzoidSct(self.b1, self.b2, self.h)
            self.bw = (self.b1 + self.b2)/2


@dataclass()
class TShapeConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    tf: float = field(default=100)
    tw: float = field(default=100)
    tf1: float|None = field(default=None)
    tw1: float|None = field(default=None)
    bw: float|None = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.TShape)
    conf_type: ConfType = field(init=False, default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type" or __name=="conf_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h" or __name=="tf" or \
            __name=="tw" or __name=="tf1" or __name=="tw1":
            self.section =  Sct.TShapeSct(self.b, self.h, self.tf, self.tw,
                                          self.tf1, self.tw1)
            self.bw = self.tw if self.tw1==None else (self.tw + self.tw1) / 2
        
        
@dataclass()
class LShapeConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    tf: float = field(default=100)
    tw: float = field(default=100)
    tf1: float|None = field(default=None)
    tw1: float|None = field(default=None)
    bw: float|None = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.LShape)
    conf_type: ConfType = field(init=False, default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type" or __name=="conf_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h" or __name=="tf" or \
            __name=="tw" or __name=="tf1" or __name=="tw1":
            self.section =  Sct.LShapeSct(self.b, self.h, self.tf, self.tw,
                                          self.tf1, self.tw1)
            self.bw = self.tw if self.tw1==None else (self.tw + self.tw1) / 2


@dataclass()
class BoxConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    th: float = field(default=50)
    bw: float = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.Box)
    conf_type: ConfType = field(init=False, default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type" or __name=="conf_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h" or __name=="th":
            self.section =  Sct.BoxSct(self.b, self.h, self.th)
            self.bw = 2 * self.th


@dataclass()
class CircConcreteSct(ConcreteSct):
    d: float = field(default=400)
    bw: float = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.Circle)
    conf_type: ConfType = field(default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="d":
            self.section =  Sct.CircleSct(self.d)
            self.bw = 0.8 * self.d


@dataclass()
class PipeConcreteSct(ConcreteSct):
    d: float = field(default=400)
    th: float = field(default=50)
    bw: float = field(init=False)
    section: Polygon = field(init=False)
    section_type: Sct.SectionType = field(init=False, default=Sct.SectionType.Pipe)
    conf_type: ConfType = field(default=ConfType.Tie)

    def __setattr__(self, __name: str, __value: Any) -> None:
        if __name=="section_type":
            raise ValueError(f"you can't change {__name}")
        super().__setattr__(__name, __value)
        if __name=="d" or __name=="th":
            self.section =  Sct.PipeSct(self.d, self.th)
            self.bw = 0.8 * 2*self.th


def showSection(concSct: ConcreteSct) -> None:
    fig = plt.figure(dpi=90)
    ax = fig.add_subplot()  # type: ignore
    plot_polygon(concSct.section, add_points=False, linewidth=1)
    if concSct.rebarCoords != None:
        plot_points([Point(rcoord.point) for rcoord in concSct.rebarCoords], color='black')
        for i in range(len(concSct.rebarCoords)):
            rcoord = concSct.rebarCoords[i]
            ax.annotate(f"{rcoord.rebar}\n{i}", xy=(rcoord.point.x, rcoord.point.y), xycoords='data',
                    xytext=(1.5, -10), textcoords='offset points')
    plt.show()
    
    
def check_As(As) -> None:
    if len(As) == 0: raise RebarCoordsError(
        "list of rebars area is empty. you must add rebar to section.")
    
    
def check_rebar_coords(coords) -> None:
    if len(coords) == 0: raise RebarCoordsError(
        "list of rebars coordination is empty. you must add rebar to section.")


def check_section(section) -> None:
    if len(list(section.exterior.coords)) < 4: raise SectionError("section must have greater than 3 edge")