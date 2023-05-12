from dataclasses import dataclass, field
from typing import Any, List, Tuple
from matplotlib import pyplot as plt

from shapely import Polygon, Point
from shapely.plotting import plot_line, plot_points, plot_polygon

from pyCivilDesign.materials import ConcreteMat, RebarMat, AIII, AII, C25
from pyCivilDesign.sections.rebarSections import RebarCoords
import pyCivilDesign.sections.section as Sct 


ListOfPoints = List[Tuple[float, float]]
C25def = lambda: C25
AIIdef = lambda: AII
AIIIdef = lambda: AIII


@dataclass()
class ConcreteSct():
    section: ListOfPoints
    sectionType: Sct.SectionType
    concMat: ConcreteMat
    lBarMat: RebarMat 
    cBarMat: RebarMat
    rebarCoords: List[RebarCoords] = field(default_factory=list)
    
    @property
    def sectionCoords(self) -> ListOfPoints:
        return list(Polygon(self.section).exterior.coords)
    
    @property
    def Coords(self) -> ListOfPoints:
        return [rcoord.point for rcoord in self.rebarCoords]
    
    @property
    def As(self) -> List[float]:
        return [rcoord.rebar.Area for rcoord in self.rebarCoords]


@dataclass()
class RectConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: ListOfPoints = field(init=False)
    sectionType: Sct.SectionType = field(init=False, default=Sct.SectionType.Rectangle)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    rebarCoords: List[RebarCoords] = field(default_factory=list)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h":
            self.section =  Sct.RectangleSct(self.b, self.h)


@dataclass()
class TrapzoidConcreteSct(ConcreteSct):
    b1: float = field(default=400)
    b2: float = field(default=300)
    h: float = field(default=600)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: ListOfPoints = field(init=False)
    sectionType: Sct.SectionType = field(init=False, default=Sct.SectionType.Rectangle)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    rebarCoords: List[RebarCoords] = field(default_factory=list)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b1" or __name=="b2" or __name=="h":
            self.section =  Sct.TrapzoidSct(self.b1, self.b2, self.h)


@dataclass()
class TShapeConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    tf: float = field(default=100)
    tw: float = field(default=100)
    tf1: float|None = field(default=None)
    tw1: float|None = field(default=None)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: ListOfPoints = field(init=False)
    sectionType: Sct.SectionType = field(init=False, default=Sct.SectionType.TShape)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    rebarCoords: List[RebarCoords] = field(default_factory=list)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h" or __name=="tf" or \
            __name=="tw" or __name=="tf1" or __name=="tw1":
            self.section =  Sct.TShapeSct(self.b, self.h, self.tf, self.tw,
                                          self.tf1, self.tw1)
        

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
