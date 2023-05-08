from dataclasses import dataclass, field
from typing import Any, List
from matplotlib import pyplot as plt

from shapely import Polygon, Point
from shapely.plotting import plot_line, plot_points, plot_polygon

from pyCivilDesign.materials import ConcreteMat, RebarMat, AIII, AII, C25
from pyCivilDesign.sections.rebarSections import RebarCoords
import pyCivilDesign.sections.section as Sct 


C25def = lambda: C25
AIIdef = lambda: AII
AIIIdef = lambda: AIII


@dataclass()
class ConcreteSct():
    section: Sct.ListOfPoints
    sectionType: Sct.SectionType
    concMat: ConcreteMat
    lBarMat: RebarMat 
    cBarMat: RebarMat
    rebarCoords: List[RebarCoords] = field(default_factory=list)
    
    @property
    def sectionCoords(self) -> Sct.ListOfPoints:
        return Sct.ListOfPoints(list(Polygon(self.section).exterior.coords))
    
    @property
    def Coords(self) -> Sct.ListOfPoints:
        return Sct.ListOfPoints([rcoord.point for rcoord in self.rebarCoords])
    
    @property
    def As(self) -> List[float]:
        return [rcoord.rebar.Area for rcoord in self.rebarCoords]


@dataclass()
class RectConcreteSct(ConcreteSct):
    b: float = field(default=400)
    h: float = field(default=600)
    concMat: ConcreteMat = field(default_factory=C25def)
    section: Sct.ListOfPoints = field(init=False)
    sectionType: Sct.SectionType = field(init=False, default=Sct.SectionType.Rectangle)
    lBarMat: RebarMat = field(default_factory=AIIIdef)
    cBarMat: RebarMat = field(default_factory=AIIdef)
    rebarCoords: List[RebarCoords] = field(default_factory=list)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if __name=="b" or __name=="h":
            self.section =  Sct.RectangleSct(self.b, self.h)
    

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
