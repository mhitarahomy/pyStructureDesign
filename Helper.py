from typing import List, Tuple
from matplotlib import pyplot as plt
from pyCivilDesign.sections.concreteSections import showSection, BoxConcreteSct, CircConcreteSct, PipeConcreteSct
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as SCT
import numpy as np
from shapely import LineString, Point, Polygon
from shapely.plotting import plot_line, plot_points, plot_polygon

d20 = Rsct.d20

# sct = TrapzoidConcreteSct(b1=500, b2=600, h=700)
# sct = TShapeConcreteSct(b=700, h=600, tf=100, tw=150, tf1=130, tw1=170)
# sct = LShapeConcreteSct(b=600, h=400, tf=100, tw=200)
# sct = BoxConcreteSct(b=500, h=600, th=100)
# sct = CircConcreteSct(d=500)
sct = PipeConcreteSct(d=500, th=100)
# sct.rebarCoords = Rsct.LinearRebars([d20 for i in range(5)], 50, 50, SCT.DistanceFrom(sct.section, 50))
# sct.rebarCoords.extend(Rsct.LinearRebars([d20 for i in range(5)], 50, 50, SCT.DistanceFrom(sct.section, 50, "bottom")))
sct.rebarCoords = Rsct.CircularBars([d20 for i in range(9)], 50, sct.d)
showSection(sct)
