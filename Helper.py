from typing import List, Tuple
from matplotlib import pyplot as plt
from pyCivilDesign.sections.concreteSections import RectConcreteSct, TShapeConcreteSct, TrapzoidConcreteSct, showSection
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as SCT
import numpy as np
from shapely import LineString, Point, Polygon
from shapely.plotting import plot_line, plot_points, plot_polygon


B300x400 = RectConcreteSct(b=300, h=400)
B300x400.rebarCoords = Rsct.RectSectRebars(B300x400.section, 4, 2, Rsct.d20, 50)
print(B300x400.bw)
showSection(B300x400)
# ListOfPoints = List[Tuple[float, float]]

# section = SCT.TrapzoidSct(b1=400, b2=600, h=700)

# print(section)
# fig = plt.figure(dpi=90)
# ax = fig.add_subplot()
# plot_polygon(Polygon(section))
# plot_line(LineString(SCT.Edge(section, "left")), add_points=False, color="red")
# # plot_points(endPoint)
# plt.show()
