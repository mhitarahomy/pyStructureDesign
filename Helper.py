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


ListOfPoints = List[Tuple[float, float]]

section = SCT.TrapzoidSct(b1=400, b2=600, h=700)

print(section)
fig = plt.figure(dpi=90)
ax = fig.add_subplot()
plot_polygon(Polygon(section))
plot_line(LineString(SCT.Edge(section, "left")), add_points=False, color="red")
# plot_points(endPoint)
plt.show()
