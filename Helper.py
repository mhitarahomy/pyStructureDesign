from matplotlib import pyplot as plt
import numpy as np
from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFDmethod.PMManalysis as PMManalysis
import pyCivilDesign.concreteDesign.LRFDmethod.Manalysis as Manalysis
import pyCivilDesign.concreteDesign.LRFD as LRFD
import pyCivilDesign.concreteDesign.designProps as props
from shapely.plotting import plot_polygon, plot_line, plot_points


d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
# showSection(sct)
# print(sct.section)

x_coords = [_coord.x for _coord in sct.Coords]
y_coords = [_coord.y for _coord in sct.Coords]
max_x = np.max(x_coords)
max_y = np.max(y_coords)
top_edge = [coord for coord in sct.Coords if coord.y == max_y]


fig, axs = plt.subplots()
plot_polygon(sct.section, add_points=False)
plot_points(top_edge)
# plot_polygon(PMManalysis.calc_neutral_region(data.section, 638, 45))
# plot_polygon(PMManalysis.calc_pressure_region(data.section, data.fc, 638, 45))
# plot_line(PMManalysis.calc_neutral_axis(data.section, 638, 45))
# plot_line(PMManalysis.calc_pressure_axis(data.section, data.fc, 638, 45))
plt.show()