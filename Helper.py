from matplotlib import pyplot as plt
import numpy as np
from pycivil.sections.concreteSections import showSection, RectConcreteSct
import pycivil.sections.rebarSections as Rsct
import pycivil.sections.section as Sct
import pycivil.ACI318_19.PMManalysis as PMManalysis
import pycivil.ACI318_19.Manalysis as Manalysis
import pycivil.ACI318_19.LRFD as LRFD
import pycivil.ACI318_19.designProps as props
from shapely.plotting import plot_polygon, plot_line, plot_points


d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
showSection(sct)



# fig, axs = plt.subplots()
# plot_polygon(sct.section, add_points=False)
# plot_points(top_edge)
# plot_polygon(PMManalysis.calc_neutral_region(data.section, 638, 45))
# plot_polygon(PMManalysis.calc_pressure_region(data.section, data.fc, 638, 45))
# plot_line(PMManalysis.calc_neutral_axis(data.section, 638, 45))
# plot_line(PMManalysis.calc_pressure_axis(data.section, data.fc, 638, 45))
# plt.show()