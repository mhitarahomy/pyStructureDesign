from matplotlib import pyplot as plt
from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFDmethod.PMManalysis as PMManalysis
import pyCivilDesign.concreteDesign.LRFDmethod.Manalysis as Manalysis
import pyCivilDesign.concreteDesign.LRFD as LRFD
import pyCivilDesign.concreteDesign.designProps as props
from shapely.plotting import plot_polygon, plot_line


d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
# showSection(sct)

data = props.setDesignDataFromSection(sct)

# Fs = PMManalysis.calc_es(data.section, data.Coords, 1750, 0)
# print(Fs)
# print(PMManalysis.calc_Mn(data, 4e6, 0))
# print(PMManalysis.calc_PM_ratio(data, 4e6, 4e8, 2e8, angle))
# print(PMManalysis.calc_percent(data, 4e6, 5e8, 4e8))
data_p = PMManalysis.set_As_percent(data, 1.691153300645965)
print(PMManalysis.calc_PM_ratio(data_p, 4e6, 4e8, 2e8))

# print(PMManalysis.calc_P(data, -226.5, 31))
# print(PMManalysis.calc_c(data, 4e6, 45))



# print(Manalysis.calc_c(data))

# print(LRFD.PMM_analyze(sct, 3e6, 5e8, 3e8))
# LRFD.show_PMM_analysis_result(sct, 2e6, 5e8, 3e8)


# fig, axs = plt.subplots()
# plot_polygon(PMManalysis.rotate_section(data.section, 45))
# plot_polygon(PMManalysis.calc_neutral_region(data.section, 638, 45))
# plot_polygon(PMManalysis.calc_pressure_region(data.section, data.fc, 638, 45))
# plot_line(PMManalysis.calc_neutral_axis(data.section, 638, 45))
# plot_line(PMManalysis.calc_pressure_axis(data.section, data.fc, 638, 45))
# plt.show()