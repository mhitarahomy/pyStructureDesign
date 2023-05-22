from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFDmethod.PMManalysis as PMManalysis
import pyCivilDesign.concreteDesign.LRFDmethod.Manalysis as Manalysis
import pyCivilDesign.concreteDesign.designProps as props


d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
# showSection(sct)

data = props.setDesignDataFromSection(sct)

print(PMManalysis.calc_As_percent(data, 3e6, 5e8, 3e8))
# print(Manalysis.calc_c(data))
