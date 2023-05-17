from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFD as LRFDsolver

d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
designData = setDesignDataFromSection(sct)
# print(LRFDsolver.CalcMn(designData, 0, 0))
# print(LRFDsolver.CalcPMRatio(designData, 0, 5e8, 0))
print(LRFDsolver.AsPercent(designData, 1e6, 2e8, 3e8))
# showSection(sct)
