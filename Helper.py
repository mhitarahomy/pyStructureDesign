from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFD as LRFDsolver
import pyCivilDesign.concreteDesign.LRFDmethod.PMMSolver as PMMsolver

d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
data = setDesignDataFromSection(sct)
print(PMMsolver.getAsPercent(data))

data = PMMsolver.setAsPercent(data, 2.03359)
result = LRFDsolver.PMRatio(data, 3e6, 6e8, 5e8)

LRFDsolver.showResult(data, result)

# print(PMMsolver.CalcAsPercent(data, 3e6, 6e8, 5e8))
