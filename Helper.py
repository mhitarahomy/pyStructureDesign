from pyCivilDesign.sections.concreteSections import showSection, RectConcreteSct
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.sections.rebarSections as Rsct
import pyCivilDesign.sections.section as Sct
import pyCivilDesign.concreteDesign.LRFD as LRFDsolver
import pyCivilDesign.concreteDesign.LRFDmethod.PMMSolver as PMMsolver

d20 = Rsct.d20

sct = RectConcreteSct(b=500, h=700)
sct.rebarCoords = Rsct.RectRebarsSct(sct.section, 5, 4, d20, 50)
print(LRFDsolver.PMM_design(sct, 3e6, 6e8, 5e8))

# print(PMMsolver.CalcAsPercent(data, 3e6, 6e8, 5e8))
