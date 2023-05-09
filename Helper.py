from pyCivilDesign.sections.concreteSections import RectConcreteSct
from pyCivilDesign.concreteDesign.designAssumptions import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct

C40x60 = RectConcreteSct(b=400, h=600)
C40x60.rebarCoords = Rsct.RectSectRebars(C40x60.section, 5, 6, Rsct.d20, 50)
# print(C40x60.rebarCoords)
data = setDesignDataFromSection(C40x60)

print(Canalysis.P0(data))