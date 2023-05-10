from pyCivilDesign.sections.concreteSections import RectConcreteSct
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct
import numpy as np

C40x60 = RectConcreteSct(b=400, h=600)
C40x60.rebarCoords = Rsct.RectSectRebars(C40x60.section, 5, 6, Rsct.d20, 50)
# print(C40x60.Coords)
data = setDesignDataFromSection(C40x60)

print(Canalysis.CalcPercent(data, 4e6, 3e8, 3e8))

# a=np.array([(1,2), (3,4), (5,6), (7,8)])
# print(a[:,1])
