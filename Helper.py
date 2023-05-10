from pyCivilDesign.sections.concreteSections import RectConcreteSct, TShapeConcreteSct, showSection
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct
import numpy as np

T600x400 = TShapeConcreteSct(b=600, h=400, tf=80, tw=200, tf1=100)
d20 = Rsct.d20
T600x400.rebarCoords = Rsct.LinearRebars([d20, d20, d20, d20], 50, 50, [(250,100),(100,78),(-100,78),(-250,100)])
showSection(T600x400)
# C40x60.rebarCoords = Rsct.RectSectRebars(C40x60.section, 5, 6, Rsct.d20, 50)
# print(C40x60.Coords)
# data = setDesignDataFromSection(C40x60)

# print(Canalysis.CalcPercent(data, 4e6, 3e8, 3e8))
# Canalysis.showResult(data, 3e6, 3e8, 2e8)

# a=np.array([(1,2), (3,4), (5,6), (7,8)])
# print(a[:,1])
