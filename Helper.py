from pyCivilDesign.sections.concreteSections import RectConcreteSct, showSection
from pyCivilDesign.concreteDesign.designAssumptions import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct

from shapely import LineString, Polygon, Point

C40x60 = RectConcreteSct(b=400, h=600)
C40x60.rebarCoords = Rsct.RectSectRebars(C40x60.section, 3, 4, Rsct.d20, 50)
data = setDesignDataFromSection(C40x60)

data = Canalysis.setAsPercent(data, 4)
# print(Canalysis.CalcPMmax(data, 0))

showSection(C40x60)
