from matplotlib import pyplot as plt
from pyCivilDesign.sections.concreteSections import RectConcreteSct, TShapeConcreteSct, showSection
from pyCivilDesign.concreteDesign.designProps import setDesignDataFromSection
import pyCivilDesign.concreteDesign.concreteAnalysis as Canalysis 
import pyCivilDesign.sections.rebarSections as Rsct
import numpy as np
from shapely import LineString, Point
from shapely.plotting import plot_line, plot_points

# T600x400 = TShapeConcreteSct(b=600, h=400, tf=80, tw=200, tf1=100)
# d20 = Rsct.d20
# T600x400.rebarCoords = Rsct.LinearRebars([d20, d20, d20, d20], 50, 50, [(250,100),(100,78),(-100,78),(-250,100)])
# showSection(T600x400)
listOfPoints = [(0,0), (1,1), (2,2), (3,3), (4,4)]
startCover = 2
endCover = 2
line = LineString(listOfPoints)
startPoint = line.interpolate(startCover)
endPoint = line.interpolate(line.length-endCover)
_points =[listOfPoints[i] for i in range(1,len(listOfPoints)-1,1) if LineString(listOfPoints[:i+1]).length<(line.length-endCover)]
_points.insert(0, listOfPoints[0])
_points.extend([(endPoint.x, endPoint.y)])
points = [_points[i+1] for i in range(len(_points)-1) if LineString(_points[:i+2]).length>startCover]
points.insert(0, (startPoint.x, startPoint.y))
print(points)
fig = plt.figure(dpi=90)
ax = fig.add_subplot()
plot_line(LineString(listOfPoints), add_points=False)
plot_line(LineString(points), add_points=False, color="red")
plot_points(startPoint)
plot_points(endPoint)
plt.show()
