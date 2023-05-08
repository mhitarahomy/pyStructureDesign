from shapely import Point, Polygon
from pyCivilDesign.sections.section import ListOfPoints

li = [(1,2), (3,4), (5,6), (7,8)]
print(ListOfPoints(Polygon(li).exterior.coords))