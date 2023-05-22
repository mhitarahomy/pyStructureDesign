from dataclasses import dataclass
from pyCivilDesign.sections.concreteSections import ConcreteSct
import numpy as np
from numpy.typing import NDArray
from shapely import Polygon, Point


@dataclass(kw_only=True)
class DesignData():
    section: Polygon
    bw: np.float32
    d: np.float32
    fy: np.float32
    fyt: np.float32
    fc: np.float32
    Coords: NDArray[Point]
    As: NDArray[np.float32]
    Es: np.float32

    @classmethod
    def fromSection(cls, section: ConcreteSct):
        return setDesignDataFromSection(section)


def setDesignDataFromSection(sct: ConcreteSct) -> DesignData:
    return DesignData(
        section=sct.section, bw=np.float32(sct.bw), d=np.float32(sct.d), fy=np.float32(sct.lBarMat.fy), 
        fyt=np.float32(sct.cBarMat.fy),fc=np.float32(sct.concMat.fc), 
        Coords=np.array(sct.Coords), As=np.array(sct.As), Es=np.float32(sct.lBarMat.Es))
