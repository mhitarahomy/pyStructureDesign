from dataclasses import dataclass
from pycivil.sections.concreteSections import ConcreteSct
import numpy as np
from numpy.typing import NDArray
from shapely import Polygon, Point

from pycivil.sections.rebarSections import ConfType


@dataclass(kw_only=True)
class DesignData():
    section: Polygon
    bw: np.float32
    fy: np.float32
    fyt: np.float32
    fc: np.float32
    Coords: NDArray[Point]
    As: NDArray[np.float32]
    Av: NDArray[np.float32]
    Es: np.float32
    cover: np.float32
    conf_type: ConfType

    @classmethod
    def fromSection(cls, section: ConcreteSct):
        return setDesignDataFromSection(section)


def setDesignDataFromSection(sct: ConcreteSct) -> DesignData:
    return DesignData(
        section=sct.section, bw=np.float32(sct.bw), fy=np.float32(sct.lBarMat.fy), 
        fyt=np.float32(sct.cBarMat.fy),fc=np.float32(sct.concMat.fc), 
        Coords=np.array(sct.Coords, dtype=np.float32), As=np.array(sct.As, dtype=np.float32), 
        Es=np.float32(sct.lBarMat.Es), cover=np.float32(sct.cover), Av=np.array(sct.Av, dtype=np.float32),
        conf_type=sct.conf_type)
