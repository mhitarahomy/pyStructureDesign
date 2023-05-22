from dataclasses import dataclass
from typing import Callable
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


beta1: Callable[[DesignData], np.float32] = lambda data: \
    np.float32(0.85) if data.fc<=28 else np.max([0.85-(0.05*(data.fc-28)/7), np.float32(0.65)])


def phif(data: DesignData, esMin: np.float32) -> np.float32:
    ety = data.fy / data.Es
    et = abs(esMin)
    return np.float32(0.65) if et <= ety else 0.65+0.25*((et-ety)/0.003)\
          if ety < et < ety+0.003 else np.float32(0.9)


@dataclass()
class Assumptions():
    beta1: Callable[[DesignData], np.float32] = beta1
    phif: Callable[[DesignData, np.float32], np.float32] = phif
    phis: np.float32 = np.float32(0.9)
    phia: np.float32 = np.float32(0.65)
    phit: np.float32 = np.float32(0.6)
    ecu: np.float32 = np.float32(0.003)


defaultAssumption = Assumptions()


def setDesignDataFromSection(sct: ConcreteSct) -> DesignData:
    return DesignData(
        section=sct.section, bw=np.float32(sct.bw), d=np.float32(sct.d), fy=np.float32(sct.lBarMat.fy), 
        fyt=np.float32(sct.cBarMat.fy),fc=np.float32(sct.concMat.fc), 
        Coords=np.array(sct.Coords), As=np.array(sct.As), Es=np.float32(sct.lBarMat.Es))
