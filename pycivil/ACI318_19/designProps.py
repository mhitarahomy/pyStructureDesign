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
    conf_dist: np.float32
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
        Coords=np.array(sct.Coords), As=np.array(sct.As, dtype=np.float32), 
        Es=np.float32(sct.lBarMat.Es), cover=np.float32(sct.cover), Av=np.array(sct.Av, dtype=np.float32),
        conf_type=sct.conf_type, conf_dist=np.float32(sct.conf_dist))


@dataclass
class PMMresults():
    P0: np.float32
    Pn_max: np.float32
    Pnt_max: np.float32
    M_max: np.float32
    c: np.float32
    angle: np.float32
    alpha: np.float32
    es: NDArray[np.float32]
    fs: NDArray[np.float32]
    Fs: NDArray[np.float32]
    Cc: np.float32
    P: np.float32
    M: np.float32
    Mx: np.float32
    My: np.float32
    AsPercent: np.float32
    msg: str
    ratio: np.float32

    def __repr__(self) -> str:
        return f'''c = {round(self.c, 0)} mm
angle = {round(self.angle, 2)} degree
alpha = {round(self.alpha, 2)} degree
Cc = {round(self.Cc, 0)} N   P = {round(self.P/1000, 0)} kN
M = {round(self.M/1000000, 0)} kN.m   Mx = {round(self.Mx/1000000, 0)} kN.m   My = {round(self.My/1000000, 0)} kN.m
AsPercent = {round(self.AsPercent, 2)} %
msg = {self.msg}
ratio = {round(self.ratio, 2)}'''