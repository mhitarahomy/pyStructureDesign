
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
import pyCivilDesign.concreteDesign.LRFDmethod.PMMSolver as PMMsolver
from pyCivilDesign.concreteDesign.designProps import Assumptions, DesignData, defaultAssumption


@dataclass
class PMMresults():
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
    ratio: np.float32|None = field(default= None)


def PMRatio(data: DesignData, P: float, Mx: float, My: float, 
                assump:Assumptions=defaultAssumption) -> PMMresults:
    _ratio = PMMsolver.CalcPMRatio(data, P, Mx, My, assump)
    _angle = float(PMMsolver.AngleFromForces(data, P, Mx, My, assump))
    _P = P / _ratio
    _M, _Mx, _My, _alpha = PMMsolver.CalcMn(data, _P, _angle, assump) # type: ignore
    _c = float(PMMsolver.C(data, _P, _angle, assump)) # type: ignore
    _es = PMMsolver.es(data, _c, _angle, assump)
    _fs = PMMsolver.fs(data, _c, _angle, assump)
    _Fs = PMMsolver.Fs(data, _c, _angle, assump)
    _Cc = PMMsolver.Cc(data, _c, _angle, assump)
    _percent = PMMsolver.getAsPercent(data)
    return PMMresults(_c, _angle, _alpha, _es, _fs, _Fs, _Cc, _P, _M, _Mx, _My, _percent, "", _ratio) # type: ignore


def Mn(data: DesignData, P: float, angle: float, 
           assump:Assumptions=defaultAssumption) -> PMMresults:
    c = float(PMMsolver.C(data, P, angle, assump))
    _M, _Mx, _My = PMMsolver.M(data, c, angle, assump)
    alpha = PMMsolver.Alpha(data, _Mx, _My, assump) # type: ignore
    _es = PMMsolver.es(data, c, angle, assump)
    _fs = PMMsolver.fs(data, c, angle, assump)
    _Fs = PMMsolver.Fs(data, c, angle, assump)
    _Cc = PMMsolver.Cc(data, c, angle, assump)
    _percent = PMMsolver.getAsPercent(data)
    return PMMresults(c, angle, alpha, _es, _fs, _Fs, _Cc, P, _M, _Mx, _My, _percent, "") # type: ignore
    

def AsPercent(data:DesignData, P: float, Mx: float, My: float,
                assump: Assumptions=defaultAssumption) -> PMMresults:
    _percent = PMMsolver.CalcAsPercent(data, P, Mx, My, assump)
    data = PMMsolver.setAsPercent(data, _percent) # type: ignore
    angle = float(PMMsolver.AngleFromForces(data,  P, Mx, My, assump))
    c = float(PMMsolver.C(data, P, angle, assump))
    _M, _Mx, _My = PMMsolver.M(data, c, angle, assump)
    alpha = PMMsolver.Alpha(data, _Mx, _My, assump) # type: ignore
    _es = PMMsolver.es(data, c, angle, assump)
    _fs = PMMsolver.fs(data, c, angle, assump)
    _Fs = PMMsolver.Fs(data, c, angle, assump)
    _Cc = PMMsolver.Cc(data, c, angle, assump)
    return PMMresults(c, angle, alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "") # type: ignore
