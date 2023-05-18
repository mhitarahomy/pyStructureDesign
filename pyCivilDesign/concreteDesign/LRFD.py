
from dataclasses import dataclass, field
from matplotlib import pyplot as plt

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
    ratio: np.float32

    def __repr__(self) -> str:
        return f'''c = {round(self.c, 0)} mm
angle = {round(self.angle, 2)} degree
alpha = {round(self.alpha, 2)} degree
es = {np.round(self.es, 5)}
fs = {np.round(self.fs, 0)} N/mm2
Fs = {np.round(self.Fs, 0)} N
Cc = {round(self.Cc, 0)} N   P = {round(self.P/1000, 0)} kN
M = {round(self.M/1000000, 0)} kN.m   Mx = {round(self.Mx/1000000, 0)} kN.m   My = {round(self.My/1000000, 0)} kN.m
AsPercent = {round(self.AsPercent, 2)} %
msg = {self.msg}
ratio = {round(self.ratio, 2)}'''


def PMRatio(data: DesignData, P: float, Mx: float, My: float, 
                assump:Assumptions=defaultAssumption) -> PMMresults:
    _ratio = PMMsolver.CalcPMRatio(data, P, Mx, My, assump)
    _M = pow(Mx**2 + My**2, 0.5)
    _angle = float(PMMsolver.AngleFromForces(data, P, Mx, My, assump))
    _alpha = PMMsolver.Alpha(data, Mx, My, assump)
    _c = float(PMMsolver.C(data, P, _angle, assump)) # type: ignore
    _es = PMMsolver.es(data, _c, _angle, assump)
    _fs = PMMsolver.fs(data, _c, _angle, assump)
    _Fs = PMMsolver.Fs(data, _c, _angle, assump)
    _Cc = PMMsolver.Cc(data, _c, _angle, assump)
    _percent = PMMsolver.getAsPercent(data)
    return PMMresults(_c, _angle, _alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "", _ratio) # type: ignore


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
    return PMMresults(c, angle, alpha, _es, _fs, _Fs, _Cc, P, _M, _Mx, _My, _percent, "", 0) # type: ignore
    

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
    return PMMresults(c, angle, alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "", 0) # type: ignore


def showResult(data: DesignData, result: PMMresults, assump: Assumptions=defaultAssumption):
    angle = PMMsolver.AngleFromForces(data, result.P, result.Mx, result.My, assump) # type: ignore
    PointNums = 20
    Paxis = np.linspace(PMMsolver.PtMax(data), PMMsolver.P0(data), PointNums, endpoint=True)
    Maxis = np.array([PMMsolver.CalcMn(data, p, angle, assump)[0] for p in Paxis]) # type: ignore

    AlphaNums = 21
    Alphas = np.linspace(0, 360, AlphaNums)
    M =  np.array([PMMsolver.CalcMn(data, result.P, alpha, assump) for alpha in Alphas]) # type: ignore
    MxAxis =  np.array([m[1] for m in M])
    MyAxis =  np.array([m[2] for m in M])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17,8))
    ax1.plot(Maxis, Paxis, linewidth=2.0)
    ax1.plot(result.M, result.P, '.', color="black", markersize=7)
    ax1.annotate(f"P={round(result.P/1000, 2)}kN \nM={round(result.M/1000000, 2)}kN.m \nratio={round(result.ratio, 2)}",
                 (result.M, result.P), 
                 textcoords="offset points", xytext=(5,0), ha="left")
    ax1.set_title(f"P-M chart for {round(result.alpha, 2)} degree")
    ax1.set_xlabel("M (N.mm)")
    ax1.set_ylabel("P (N)")
    ax1.axhline(y=0, color='b', linestyle='-')
    ax1.axvline(x=0, color='b', linestyle='-')
    ax1.grid(True)
    
    ax2.plot(MxAxis, MyAxis, linewidth=2.0)
    ax2.plot(result.Mx, result.My, '.', color="black", markersize=7)
    ax2.annotate(f"Mx={round(result.Mx/1000000, 2)}kN.m \nMy={round(result.My/1000000, 2)}kN.m", (result.Mx, result.My), 
                 textcoords="offset points", xytext=(5,0), ha="left")
    ax2.set_title(f"Mx-My chart on P={round(result.P/1000)} kN")
    ax2.set_xlabel("Mx (N.mm)")
    ax2.set_ylabel("My (N.mm)")
    ax2.axhline(y=0, color='b', linestyle='-')
    ax2.axvline(x=0, color='b', linestyle='-')
    ax2.axis("equal")
    ax2.grid(True)
    
    plt.show()
