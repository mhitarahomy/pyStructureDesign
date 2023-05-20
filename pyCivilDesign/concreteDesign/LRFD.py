
from dataclasses import dataclass, field
from matplotlib import pyplot as plt

import numpy as np
from numpy.typing import NDArray
import pyCivilDesign.concreteDesign.LRFDmethod.PMMSolver as PMMsolver
from pyCivilDesign.concreteDesign.LRFDmethod.PMMSolver import Assumptions, DesignData, defaultAssumption
from pyCivilDesign.sections.concreteSections import ConcreteSct


@dataclass
class PMMresults():
    P0: np.float32
    Pnmax: np.float32
    Ptmax: np.float32
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


def PMM_analyze(section: ConcreteSct, P: float, Mx: float, My: float, 
                   assump: Assumptions=defaultAssumption):
    data = DesignData.fromSection(section)
    _ratio = PMMsolver.CalcPMRatio(data, P, Mx, My, assump)
    _M = pow(Mx**2 + My**2, 0.5)
    _angle = float(PMMsolver.AngleFromForces(data, P, Mx, My, assump))
    _alpha = PMMsolver.Alpha(data, Mx, My, assump)
    _c = float(PMMsolver.C(data, P, _angle, assump)) # type: ignore
    _es = PMMsolver.es(data, _c, _angle, assump)
    _fs = PMMsolver.fs(data, _c, _angle, assump)
    _Fs = PMMsolver.Fs(data, _c, _angle, assump)
    _Cc = PMMsolver.Cc(data, _c, _angle, assump)
    _percent = PMMsolver.AsPercent(data)
    return PMMresults(_c, _angle, _alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "", _ratio) # type: ignore
        

def PMM_design(section: ConcreteSct, P: float, Mx: float, My: float,
                assump: Assumptions=defaultAssumption):
    data = DesignData.fromSection(section)
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


def show_PMM_analysis_result(section: ConcreteSct, P: float, Mx: float, My: float,
                              assump: Assumptions=defaultAssumption):
    data = DesignData.fromSection(section)
    M = pow(Mx**2 + My**2, 0.5)
    ratio = PMMsolver.CalcPMRatio(data, P, Mx, My, assump)

    angle = PMMsolver.AngleFromForces(data, P, Mx, My, assump) # type: ignore
    alpha = PMMsolver.Alpha(data, Mx, My, assump)
    PointNums = 20
    Paxis = np.linspace(PMMsolver.PtMax(data), PMMsolver.P0(data), PointNums, endpoint=True)
    Maxis = np.array([PMMsolver.CalcMn(data, p, angle, assump)[0] for p in Paxis]) # type: ignore

    AlphaNums = 21
    Alphas = np.linspace(0, 360, AlphaNums)
    _M =  np.array([PMMsolver.CalcMn(data, P, alpha, assump) for alpha in Alphas]) # type: ignore
    MxAxis =  np.array([m[1] for m in _M])
    MyAxis =  np.array([m[2] for m in _M])

    sectionCoords = np.array(section.section.exterior.coords)
    sectionXcoords = sectionCoords[:,0]
    sectionYcoords = sectionCoords[:,1]
    
    fig, axs = plt.subplots(2, 2, figsize=(10,8))
    axs[0, 0].plot(sectionXcoords, sectionYcoords)
    axs[0, 0].set_aspect('equal')
    axs[0, 0].fill(sectionXcoords, sectionYcoords)
    axs[0, 0].axhline(y=0, color='k', linestyle='-')
    axs[0, 0].axvline(x=0, color='k', linestyle='-')
    axs[0, 0].grid(True)
    
    axs[0, 1].plot(sectionXcoords, sectionYcoords)
    axs[0, 1].set_aspect('equal')
    
    axs[1, 0].plot(Maxis, Paxis, linewidth=2.0)
    axs[1, 0].plot(M, P, '.', color="black", markersize=7)
    axs[1, 0].annotate(f"P={round(P/1000, 2)}kN \nM={round(M/1000000, 2)}kN.m \nratio={round(ratio, 2)}",
                 (M, P), 
                 textcoords="offset points", xytext=(5,0), ha="left")
    axs[1, 0].set_title(f"P-M chart for {round(alpha, 2)} degree")
    axs[1, 0].set_xlabel("M (N.mm)")
    axs[1, 0].set_ylabel("P (N)")
    axs[1, 0].axhline(y=0, color='b', linestyle='-')
    axs[1, 0].axvline(x=0, color='b', linestyle='-')
    axs[1, 0].grid(True)
    
    axs[1, 1].plot(MxAxis, MyAxis, linewidth=2.0)
    axs[1, 1].plot(Mx, My, '.', color="black", markersize=7)
    axs[1, 1].annotate(f"Mx={round(Mx/1000000, 2)}kN.m \nMy={round(My/1000000, 2)}kN.m", (Mx, My), 
                 textcoords="offset points", xytext=(5,0), ha="left")
    axs[1, 1].set_title(f"Mx-My chart on P={round(P/1000)} kN")
    axs[1, 1].set_xlabel("Mx (N.mm)")
    axs[1, 1].set_ylabel("My (N.mm)")
    axs[1, 1].axhline(y=0, color='b', linestyle='-')
    axs[1, 1].axvline(x=0, color='b', linestyle='-')
    axs[1, 1].axis("equal")
    axs[1, 1].grid(True)
    
    plt.show()
