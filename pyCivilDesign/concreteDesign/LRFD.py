
from dataclasses import dataclass, field
from matplotlib import pyplot as plt

import numpy as np
from numpy.typing import NDArray
import pyCivilDesign.concreteDesign.LRFDmethod.PMManalysis as PMMsolver
from pyCivilDesign.concreteDesign.designProps import DesignData, setDesignDataFromSection
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


def PMM_analyze(section: ConcreteSct, P: float, Mx: float, My: float):
    data = DesignData.fromSection(section)
    _ratio = PMMsolver.calc_PM_ratio(data, P, Mx, My)
    _M = pow(Mx**2 + My**2, 0.5)
    _angle = float(PMMsolver.calc_angle_from_forces(data, P, Mx, My))
    _alpha = PMMsolver.calc_alpha(Mx, My)
    _c = float(PMMsolver.calc_c(data, P, _angle)) # type: ignore
    _es = PMMsolver.calc_es(data.section, data.Coords, _c, _angle)
    _fs = PMMsolver.calc_fs(data, _c, _angle)
    _Fs = PMMsolver.calc_Fs(data, _c, _angle)
    _Cc = PMMsolver.calc_Cc(data, _c, _angle)
    _percent = PMMsolver.As_percent(data)
    return PMMresults(_c, _angle, _alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "", _ratio) # type: ignore
        

def PMM_design(section: ConcreteSct, P: float, Mx: float, My: float):
    data = DesignData.fromSection(section)
    _percent = PMMsolver.calc_As_percent(data, P, Mx, My)
    data = PMMsolver.set_As_percent(data, _percent) # type: ignore
    angle = float(PMMsolver.calc_angle_from_forces(data,  P, Mx, My))
    c = float(PMMsolver.calc_c(data, P, angle))
    _M, _Mx, _My = PMMsolver.calc_M(data, c, angle)
    alpha = PMMsolver.calc_alpha(data, _Mx, _My) # type: ignore
    _es = PMMsolver.calc_es(data.section, data.Coords, c, angle)
    _fs = PMMsolver.calc_fs(data, c, angle)
    _Fs = PMMsolver.calc_Fs(data, c, angle)
    _Cc = PMMsolver.calc_Cc(data, c, angle)
    return PMMresults(c, angle, alpha, _es, _fs, _Fs, _Cc, P, _M, Mx, My, _percent, "", 0) # type: ignore


def show_PMM_analysis_result(section: ConcreteSct, P: float, Mx: float, My: float):
    data = setDesignDataFromSection(section)
    print(data.section)
    M = pow(Mx**2 + My**2, 0.5)
    ratio = PMMsolver.calc_PM_ratio(data, P, Mx, My)

    angle = PMMsolver.calc_angle_from_forces(data, P, Mx, My) # type: ignore
    alpha = PMMsolver.calc_alpha(Mx, My)
    PointNums = 20
    Paxis = np.linspace(PMMsolver.calc_Pt_max(data), PMMsolver.calc_P0(data), PointNums, endpoint=True)
    Maxis = np.array([PMMsolver.calc_Mn(data, p, angle)[0] for p in Paxis]) # type: ignore

    AlphaNums = 21
    Alphas = np.linspace(0, 360, AlphaNums)
    _M =  np.array([PMMsolver.calc_Mn(data, P, alpha) for alpha in Alphas]) # type: ignore
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
