from math import sqrt

import numpy as np
from pycivil.ACI318_19.designProps import DesignData
import pycivil.ACI318_19.assumptions as assump
from pycivil.sections.rebarSections import Rebar, GRebars


def calc_d(data: DesignData) -> float:
    _, _, _, maxy = data.section.bounds
    return max(*[maxy - point.y for point in data.Coords]) if len(data.Coords)!=0 else 0


def calc_Av_S_min(data: DesignData) -> float: 
    return max(0.062 * np.sqrt(data.fc) * (data.bw/data.fyt),
               0.35 * (data.bw/data.fyt))


def calc_Vc_max(data: DesignData, *, d: float|None = None) -> float:
    _d = calc_d(data) if d == None else d 
    return float(0.42*sqrt(data.fc) * data.bw * _d)


def calc_Vc(data: DesignData, *, Nu: float|None=None, rho: float|None=None, d: float|None = None) -> float:
    _d = calc_d(data) if d == None else d 
    _Nu = 0 if Nu==None else Nu
    Vc: float = 0
    if rho == None:
        Vc = float((0.17*sqrt(data.fc) + min(_Nu/(6*data.section.area), float(0.05*data.fc))) * data.bw * _d)
    else:
        Vc = float((0.66*pow(rho, (1/3))*sqrt(data.fc) + min(_Nu/(6*data.section.area), float(0.05*data.fc))) * data.bw * _d)
    Vc = max(Vc, 0)
    Vc = min(Vc, calc_Vc_max(data, d=d))
    return Vc


def calc_Vs(data: DesignData, Vu: float, *, Nu: float|None=None, rho: float|None=None, d: float|None = None) -> float:
    Vc = calc_Vc(data, Nu=Nu, rho=rho, d=d)
    return 0 if Vu <= assump.phis*Vc else float((Vu/assump.phis) - Vc)


def calc_Av_S(data: DesignData, Vu: float, *, Nu: float|None=None, rho: float|None=None, d: float|None = None) -> float:
    _d = calc_d(data) if d == None else d 
    Vc = calc_Vc(data, Nu=Nu, rho=rho, d=d)
    return 0 if Vu <= 0.5*assump.phis*Vc else \
        calc_Av_S_min(data) if assump.phis*Vc>Vu>0.5*assump.phis*Vc else \
            max(float(calc_Vs(data, Vu, Nu=Nu, rho=rho, d=d)/(data.fyt * _d)), calc_Av_S_min(data))


def calc_conf_rebars_dist(data: DesignData, Vu: float, rebar:Rebar|GRebars, *, 
                          Nu: float|None=None, rho: float|None=None, d: float|None = None) -> float:
    Vc = calc_Vc(data, Nu=Nu, rho=rho, d=d)
    return 0 if Vu<=0.5*assump.phis*Vc else round(rebar.Area/calc_Av_S(data, Vu, Nu=Nu, rho=rho, d=d))


def calc_Acp(data: DesignData) -> float:
    return data.section.area


def calc_Pcp(data: DesignData) -> float:
    return data.section.length


def calc_Aoh(data: DesignData) -> float:
    return data.section.buffer(-1*data.cover).area


def calc_Ph(data: DesignData) -> float:
    return data.section.buffer(-1*data.cover).length


def calc_Tth(data: DesignData) -> float:
    return 0.083 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def calc_Tcr(data: DesignData) -> float:
    return 0.33 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def Calc_At_S(data: DesignData, Tu: float, *, reDistribution=False, Aoh: float|None=None):
    _Aoh = calc_Aoh(data) if Aoh==None else Aoh
    if Tu < assump.phit*calc_Tth(data):
        At_S = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data))
        At_S=round((Tu/assump.phit)/(2 * 0.85 * _Aoh * data.fyt),2)
    return At_S


def Calc_Al(data: DesignData, Tu: float, *, reDistribution=False, Aoh: float|None=None, Ph: float|None=None):
    _Aoh = calc_Aoh(data) if Aoh==None else Aoh
    _Ph = calc_Ph(data) if Ph==None else Ph
    if Tu < assump.phit*calc_Tth(data):
        Al = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data))
        Al=round((Tu/assump.phit)*_Ph/(2*0.85*_Aoh*data.fy))
    return Al
