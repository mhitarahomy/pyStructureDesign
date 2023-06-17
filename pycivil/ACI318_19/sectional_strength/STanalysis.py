from math import sqrt

import numpy as np

from pycivil.ACI318_19.designProps import DesignData
import pycivil.ACI318_19.assumptions as assump
from pycivil.sections.rebarSections import Rebar, GRebars


def calc_d(data: DesignData) -> np.float32:
    _, miny, _, maxy = data.section.bounds
    return maxy-miny-data.clear_cover-data.max_conf_rebar_size-(data.max_rebar_size/2)


def calc_Av_S_min(data: DesignData) -> np.float32: 
    return max(0.062 * np.sqrt(data.fc) * (data.bw/data.fyt),
               0.35 * (data.bw/data.fyt))


def calc_Vc_max(data: DesignData, _d: np.float32|None = None) -> np.float32:
    d = calc_d(data) if _d == None else _d 
    return 0.42*sqrt(data.fc) * data.bw * d


def calc_Vc(data: DesignData, _Nu: np.float32|None=None, _rho: np.float32|None=None, _d: np.float32|None = None) -> np.float32:
    d = calc_d(data) if _d == None else _d 
    Nu = 0 if _Nu==None else _Nu
    Vc: float = 0
    if _rho == None:
        Vc = (0.17*sqrt(data.fc) + min(Nu/(6*data.section.area), 0.05*data.fc)) * data.bw * d # type: ignore
    else:
        Vc = (0.66*pow(_rho, (1/3))*sqrt(data.fc) + min(Nu/(6*data.section.area), 0.05*data.fc)) * data.bw * d # type: ignore
    Vc = max(Vc, 0)
    Vc = min(Vc, calc_Vc_max(data, _d=d))
    return np.float32(Vc)


def calc_Vs(data: DesignData, Vu: np.float32, _Nu: np.float32|None=None, _rho: np.float32|None=None, _d: np.float32|None = None) -> np.float32:
    Vc = calc_Vc(data, _Nu=_Nu, _rho=_rho, _d=_d)
    return np.float32(0) if Vu <= assump.PHI_SHEAR*Vc else (Vu/assump.PHI_SHEAR) - Vc


def calc_Av_S(data: DesignData, Vu: np.float32, _Nu: np.float32|None=None, _rho: np.float32|None=None, _d: np.float32|None = None) -> np.float32:
    d = calc_d(data) if _d == None else _d 
    Vc = calc_Vc(data, _Nu=_Nu, _rho=_rho, _d=_d)
    return np.float32(0) if Vu <= 0.5*assump.PHI_SHEAR*Vc else \
        calc_Av_S_min(data) if assump.PHI_SHEAR*Vc>Vu>0.5*assump.PHI_SHEAR*Vc else \
            max(calc_Vs(data, Vu, _Nu=_Nu, _rho=_rho, _d=_d)/(data.fyt * d), calc_Av_S_min(data)) # type: ignore


def calc_conf_rebars_dist(data: DesignData, Vu: np.float32, rebar:Rebar|GRebars, *, 
                          _Nu: np.float32|None=None, _rho: np.float32|None=None, _d: np.float32|None = None) -> np.float32:
    Vc = calc_Vc(data, _Nu=_Nu, _rho=_rho, _d=_d)
    return np.float32(0) if Vu<=0.5*assump.PHI_SHEAR*Vc else round(rebar.area/calc_Av_S(data, Vu, _Nu=_Nu, _rho=_rho, _d=_d))


def calc_Acp(data: DesignData) -> np.float32:
    return np.float32(data.section.area)


def calc_Pcp(data: DesignData) -> np.float32:
    return np.float32(data.section.length)


def calc_Aoh(data: DesignData) -> np.float32:
    return np.float32(data.section.buffer(-1*data.clear_cover).area)


def calc_Ph(data: DesignData) -> np.float32:
    return np.float32(data.section.buffer(-1*data.clear_cover).length)


def calc_Tth(data: DesignData) -> np.float32:
    return 0.083 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def calc_Tcr(data: DesignData) -> np.float32:
    return 0.33 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def Calc_At_S(data: DesignData, Tu: np.float32, reDistribution=False, Aoh: np.float32|None=None):
    _Aoh = calc_Aoh(data) if Aoh==None else Aoh
    if Tu < assump.PHI_TORSION*calc_Tth(data):
        At_S = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data)) # type: ignore
        At_S=round((Tu/assump.PHI_TORSION)/(2 * 0.85 * _Aoh * data.fyt),2)
    return At_S


def Calc_Al(data: DesignData, Tu: np.float32, reDistribution=False, Aoh: np.float32|None=None, Ph: np.float32|None=None):
    _Aoh = calc_Aoh(data) if Aoh==None else Aoh
    _Ph = calc_Ph(data) if Ph==None else Ph
    if Tu < assump.PHI_TORSION*calc_Tth(data):
        Al = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data)) # type: ignore
        Al=round((Tu/assump.PHI_TORSION)*_Ph/(2*0.85*_Aoh*data.fy))
    return Al
