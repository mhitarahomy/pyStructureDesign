from math import sqrt

import numpy as np
from pyCivilDesign.concreteDesign.designProps import DesignData
import pyCivilDesign.concreteDesign.LRFDmethod.assumptions as assump
from pyCivilDesign.sections.rebarSections import Rebar, GRebars


def calc_Av_S_min(data: DesignData): 
    return max(0.062*np.sqrt(data.fc)*(data.bw/data.fyt), 0.35*(data.bw/data.fyt))


def calc_Vc_max(data: DesignData): 
    return 0.42*sqrt(data.fc) * data.bw * data.d


def calc_Vc(data: DesignData, Nu=0, rho=0) -> float:
    Vc1 = (0.17*sqrt(data.fc) + min(Nu/(6*data.section.area), 0.05*data.fc)) * data.bw * data.d # type: ignore
    Vc2 = (0.66*pow(rho, (1/3))*sqrt(data.fc) + min(Nu/(6*data.section.area), 0.05*data.fc)) * data.bw * data.d # type: ignore
    Vc = Vc1 if rho == 0 else Vc2
    if Vc < 0:
        Vc = 0
    Vc = min(Vc, calc_Vc_max(data)) # type: ignore
    return Vc


def calc_Vs(data: DesignData, Vu: float, Nu=0, rho=0):
    Vc = calc_Vc(data, Nu, rho)
    return 0 if Vu <= assump.phis*Vc else (Vu/assump.phis) - Vc


def calc_Av_S(data: DesignData, Vu: float, Nu=0, rho=0):
    Vc = calc_Vc(data, Nu, rho)
    return 0 if Vu <= 0.5*assump.phis*Vc else \
        calc_Av_S_min(data) if assump.phis*Vc>Vu>0.5*assump.phis*Vc else \
            round(max(calc_Vs(data, Vu, Nu, rho)/(data.fyt * data.d), calc_Av_S_min(data)),2)


def calc_conf_rebars_dist(data: DesignData, Vu: float, rebar:Rebar|GRebars, Nu=0, rho=0):
    Vc = calc_Vc(data, Nu, rho)
    return 0 if Vu<=0.5*assump.phis*Vc else round(rebar.Area/calc_Av_S(data, Vu, Nu, rho))


def calc_Acp(data: DesignData) -> float:
    return data.section.area


def calc_Pcp(data: DesignData) -> float:
    return data.section.length


def calc_Aoh(data: DesignData) -> float:
    pass


def calc_Ph(data: DesignData) -> float:
    pass


def calc_Tth(data: DesignData) -> float:
    return 0.083 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def calc_Tcr(data: DesignData) -> float:
    return 0.33 * sqrt(data.fc) * calc_Acp(data)**2 / calc_Pcp(data)


def Calc_At_S(data: DesignData, Tu: float, reDistribution=False):
    if Tu < assump.phit*calc_Tth(data):
        At_S = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data))
        At_S=round((Tu/assump.phit)/(2 * 0.85 * calc_Aoh(data) * data.fyt),2)
    return At_S


def Calc_Al(data: DesignData, Tu: float, reDistribution=False):
    if Tu < assump.phit*calc_Tth(data):
        Al = 0
    else:
        if reDistribution:
            Tu=min(Tu, calc_Tcr(data))
        Al=round((Tu/assump.phit)*calc_Ph(data)/(2*0.85*calc_Aoh(data)*data.fy))
    return Al
