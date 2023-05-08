from dataclasses import dataclass
from typing import Callable

from pyCivilDesign.concreteBuilding.concreteAnalysis import DesignData


beta1: Callable[[DesignData], float] = lambda data: 0.85 if data.fc<=28 else max(0.85-(0.05*(data.fc-28)/7), 0.65)

def phif(data: DesignData, esMin: float) -> float:
    ety = data.fy / data.Es
    et = abs(esMin)
    return 0.65 if et <= ety else 0.65+0.25*((et-ety)/0.003) if ety < et < ety+0.003 else 0.9


@dataclass()
class Assumptions():
    beta1: Callable[[DesignData], float] = beta1
    phif: Callable[[DesignData, float], float] = phif
    phis: float = 0.9
    phia: float = 0.65
    phit: float = 0.6
    ecu = 0.003

defaultAssumption = Assumptions()
