import numpy as np


PHI_MOMENT_AXIAL_MIN = 0.65
PHI_MOMENT_AXIAL_MAX = 0.9
PHI_SHEAR = 0.75
PHI_TORSION = 0.75
PHI_BEARING = 0.65

ECU = 0.003

def ETY(fy: np.float32, Es: np.float32) -> np.float32:
    return fy / Es


def BETA1(fc: np.float32) -> np.float32:
    return np.float32(0.85 if fc<=28 else max(0.85-(0.05*(fc-28)/7), 0.65))
    

def PHI_MOMENRT_AXIAL(fy: np.float32, Es: np.float32, esMin: np.float32) -> np.float32:
    ety = fy / Es
    et = abs(esMin)
    return np.float32(PHI_MOMENT_AXIAL_MIN if et <= ety else PHI_MOMENT_AXIAL_MIN+0.25*((et-ety)/0.003)\
          if ety < et < ety+0.003 else PHI_MOMENT_AXIAL_MAX)


