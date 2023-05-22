import numpy as np


phis = np.float32(0.9)
phia = np.float32(0.65)
phit = np.float32(0.6)
ecu = np.float32(0.003)


def beta1(fc: float):
    return np.float32(0.85 if fc<=28 else max(0.85-(0.05*(fc-28)/7), 0.65))
    

def phif(fy: float, Es: float, esMin: np.float32) -> np.float32:
    ety = fy / Es
    et = abs(esMin)
    return np.float32(0.65 if et <= ety else 0.65+0.25*((et-ety)/0.003)\
          if ety < et < ety+0.003 else 0.9)


