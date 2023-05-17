from math import sqrt

import numpy as np
from pyCivilDesign.concreteDesign.designProps import DesignData
# TODO: Shear & Torsion functions must be here


def Av_S_Min(data: DesignData): 
    return max(0.062*np.sqrt(data.fc)*(data.bw/data.fyt), 0.35*(data.bw/data.fyt))


def VcMax(data: DesignData): 
    return 0.42*sqrt(data.fc) * data.bw * data.d