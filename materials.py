from dataclasses import dataclass, field
from typing import Callable, Any


@dataclass(kw_only=True)
class Material():
    density: float


@dataclass(kw_only=True)
class StructuralMaterial(Material):
    E: float
    nu: float
    A: float


@dataclass(kw_only=True)
class ConcreteMat(StructuralMaterial):
    fc: float = field(default=25)
    calcEc: bool = field(default=True)
    density: float = field(default=2300)
    E: float = field(init=False)
    nu: float = field(init=False, default=0.2)
    A: float = field(init=False, default=10e-6)

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        if (__name=="fc" or __name=="density") and self.calcEc:
            self.E =  CalcEc(self.fc, self.density)
            
    wc = property(lambda self: self.density)
    Ec = property(lambda self: self.E)
    

@dataclass(kw_only=True)
class SteelMat(StructuralMaterial):
    fy: float
    fu: float
    alpha: float = field(default=1.2)
    density: float = field(init=False, default=7849.05)
    E: float = field(init=False, default=2000000)
    nu: float = field(init=False, default=0.3)
    A: float = field(init=False, default=12e-6)
    
    Es = property(lambda self: self.E)
    fye = property(lambda self: self.alpha * self.fy)
    fue = property(lambda self: self.alpha * self.fu)


@dataclass(kw_only=True)
class RebarMat(StructuralMaterial):
    fy: float
    fu: float
    E: float = field(init=False, default=2000000)
    density: float = field(init=False, default=7849.05)
    nu: float = field(init=False, default=0.3)
    A: float = field(init=False, default=12e-6)

    Es = property(lambda self: self.E)
    es = property(lambda self: self.fy/self.E)


CalcEc: Callable[[float, float], float] = lambda fc, wc=2300: round(0.043 * pow(wc, 1.5) * pow(fc, 0.5))

C10 = ConcreteMat(fc=10)
C12 = ConcreteMat(fc=12)
C16 = ConcreteMat(fc=16)
C20 = ConcreteMat(fc=20)
C25 = ConcreteMat(fc=25)
C30 = ConcreteMat(fc=30)
C35 = ConcreteMat(fc=35)
C40 = ConcreteMat(fc=40)
C45 = ConcreteMat(fc=45)
C50 = ConcreteMat(fc=50)
C55 = ConcreteMat(fc=55)
C60 = ConcreteMat(fc=60)
C65 = ConcreteMat(fc=65)
C70 = ConcreteMat(fc=70)

RolledST52 = SteelMat(fy=3600, fu=5200, alpha=1.2)
RolledST37 = SteelMat(fy=2400, fu=3700, alpha=1.2)
PlateST37 = SteelMat(fy=2400, fu=3700, alpha=1.15)
PlateST52 = SteelMat(fy=3600, fu=5200, alpha=1.15)

S240 = RebarMat(fy=240, fu=360)
S340 = RebarMat(fy=340, fu=500)
S350 = RebarMat(fy=350, fu=500)
S400 = RebarMat(fy=400, fu=600)
S420 = RebarMat(fy=420, fu=600)
S500 = RebarMat(fy=500, fu=650)
S500C = RebarMat(fy=500, fu=550)
S520 = RebarMat(fy=520, fu=690)
AI = S240
AII = S340
AIII = S400
AIV = S500
