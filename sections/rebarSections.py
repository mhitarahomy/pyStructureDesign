from math import pi, ceil
from dataclasses import dataclass
from typing import Callable


@dataclass
class Rebar():
    d: float
    Area = property(lambda self: CalcArea(self.d, 1))

    def __repr__(self) -> str:
        return f"\u03C6{self.d}"
    
    def __mul__(self, num:int):
        return GRebars({f"{self.d}": num})
    
    def __rmul__(self, num:int):
        return GRebars({f"{self.d}": num})
    
    def __add__(self, rebar):
        return GRebars({f"{self.d}": 2}) if rebar.d==self.d else GRebars({f"{self.d}": 1, f"{rebar.d}": 1})
        

@dataclass
class GRebars():
    listOfRebars: dict[str, float]
    Area = property(lambda self: CalcRebarsArea(**self.listOfRebars))

    def __repr__(self) -> str:
        rebarStr = ""
        for k, v in list(self.listOfRebars.items()):
            rebarStr += f"{v}\u03C6{k}+"
        return rebarStr[:-1]
    
    def __add__(self, rebar):
        output = GRebars(self.listOfRebars)
        if type(rebar) == Rebar:
            if str(rebar.d) in self.listOfRebars: 
                output.listOfRebars[f"{rebar.d}"] += 1
            else:
                output.listOfRebars.update({f"{rebar.d}": 1})
        elif type(rebar) == GRebars:
            for k, v in list(rebar.listOfRebars.items()):
                if k in output.listOfRebars:
                    output.listOfRebars[k] += v
                else:
                    output.listOfRebars.update({k: v})
        output.listOfRebars = {k: output.listOfRebars[k] for k in sorted(output.listOfRebars)}
        return output

    
CalcArea: Callable[[float, int], float] = lambda d, num=1: num * (pi*d**2)/4

def CalcRebarsArea(**listOfRebars) -> float:
    area: float = 0
    for k,v in list(listOfRebars.items()):
        area += CalcArea(float(k),v)
    return round(area,2)

def CalcNumOfBars(area:float, bar:float) -> int:
    return ceil(area/bar)

d8 =  Rebar(8)
d10 = Rebar(10)
d12 = Rebar(12)
d14 = Rebar(14)
d16 = Rebar(16)
d18 = Rebar(18)
d20 = Rebar(20)
d22 = Rebar(22)
d25 = Rebar(25)
d28 = Rebar(28)
d32 = Rebar(32)
d36 = Rebar(36)
