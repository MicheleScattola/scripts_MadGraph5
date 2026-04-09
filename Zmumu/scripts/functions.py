import math

class restrictions:

    def __init__(self,pt_min,eta_max,iso_max):
        self.pt_min = pt_min
        self.eta_max = eta_max
        self.iso_max = iso_max


    def isSel(self,mu):
        if (
            mu.PT > self.pt_min
            and abs(mu.Eta) < self.eta_max
            and abs(mu.IsolationVar) < self.iso_max
        ):
            return True
        return False
    
    def Zsel(self,muons):
        # check if both mu are tight (self selection)
        if (
            len(muons)==2
            and self.isSel(muons[0])
            and self.isSel(muons[1])
            and muons[0].Charge * muons[1].Charge < 0
            and max([muons[0].PT,muons[1].PT]) > 25.0
            and min([muons[0].PT,muons[1].PT]) > 20.0
        ):
            return True
        return False

# DELTA R
def delta_r_eta_phi(eta1, phi1, eta2, phi2):
    dphi = math.atan2(math.sin(phi1 - phi2), math.cos(phi1 - phi2))
    deta = eta1 - eta2
    return math.sqrt(deta * deta + dphi * dphi)

def delta_r(obj1, obj2):
    return delta_r_eta_phi(obj1.Eta, obj1.Phi, obj2.Eta, obj2.Phi)

# Drawig

