import math
import ROOT

ELECTRON_MASS = 0.511e-3

class restrictions:

    def __init__(self,pt_min,eta_max,tr_low,tr_high,d0_res_max,z0_max):
        self.pt_min = pt_min
        self.eta_max = eta_max
        self.transition_low=tr_low
        self.transition_high=tr_high
        self.d0_res_max=d0_res_max
        self.z0_max=z0_max

    def isSel(self,el):
        if (el.PT > self.pt_min
            and abs(el.Eta) < self.eta_max
            and (abs(el.Eta) < self.transition_low or abs(el.Eta) > self.transition_high)
            #and abs(el.D0)/(el.ErrorD0) < self.d0_res_max
                ):
            # build P4 to get theta
            # compute other restrictions later; if they are not necessary they save computation time
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM(el.PT,el.Eta,el.Phi,ELECTRON_MASS)
            theta = p4.Theta()

            if (abs(el.DZ * math.sin(theta)) < self.z0_max):
                return True
        return False
    
    
    

## OTHER FUNCTIONS
# VARCONE
def delta_r_varcone(p,dr_max):
    return min(10.0/p.Et(),dr_max)

# isolation with varcone
def isIso(p,tracks,dr_max,iso_max):
    # look for particles in varcone
    p4 = ROOT.TLorentzVector()
    p4.SetPtEtaPhiM(p.PT,p.Eta,p.Phi,ELECTRON_MASS)
    dr = delta_r_varcone(p4,dr_max)
    sum_pt = 0.0
    for trk in tracks:
        # remove the electron's pt (0.0001 safety cone)
        distance = delta_r(trk,p)
        if distance < dr and distance > 0.0001:
            sum_pt += trk.PT
        
    iso = sum_pt/p.PT
    return (iso < iso_max)
  
# return index for electrons coming from Z > e e 
def check_z(particles):

    n = particles.GetEntries()
    for i in range(n):
        p = particles.At(i)
        id1 = p.D1
        id2 = p.D2
        if id1 != -1 and id2 != -1:
            d1 = particles.At(id1)
            d2 = particles.At(id2)

            if p.PID == 23 and abs(d1.PID) == 11 and abs(d2.PID) == 11:
                print(p.Status)
            
# return daughter ids given status and pdg            
def getDaughters(particles,status,pdg):

    n_entries = particles.GetEntries()
    for i in range(n_entries):
        p = particles.At(i)
        if p.Status == status and abs(p.PID) == pdg:
            return i, p.D1, p.D2
        else:
            continue

    return -1,-1,-1
            

# DELTA R
def delta_r_eta_phi(eta1, phi1, eta2, phi2):
    dphi = math.atan2(math.sin(phi1 - phi2), math.cos(phi1 - phi2))
    deta = eta1 - eta2
    return math.sqrt(deta * deta + dphi * dphi)

def delta_r(obj1, obj2):

    return delta_r_eta_phi(obj1.Eta, obj1.Phi, obj2.Eta, obj2.Phi)


