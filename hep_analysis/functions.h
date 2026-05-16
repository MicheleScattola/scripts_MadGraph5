#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <tuple>

using ROOT::RVec;
using ROOT::Math::PtEtaPhiMVector; 

const float Z_MASS = 91.118;
const float H_MASS = 125.0;
const float E_MASS = 511.0e-6;
const float MU_MASS = 106.0e-3;

// Forward declaration
struct Particle; 
float getMass(const Particle &p);

// struct for particle and jet selection
struct Particle {

  float pt;
  float eta;
  float phi;
  float d0;
  float sigma_d0;
  int charge;
  int type;

  Particle() : pt(0), eta(0), phi(0), d0(0), sigma_d0(0), charge(0), type(0) {}

  Particle(float pt_, float eta_, float phi_, float d0_, float sigma_d0_, int charge_, int type_)
    : pt(pt_/1000.), eta(eta_), phi(phi_), d0(d0_), sigma_d0(sigma_d0_), charge(charge_), type(type_) {}

  PtEtaPhiMVector p4() const {
    return PtEtaPhiMVector(pt, eta, phi, 0.0);
  }

  bool isGoodMu(float pt_min=5.0, float eta_max=2.7, float d0_err=3.) const {
    return (std::abs(type)==13 && pt > pt_min && std::abs(eta) < eta_max && std::abs(d0)/sigma_d0 < d0_err);
  }

  bool isGoodEl(float pt_min=7.0, float eta_max=2.47, float d0_err=5.) const {
    return (std::abs(type)==11 && pt > pt_min && std::abs(eta) < eta_max && std::abs(d0)/sigma_d0 < d0_err);
  }

  bool isGoodLep(float el_pt_min=7.0, float mu_pt_min=5.0, float el_eta_max=2.47, float mu_eta_max=2.7,float el_d0_err=5., float mu_d0_err=3.) const {
    return (isGoodEl(el_pt_min, el_eta_max, el_d0_err) || isGoodMu(mu_pt_min, mu_eta_max, mu_d0_err));
  }

  bool isGoodJet(float pt_min=30.0, float eta_max=4.5) const {
    return (pt > pt_min && std::abs(eta) < eta_max);
  }

  PtEtaPhiMVector constructP4() const {
    return PtEtaPhiMVector(pt, eta, phi, getMass(*this));
  } 
};

// Fixed anti-particle mass bug
float getMass(const Particle &p){
  if(std::abs(p.type) == 11) return E_MASS;
  if(std::abs(p.type) == 13) return MU_MASS;
  return 0.0;
}

// Pass by const reference to avoid copying
RVec<Particle> selectGoodLeptons(
  const RVec<Particle>& all_leptons,
  float el_pt_min=7.0,
  float mu_pt_min=5.0,
  float el_eta_max=2.47,
  float mu_eta_max=2.7,
  float el_d0_err=5.0,
  float mu_d0_err=3.0)
{
  RVec<Particle> out;
  out.reserve(all_leptons.size());
  // Range-based for loop is cleaner
  for(const auto& lep : all_leptons){
    if(lep.isGoodLep(el_pt_min, mu_pt_min, el_eta_max, mu_eta_max, el_d0_err, mu_d0_err)){
      out.push_back(lep);
    }
  }
  return out;
}

// Pass by const reference
bool isGoodQuadruplet(const RVec<Particle>& lep){
  if(lep.size() != 4) return false;
  int sum = 0;
  for(const auto& l : lep){
    sum += l.charge * l.type;
  }
  return sum == 0;
}

RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> leptonPairing(const RVec<Particle>& lep){
  RVec<Particle> sorted = lep;
  std::sort(sorted.begin(), sorted.end(), [](const Particle& a, const Particle& b){
      return a.pt > b.pt;
  });

  if(sorted[0].pt < 20.0 || sorted[1].pt < 15.0 || sorted[2].pt < 10.0){
    return RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>>();
  } 

  std::tuple<PtEtaPhiMVector,PtEtaPhiMVector> leading;
  std::tuple<PtEtaPhiMVector,PtEtaPhiMVector> sub_leading;
  float best_delta_mass = 9999.0;
  bool isBestFound = false;
  float m_leading = 0.;
  float m_subleading = 0.;

  for(int i=0; i<sorted.size()-3; i++){
    for(int j=i+1; j<sorted.size()-2; j++){
      for(int k=j+1; k<sorted.size()-1; k++){
        for(int l=k+1; l<sorted.size(); l++){

          if(sorted[i].type*sorted[i].charge + sorted[j].type*sorted[j].charge +
             sorted[k].type*sorted[k].charge + sorted[l].type*sorted[l].charge != 0) continue;
             
          Particle quad[4] = {sorted[i], sorted[j], sorted[k], sorted[l]};
          int pairings[3][4] = {{0,1,2,3}, {0,2,1,3}, {0,3,1,2}};
          
          // Fixed out-of-bounds loop (s < 3)
          for(int s=0; s<3; s++){
            Particle lep1 = quad[pairings[s][0]];
            Particle lep2 = quad[pairings[s][1]];
            Particle lep3 = quad[pairings[s][2]];
            Particle lep4 = quad[pairings[s][3]];
            
            if (std::abs(lep1.type) != std::abs(lep2.type) || lep1.charge * lep2.charge > 0) continue;
            if (std::abs(lep3.type) != std::abs(lep4.type) || lep3.charge * lep4.charge > 0) continue;
            
            auto p1 = lep1.constructP4();
            auto p2 = lep2.constructP4();
            auto p3 = lep3.constructP4();
            auto p4 = lep4.constructP4(); 

            auto p_lead = p1 + p2;
            auto p_sublead = p3 + p4;
            m_leading = p_lead.M();
            m_subleading = p_sublead.M();
            
            float delta_lead = std::abs(m_leading - Z_MASS);
            float delta_sublead = std::abs(m_subleading - Z_MASS);
            
            // if inverted swap
            if(delta_sublead < delta_lead){
              std::swap(lep1, lep3);
              std::swap(lep2, lep4);
              std::swap(delta_lead, delta_sublead);
              std::swap(p_lead, p_sublead);
              std::swap(m_leading, m_subleading);
              std::swap(p1, p3);
              std::swap(p2, p4);
            }
            /*
            // mass requirements
            if(m_leading < 50. || m_leading > 106.) continue;
            if(m_subleading < 12. || m_subleading > 115.) continue;
            float tot_mass = (p_lead + p_sublead).M();
            if(tot_mass < 105. || tot_mass > 160.) continue;
            
            
            // lepton separation
            if(ROOT::Math::VectorUtil::DeltaR(p1, p2) < 0.1) continue;
            if(ROOT::Math::VectorUtil::DeltaR(p3, p4) < 0.1) continue;*/
            
            // check if it's the best pair found
            if(delta_lead < best_delta_mass){
              best_delta_mass = delta_lead;
              leading = std::make_tuple(p1,p2);
              sub_leading = std::make_tuple(p3,p4);
              isBestFound = true;
            }
          }
        }
      }
    }
  }

  RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> out;
  if(isBestFound){
    //std::cout << "Best found\n";
    out.push_back(leading);
    out.push_back(sub_leading);
  }
  
  return out;
}

float fourLeptonMass(const RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> &tuple){
  if(tuple.size()!=2){
    throw std::runtime_error("tuple is not of size 2: leading & sub-leading pair");
    return 0.;
  }
  PtEtaPhiMVector p1 = std::get<0>(tuple[0]);
  PtEtaPhiMVector p2 = std::get<1>(tuple[0]);
  PtEtaPhiMVector p3 = std::get<0>(tuple[1]);
  PtEtaPhiMVector p4 = std::get<1>(tuple[1]);

  float tot_mass = (p1+p2+p3+p4).M();
  return tot_mass;
}
float diLeptonMass(const RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> &tuple, const bool &leading){
  if(tuple.size()!=2){
    throw std::runtime_error("tuple is not of size 2: leading & sub-leading pair");
    return 0.;
  }
  PtEtaPhiMVector p1 = std::get<0>(tuple[0]);
  PtEtaPhiMVector p2 = std::get<1>(tuple[0]);
  PtEtaPhiMVector p3 = std::get<0>(tuple[1]);
  PtEtaPhiMVector p4 = std::get<1>(tuple[1]);

  if(leading){
    return (p1+p2).M();
  } else {
    return (p3+p4).M();
  }
}

float diLeptonDeltaR(const RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> &tuple, const bool &leading){
  if(tuple.size()!=2){
    throw std::runtime_error("tuple is not of size 2: leading & sub-leading pair");
    return 0.;
  }
  PtEtaPhiMVector p1 = std::get<0>(tuple[0]);
  PtEtaPhiMVector p2 = std::get<1>(tuple[0]);
  PtEtaPhiMVector p3 = std::get<0>(tuple[1]);
  PtEtaPhiMVector p4 = std::get<1>(tuple[1]);

  if(leading){
    return ROOT::Math::VectorUtil::DeltaR(p1,p2);
  } else {
    return ROOT::Math::VectorUtil::DeltaR(p3,p4);
  }
}

RVec<float> diLeptonEta(const RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> &tuple, const bool &leading){
  
  RVec<float> out;
  if(tuple.size()!=2){
    throw std::runtime_error("tuple is not of size 2: leading & sub-leading pair");
    return out;
  }
  PtEtaPhiMVector p1 = std::get<0>(tuple[0]);
  PtEtaPhiMVector p2 = std::get<1>(tuple[0]);
  PtEtaPhiMVector p3 = std::get<0>(tuple[1]);
  PtEtaPhiMVector p4 = std::get<1>(tuple[1]);

  out.reserve(tuple.size());

  if(leading){
    out.push_back(p1.Eta());
    out.push_back(p2.Eta());
  } else {
    out.push_back(p3.Eta());
    out.push_back(p4.Eta());
  }
  return out;
}

RVec<float> diLeptonPhi(const RVec<std::tuple<PtEtaPhiMVector,PtEtaPhiMVector>> &tuple, const bool &leading){
  
  RVec<float> out;
  if(tuple.size()!=2){
    throw std::runtime_error("tuple is not of size 2: leading & sub-leading pair");
    return out;
  }
  PtEtaPhiMVector p1 = std::get<0>(tuple[0]);
  PtEtaPhiMVector p2 = std::get<1>(tuple[0]);
  PtEtaPhiMVector p3 = std::get<0>(tuple[1]);
  PtEtaPhiMVector p4 = std::get<1>(tuple[1]);

  out.reserve(tuple.size());

  if(leading){
    out.push_back(p1.Phi());
    out.push_back(p2.Phi());
  } else {
    out.push_back(p3.Phi());
    out.push_back(p4.Phi());
  }
  return out;
}

#endif