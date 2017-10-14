#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "MyTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
//#include "TLorentzVector.h"
#include "TString.h"
#include "TVector2.h"

using namespace std;
using fastjet::PseudoJet;


// Constructor 
MyTools::MyTools(){
    m_test = 0;
}

double MyTools::JetPull(fastjet::PseudoJet jet,fastjet::PseudoJet jet2){

  double pi = 4.*atan(1.);
  TVector2* Pull = new TVector2(0.,0.);
  double PT = jet.perp();
  TVector2 J = TVector2(jet.rapidity(),jet.phi());

  for (unsigned int i=0; i<jet.constituents().size(); i++){
    TVector2 ci = TVector2(jet.constituents()[i].rapidity(),jet.constituents()[i].phi());
    TVector2 ri = ci-J;
    while (ri.Y() > pi){
      ri.Set(ri.X(),ri.Y()-2*pi);
    }
    while (ri.Y() < -pi){
      ri.Set(ri.X(),ri.Y()+2*pi);
    }
    Pull->Set(Pull->X()+ri.Mod()*ri.X()*jet.constituents()[i].perp()/PT,Pull->Y()+ri.Mod()*ri.Y()*jet.constituents()[i].perp()/PT);
  }
  
  TVector2 J2 = TVector2(jet2.rapidity() - jet.rapidity(),jet2.phi() - jet.phi());

  /*
  TVector2 J3 = TVector2(jet2.rapidity() - jet.rapidity(),jet2.phi() - jet.phi());

  while (J3.Y() > pi){
    J3.Set(J3.X(),J3.Y()-2*pi);
  }
  while (J3.Y() < -pi){
    J3.Set(J3.X(),J3.Y()+2*pi);
  }

  std::cout << Pull->DeltaPhi(J2) <<"  "  << Pull->DeltaPhi(J2) << std::endl;
  */

  return Pull->DeltaPhi(J2);

}

int MyTools::Match(fastjet::PseudoJet jet,vector<fastjet::PseudoJet> jets){

  double close=500;
  int found=-1;
  double closePt = 0.;
  for (unsigned int i=0; i<jets.size(); i++){
    double myR = jet.delta_R(jets[i]);
    if (jet.pt() > closePt && myR < 0.7){
      found=i;
      myR=close;
      closePt = jet.pt();
    }
  }

  return found;

}

double MyTools::JetCharge(fastjet::PseudoJet jet,double kappa){
  //Returns the jet charge with weighting factor kappa
  double charge=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
    charge+=jet.constituents()[i].user_info<MyToolsUserInfo>().charge()*pow(jet.constituents()[i].pt(),kappa);
  }
  return charge/pow(jet.pt(),kappa);
}

bool MyTools::IsBHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=500 && abs_pdgId_mod10k<600) /*mesons*/  ||
      (abs_pdgId>=5000      && abs_pdgId<6000)      /*baryons*/   )
    return true;

  return false;
}

bool MyTools::IsCHadron(int pdgId){
  int abs_pdgId = abs(pdgId);
  int abs_pdgId_mod10k = (abs(pdgId)%10000);
  if( (abs_pdgId_mod10k>=400 && abs_pdgId_mod10k<500) /*mesons*/  ||
      (abs_pdgId>=4000      && abs_pdgId<5000)      /*baryons*/   )
    return true;

  return false;
}

bool MyTools::Btag(fastjet::PseudoJet jet,vector<fastjet::PseudoJet> bhadrons,vector<fastjet::PseudoJet> chadrons,double jetrad,double b,double c,double uds){

  TRandom3 *rand = new TRandom3(0);

  int foundb=0;
  int foundc=0;
  
  for (unsigned int i=0; i<bhadrons.size(); i++){
    if (bhadrons[i].delta_R(jet)<jetrad){
      foundb=1;
    }
  }

  for (unsigned int i=0; i<chadrons.size(); i++){
    if (chadrons[i].delta_R(jet)<jetrad){
      foundc=1;
    }
  }

  if (foundb==1){
    double flip = rand->Uniform(0.,1.);
    if (flip < b){
      delete rand;
      return true;
    }
  }
  if (foundc==1){
    double flip= rand->Uniform(0.,1.);
    if (flip < 1./c){
      delete rand;
      return true;
    }
  }
  double flip= rand->Uniform(0.,1.);
  if (flip < 1./uds){
    delete rand;
    return true;
  }
  
  delete rand;
  return false;
}

bool MyTools::BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, double jetrad, int BosonID ){

  for (unsigned int i=0; i<Bosons.size(); i++){
      if (Bosons[i].user_info<MyToolsUserInfo>().pdg_id() != BosonID) continue;
      if (Bosons[i].delta_R(jet)<jetrad){
        return true;
      }
  }
  return false;
}

bool MyTools::IsIsolated(Pythia8::Particle* particle, Pythia8::Pythia* pythia8, float rel_iso, float ConeSize){
    float sumpT=0;
    fastjet::PseudoJet part(particle->px(), particle->py(), particle->pz(),particle->e() );
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){
        if (!pythia8->event[ip].isFinal() )      continue;
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;
        if (pythia8->event[ip].pT()       < 0.5) continue;
        if (&pythia8->event[ip] == particle    ) continue; //same particle
        fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
        if(p.delta_R(part)>ConeSize)             continue;
        sumpT+=p.pt();
    }
    if(sumpT/part.pt()>rel_iso) return false;
    else return true;
}
