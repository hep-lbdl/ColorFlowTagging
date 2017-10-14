#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TGenPhaseSpace.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
//#include "TRoot.h"

TString addspace(double x){
  if (x<0){
    return "";
  }
  else
    return " ";
}

TString numb_convert(float val){

  if (val==0){
    return " 0.00000E+00";
  }

  TString addsign = addspace(val);
  int mysign = 1;
  if (val < 0) mysign = -1;
  val = fabs(val);
  double base_10 = log(val)/log(10);

  int expo = abs((int)TMath::Floor(base_10+0.5));
  double rest = mysign*val/pow(10,TMath::Floor(base_10+0.5));

  if (expo < 10){
    if (base_10 > 0) return addsign+TString::Format("%1.7f",rest)+TString::Format("E+0%i",expo);
    else return addsign+TString::Format("%1.7f",rest)+TString::Format("E-0%i",expo);
  }
  else{
    if (base_10 > 0) return addsign+TString::Format("%1.7f",rest)+TString::Format("E+%i",expo);
    else return addsign+TString::Format("%1.7f",rest)+TString::Format("E-%i",expo);
  }

}

void mymain(){

  int nevents = 5000000;
  double MotherHiggsMass = 1500; //GeV
  double DaughterHiggsMass = 125; //GeV
  TString ColorRep = "Singlet"; //Other options: "Triplet, Octet"
  
  ofstream myfile;
  myfile.open("eventsSinglet.lhe");
  
  ofstream myfile2;
  myfile2.open("eventsOctet.lhe");

  //First, write the header for the LHE file.
  myfile << "<LesHouchesEvents version=\"1.0\">" << std::endl;
  myfile << "<header>" << std::endl;

  myfile << "<!--" << std::endl;
  myfile << "#*********************************************************************" << std::endl;
  myfile << "#                                                                    *" << std::endl;
  myfile << "#                               HiggsGun                             *" << std::endl;
  myfile << "#                                                                    *" << std::endl;  
  myfile << "#    This is an event generator that simulates H'-> HA, H->bb,       *" << std::endl;
  myfile << "#    A->invisible.  The color representation of H is configurable.   *" << std::endl;
  myfile << "#    The output is an LHE file that can be showered with e.g.        *" << std::endl;
  myfile << "#    Pythia/Herwig.                                                  *" << std::endl;
  myfile << "#                                                                    *" << std::endl;
  myfile << "#           Author: Benjamin Nachman, bnachman@cern.ch               *" << std::endl;
  myfile << "#           Last Revision: Jan. 31 2016                              *" << std::endl;
  myfile << "#                                                                    *" << std::endl;
  myfile << "#*********************************************************************" << std::endl;
  myfile << "--> " << std::endl;

  myfile2 << "<LesHouchesEvents version=\"1.0\">" << std::endl;
  myfile2 << "<header>" << std::endl;

  myfile2 << "<!--" << std::endl;
  myfile2 << "#*********************************************************************" << std::endl;
  myfile2 << "#                                                                    *" << std::endl;
  myfile2 << "#                               HiggsGun                             *" << std::endl;
  myfile2 << "#                                                                    *" << std::endl;
  myfile2 << "#    This is an event generator that simulates H'-> HA, H->bb,       *" << std::endl;
  myfile2 << "#    A->invisible.  The color representation of H is configurable.   *" << std::endl;
  myfile2 << "#    The output is an LHE file that can be showered with e.g.        *" << std::endl;
  myfile2 << "#    Pythia/Herwig.                                                  *" << std::endl;
  myfile2 << "#                                                                    *" << std::endl;
  myfile2 << "#           Author: Benjamin Nachman, bnachman@cern.ch               *" << std::endl;
  myfile2 << "#           Last Revision: Jan. 31 2016                              *" << std::endl;
  myfile2 << "#                                                                    *" << std::endl;
  myfile2 << "#*********************************************************************" << std::endl;
  myfile2 << "--> " << std::endl;

  myfile << "<GenerationInfo>" << std::endl;
  myfile << "#  Number of Events        : " << nevents << std::endl;
  myfile << "#  Color Representation  :  " << "Singlet" << std::endl;
  myfile << "</GenerationInfo>" << std::endl;
  myfile << "</header>" << std::endl;
  myfile << "<init>" << std::endl;
  myfile << "2212     2212  0.65000000000E+04  0.65000000000E+04 0 0   10042   10042 2   1" << std::endl;
  myfile << "0.32537690924E-01  0.55367175698E-04  0.43599910000E-05   0" << std::endl;
  myfile << "</init>" << std::endl;

  myfile2 << "<GenerationInfo>" << std::endl;
  myfile2 << "#  Number of Events        : " << nevents << std::endl;
  myfile2 << "#  Color Representation  :  " << "Octet" << std::endl;
  myfile2 << "</GenerationInfo>" << std::endl;
  myfile2 << "</header>" << std::endl;
  myfile2 << "<init>" << std::endl;
  myfile2 << "2212     2212  0.65000000000E+04  0.65000000000E+04 0 0   10042   10042 2   1" << std::endl;
  myfile2 << "0.32537690924E-01  0.55367175698E-04  0.43599910000E-05   0" << std::endl;
  myfile2 << "</init>" << std::endl;
  
  TLorentzVector MotherHiggs = TLorentzVector(0,0,0,0);
  MotherHiggs.SetPtEtaPhiM(0.,0.,0.,MotherHiggsMass);
  vector<double> masses;
  masses.push_back(DaughterHiggsMass);
  masses.push_back(15); //A-mass, doesn't really matter as it decays invisibally.
  TGenPhaseSpace event;
  event.SetDecay(MotherHiggs, masses.size(), &masses[0]); 

  vector<double> daughtermasses;
  daughtermasses.push_back(5);
  daughtermasses.push_back(5);

  int accepted = 0;
  for (int i=0; i<nevents; i++){

    event.Generate();
    TLorentzVector *Higgs=event.GetDecay(0); 
    TLorentzVector *A=event.GetDecay(1); 

    if (Higgs->Pt() < 400 || Higgs->Pt() > 450){
      continue;
    }
    accepted++;

    myfile << "<event>" << std::endl;
    myfile2 << "<event>" << std::endl;

    TGenPhaseSpace HiggsDecay;
    HiggsDecay.SetDecay(*Higgs, daughtermasses.size(), &daughtermasses[0]); 
    HiggsDecay.Generate();

    TLorentzVector *b1=HiggsDecay.GetDecay(0); 
    TLorentzVector *b2=HiggsDecay.GetDecay(1); 

    //std::cout << b1->DeltaR(*b2) << std::endl;

    TGenPhaseSpace ADecay;
    ADecay.SetDecay(*A, daughtermasses.size(), &daughtermasses[0]); 
    ADecay.Generate();

    TLorentzVector *nu1=ADecay.GetDecay(0); 
    TLorentzVector *nu2=ADecay.GetDecay(1); 

    double pZ = MotherHiggs.M()/2.;

    TLorentzVector b1_truncated = TLorentzVector(atof(numb_convert(b1->Px())),atof(numb_convert(b1->Py())),atof(numb_convert(b1->Pz())),atof(numb_convert(b1->E())));
    TLorentzVector b2_truncated = TLorentzVector(atof(numb_convert(b2->Px())),atof(numb_convert(b2->Py())),atof(numb_convert(b2->Pz())),atof(numb_convert(b2->E())));
    TLorentzVector nu1_truncated = TLorentzVector(atof(numb_convert(nu1->Px())),atof(numb_convert(nu1->Py())),atof(numb_convert(nu1->Pz())),atof(numb_convert(nu1->E())));
    TLorentzVector nu2_truncated = TLorentzVector(atof(numb_convert(nu2->Px())),atof(numb_convert(nu2->Py())),atof(numb_convert(nu2->Pz())),atof(numb_convert(nu2->E())));

    if (1==1){

      //6 0 weight scale alpha_em alpha_strong
      myfile << "6   0  0.1587097E-05  0.5552668E+03  0.7546771E-02  0.1010447E+00" << std::endl;
      myfile << "   1   -1    0    0  501    0  " <<  numb_convert(0.) << " " << numb_convert(0.) << " " << numb_convert(pZ) << " " << numb_convert(pZ) << " " << numb_convert(0.) << "  0. -1. " << std::endl;
      myfile << "  -1   -1    0    0    0  501  " <<  numb_convert(0.) << " " << numb_convert(0.) << " " << numb_convert(-pZ) <<" " << numb_convert(pZ) << " " << numb_convert(0.)  << "  0. -1. " << std::endl;
      myfile << "   2    1    1    2  502    0  " <<  numb_convert(b1->Px()) << " " << numb_convert(b1->Py()) << " " << numb_convert(b1->Pz()) << " " << numb_convert(b1->E()) << " " << numb_convert(b1_truncated.M())  << "  0. -1. " << std::endl;
      myfile << "  -2    1    1    2    0  502  " <<  numb_convert(b2->Px()) << " " << numb_convert(b2->Py()) << " " << numb_convert(b2->Pz()) <<" " << numb_convert(b2->E()) << " " << numb_convert(b2_truncated.M())  << "  0. -1. " <<std::endl;
      myfile << "  14    1    1    2    0    0  " <<  numb_convert(nu1->Px()) << " " << numb_convert(nu1->Py()) << " " << numb_convert(nu1->Pz()) << " " << numb_convert(nu1->E()) << " " << numb_convert(nu1_truncated.M())  << "  0. -1. " <<std::endl;
      myfile << " -14    1    1    2    0    0  " <<  numb_convert(nu2->Px()) << " " << numb_convert(nu2->Py()) << " " << numb_convert(nu2->Pz()) <<" " << numb_convert(nu2->E()) << " " << numb_convert(nu2_truncated.M())  << "  0. -1. " << std::endl;

    }

    if (1==1){

      //6 0 weight scale alpha_em alpha_strong
      myfile2 << "6   0  0.1587097E-05  0.5552668E+03  0.7546771E-02  0.1010447E+00" << std::endl;
      myfile2 << "  21   -1    0    0  501  502  " <<  numb_convert(0.) << " " << numb_convert(0.) << " " << numb_convert(pZ) << " " << numb_convert(pZ) << " " << numb_convert(0.) << "  0. -1. " << std::endl;
      myfile2 << "  21   -1    0    0  503  501  " <<  numb_convert(0.) << " " << numb_convert(0.) << " " << numb_convert(-pZ) <<" " << numb_convert(pZ) << " " << numb_convert(0.)  << "  0. -1. " << std::endl;
      myfile2 << "   2    1    1    2  503    0  " <<  numb_convert(b1->Px()) << " " << numb_convert(b1->Py()) << " " << numb_convert(b1->Pz()) << " " << numb_convert(b1->E()) << " " << numb_convert(b1_truncated.M())  << "  0. -1. " << std::endl;
      myfile2 << "  -2    1    1    2    0  502  " <<  numb_convert(b2->Px()) << " " << numb_convert(b2->Py()) << " " << numb_convert(b2->Pz()) <<" " << numb_convert(b2->E()) << " " << numb_convert(b2_truncated.M())  << "  0. -1. " <<std::endl;
      myfile2 << "  14    1    1    2    0    0  " <<  numb_convert(nu1->Px()) << " " << numb_convert(nu1->Py()) << " " << numb_convert(nu1->Pz()) << " " << numb_convert(nu1->E()) << " " << numb_convert(nu1_truncated.M())  << "  0. -1. " <<std::endl;
      myfile2 << " -14    1    1    2    0    0  " <<  numb_convert(nu2->Px()) << " " << numb_convert(nu2->Py()) << " " << numb_convert(nu2->Pz()) <<" " << numb_convert(nu2->E()) << " " << numb_convert(nu2_truncated.M())  << "  0. -1. " << std::endl;

    }

    else if (ColorRep=="Triplet"){

      std::cout << "I'm sorry, that color representation is not implemented.  Please use Octet or Singlet.";
      exit(1);

    }

    else{
      std::cout << "I'm sorry, that color representation is not implemented.  Please use Octet or Singlet.";
      exit(1);
    }

    myfile << "</event>" << std::endl;
    myfile2 << "</event>" << std::endl;
  }

  myfile << "</LesHouchesEvents>" << std::endl;
  myfile2 << "</LesHouchesEvents>" << std::endl;
  std::cout << accepted << std::endl;

}
