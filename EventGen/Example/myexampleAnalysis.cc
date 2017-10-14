#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "myexampleAnalysis.h"
#include "myTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// Constructor 
myexampleAnalysis::myexampleAnalysis(){
    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new myTools();

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_largeR_ALTAS  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis End " << endl;
}

// Destructor 
myexampleAnalysis::~myexampleAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_largeR_ALTAS;
}

// Begin method
void myexampleAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for myexample");
    
   DeclareBranches();
   ResetBranches();
   

   return;
}

// End
void myexampleAnalysis::End(){
    
    tT->Write();
    tF->Close();
    return;
}

// Analyze
void myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    //pythia8->event.list();

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    std::vector <fastjet::PseudoJet>           bhadrons;
    std::vector <fastjet::PseudoJet>           chadrons;
    std::vector <fastjet::PseudoJet>           HS;
    
    // Particle loop -----------------------------------------------------------
    for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

      fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
      p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));
      
      if (tool->IsBHadron(pythia8->event[ip].id())){
        bhadrons.push_back(p);
      }

      // particles for jets --------------
      if (!pythia8->event[ip].isFinal() )      continue;
      if (fabs(pythia8->event[ip].id())  ==12) continue;
      if (fabs(pythia8->event[ip].id())  ==13) continue;
      if (fabs(pythia8->event[ip].id())  ==14) continue;
      if (fabs(pythia8->event[ip].id())  ==16) continue;

      particlesForJets.push_back(p);
      
    } // end particle loop -----------------------------------------------

    // small R jets: ATLAS Style ------------------------------------------
    fastjet::ClusterSequence cs(particlesForJets, *m_jet_def);
    vector<fastjet::PseudoJet> myJets = fastjet::sorted_by_pt(cs.inclusive_jets(25.0)); //was 10
    int nbjets = 0;
    TRandom3 *rand = new TRandom3(0);
    for (unsigned int ij = 0; ij < myJets.size(); ij++){
        if(fTNJetsSmallRFilled == MaxNJetSmallR) {cout << "Warning: More than " << MaxNJetSmallR << " small R jets" << endl; continue;}
        if(myJets[ij].pt() < 15) continue;
        fTJsmallPt       [fTNJetsSmallRFilled] = myJets[ij].pt();
        fTJsmallEta      [fTNJetsSmallRFilled] = myJets[ij].eta();
        fTJsmallPhi      [fTNJetsSmallRFilled] = myJets[ij].phi();
        fTJsmallM        [fTNJetsSmallRFilled] = myJets[ij].m();
	fTJsmallBtag     [fTNJetsSmallRFilled] = tool->Btag(myJets[ij],bhadrons,chadrons,0.4,1.,1000000000,1000000000);
	fTJsmallTrackMass[fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],1);
	fTJsmallTrackpT  [fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],2);
	fTJsmallTrackMassR[fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],-1);
        fTJsmallTrackpTR  [fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],-2);
	fTsmallntrack    [fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],3);
	fTsmallntrackR    [fTNJetsSmallRFilled] = tool->JetTrackMass(myJets[ij],-3);
	nbjets+=fTJsmallBtag     [fTNJetsSmallRFilled];
	fTNJetsSmallRFilled++;
    }
    tT->Fill();
      
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void myexampleAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   
   // smallR jets
   tT->Branch("NJetsFilledSmallR",         &fTNJetsSmallRFilled,       "NJetsFilledSmallR/I");
   tT->Branch("JsmallPt",                  &fTJsmallPt,                "JsmallPt[NJetsFilledSmallR]/F");
   tT->Branch("JsmallEta",                 &fTJsmallEta,               "JsmallEta[NJetsFilledSmallR]/F");
   tT->Branch("JsmallPhi",                 &fTJsmallPhi,               "JsmallPhi[NJetsFilledSmallR]/F");
   tT->Branch("JsmallM",                   &fTJsmallM,                 "JsmallM[NJetsFilledSmallR]/F");
   tT->Branch("JsmallBtag", &fTJsmallBtag, "JsmallBtag[NJetsFilledSmallR]/I");
   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void myexampleAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;

      fTNJetsSmallRFilled=0;
      for (int iP=0; iP < MaxNJetSmallR; ++iP){
          fTJsmallPt      [iP]= -999;
          fTJsmallPhi     [iP]= -999;
          fTJsmallEta     [iP]= -999;
          fTJsmallM       [iP]= -999;
	  fTJsmallBtag[iP]=-999;
	  fTsmallntrack[iP]=-999;
	  fTJsmallTrackpT[iP]=-999;
	  fTJsmallTrackMass[iP]=-999;
	  fTsmallntrackR[iP]=-999;
          fTJsmallTrackpTR[iP]=-999;
          fTJsmallTrackMassR[iP]=-999;
      }
}
