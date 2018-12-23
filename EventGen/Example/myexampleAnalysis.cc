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
#include "TH2F.h"

#include "myexampleAnalysis.h"
#include "myTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"

#include "Pythia8/Pythia.h"

// #include "Nsubjettiness.h"
//#include "Njettiness.hh"
//#include "Nsubjettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

typedef pair<double, double> point;

double euclidean_distance(const point &x, const point &y)
{
    double d1 = x.first - y.first;
    double d2 = x.second - y.second;

    return sqrt(d1 * d1 + d2 * d2);
}

// Constructor 
myexampleAnalysis::myexampleAnalysis(int imagesize)
{
    int radius = imagesize;
    imagesize *= imagesize;
    MaxN = imagesize;
    fTIntensity = new float[imagesize];
    fTIntensity_pT = new float[imagesize];
    // fTRotatedIntensity = new float[imagesize];
    // fTLocalDensity = new float[imagesize];
    // fTGlobalDensity = new float[imagesize];


    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "octet_test.root";
    tool = new myTools();

    //model the detector as a 2D histogram   
    //                         xbins       y bins
    detector = new TH2D("", "", radius*5, -6.25, 6.25, radius*9, -11.25, 11.25);
    for(int i = 1; i <= radius*5; i++)
    {
        for (int j = 1; j <= radius*9; j++)
        {
            detector->SetBinContent(i,j,0);
        }
    }

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis End " << endl;
}

// Destructor 
myexampleAnalysis::~myexampleAnalysis()
{
    delete tool;

    delete[] fTIntensity;
    delete[] fTIntensity_pT;
    // delete[] fTRotatedIntensity;
    // delete[] fTGlobalDensity;
    // delete[] fTLocalDensity;
}

// Begin method
void myexampleAnalysis::Begin()
{
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for myexample");
   
   // for shit you want to do by hand
   DeclareBranches();
   ResetBranches();
   
   return;
}

// End
void myexampleAnalysis::End()
{
    tT->Write();
    tF->Close();
    return;
}

// Analyze
void myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV,
    int pixels, float range, bool onlyCharged)
{

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches 
    ResetBranches();
    
    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet> particlesForJets;
    std::vector <fastjet::PseudoJet> particlesForJets_nopixel;

    detector->Reset();
    TLorentzVector z1 = TLorentzVector();
    TLorentzVector z2 =TLorentzVector();

    //pythia8->event.list(); 

    // Particle loop ----------------------------------------------------------
    for (int ip=0; ip<pythia8->event.size(); ++ip){

      if (pythia8->event[ip].id()==14) z1.SetPxPyPzE(pythia8->event[ip].px(),pythia8->event[ip].py(),pythia8->event[ip].pz(),pythia8->event[ip].e());
      if (pythia8->event[ip].id()==-14) z1.SetPxPyPzE(pythia8->event[ip].px(),pythia8->event[ip].py(),pythia8->event[ip].pz(),pythia8->event[ip].e());

        fastjet::PseudoJet p(pythia8->event[ip].px(),
                             pythia8->event[ip].py(), 
                             pythia8->event[ip].pz(),
                             pythia8->event[ip].e() );

        // particles for jets --------------
        if (!pythia8->event[ip].isFinal()) continue;

	// only charged if wanted
	if (onlyCharged && !pythia8->event[ip].isCharged()) continue;

        //Skip neutrinos, PDGid = 12, 14, 16
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;


        // find the particles rapidity and phi, then get the detector bins
        int ybin = detector->GetXaxis()->FindBin(p.rapidity());
        int phibin = detector->GetYaxis()->FindBin(p.phi());

        // do bin += value in the associated detector bin
        detector->SetBinContent(ybin, phibin, 
                                detector->GetBinContent(ybin, phibin) + p.e());
	fastjet::PseudoJet p_nopix(p.px(),p.py(),p.pz(),p.e());
	particlesForJets_nopixel.push_back(p_nopix);
    }  
    // end particle loop -----------------------------------------------  

    //Now, we extract the energy from the calorimeter for processing by fastjet
    for (int i = 1; i <= detector->GetNbinsX(); i++)
    {
        for (int j = 1; j <= detector->GetNbinsY(); j++)
        {
            if (detector->GetBinContent(i, j) > 0)
            {
                double phi = detector->GetYaxis()->GetBinCenter(j);
                double eta = detector->GetXaxis()->GetBinCenter(i);
                double E = detector->GetBinContent(i, j);
                fastjet::PseudoJet p(0., 0., 0., 0.);

                //We measure E (not pT)!  And treat 'clusters' as massless.
                p.reset_PtYPhiM(E/cosh(eta), eta, phi, 0.); 
                particlesForJets.push_back(p);
            }
        }
    }

    fastjet::JetDefinition *m_jet_def = new fastjet::JetDefinition(
        fastjet::antikt_algorithm, 1.);

    fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
        fastjet::SelectorPtFractionMin(0.05));

    fastjet::ClusterSequence csLargeR(particlesForJets, *m_jet_def);
    fastjet::ClusterSequence csLargeR_nopix(particlesForJets_nopixel, *m_jet_def);

    vector<fastjet::PseudoJet> considered_jets = fastjet::sorted_by_pt(
								       csLargeR.inclusive_jets(10.0));
    vector<fastjet::PseudoJet> considered_jets_nopix = fastjet::sorted_by_pt(
									     csLargeR_nopix.inclusive_jets(10.0));
    fastjet::PseudoJet leading_jet = trimmer(considered_jets[0]);
    fastjet::PseudoJet leading_jet_nopix = trimmer(considered_jets_nopix[0]);

    fTLeadingEta = leading_jet.eta();
    fTLeadingPhi = leading_jet.phi();
    fTLeadingPt = leading_jet.perp();
    fTLeadingM = leading_jet.m();
    fTLeadingEta_nopix = leading_jet_nopix.eta();
    fTLeadingPhi_nopix = leading_jet_nopix.phi();
    fTLeadingPt_nopix = leading_jet_nopix.perp();
    fTLeadingM_nopix = leading_jet_nopix.m();

    vector<fastjet::PseudoJet> subjets_nopix = leading_jet_nopix.pieces();
    
    fTpull1 = -1;
    fTpull2 = -1;

    fTpull1_nopix = -1;
    fTpull2_nopix = -1;

    if (subjets_nopix.size() > 1){
      fTpull1_nopix = tool->JetPull(subjets_nopix[0],subjets_nopix[1]);
      fTpull2_nopix = tool->JetPull(subjets_nopix[1],subjets_nopix[0]);
    }
   
    fTdeltaR = 0.;
    if (leading_jet.pieces().size() > 1){
      vector<fastjet::PseudoJet> subjets = leading_jet.pieces();
      TLorentzVector l(subjets[0].px(),subjets[0].py(),subjets[0].pz(),subjets[0].E());
      TLorentzVector sl(subjets[1].px(),subjets[1].py(),subjets[1].pz(),subjets[1].E());
      fTdeltaR = l.DeltaR(sl); 
      fTSubLeadingEta = sl.Eta()-l.Eta();
      fTSubLeadingPhi = subjets[1].delta_phi_to(subjets[0]);

      //Let's compute the jet pull.
      fTpull1 = tool->JetPull(subjets[0],subjets[1]);
      fTpull2 = tool->JetPull(subjets[1],subjets[0]);
    }
    
    vector<pair<double, double>  > consts_image;
    vector<fastjet::PseudoJet> sorted_consts = sorted_by_pt(leading_jet.constituents());

    for(int i = 0; i < sorted_consts.size(); i++)
    {
        pair<double, double> const_hold;
        const_hold.first = sorted_consts[i].eta();
        const_hold.second = sorted_consts[i].phi();
        consts_image.push_back(const_hold);
    }

    vector<fastjet::PseudoJet> subjets = leading_jet.pieces();

    //Step 1: Center on the jet axis.
    fastjet::PseudoJet shiftjet;
    shiftjet.reset_momentum_PtYPhiM(1.,detector->GetXaxis()->GetBinCenter(detector->GetXaxis()->FindBin(subjets[0].eta())),detector->GetYaxis()->GetBinCenter(detector->GetYaxis()->FindBin(subjets[0].phi())),0.);
    for (int i =0; i < sorted_consts.size(); i++)
    {
      //std::cout << "check this " << i << " " << 1000*detector->GetXaxis()->GetBinCenter(detector->GetXaxis()->FindBin(consts_image[i].first-shiftjet.eta())) << " " << 1000*(consts_image[i].first-shiftjet.eta()) << " " << 1000*detector->GetXaxis()->GetBinCenter(detector->GetXaxis()->FindBin(consts_image[i].first-shiftjet.eta()))/(consts_image[i].first-shiftjet.eta()) << " " << 1000*(consts_image[i].first-shiftjet.eta()) << " " << 1000*consts_image[i].first << " " << 1000*shiftjet.eta() << std::endl;
      consts_image[i].first = consts_image[i].first-shiftjet.eta();
      consts_image[i].second = sorted_consts[i].delta_phi_to(shiftjet);
      //There is a slight bias due to the occations when phi is re-shifted.
      //std::cout << sorted_consts[i].delta_phi_to(shiftjet) << " " << -sorted_consts[i].phi()+shiftjet.phi() << " " << -sorted_consts[i].phi()+shiftjet.phi()-8*atan(1) << " " << -sorted_consts[i].phi()+shiftjet.phi()+8*atan(1) << std::endl;
      /*
      double mydphi = -sorted_consts[i].phi()+shiftjet.phi();
      if (mydphi < -4*atan(1)) mydphi+8*atan(1);
      if (mydphi > 4*atan(1)) mydphi-=8*atan(1);
      std::cout << sorted_consts[i].delta_phi_to(shiftjet) << " " << mydphi << " " << -sorted_consts[i].phi()+shiftjet.phi() << std::endl;
      if (fabs(1000*(detector->GetYaxis()->GetBinCenter(detector->GetYaxis()->FindBin(sorted_consts[i].delta_phi_to(shiftjet))))/consts_image[i].second - 1000) > 25){
	std::cout << "hello " << 1000*(sorted_consts[i].phi()-shiftjet.phi()) << " " << 1000*sorted_consts[i].phi() << " " << 1000*shiftjet.phi() << " " << 1000*(detector->GetYaxis()->GetBinCenter(detector->GetYaxis()->FindBin(sorted_consts[i].delta_phi_to(shiftjet)))) << " " << 1000*consts_image[i].second << " " << 1000*(detector->GetYaxis()->GetBinCenter(detector->GetYaxis()->FindBin(sorted_consts[i].delta_phi_to(shiftjet))))/consts_image[i].second << std::endl;
	exit(1);
      } 
      */
      //consts_image[i].first = consts_image[i].first-subjets[0].eta();
      //consts_image[i].second = sorted_consts[i].delta_phi_to(subjets[0]); //use delta phi to take care of the dis-continuity in phi
    }
    
    //Quickly run PCA for the rotation.
    double xbar = 0.;
    double ybar = 0.;
    double x2bar = 0.;
    double y2bar = 0.;
    double xybar = 0.;
    double n = 0;

    for(int i = 0; i < leading_jet.constituents().size(); i++)
      {
        double x = consts_image[i].first;
        double y = consts_image[i].second;
	double E = sorted_consts[i].e();
        n+=E;
        xbar+=x*E;
        ybar+=y*E;
      }

    double mux = xbar / n;
    double muy = ybar / n;

    xbar = 0.;
    ybar = 0.;
    n = 0.;

    for(int i = 0; i < leading_jet.constituents().size(); i++)
      {
        double x = consts_image[i].first - mux;
        double y = consts_image[i].second - muy;
        double E = sorted_consts[i].e();
        n+=E;
        xbar+=x*E;
        ybar+=y*E;
        x2bar+=x*x*E;
        y2bar+=y*y*E;
        xybar+=x*y*E;
      }

    double sigmax2 = x2bar / n - mux*mux;
    double sigmay2 = y2bar / n - muy*muy;
    double sigmaxy = xybar / n - mux*muy;
    double lamb_min = 0.5* ( sigmax2 + sigmay2 - sqrt( (sigmax2-sigmay2)*(sigmax2-sigmay2) + 4*sigmaxy*sigmaxy) );
    double lamb_max = 0.5* ( sigmax2 + sigmay2 + sqrt( (sigmax2-sigmay2)*(sigmax2-sigmay2) + 4*sigmaxy*sigmaxy) );

    double dir_x = sigmax2+sigmaxy-lamb_min;
    double dir_y = sigmay2+sigmaxy-lamb_min;

    //The first PC is only defined up to a sign.  Let's have it point toward the side of the jet with the most energy.

    double Eup = 0.;
    double Edn = 0.;

    for(int i = 0; i < leading_jet.constituents().size(); i++)
      {
	double x = consts_image[i].first - mux;
        double y = consts_image[i].second - muy;
	double E = sorted_consts[i].e();
	double dotprod = dir_x*x+dir_y*y;
	if (dotprod > 0) Eup+=E;
	else Edn+=E;
      }
    
    if (Edn < Eup){
      dir_x = -dir_x;
      dir_y = -dir_y;
    }

    fTPCEta = dir_x;
    fTPCPhi = dir_y;

    //Doing a little check to see how often the PC points in the direction of the subleading subjet if it exists.
    //std::cout << "new event " << 100*(dir_x) << " " << 100*(dir_y) << " " << leading_jet.m() << " " << leading_jet.perp() << " " << sigmax2*100 << " " << sigmay2*100 << " " << sigmaxy*100 << " " << 100*lamb_min << " " << 100*lamb_max << " " << 100*(sigmax2+sigmaxy-lamb_max) << " " << 100*(sigmay2+sigmaxy-lamb_max) << std::endl; 
    //for (int i=0; i<leading_jet.pieces().size(); i++){
    //  std::cout << i << " " << (leading_jet.pieces()[i].eta() - subjets[0].eta())*100 << " " << leading_jet.pieces()[i].delta_phi_to(subjets[0])*100 << " " << leading_jet.pieces()[i].e() << std::endl; 
    //}
      
    //return;

    //for (int i =1; i >= 0; i--)
    //{
    //    subjets_image[i].first = subjets_image[i].first - subjets_image[0].first;
    //    subjets_image[i].second = subjets_image[i].second- subjets_image[0].second;
    //}

    //Step 2: Fill in the unrotated image
    //-------------------------------------------------------------------------   
    range = 1.25;
    TH2D* orig_im = new TH2D("", "", pixels, -range, range, pixels, -range, range);
    /*
    for (int i=1; i<=orig_im->GetNbinsX(); i++){
      for (int j=1; j<=orig_im->GetNbinsY(); j++){
	orig_im->SetBinContent(i,j,0.);
      }
    }
    */
    //std::cout << 100*(orig_im->GetXaxis()->GetBinCenter(1)-0.5*orig_im->GetXaxis()->GetBinWidth(2)) << " " << 100*orig_im->GetXaxis()->GetBinWidth(2) << " " << pixels << " " << range << std::endl;

    TLorentzVector image_mass = TLorentzVector();
    for (int i = 0; i < sorted_consts.size(); i++)
    {
      TLorentzVector hold = TLorentzVector();
      hold.SetPtEtaPhiM(sorted_consts[i].perp(),consts_image[i].first,consts_image[i].second,0.);
      //image_mass+=hold;
      //int myi = orig_im->GetXaxis()->FindBin(consts_image[i].first);
      //int myj =orig_im->GetYaxis()->FindBin(consts_image[i].second);
      //orig_im->SetBinContent(myi,myj,hold.E());
      orig_im->Fill(consts_image[i].first,consts_image[i].second,hold.E());
      //double myeta = orig_im->GetXaxis()->GetBinCenter(myi);
      //double myphi = orig_im->GetYaxis()->GetBinCenter(myj);
      //std::cout << "check this " << myi << " " << myj << " " << 1000*myeta << " " << 1000*consts_image[i].first << " " << 1000*myeta/consts_image[i].first << " " << hold.Pt()*1000 << " " << 1000*hold.Eta() << " " << hold.Phi()*1000 << " " << hold.M()*1000 << std::endl;
      //std::cout << "     " << hold.Eta() << " " << hold.Phi() << " " << hold.Pt() << " " << orig_im->GetBinContent(myi,myj)/cosh(myeta) << " " << orig_im->GetBinContent(myi,myj)/cosh(hold.Eta()) << " " << 1000*myeta/hold.Eta() << " " << 1000*orig_im->GetBinContent(myi,myj)/hold.E() << std::endl;
    }

    /*
    TLorentzVector image_mass2 = TLorentzVector(0.,0.,0.,0.);
    for (int i=1; i<=orig_im->GetNbinsX(); i++){
      for (int j=1; j<=orig_im->GetNbinsY(); j++){
	TLorentzVector myhold = TLorentzVector();
	double myeta = orig_im->GetXaxis()->GetBinCenter(i);
	double myphi = orig_im->GetYaxis()->GetBinCenter(j);
	myhold.SetPtEtaPhiM(orig_im->GetBinContent(i,j)/cosh(myeta),myeta,myphi,0.);
	image_mass2+=myhold;
	//if (orig_im->GetBinContent(i,j) > 0) std::cout << "     bbb " << " " << i << " " << j << " " << myeta*1000 << " " << myphi*1000 << " " << myhold.Pt()*1000 << std::endl;  
      }
    }
    */

    //std::cout << fTLeadingM << " " << image_mass2.M() << std::endl;
    
    /*
    if (fabs(fTLeadingM-image_mass2.M()) > 0.1){
      std::cout << fTLeadingM << " " << leading_jet.perp() << " " << leading_jet.e() << " " << leading_jet.eta() << " " << image_mass.M() << " " << image_mass.Pt() << " " << image_mass.E() << " " << image_mass.Eta() << " " << image_mass2.M() << " " << image_mass2.Pt() << " " << image_mass2.E() << " " << image_mass2.Eta() << std::endl;
      for (int i = 0; i < sorted_consts.size(); i++)
	{
	  TLorentzVector hold = TLorentzVector();
	  hold.SetPtEtaPhiM(sorted_consts[i].perp(),consts_image[i].first,consts_image[i].second,0.);
	  std::cout << "     " << i << " " <<  hold.Pt()*1000 << " " << hold.Eta()*1000 << " " << hold.Phi()*1000 << std::endl;   
	}
      std::cout << "  and now for the image " << std::endl;
      for (int i=1; i<=orig_im->GetNbinsX(); i++){
	for (int j=1; j<=orig_im->GetNbinsY(); j++){
	  TLorentzVector myhold = TLorentzVector();
	  double myeta = orig_im->GetXaxis()->GetBinCenter(i);
	  double myphi = orig_im->GetYaxis()->GetBinCenter(j);
	  myhold.SetPtEtaPhiM(orig_im->GetBinContent(i,j)/cosh(myeta),myeta,myphi,0.);
	  if (orig_im->GetBinContent(i,j) > 0) std::cout << "     " << i << " " <<  j << " " << myhold.Pt()*1000 << " "<< myhold.Eta()*1000<< " " << myhold.Phi()*1000 << std::endl;
	}
      }
    }
    */

    //Step 2b): fill in the density
    //-------------------------------------------------------------------------
    // double Rlocal = 0.5;
    // double Rglobal = 1.0;
    // TH2F* localdensity = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // TH2F* globaldensity = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // for (int i = 0; i < sorted_consts.size(); ++i)
    // {
    //     auto reference_px = consts_image[i];
    //     auto ref_pt = sorted_consts[i].e();

    //     double global_acc = 0;
    //     double local_acc = 0;

    //     for (int j = 0; j < sorted_consts.size(); ++j)
    //     {   
    //         if (j != i)
    //         {
    //             auto dR = euclidean_distance(reference_px, consts_image[j]);
    //             if (dR <= Rglobal)
    //             {
    //                 global_acc += (ref_pt / dR);
    //                 if (dR <= Rlocal)
    //                 {
    //                     local_acc += (ref_pt / dR);
    //                 }
    //             }
    //         }
    //     }
    //     localdensity->Fill(reference_px.first, reference_px.second, local_acc);
    //     globaldensity->Fill(reference_px.first, reference_px.second, global_acc);
    // }


    

    //Step 3: Rotate so the subleading subjet is at -pi/2
    //-------------------------------------------------------------------------
    /*
    fTSubLeadingEta = subjets_image[1].first;
    fTSubLeadingPhi = subjets_image[1].second;
    double theta = atan(subjets_image[1].second/subjets_image[1].first)+2.*atan(1.); //atan(1) = pi/4
    
    if (-sin(theta)*subjets_image[1].first+cos(theta)*subjets_image[1].second > 0)
    {
        theta+=-4.*atan(1.);
    }

    fTRotationAngle = theta;

    for (int i=0; i<sorted_consts.size(); i++)
    {
        double x = consts_image[i].first;
        double y = consts_image[i].second;
        consts_image[i].first = cos(theta)*x + sin(theta)*y;
        consts_image[i].second = -sin(theta)*x + cos(theta)*y;
    }

    for (int i=0; i<2; i++)
    {
        double x = subjets_image[i].first;
        double y = subjets_image[i].second;
        subjets_image[i].first = cos(theta)*x + sin(theta)*y;
        subjets_image[i].second =  -sin(theta)*x + cos(theta)*y;
    }

    theta = atan(subjets_image[1].second/subjets_image[1].first);
    */

    // //Step 4b): fill in the rotated image
    // //-------------------------------------------------------------------------
    // TH2F* rotatedimage = new TH2F("", "", pixels, -range, range, pixels, -range, range);
    // for (int i = 0; i < sorted_consts.size(); i++)
    // {
    //     rotatedimage->Fill(consts_image[i].first,consts_image[i].second,sorted_consts[i].e());
    // }
    

    //Step 5: Dump the images in the tree!
    //-------------------------------------------------------------------------
    int counter=0;
    for (int i=1; i<=orig_im->GetNbinsX(); i++)
    {
        for (int j=1; j<=orig_im->GetNbinsY(); j++)
        {
            // fTRotatedIntensity[counter] = rotatedimage->GetBinContent(i,j);
            fTIntensity_pT[counter] = orig_im->GetBinContent(i,j);
	    double myeta = orig_im->GetXaxis()->GetBinCenter(i);
	    fTIntensity[counter] = orig_im->GetBinContent(i,j)/cosh(myeta);
            // fTLocalDensity[counter] = localdensity->GetBinContent(i, j);
            // fTGlobalDensity[counter] = globaldensity->GetBinContent(i, j);

            counter++;
        }
    }

    // Step 6: Fill in nsubjettiness (new)
    //----------------------------------------------------------------------------
    OnePass_WTA_KT_Axes axis_spec;
    NormalizedMeasure parameters(1.0, 1.0);

    // NormalizedMeasure parameters(1.0, 1.0);
    Nsubjettiness subjettiness_1(1, axis_spec, parameters);
    Nsubjettiness subjettiness_2(2, axis_spec, parameters);
    Nsubjettiness subjettiness_3(3, axis_spec, parameters);

    fTTau1 = (float) subjettiness_1.result(leading_jet);
    fTTau2 = (float) subjettiness_2.result(leading_jet);
    fTTau3 = (float) subjettiness_3.result(leading_jet);

    fTTau32 = (abs(fTTau2) < 1e-4 ? -10 : fTTau3 / fTTau2);
    fTTau21 = (abs(fTTau1) < 1e-4 ? -10 : fTTau2 / fTTau1);

    fTTau1_nopix = (float) subjettiness_1.result(leading_jet_nopix);
    fTTau2_nopix = (float) subjettiness_2.result(leading_jet_nopix);
    fTTau3_nopix = (float) subjettiness_3.result(leading_jet_nopix);

    fTTau32_nopix = (abs(fTTau2_nopix) < 1e-4 ? -10 : fTTau3_nopix / fTTau2_nopix);
    fTTau21_nopix = (abs(fTTau1_nopix) < 1e-4 ? -10 : fTTau2_nopix / fTTau1_nopix);

    // // Step 7: Fill in nsubjettiness (old)
    // //----------------------------------------------------------------------------
    // OnePass_KT_Axes axis_spec_old;
    // NormalizedMeasure parameters_old(1.0, 1.0);

    // Nsubjettiness old_subjettiness_1(1, axis_spec_old, parameters_old);
    // Nsubjettiness old_subjettiness_2(2, axis_spec_old, parameters_old);
    // Nsubjettiness old_subjettiness_3(3, axis_spec_old, parameters_old);

    // fTTau1 = (float) old_subjettiness_1.result(leading_jet);
    // fTTau2 = (float) old_subjettiness_2.result(leading_jet);
    // fTTau3 = (float) old_subjettiness_3.result(leading_jet);

    // fTTau32old = (abs(fTTau2) < 1e-4 ? -10 : fTTau3 / fTTau2);
    // fTTau21old = (abs(fTTau1) < 1e-4 ? -10 : fTTau2 / fTTau1);

    // if (fTLeadingM > 65 && fTLeadingM < 95) {
    // if (fTTau21 > 0.19 && fTTau21 < 0.21){
	tT->Fill();
	///}
	// }

    return;
}

// declate branches
void myexampleAnalysis::DeclareBranches()
{

    // Event Properties 
    tT->Branch("NFilled", &fTNFilled, "NFilled/I");

    tT->Branch("Intensity", *&fTIntensity, "Intensity[NFilled]/F");
    tT->Branch("Intensity_pT", *&fTIntensity_pT, "Intensity_pT[NFilled]/F");

    // tT->Branch("LocalDensity", *&fTLocalDensity, "LocalDensity[NFilled]/F");
    // tT->Branch("GlobalDensity", *&fTGlobalDensity, "GlobalDensity[NFilled]/F");

    // tT->Branch("RotatedIntensity", 
    //     *&fTRotatedIntensity, "RotatedIntensity[NFilled]/F");

    tT->Branch("SubLeadingEta", &fTSubLeadingEta, "SubLeadingEta/F");
    tT->Branch("SubLeadingPhi", &fTSubLeadingPhi, "SubLeadingPhi/F");

    tT->Branch("PCEta", &fTPCEta, "PCEta/F");
    tT->Branch("PCPhi", &fTPCPhi, "PCPhi/F");

    tT->Branch("LeadingEta", &fTLeadingEta, "LeadingEta/F");
    tT->Branch("LeadingPhi", &fTLeadingPhi, "LeadingPhi/F");
    tT->Branch("LeadingPt", &fTLeadingPt, "LeadingPt/F");
    tT->Branch("LeadingM", &fTLeadingM, "LeadingM/F");
    //tT->Branch("RotationAngle", &fTRotationAngle, "RotationAngle/F");

    tT->Branch("LeadingEta_nopix", &fTLeadingEta_nopix, "LeadingEta_nopix/F");
    tT->Branch("LeadingPhi_nopix", &fTLeadingPhi_nopix, "LeadingPhi_nopix/F");
    tT->Branch("LeadingPt_nopix", &fTLeadingPt_nopix, "LeadingPt_nopix/F");
    tT->Branch("LeadingM_nopix", &fTLeadingM_nopix, "LeadingM_nopix/F");

    tT->Branch("Tau1", &fTTau1, "Tau1/F");
    tT->Branch("Tau2", &fTTau2, "Tau2/F");
    tT->Branch("Tau3", &fTTau3, "Tau3/F");

    tT->Branch("Tau1_nopix", &fTTau1_nopix, "Tau1_nopix/F");
    tT->Branch("Tau2_nopix", &fTTau2_nopix, "Tau2_nopix/F");
    tT->Branch("Tau3_nopix", &fTTau3_nopix, "Tau3_nopix/F");

    tT->Branch("DeltaR", &fTdeltaR, "DeltaR/F");

    tT->Branch("Tau32", &fTTau32, "Tau32/F");
    tT->Branch("Tau21", &fTTau21, "Tau21/F");
    
    tT->Branch("Tau32_nopix", &fTTau32_nopix, "Tau32_nopix/F");
    tT->Branch("Tau21_nopix", &fTTau21_nopix, "Tau21_nopix/F");

    tT->Branch("pull1", &fTpull1, "pull1/F");
    tT->Branch("pull2", &fTpull2, "pull2/F");

    tT->Branch("pull1_nopix", &fTpull1_nopix, "pull1_nopix/F");
    tT->Branch("pull2_nopix", &fTpull2_nopix, "pull2_nopix/F");

    //Need to add the jet pull and a few other things...
    
    // tT->Branch("Tau32old", &fTTau32old, "Tau32old/F");
    // tT->Branch("Tau21old", &fTTau21old, "Tau21old/F");
    return;
}

// resets vars
void myexampleAnalysis::ResetBranches(){
    // reset branches 
    fTNFilled = MaxN;
    fTSubLeadingPhi = -999;
    fTSubLeadingEta = -999;
    fTPCPhi = -999;
    fTPCEta = -999;
    //fTRotationAngle = -999;
  
    fTTau32 = -999;
    fTTau21 = -999;

    fTTau1 = -999;
    fTTau2 = -999;
    fTTau3 = -999;

    fTTau32_nopix = -999;
    fTTau21_nopix = -999;

    fTTau1_nopix = -999;
    fTTau2_nopix = -999;
    fTTau3_nopix = -999;

    // fTTau32old = -999;
    // fTTau21old = -999;

    fTLeadingEta = -999;
    fTLeadingPhi = -999;
    fTLeadingPt = -999;
    fTLeadingM = -999;

    fTLeadingEta_nopix = -999;
    fTLeadingPhi_nopix = -999;
    fTLeadingPt_nopix = -999;
    fTLeadingM_nopix = -999;

    for (int iP=0; iP < MaxN; ++iP)
    {
        fTIntensity[iP]= -999;
	fTIntensity_pT[iP]= -999;
        // fTRotatedIntensity[iP]= -999;
        // fTLocalDensity[iP]= -999;
        // fTGlobalDensity[iP]= -999;
    }
}
