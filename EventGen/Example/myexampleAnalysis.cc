#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <tuple>

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
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "Pythia8/Pythia.h"

#include "fastjet/contrib/Nsubjettiness.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

typedef pair<double, double> point;


typedef vector<double> Row; // One row of the matrix
typedef vector<Row> Matrix; // Matrix: a vector of rows
typedef vector<Matrix> Mat3d; // Matrix3D: a vector of Matrices

double euclidean_distance(const point &x, const point &y)
{
    double d1 = x.first - y.first;
    double d2 = x.second - y.second;

    return sqrt(d1 * d1 + d2 * d2);
}

std::tuple<double, double> calculate_pull(vector<fastjet::PseudoJet> subjets, myTools* tool) {
    if (subjets.size() > 1){
      double p1 = tool->JetPull(subjets[0],subjets[1]);
      double p2 = tool->JetPull(subjets[1],subjets[0]);
      return std::make_tuple(p1, p2);
    }
    else {
	throw 911;
        return std::make_tuple(NULL, NULL);
    }
}

// Constructor
myexampleAnalysis::myexampleAnalysis(int imagesize)
{
    int radius = imagesize;
    imagesize *= imagesize;
    MaxN = imagesize;
    fTIntensity_standard = new float[imagesize];
    fTIntensity_pT_standard = new float[imagesize];

    fTIntensity_charged = new float[imagesize];
    fTIntensity_pT_charged = new float[imagesize];

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    this->SetOutName("octet_test");
    tool = new myTools();

    //model the detector as a 2D histogram
    //                         xbins       y bins
    detector_standard = new TH2D("", "", radius*5, -6.25, 6.25, radius*9, -11.25, 11.25);
    for(int i = 1; i <= radius*5; i++)
    {
        for (int j = 1; j <= radius*9; j++)
        {
            detector_standard->SetBinContent(i,j,0);
        }
    }

    detector_charged = new TH2D("", "", radius*5, -6.25, 6.25, radius*9, -11.25, 11.25);
    for(int i = 1; i <= radius*5; i++)
    {
        for (int j = 1; j <= radius*9; j++)
        {
            detector_charged->SetBinContent(i,j,0);
        }
    }

    if(fDebug) cout << "myexampleAnalysis::myexampleAnalysis End " << endl;
}

// Destructor
myexampleAnalysis::~myexampleAnalysis()
{
    delete tool;

    delete[] fTIntensity_standard;
    delete[] fTIntensity_pT_standard;

    delete[] fTIntensity_charged;
    delete[] fTIntensity_pT_charged;
}

// Begin method
void myexampleAnalysis::Begin()
{
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");

   tT_standard = new TTree("StandardEventTree", "Event Tree for myexample");
   tT_charged = new TTree("ChargedEventTree", "Event Tree for myexample");

   // for stuff you want to do by hand
   DeclareBranches();
   ResetBranches();

   return;
}

// End
void myexampleAnalysis::End()
{
    tT_standard->Write();
    tT_charged->Write();

    tF->Close();
    return;
}

// Analyze
void myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV,
    int pixels, float range, float ptjMin, float ptjMax, float etaMax, float massMin, float massMax,
    bool untrim, bool cambridge, bool reproduce)
{

    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if (!pythia8->next()) return;
    if(fDebug) cout << "myexampleAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches
    ResetBranches();

    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet> particlesForJets_standard;
    std::vector <fastjet::PseudoJet> particlesForJets_nopixel_standard;


    std::vector <fastjet::PseudoJet> particlesForJets_charged;
    std::vector <fastjet::PseudoJet> particlesForJets_nopixel_charged;

    detector_standard->Reset();
    detector_charged->Reset();

    TLorentzVector z1 = TLorentzVector();
    TLorentzVector z2 =TLorentzVector();

    for (int ip=0; ip<pythia8->event.size(); ++ip){

        if (pythia8->event[ip].id()==14) z1.SetPxPyPzE(pythia8->event[ip].px(),pythia8->event[ip].py(),pythia8->event[ip].pz(),pythia8->event[ip].e());
        if (pythia8->event[ip].id()==-14) z1.SetPxPyPzE(pythia8->event[ip].px(),pythia8->event[ip].py(),pythia8->event[ip].pz(),pythia8->event[ip].e());

        fastjet::PseudoJet p(pythia8->event[ip].px(),
                             pythia8->event[ip].py(),
                             pythia8->event[ip].pz(),
                             pythia8->event[ip].e());

        // particles for jets --------------
        if (!pythia8->event[ip].isFinal()) continue;

        //Skip neutrinos, PDGid = 12, 14, 16
        if (fabs(pythia8->event[ip].id())  ==12) continue;
        if (fabs(pythia8->event[ip].id())  ==14) continue;
        if (fabs(pythia8->event[ip].id())  ==16) continue;
        if (reproduce && !isjetc(&pythia8->event[ip])) continue;

        // find the particles rapidity and phi, then get the detector bins
        int ybin = detector_standard->GetXaxis()->FindBin(p.rapidity());
        int phibin = detector_standard->GetYaxis()->FindBin(p.phi());

        // do bin += value in the associated detector bin
        detector_standard->SetBinContent(ybin, phibin,
            detector_standard->GetBinContent(ybin, phibin) + p.e());

        fastjet::PseudoJet p_nopix(p.px(),p.py(),p.pz(),p.e());
	    particlesForJets_nopixel_standard.push_back(p_nopix);

        if (pythia8->event[ip].isCharged()) {
            int ybin = detector_charged->GetXaxis()->FindBin(p.rapidity());
            int phibin = detector_charged->GetYaxis()->FindBin(p.phi());

            detector_charged->SetBinContent(ybin, phibin,
                detector_charged->GetBinContent(ybin, phibin) + p.e());

            particlesForJets_nopixel_charged.push_back(p_nopix);
        }
    }
    // end particle loop -----------------------------------------------

    // Calculate the nopixalated leading standard jet (to go to the cutoff as fast as possible).
    fastjet::JetDefinition *m_jet_def;
    if (cambridge) {
        m_jet_def = new fastjet::JetDefinition(
                fastjet::cambridge_algorithm, 1.);
    } else {
        m_jet_def = new fastjet::JetDefinition(
                fastjet::antikt_algorithm, 1.);
    }

    // Trimming
    fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
        fastjet::SelectorPtFractionMin(0.05));


    fastjet::ClusterSequence csLargeR_nopix_standard(particlesForJets_nopixel_standard, *m_jet_def);

    vector<fastjet::PseudoJet> considered_jets_nopix_standard = fastjet::sorted_by_pt(
        csLargeR_nopix_standard.inclusive_jets(10.0));


    fastjet::PseudoJet leading_jet_nopix_standard = trimmer(considered_jets_nopix_standard[0]);

    fTLeadingEta_nopix_standard = leading_jet_nopix_standard.eta();
    fTLeadingPt_nopix_standard = leading_jet_nopix_standard.perp();
    fTLeadingM_nopix_standard = leading_jet_nopix_standard.m();

    // Cut condition for the trimmed jets
    if (!untrim && !(fTLeadingPt_nopix_standard < ptjMax && fTLeadingPt_nopix_standard > ptjMin &&
            fTLeadingEta_nopix_standard < etaMax && fTLeadingM_nopix_standard > massMin && 
            fTLeadingM_nopix_standard < massMax)) {
        // One of the conditions did not pass.
        // TODO: ADD GARBAGE COLLECTION/FREE MEMORY UP
        return;
    } else if (untrim && !(fTLeadingM_nopix_standard > massMin && fTLeadingM_nopix_standard < massMax)) {
        return;
    }

    if (untrim) {
        leading_jet_nopix_standard = considered_jets_nopix_standard[0];
        fTLeadingEta_nopix_standard = leading_jet_nopix_standard.eta();
        fTLeadingPt_nopix_standard = leading_jet_nopix_standard.perp();
        fTLeadingM_nopix_standard = leading_jet_nopix_standard.m();

        // Cut condition for the untrimmed jets
        if (fTLeadingPt_nopix_standard < ptjMin) {
            return;
        }
    }

    // Repeat the above for charged, without cutoff.
    fTLeadingPhi_nopix_standard = leading_jet_nopix_standard.phi();

    fastjet::ClusterSequence csLargeR_nopix_charged(particlesForJets_nopixel_charged, *m_jet_def);
    vector<fastjet::PseudoJet> considered_jets_nopix_charged = fastjet::sorted_by_pt(
        csLargeR_nopix_charged.inclusive_jets(10.0));

    fastjet::PseudoJet leading_jet_nopix_charged;
    if (untrim) {
        leading_jet_nopix_charged = considered_jets_nopix_charged[0];
    } else {
        leading_jet_nopix_charged = trimmer(considered_jets_nopix_charged[0]);
    }

    fTLeadingEta_nopix_charged = leading_jet_nopix_charged.eta();
    fTLeadingPt_nopix_charged = leading_jet_nopix_charged.perp();
    fTLeadingM_nopix_charged = leading_jet_nopix_charged.m();
    fTLeadingPhi_nopix_charged = leading_jet_nopix_charged.phi();

    //Now, we extract the energy from the calorimeter for processing by fastjet
    for (int i = 1; i <= detector_standard->GetNbinsX(); i++)
    {
        for (int j = 1; j <= detector_standard->GetNbinsY(); j++)
        {
            if (detector_standard->GetBinContent(i, j) > 0)
            {
                double phi = detector_standard->GetYaxis()->GetBinCenter(j);
                double eta = detector_standard->GetXaxis()->GetBinCenter(i);
                double E = detector_standard->GetBinContent(i, j);
                fastjet::PseudoJet p(0., 0., 0., 0.);

                //We measure E (not pT)!  And treat 'clusters' as massless.
                p.reset_PtYPhiM(E/cosh(eta), eta, phi, 0.);
                particlesForJets_standard.push_back(p);
            }
        }
    }

    for (int i = 1; i <= detector_charged->GetNbinsX(); i++)
    {
        for (int j = 1; j <= detector_charged->GetNbinsY(); j++)
        {
            if (detector_charged->GetBinContent(i, j) > 0)
            {
                double phi = detector_charged->GetYaxis()->GetBinCenter(j);
                double eta = detector_charged->GetXaxis()->GetBinCenter(i);
                double E = detector_charged->GetBinContent(i, j);
                fastjet::PseudoJet p(0., 0., 0., 0.);

                //We measure E (not pT)!  And treat 'clusters' as massless.
                p.reset_PtYPhiM(E/cosh(eta), eta, phi, 0.); 
                particlesForJets_charged.push_back(p);
            }
        }
    }

    fastjet::ClusterSequence csLargeR_standard(particlesForJets_standard, *m_jet_def);
    fastjet::ClusterSequence csLargeR_charged(particlesForJets_charged, *m_jet_def);

    vector<fastjet::PseudoJet> considered_jets_standard = fastjet::sorted_by_pt(
								       csLargeR_standard.inclusive_jets(10.0));

    vector<fastjet::PseudoJet> considered_jets_charged = fastjet::sorted_by_pt(
								       csLargeR_charged.inclusive_jets(10.0));

    fastjet::PseudoJet leading_jet_standard;
    fastjet::PseudoJet leading_jet_charged;
    if (untrim) {
        leading_jet_standard = considered_jets_standard[0];
        leading_jet_charged = considered_jets_charged[0];
    } else {
        leading_jet_standard = trimmer(considered_jets_standard[0]);
        leading_jet_charged = trimmer(considered_jets_charged[0]);
    }

    // Standard, Pixelated
    vector<float> ec_standard = Corelators(considered_jets_standard, leading_jet_standard);
    fTLeadingEta_standard = leading_jet_standard.eta();
    fTLeadingM_standard = leading_jet_standard.m();
    fTLeadingPhi_standard = leading_jet_standard.phi();
    fTLeadingPt_standard = leading_jet_standard.perp();
    ec1_standard = ec_standard[0];
    cout << "ec1 standard: " << ec1_standard << endl;
    ec2_standard = ec_standard[1];
    ec3_standard = ec_standard[2];

    // Standard, Not Pixelated
    vector<float> ec_nopix_standard = Corelators(considered_jets_nopix_standard, leading_jet_nopix_standard);
    fTLeadingEta_nopix_standard = leading_jet_nopix_standard.eta();
    fTLeadingPhi_nopix_standard = leading_jet_nopix_standard.phi();
    fTLeadingPt_nopix_standard = leading_jet_nopix_standard.perp();
    fTLeadingM_nopix_standard = leading_jet_nopix_standard.m();
    ec1_nopix_standard = ec_nopix_standard[0];
    ec2_nopix_standard = ec_nopix_standard[1];
    ec3_nopix_standard = ec_nopix_standard[2];

    // Charged, Pixelated
    vector<float> ec_charged = Corelators(considered_jets_charged, leading_jet_charged);
    fTLeadingEta_charged = leading_jet_charged.eta();
    fTLeadingM_charged = leading_jet_charged.m();
    fTLeadingPhi_charged = leading_jet_charged.phi();
    fTLeadingPt_charged = leading_jet_charged.perp();
    ec1_charged = ec_charged[0];
    ec2_charged = ec_charged[1];
    ec3_charged = ec_charged[2];

    // Charged, Not Pixelated
    vector<float> ec_nopix_charged = Corelators(considered_jets_nopix_charged, leading_jet_nopix_charged);
    fTLeadingEta_nopix_charged = leading_jet_nopix_charged.eta();
    fTLeadingPhi_nopix_charged = leading_jet_nopix_charged.phi();
    fTLeadingPt_nopix_charged = leading_jet_nopix_charged.perp();
    fTLeadingM_nopix_charged = leading_jet_nopix_charged.m();
    ec1_nopix_charged = ec_nopix_charged[0];
    ec2_nopix_charged = ec_nopix_charged[1];
    ec3_nopix_charged = ec_nopix_charged[2];

    vector<fastjet::PseudoJet> subjets_nopix_standard = leading_jet_nopix_standard.pieces();
    vector<fastjet::PseudoJet> subjets_nopix_charged = leading_jet_nopix_charged.pieces();

    fTpull1_standard = -1;
    fTpull2_standard = -1;

    fTpull1_charged = -1;
    fTpull2_charged = -1;

    fTpull1_nopix_standard = -1;
    fTpull2_nopix_standard = -1;

    fTpull1_nopix_charged = -1;
    fTpull2_nopix_charged = -1;
    int e = 911;
    try {
        std::tuple<double, double> no_pix_pulls_standard = calculate_pull(subjets_nopix_standard, tool);
        fTpull1_nopix_standard = get<0>(no_pix_pulls_standard);
        fTpull2_nopix_standard = get<1>(no_pix_pulls_standard);

        std::tuple<double, double> no_pix_pulls_charged = calculate_pull(subjets_nopix_charged, tool);
        fTpull1_nopix_charged = get<0>(no_pix_pulls_charged);
        fTpull2_nopix_charged = get<1>(no_pix_pulls_charged);
    } catch(int e) {
        return;
    }

    fTdeltaR_standard = 0.;
    fTdeltaR_charged = 0.;

    if (leading_jet_standard.pieces().size() > 1){
      vector<fastjet::PseudoJet> subjets_standard = leading_jet_standard.pieces();
      TLorentzVector l(subjets_standard[0].px(),subjets_standard[0].py(),subjets_standard[0].pz(),subjets_standard[0].E());
      TLorentzVector sl(subjets_standard[1].px(),subjets_standard[1].py(),subjets_standard[1].pz(),subjets_standard[1].E());
      fTdeltaR_standard = l.DeltaR(sl);
      fTSubLeadingEta_standard = sl.Eta()-l.Eta();
      fTSubLeadingPhi_standard = subjets_standard[1].delta_phi_to(subjets_standard[0]);

      //Let's compute the jet pull.
      fTpull1_standard = tool->JetPull(subjets_standard[0],subjets_standard[1]);
      fTpull2_standard = tool->JetPull(subjets_standard[1],subjets_standard[0]);
    }

    if (leading_jet_charged.pieces().size() > 1){
      vector<fastjet::PseudoJet> subjets_charged = leading_jet_charged.pieces();
      TLorentzVector l(subjets_charged[0].px(),subjets_charged[0].py(),subjets_charged[0].pz(),subjets_charged[0].E());
      TLorentzVector sl(subjets_charged[1].px(),subjets_charged[1].py(),subjets_charged[1].pz(),subjets_charged[1].E());
      fTdeltaR_charged = l.DeltaR(sl);
      fTSubLeadingEta_charged = sl.Eta()-l.Eta();
      fTSubLeadingPhi_charged = subjets_charged[1].delta_phi_to(subjets_charged[0]);

      //Let's compute the jet pull.
      fTpull1_charged = tool->JetPull(subjets_charged[0],subjets_charged[1]);
      fTpull2_charged = tool->JetPull(subjets_charged[1],subjets_charged[0]);
    }

    vector<pair<double, double>  > consts_image_standard;
    vector<fastjet::PseudoJet> sorted_consts_standard = sorted_by_pt(leading_jet_standard.constituents());

    vector<pair<double, double>  > consts_image_charged;
    vector<fastjet::PseudoJet> sorted_consts_charged = sorted_by_pt(leading_jet_charged.constituents());

    for(int i = 0; i < sorted_consts_standard.size(); i++)
    {
        pair<double, double> const_hold;
        const_hold.first = sorted_consts_standard[i].eta();
        const_hold.second = sorted_consts_standard[i].phi();
        consts_image_standard.push_back(const_hold);
    }

    for(int i = 0; i < sorted_consts_charged.size(); i++)
    {
        pair<double, double> const_hold;
        const_hold.first = sorted_consts_charged[i].eta();
        const_hold.second = sorted_consts_charged[i].phi();
        consts_image_charged.push_back(const_hold);
    }


    vector<fastjet::PseudoJet> subjets_standard = leading_jet_standard.pieces();
    vector<fastjet::PseudoJet> subjets_charged = leading_jet_charged.pieces();

    //Step 1: Center on the jet axis.
    fastjet::PseudoJet shiftjet_charged;
    shiftjet_charged.reset_momentum_PtYPhiM(1.,detector_charged->GetXaxis()->GetBinCenter(detector_charged->GetXaxis()->FindBin(subjets_charged[0].eta())),detector_charged->GetYaxis()->GetBinCenter(detector_charged->GetYaxis()->FindBin(subjets_charged[0].phi())),0.);
    for (int i =0; i < sorted_consts_charged.size(); i++)
    {
      consts_image_charged[i].first = consts_image_charged[i].first-shiftjet_charged.eta();
      consts_image_charged[i].second = sorted_consts_charged[i].delta_phi_to(shiftjet_charged);
    }

    fastjet::PseudoJet shiftjet_standard;
    shiftjet_standard.reset_momentum_PtYPhiM(1.,detector_standard->GetXaxis()->GetBinCenter(detector_standard->GetXaxis()->FindBin(subjets_standard[0].eta())),detector_standard->GetYaxis()->GetBinCenter(detector_standard->GetYaxis()->FindBin(subjets_standard[0].phi())),0.);
    for (int i =0; i < sorted_consts_standard.size(); i++)
    {
      consts_image_standard[i].first = consts_image_standard[i].first-shiftjet_standard.eta();
      consts_image_standard[i].second = sorted_consts_standard[i].delta_phi_to(shiftjet_standard);
    }

    //Quickly run PCA for the rotation.
    double xbar_standard = 0.;
    double ybar_standard = 0.;
    double x2bar_standard = 0.;
    double y2bar_standard = 0.;
    double xybar_standard = 0.;
    double n_standard = 0;

    double xbar_charged = 0.;
    double ybar_charged = 0.;
    double x2bar_charged = 0.;
    double y2bar_charged = 0.;
    double xybar_charged = 0.;
    double n_charged = 0;

    for(int i = 0; i < leading_jet_standard.constituents().size(); i++)
    {
        double x = consts_image_standard[i].first;
        double y = consts_image_standard[i].second;
	    double E = sorted_consts_standard[i].e();
        n_standard+=E;
        xbar_standard+=x*E;
        ybar_standard+=y*E;
    }

    for(int i = 0; i < leading_jet_charged.constituents().size(); i++)
    {
        double x = consts_image_charged[i].first;
        double y = consts_image_charged[i].second;
	    double E = sorted_consts_charged[i].e();
        n_charged+=E;
        xbar_charged+=x*E;
        ybar_charged+=y*E;
    }

    double mux_standard = xbar_standard / n_standard;
    double muy_standard = ybar_standard / n_standard;

    double mux_charged = xbar_charged / n_charged;
    double muy_charged = ybar_charged / n_charged;


    xbar_standard = 0.;
    ybar_standard = 0.;
    n_standard = 0.;

    xbar_charged = 0.;
    ybar_charged = 0.;
    n_charged = 0.;

    for(int i = 0; i < leading_jet_standard.constituents().size(); i++)
    {
        double x = consts_image_standard[i].first - mux_standard;
        double y = consts_image_standard[i].second - muy_standard;
        double E = sorted_consts_standard[i].e();
        n_standard+=E;
        xbar_standard+=x*E;
        ybar_standard+=y*E;
        x2bar_standard+=x*x*E;
        y2bar_standard+=y*y*E;
        xybar_standard+=x*y*E;
    }

    for(int i = 0; i < leading_jet_charged.constituents().size(); i++)
    {
        double x = consts_image_charged[i].first - mux_charged;
        double y = consts_image_charged[i].second - muy_charged;
        double E = sorted_consts_charged[i].e();
        n_charged+=E;
        xbar_charged+=x*E;
        ybar_charged+=y*E;
        x2bar_charged+=x*x*E;
        y2bar_charged+=y*y*E;
        xybar_charged+=x*y*E;
    }

    double sigmax2_standard = x2bar_standard / n_standard - mux_standard*mux_standard;
    double sigmay2_standard = y2bar_standard / n_standard - muy_standard*muy_standard;
    double sigmaxy_standard = xybar_standard / n_standard - mux_standard*muy_standard;
    double lamb_min_standard = 0.5* ( sigmax2_standard + sigmay2_standard - sqrt( (sigmax2_standard-sigmay2_standard)*(sigmax2_standard-sigmay2_standard) + 4*sigmaxy_standard*sigmaxy_standard) );
    double lamb_max_standard = 0.5* ( sigmax2_standard + sigmay2_standard + sqrt( (sigmax2_standard-sigmay2_standard)*(sigmax2_standard-sigmay2_standard) + 4*sigmaxy_standard*sigmaxy_standard) );

    double sigmax2_charged = x2bar_charged / n_charged - mux_charged*mux_charged;
    double sigmay2_charged = y2bar_charged / n_charged - muy_charged*muy_charged;
    double sigmaxy_charged = xybar_charged / n_charged - mux_charged*muy_charged;
    double lamb_min_charged = 0.5* ( sigmax2_charged + sigmay2_charged - sqrt( (sigmax2_charged-sigmay2_charged)*(sigmax2_charged-sigmay2_charged) + 4*sigmaxy_charged*sigmaxy_charged) );
    double lamb_max_charged = 0.5* ( sigmax2_charged + sigmay2_charged + sqrt( (sigmax2_charged-sigmay2_charged)*(sigmax2_charged-sigmay2_charged) + 4*sigmaxy_charged*sigmaxy_charged) );


    double dir_x_standard = sigmax2_standard+sigmaxy_standard-lamb_min_standard;
    double dir_y_standard = sigmay2_standard+sigmaxy_standard-lamb_min_standard;

    double dir_x_charged = sigmax2_charged+sigmaxy_charged-lamb_min_charged;
    double dir_y_charged = sigmay2_charged+sigmaxy_charged-lamb_min_charged;

    //The first PC is only defined up to a sign.  Let's have it point toward the side of the jet with the most energy.

    double Eup_standard = 0.;
    double Edn_standard = 0.;

    double Eup_charged = 0.;
    double Edn_charged = 0.;

    for(int i = 0; i < leading_jet_standard.constituents().size(); i++)
    {
        double x = consts_image_standard[i].first - mux_standard;
        double y = consts_image_standard[i].second - muy_standard;
        double E = sorted_consts_standard[i].e();
        double dotprod = dir_x_standard*x+dir_y_standard*y;
        if (dotprod > 0) {
            Eup_standard+=E;
        } else {
            Edn_standard+=E;
        }
    }

    for(int i = 0; i < leading_jet_charged.constituents().size(); i++)
    {
	double x = consts_image_charged[i].first - mux_charged;
        double y = consts_image_charged[i].second - muy_charged;
        double E = sorted_consts_charged[i].e();
        double dotprod = dir_x_charged*x+dir_y_charged*y;
        if (dotprod > 0) {
            Eup_charged+=E;
        } else {
            Edn_charged+=E;
        }
    }

    if (Edn_standard < Eup_standard){
        dir_x_standard = -dir_x_standard;
        dir_y_standard = -dir_y_standard;
    }

    if (Edn_charged < Eup_charged){
        dir_x_charged = -dir_x_charged;
        dir_y_charged = -dir_y_charged;
    }

    fTPCEta_standard = dir_x_standard;
    fTPCPhi_standard = dir_y_standard;

    fTPCEta_charged = dir_x_charged;
    fTPCPhi_charged = dir_y_charged;

    //Step 2: Fill in the unrotated image
    //-------------------------------------------------------------------------
    range = 1.25;
    TH2D* orig_im_standard = new TH2D("", "", pixels, -range, range, pixels, -range, range);
    TH2D* orig_im_charged = new TH2D("", "", pixels, -range, range, pixels, -range, range);

    // TLorentzVector image_mass_standard = TLorentzVector();
    for (int i = 0; i < sorted_consts_standard.size(); i++)
    {
      TLorentzVector hold = TLorentzVector();
      hold.SetPtEtaPhiM(sorted_consts_standard[i].perp(),consts_image_standard[i].first,consts_image_standard[i].second,0.);
      orig_im_standard->Fill(consts_image_standard[i].first,consts_image_standard[i].second,hold.E());
    }

    // TLorentzVector image_mass_charged = TLorentzVector();
    for (int i = 0; i < sorted_consts_charged.size(); i++)
    {
      TLorentzVector hold = TLorentzVector();
      hold.SetPtEtaPhiM(sorted_consts_charged[i].perp(),consts_image_charged[i].first,consts_image_charged[i].second,0.);
      orig_im_charged->Fill(consts_image_charged[i].first,consts_image_charged[i].second,hold.E());
    }

    //Step 5: Dump the images in the tree!
    //-------------------------------------------------------------------------
    int counter=0;
    for (int i=1; i<=orig_im_standard->GetNbinsX(); i++)
    {
        for (int j=1; j<=orig_im_standard->GetNbinsY(); j++)
        {
            fTIntensity_pT_standard[counter] = orig_im_standard->GetBinContent(i,j);
	        double myeta = orig_im_standard->GetXaxis()->GetBinCenter(i);
	        fTIntensity_standard[counter] = orig_im_standard->GetBinContent(i,j)/cosh(myeta);

            counter++;
        }
    }
    counter=0;
    for (int i=1; i<=orig_im_charged->GetNbinsX(); i++)
    {
        for (int j=1; j<=orig_im_charged->GetNbinsY(); j++)
        {
            fTIntensity_pT_charged[counter] = orig_im_charged->GetBinContent(i,j);
	        double myeta = orig_im_charged->GetXaxis()->GetBinCenter(i);
	        fTIntensity_charged[counter] = orig_im_charged->GetBinContent(i,j)/cosh(myeta);

            counter++;
        }
    }

    // Step 6: Fill in nsubjettiness (new)
    //----------------------------------------------------------------------------
    OnePass_WTA_KT_Axes axis_spec_standard;
    NormalizedMeasure parameters_standard(1.0, 1.0);

    // NormalizedMeasure parameters(1.0, 1.0);
    Nsubjettiness subjettiness_1_standard(1, axis_spec_standard, parameters_standard);
    Nsubjettiness subjettiness_2_standard(2, axis_spec_standard, parameters_standard);
    Nsubjettiness subjettiness_3_standard(3, axis_spec_standard, parameters_standard);

    OnePass_WTA_KT_Axes axis_spec_charged;
    NormalizedMeasure parameters_charged(1.0, 1.0);

    // NormalizedMeasure parameters(1.0, 1.0);
    Nsubjettiness subjettiness_1_charged(1, axis_spec_charged, parameters_charged);
    Nsubjettiness subjettiness_2_charged(2, axis_spec_charged, parameters_charged);
    Nsubjettiness subjettiness_3_charged(3, axis_spec_charged, parameters_charged);


    fTTau1_standard = (float) subjettiness_1_standard.result(leading_jet_standard);
    fTTau2_standard = (float) subjettiness_2_standard.result(leading_jet_standard);
    fTTau3_standard = (float) subjettiness_3_standard.result(leading_jet_standard);

    fTTau1_charged = (float) subjettiness_1_charged.result(leading_jet_charged);
    fTTau2_charged = (float) subjettiness_2_charged.result(leading_jet_charged);
    fTTau3_charged = (float) subjettiness_3_charged.result(leading_jet_charged);

    fTTau32_standard = (abs(fTTau2_standard) < 1e-4 ? -10 : fTTau3_standard / fTTau2_standard);
    fTTau21_standard = (abs(fTTau1_standard) < 1e-4 ? -10 : fTTau2_standard / fTTau1_standard);

    fTTau32_charged = (abs(fTTau2_charged) < 1e-4 ? -10 : fTTau3_charged / fTTau2_charged);
    fTTau21_charged = (abs(fTTau1_charged) < 1e-4 ? -10 : fTTau2_charged / fTTau1_charged);

    fTTau1_nopix_standard = (float) subjettiness_1_standard.result(leading_jet_nopix_standard);
    fTTau2_nopix_standard = (float) subjettiness_2_standard.result(leading_jet_nopix_standard);
    fTTau3_nopix_standard = (float) subjettiness_3_standard.result(leading_jet_nopix_standard);

    fTTau1_nopix_charged = (float) subjettiness_1_charged.result(leading_jet_nopix_charged);
    fTTau2_nopix_charged = (float) subjettiness_2_charged.result(leading_jet_nopix_charged);
    fTTau3_nopix_charged = (float) subjettiness_3_charged.result(leading_jet_nopix_charged);

    fTTau32_nopix_standard = (abs(fTTau2_nopix_standard) < 1e-4 ? -10 : fTTau3_nopix_standard / fTTau2_nopix_standard);
    fTTau21_nopix_standard = (abs(fTTau1_nopix_standard) < 1e-4 ? -10 : fTTau2_nopix_standard / fTTau1_nopix_standard);

    fTTau32_nopix_charged = (abs(fTTau2_nopix_charged) < 1e-4 ? -10 : fTTau3_nopix_charged / fTTau2_nopix_charged);
    fTTau21_nopix_charged = (abs(fTTau1_nopix_charged) < 1e-4 ? -10 : fTTau2_nopix_charged / fTTau1_nopix_charged);

    tT_standard->Fill();

    tT_charged->Fill();

    return;
}

// declare branches
void myexampleAnalysis::DeclareBranches()
{
    // Event Standard Properties
    tT_standard->Branch("NFilled", &fTNFilled_standard, "NFilled/I");

    tT_standard->Branch("Intensity", *&fTIntensity_standard, "Intensity[NFilled]/F");
    tT_standard->Branch("Intensity_pT", *&fTIntensity_pT_standard, "Intensity_pT[NFilled]/F");

    tT_standard->Branch("SubLeadingEta", &fTSubLeadingEta_standard, "SubLeadingEta/F");
    tT_standard->Branch("SubLeadingPhi", &fTSubLeadingPhi_standard, "SubLeadingPhi/F");

    tT_standard->Branch("PCEta", &fTPCEta_standard, "PCEta/F");
    tT_standard->Branch("PCPhi", &fTPCPhi_standard, "PCPhi/F");

    tT_standard->Branch("LeadingEta", &fTLeadingEta_standard, "LeadingEta/F");
    tT_standard->Branch("LeadingPhi", &fTLeadingPhi_standard, "LeadingPhi/F");
    tT_standard->Branch("LeadingPt", &fTLeadingPt_standard, "LeadingPt/F");
    tT_standard->Branch("LeadingM", &fTLeadingM_standard, "LeadingM/F");

    tT_standard->Branch("LeadingEta_nopix", &fTLeadingEta_nopix_standard, "LeadingEta_nopix/F");
    tT_standard->Branch("LeadingPhi_nopix", &fTLeadingPhi_nopix_standard, "LeadingPhi_nopix/F");
    tT_standard->Branch("LeadingPt_nopix", &fTLeadingPt_nopix_standard, "LeadingPt_nopix/F");
    tT_standard->Branch("LeadingM_nopix", &fTLeadingM_nopix_standard, "LeadingM_nopix/F");

    tT_standard->Branch("Tau1", &fTTau1_standard, "Tau1/F");
    tT_standard->Branch("Tau2", &fTTau2_standard, "Tau2/F");
    tT_standard->Branch("Tau3", &fTTau3_standard, "Tau3/F");

    tT_standard->Branch("Tau1_nopix", &fTTau1_nopix_standard, "Tau1_nopix/F");
    tT_standard->Branch("Tau2_nopix", &fTTau2_nopix_standard, "Tau2_nopix/F");
    tT_standard->Branch("Tau3_nopix", &fTTau3_nopix_standard, "Tau3_nopix/F");

    tT_standard->Branch("DeltaR", &fTdeltaR_standard, "DeltaR/F");

    tT_standard->Branch("Tau32", &fTTau32_standard, "Tau32/F");
    tT_standard->Branch("Tau21", &fTTau21_standard, "Tau21/F");

    tT_standard->Branch("Tau32_nopix", &fTTau32_nopix_standard, "Tau32_nopix/F");
    tT_standard->Branch("Tau21_nopix", &fTTau21_nopix_standard, "Tau21_nopix/F");

    tT_standard->Branch("pull1", &fTpull1_standard, "pull1/F");
    tT_standard->Branch("pull2", &fTpull2_standard, "pull2/F");

    tT_standard->Branch("pull1_nopix", &fTpull1_nopix_standard, "pull1_nopix/F");
    tT_standard->Branch("pull2_nopix", &fTpull2_nopix_standard, "pull2_nopix/F");

    tT_standard->Branch("ec_1", &ec1_standard, "ec_1/F");
    tT_standard->Branch("ec_2", &ec2_standard, "ec_2/F");
    tT_standard->Branch("ec_3", &ec3_standard, "ec_3/F");

    tT_standard->Branch("ec_1_nopix", &ec1_nopix_standard, "ec_1_nopix/F");
    tT_standard->Branch("ec_2_nopix", &ec2_nopix_standard, "ec_2_nopix/F");
    tT_standard->Branch("ec_3_nopix", &ec3_nopix_standard, "ec_3_nopix/F");

    // Event Charged Properties
    tT_charged->Branch("NFilled", &fTNFilled_charged, "NFilled/I");

    tT_charged->Branch("Intensity", *&fTIntensity_charged, "Intensity[NFilled]/F");
    tT_charged->Branch("Intensity_pT", *&fTIntensity_pT_charged, "Intensity_pT[NFilled]/F");

    tT_charged->Branch("SubLeadingEta", &fTSubLeadingEta_charged, "SubLeadingEta/F");
    tT_charged->Branch("SubLeadingPhi", &fTSubLeadingPhi_charged, "SubLeadingPhi/F");

    tT_charged->Branch("PCEta", &fTPCEta_charged, "PCEta/F");
    tT_charged->Branch("PCPhi", &fTPCPhi_charged, "PCPhi/F");

    tT_charged->Branch("LeadingEta", &fTLeadingEta_charged, "LeadingEta/F");
    tT_charged->Branch("LeadingPhi", &fTLeadingPhi_charged, "LeadingPhi/F");
    tT_charged->Branch("LeadingPt", &fTLeadingPt_charged, "LeadingPt/F");
    tT_charged->Branch("LeadingM", &fTLeadingM_charged, "LeadingM/F");

    tT_charged->Branch("LeadingEta_nopix", &fTLeadingEta_nopix_charged, "LeadingEta_nopix/F");
    tT_charged->Branch("LeadingPhi_nopix", &fTLeadingPhi_nopix_charged, "LeadingPhi_nopix/F");
    tT_charged->Branch("LeadingPt_nopix", &fTLeadingPt_nopix_charged, "LeadingPt_nopix/F");
    tT_charged->Branch("LeadingM_nopix", &fTLeadingM_nopix_charged, "LeadingM_nopix/F");

    tT_charged->Branch("Tau1", &fTTau1_charged, "Tau1/F");
    tT_charged->Branch("Tau2", &fTTau2_charged, "Tau2/F");
    tT_charged->Branch("Tau3", &fTTau3_charged, "Tau3/F");

    tT_charged->Branch("Tau1_nopix", &fTTau1_nopix_charged, "Tau1_nopix/F");
    tT_charged->Branch("Tau2_nopix", &fTTau2_nopix_charged, "Tau2_nopix/F");
    tT_charged->Branch("Tau3_nopix", &fTTau3_nopix_charged, "Tau3_nopix/F");

    tT_charged->Branch("DeltaR", &fTdeltaR_charged, "DeltaR/F");

    tT_charged->Branch("Tau32", &fTTau32_charged, "Tau32/F");
    tT_charged->Branch("Tau21", &fTTau21_charged, "Tau21/F");

    tT_charged->Branch("Tau32_nopix", &fTTau32_nopix_charged, "Tau32_nopix/F");
    tT_charged->Branch("Tau21_nopix", &fTTau21_nopix_charged, "Tau21_nopix/F");

    tT_charged->Branch("pull1", &fTpull1_charged, "pull1/F");
    tT_charged->Branch("pull2", &fTpull2_charged, "pull2/F");

    tT_charged->Branch("pull1_nopix", &fTpull1_nopix_charged, "pull1_nopix/F");
    tT_charged->Branch("pull2_nopix", &fTpull2_nopix_charged, "pull2_nopix/F");

    tT_charged->Branch("ec_1", &ec1_charged, "ec_1/F");
    tT_charged->Branch("ec_2", &ec2_charged, "ec_2/F");
    tT_charged->Branch("ec_3", &ec3_charged, "ec_3/F");

    tT_charged->Branch("ec_1_nopix", &ec1_nopix_charged, "ec_1_nopix/F");
    tT_charged->Branch("ec_2_nopix", &ec2_nopix_charged, "ec_2_nopix/F");
    tT_charged->Branch("ec_3_nopix", &ec3_nopix_charged, "ec_3_nopix/F");

    return;
}

// resets vars
void myexampleAnalysis::ResetBranches(){
    // reset branches

    // Standard Jets
    fTNFilled_standard = MaxN;
    fTSubLeadingPhi_standard = -999;
    fTSubLeadingEta_standard = -999;
    fTPCPhi_standard = -999;
    fTPCEta_standard = -999;

    fTTau32_standard = -999;
    fTTau21_standard = -999;

    fTTau1_standard = -999;
    fTTau2_standard = -999;
    fTTau3_standard = -999;

    fTTau32_nopix_standard = -999;
    fTTau21_nopix_standard = -999;

    fTTau1_nopix_standard = -999;
    fTTau2_nopix_standard = -999;
    fTTau3_nopix_standard = -999;

    fTLeadingEta_standard = -999;
    fTLeadingPhi_standard = -999;
    fTLeadingPt_standard = -999;
    fTLeadingM_standard = -999;

    fTLeadingEta_nopix_standard = -999;
    fTLeadingPhi_nopix_standard = -999;
    fTLeadingPt_nopix_standard = -999;
    fTLeadingM_nopix_standard = -999;

    ec1_standard = -999;
    ec2_standard = -999;
    ec3_standard = -999;

    ec1_nopix_standard = -999;
    ec2_nopix_standard = -999;
    ec3_nopix_standard = -999;

    for (int iP=0; iP < MaxN; ++iP)
    {
        fTIntensity_standard[iP]= -999;
	fTIntensity_pT_standard[iP]= -999;
    }

    // Charged Jets
    fTNFilled_charged = MaxN;
    fTSubLeadingPhi_charged = -999;
    fTSubLeadingEta_charged = -999;
    fTPCPhi_charged = -999;
    fTPCEta_charged = -999;

    fTTau32_charged = -999;
    fTTau21_charged = -999;

    fTTau1_charged = -999;
    fTTau2_charged = -999;
    fTTau3_charged = -999;

    fTTau32_nopix_charged = -999;
    fTTau21_nopix_charged = -999;

    fTTau1_nopix_charged = -999;
    fTTau2_nopix_charged = -999;
    fTTau3_nopix_charged = -999;

    fTLeadingEta_charged = -999;
    fTLeadingPhi_charged = -999;
    fTLeadingPt_charged = -999;
    fTLeadingM_charged = -999;

    fTLeadingEta_nopix_charged = -999;
    fTLeadingPhi_nopix_charged = -999;
    fTLeadingPt_nopix_charged = -999;
    fTLeadingM_nopix_charged = -999;

    ec1_charged = -999;
    ec2_charged = -999;
    ec3_charged = -999;

    ec1_nopix_charged = -999;
    ec2_nopix_charged = -999;
    ec3_nopix_charged = -999;

    for (int iP=0; iP < MaxN; ++iP)
    {
        fTIntensity_charged[iP]= -999;
        fTIntensity_pT_charged[iP]= -999;
    }
}

// Get Corelator vars
vector<float> myexampleAnalysis::Corelators(const vector<PseudoJet> & input_particles,  PseudoJet & resonance) {
    Mat3d MCorels = myexampleAnalysis::Ecorel(input_particles, resonance);

    vector<float> result;
    result.push_back(MCorels[0][5][1]);
    result.push_back(MCorels[0][5][2]);
    result.push_back(MCorels[0][5][3]);
    return result;
}

bool myexampleAnalysis::isjetc( const Pythia8::Particle* p ) {
    if(p->isFinal() && p->pT() > 0.5 && fabs(p->eta()) < 4.0 )
    {
        int pid = p->id();

        int lepArray[5] ={11,13};
        int neuArray[5] ={12,14,16};
        bool isNeu = std::any_of(
            std::begin(neuArray),
            std::end(neuArray),
            [&](int i) {
                return i == abs(pid);
            }
        );
        bool isLep = std::any_of(
            std::begin(lepArray),
            std::end(lepArray),
            [&](int i) {
                return i == abs(pid);
            }
        );
        if((isLep && p->pT() > 5.) || isNeu) return 0;
        return 1;
    }
    return 0;
}

// Kirtimaan's code
//********************************************************************//
//========For calculatng energy correlators
//********************************************************************//
Mat3d myexampleAnalysis::Ecorel( const vector<PseudoJet> & input_particles,  PseudoJet & resonance) {
    JetAlgorithm algorithm = cambridge_algorithm;
    double jet_rad = 1.0;
    JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
    ClusterSequence clust_seq(input_particles, jetDef);
    vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());

    // cout << "Resonance px: " << resonance.px();
    // cout << " py: " << resonance.py();
    // cout << " pz: " << resonance.pz();
    // cout << " e: " << resonance.e();
    // cout << " size of constituents: " << sizeof(resonance.constituents());
    // cout << "\n";
    // cout << "Size of input particles: " << sizeof(input_particles) << "\n";
    // cout << "input resonance: " << input_particles;
    //====== EnergyCorrelator ====================//
    //======= various values of beta ==============//
    vector<double> betalist;
    betalist.push_back(0.1);
    betalist.push_back(0.2);
    betalist.push_back(0.5);
    betalist.push_back(1.0);
    betalist.push_back(1.5);
    betalist.push_back(2.0);

    //==== checking the two energy/angle modes=======//
    vector<EnergyCorrelator::Measure> measurelist;
    measurelist.push_back(EnergyCorrelator::pt_R);
    measurelist.push_back(EnergyCorrelator::E_theta);
    //measurelist.push_back(EnergyCorrelator::E_inv);

    //========Store correlators==================//
    Mat3d mat_corels(measurelist.size(), Matrix(betalist.size(), Row(5)));

    //====Decalre objects on which correlator is run
    PseudoJet myJet, myJeti;

	double dR=9999.0;
	for (int j = 0; j < antikt_jets.size(); j++) {
        double dR2= resonance.delta_R(antikt_jets[j]);
		if( dR2 < dR) {
    		myJeti=antikt_jets[j]; dR= dR2;
        }
    }

    fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm,0.3), fastjet::SelectorPtFractionMin(0.00));
    PseudoJet trimmed = trimmer(myJeti);
    myJet=trimmed;

    cout << "My Jet px: " << myJet.px();
    cout << " py: " << myJet.py();
    cout << " pz: " << myJet.pz();
    cout << " e: " << myJet.e();
    cout << " size of constituents: " << sizeof(myJet.constituents());
    cout << "\n";

    vector<string> modename;
    modename.push_back("pt_R");
    modename.push_back("E_theta");

    for (unsigned int M = 0; M < measurelist.size(); M++) {
        for (unsigned int B = 0; B < betalist.size(); B++) {
            double beta = betalist[B];

            EnergyCorrelatorDoubleRatio C1(1,beta,measurelist[M]);
            EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M]);
            EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M]);

            mat_corels[M][B][0]=beta;
            mat_corels[M][B][1]=C1(myJet);
            // cout << "C1: " << C1(myJet) << endl;
            // cout << "C2: " << C2(myJet) << endl;
            // cout << "C3: " << C3(myJet) << endl;
            mat_corels[M][B][2]=C2(myJet);
            mat_corels[M][B][3]=C3(myJet);
            mat_corels[M][B][4]=0;
        }
    }

	return mat_corels;
}
