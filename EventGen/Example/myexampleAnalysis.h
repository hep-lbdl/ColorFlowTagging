#ifndef  myexampleAnalysis_H
#define  myexampleAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "myTools.h"
#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

#include "TH2F.h"

using namespace std;
using namespace fastjet;

typedef vector<double> Row; // One row of the matrix
typedef vector<Row> Matrix; // Matrix: a vector of rows
typedef vector<Matrix> Mat3d; // Matrix3D: a vector of Matrices

class myexampleAnalysis
{
    public:
        myexampleAnalysis(int imagesize);
        ~myexampleAnalysis();

        void Begin();
        void AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB,
                int NPV, int pixels, float range, float ptjMin, float ptjMax, float etaMax, float massMin, float massMax,
                bool cambridge, bool reproduce);



        void End();
        void DeclareBranches();
        void ResetBranches();

        void Debug(int debug)
        {
            fDebug = debug;
        }


        void SetOutName(const string &outname)
        {
            fOutName = outname;
        }

        vector<float> Corelators(const vector<PseudoJet> & input_particles, PseudoJet & resonance);
    private:
        int  ftest;
        int  fDebug;
        string fOutName;

        TFile *tF;
        TTree *tT_standard;
        TTree *tT_charged;

        myTools *tool;

        // Tree Vars ---------------------------------------
        int fTEventNumber;
        int fTNPV;

        void SetupInt(int & val, TString name);
        void SetupFloat(float & val, TString name);

        bool isjetc(const Pythia8::Particle* p);

        vector<TString> names;
        vector<float> pts;
        vector<float> ms;
        vector<float> etas;
        vector<float> nsub21s;
        vector<float> nsub32s;
        vector<int>   nsubs;

        TH2D* detector_standard;
        TH2D* detector_charged;

        int MaxN;

        // Standard jet parameters
        int fTNFilled_standard;

        float fTLeadingEta_standard;
        float fTLeadingPhi_standard;
        float fTLeadingPt_standard;
        float fTLeadingM_standard;

	float fTLeadingEta_nopix_standard;
        float fTLeadingPhi_nopix_standard;
        float fTLeadingPt_nopix_standard;
        float fTLeadingM_nopix_standard;

	float fTpull1_standard;
	float fTpull2_standard;
	float fTpull1_nopix_standard;
        float fTpull2_nopix_standard;

        float fTSubLeadingEta_standard;
        float fTSubLeadingPhi_standard;

	float fTPCEta_standard;
	float fTPCPhi_standard;

        float fTTau1_standard;
        float fTTau2_standard;
        float fTTau3_standard;

        float fTTau21_standard;
        float fTTau32_standard;

	float fTTau1_nopix_standard;
        float fTTau2_nopix_standard;
        float fTTau3_nopix_standard;

	float fTTau21_nopix_standard;
	float fTTau32_nopix_standard;

        float ec1_standard;
        float ec2_standard;
        float ec3_standard;

        float ec1_nopix_standard;
        float ec2_nopix_standard;
        float ec3_nopix_standard;

        float fTdeltaR_standard;

        float *fTIntensity_standard;
	float *fTIntensity_pT_standard;

        // Charged jet parameters
        int fTNFilled_charged;

        float fTLeadingEta_charged;
        float fTLeadingPhi_charged;
        float fTLeadingPt_charged;
        float fTLeadingM_charged;

	float fTLeadingEta_nopix_charged;
        float fTLeadingPhi_nopix_charged;
        float fTLeadingPt_nopix_charged;
        float fTLeadingM_nopix_charged;

	float fTpull1_charged;
	float fTpull2_charged;
	float fTpull1_nopix_charged;
        float fTpull2_nopix_charged;

        float fTSubLeadingEta_charged;
        float fTSubLeadingPhi_charged;

	float fTPCEta_charged;
	float fTPCPhi_charged;

        float fTTau1_charged;
        float fTTau2_charged;
        float fTTau3_charged;

        float fTTau21_charged;
        float fTTau32_charged;

	float fTTau1_nopix_charged;
        float fTTau2_nopix_charged;
        float fTTau3_nopix_charged;

	float fTTau21_nopix_charged;
	float fTTau32_nopix_charged;

        float ec1_charged;
        float ec2_charged;
        float ec3_charged;

        float ec1_nopix_charged;
        float ec2_nopix_charged;
        float ec3_nopix_charged;

        float fTdeltaR_charged;

        float *fTIntensity_charged;
	float *fTIntensity_pT_charged;

        Mat3d Ecorel(const vector<PseudoJet> & input_particles, PseudoJet & resonance);
};

#endif
