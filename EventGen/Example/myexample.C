#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"
#include "TSystem.h"
#include "TError.h"
//#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "myTools.h"
#include "myexampleAnalysis.h"

#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;
namespace po = boost::program_options;

int getSeed(int seed){                                                                                                                                               
  if (seed > -1) return seed;                                                                                                                                      \
  int timeSeed = time(NULL);                                                                                                                                       \
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);                                                                                                         \
}

int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments 
    int nEvents = 0;
    int fDebug  = 0;
    string outName = "test.root";
    string inName = "test_input.root";

    int pixels = 25;
    double image_range = 1.25;
    int pileup =0; //number of extra collisions that happen on top of the main one.  For now, set this to 0.
    bool isCharged = false;
    float pTmin, pTmax, etamax, massmin, massmax;

    po::options_description desc("Allowed options");
    desc.add_options() //100000
      ("help", "produce help message")
      ("NEvents", po::value<int>(&nEvents)->default_value(10000) ,    "Number of Events ")
      ("NPixels", po::value<int>(&pixels)->default_value(25) ,    "Number of Pixels ")
      ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile", po::value<string>(&outName)->default_value("test.root"), "output file name")
      ("InFile", po::value<string>(&inName)->default_value("test_input.root"), "input file name")
      ("IsCharged", po::value<bool>(&isCharged)->default_value(false), "Take only charged jets")
      ("pTMin", po::value<float>(&pTmin)->default_value(300), "Upper bound of cut on pT of leading jet")
      ("pTMax", po::value<float>(&pTmax)->default_value(600), "Lower bound of cut on pT of leading jet")
      ("etaMax", po::value<float>(&etamax)->default_value(2), "Upper bound of cut on eta of leading jet")
      ("massMin", po::value<float>(&massmin)->default_value(100), "Lower bound of cut on mass of leading jet")
      ("massMax", po::value<float>(&massmax)->default_value(150), "Upper bound of cut on mass of leading jet")
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }
    //------

   Pythia8::Pythia* pythia8b = new Pythia8::Pythia();

   int    seed      =-1;
   seed = getSeed(seed);

   // Configure and initialize pythia                                                                                                                                
   pythia8b->readString("Random:setSeed = on");
   std::stringstream ss; ss << "Random:seed = " << seed;
   cout << ss.str() << endl;
   pythia8b->readString(ss.str());

   myexampleAnalysis * analysis1 = new myexampleAnalysis(pixels);

   pythia8b->readString("Beams:frameType = 4");
   std::stringstream in; in << "Beams:LHEF = lhe/" << inName << ".lhe";
   pythia8b->readString(in.str());
   pythia8b->init();

   //Setup the pileup
   Pythia8::Pythia* pythia_MB = new Pythia8::Pythia();
   pythia_MB->readString("Random:setSeed = on");
   ss.clear(); ss.str(""); ss << "Random:seed = " << seed+1;
   cout << ss.str() << endl;
   pythia_MB->readString(ss.str());
   pythia_MB->readString("SoftQCD:nonDiffractive = on");
   pythia_MB->readString("HardQCD:all = off");
   pythia_MB->readString("PhaseSpace:pTHatMin  = .1");
   pythia_MB->readString("PhaseSpace:pTHatMax  = 20000");
   pythia_MB->init();

   analysis1->SetOutName(outName);
   analysis1->Begin();
   cout << "running on " << nEvents << " events " << endl;
   for (Int_t iev = 0; iev < nEvents; iev++) {
     if (iev%100==0) cout << iev << " " << nEvents << endl;
     analysis1->AnalyzeEvent(iev, pythia8b, pythia_MB, pileup, pixels, image_range, pTmin, pTmax, etamax, massmin, massmax, isCharged);
   }
   analysis1->End();
   pythia8b->stat();
   delete analysis1;

   // that was it
   delete pythia8b;
   return 0;
}
