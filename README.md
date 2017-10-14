# ColorFlowTagging

Event Generation

(1) The first step is to generate .lhe events with the Higgs boson decay products.  This is then fed to Pythia in (2).  There is a script in EventGen called HiggsGun.C.  This generates both the signal ("singlet") and background ("octet") configurations.  The only parameters in the code you may want to change are MotherHiggsMass and nevents.  The former sets the energy scale of the Higgses. The rule of thumb is that the momentum of the Higgs will be ~mMother/2 and the angular separation between the Higgs daughter quarks will be ~2*(Higgs mass) / (Higgs momentum), where Higgs mass = 125 GeV (proton mass ~ 1 GeV).

(2) Run Pythia using the above LHE events as input. You can do this with commands like

        pythia8->readString("Beams:frameType = 4");
        pythia8->readString("Beams:LHEF = lhe/eventsSinglet.lhe");
        pythia8->init();

The EventGen code directory has my personal setup - it depends on ROOT and fastjet.