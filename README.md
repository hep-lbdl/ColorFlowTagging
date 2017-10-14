# ColorFlowTagging

Event Generation

(1) The first step is to generate .lhe events with the Higgs boson decay products.  This is then fed to Pythia in (2).  There is a script in EventGen called HiggsGun.C.  This generates both the signal ("singlet") and background ("octet") configurations.  The only parameters in the code you may want to change are MotherHiggsMass and nevents.  The former sets the energy scale of the Higgses. The rule of thumb is that the momentum of the Higgs will be ~mMother/2 and the angular separation between the Higgs daughter quarks will be ~2*(Higgs mass) / (Higgs momentum), where Higgs mass = 125 GeV (proton mass ~ 1 GeV).

N.B. the rest of the steps are based on the same ones used in https://github.com/ml-slac/jet-simulations (written largely by Luke de Oliveira).

(2) Run Pythia using the above LHE events as input. You can do this with commands like

        pythia8->readString("Beams:frameType = 4");
        pythia8->readString("Beams:LHEF = lhe/eventsSinglet.lhe");
        pythia8->init();

The EventGen code directory has my personal setup - it depends on ROOT and fastjet.

(3) Next, is pre-processing.  We may want to not do all of the steps here, but we can discuss.  The script is jetconverter_withpull.py in the PreProcessing directory.  Anyway, you should be able to do something like

./jetconverter_withpull.py --save save_dir --dump Octet_withpull input_root_file --chunk 1

(4) This step is optional - probably better to have things setup so that the previous step gives you the numpy array or hdf5 file you need.  However, as a carry-over from the past, there is a script to dump a text file.  It is also in the PreProcessing dir.