# ColorFlowTagging

## Physics

Searching for highly Lorentz-boosted Higgs boson is a very important and hot topic.
The biggest background for most of these searches is gluon splitting to bottom quark pairs, as the Higgs mostly also decays to bottom quark pairs.
There are many tools for two-prong tagging, but in this case, the g->bb background has two prong and two b-quark jets.
So the question is: how well can we do with just colorflow information?
The state-of-the-art colorflow tagging observable is jet pull (Phys.Rev.Lett.105:022001,2010).
Can we do (much) better than that?  What is there to learn?

## Event Generation

1. The first step is to generate .lhe events with the Higgs boson decay products. 

    This is then fed to Pythia in (2).  There is a script in EventGen called HiggsGun.C.
    This generates both the signal ("singlet") and background ("octet") configurations.
    The only parameters in the code you may want to change are MotherHiggsMass and nevents.
    The former sets the energy scale of the Higgses.
    The rule of thumb is that the momentum of the Higgs will be ~mMother/2 and 
    the angular separation between the Higgs daughter quarks will be 
    ~2*(Higgs mass) / (Higgs momentum), where Higgs mass = 125 GeV (proton mass ~ 1 GeV).
    
    N.B. the rest of the steps are based on the same ones used in this [github repo (written largely by Luke de Oliveira)](https://github.com/ml-slac/jet-simulations).

2. Run Pythia using the above LHE events as input. You can do this with commands like
    ```
    pythia8->readString("Beams:frameType = 4");
    pythia8->readString("Beams:LHEF = lhe/eventsSinglet.lhe");
    pythia8->init();
    ```
    The EventGen code directory has my personal setup - it depends on ROOT and fastjet.
    If you have these two and change the setup file,
    it *should work* but there may be some small tinkering required.

3. Next, is pre-processing.
    We may want to not do all of the steps here, but we can discuss.
    The script is jetconverter_withpull.py in the PreProcessing directory.
    Anyway, you should be able to do something like:
    
    ```
    ./jetconverter_withpull.py --save save_dir --dump Octet_withpull input_root_file --chunk 1
    ```
4. This step is optional - probably better to have things setup so that the previous step gives you the numpy array or hdf5 file you need.
However, as a carry-over from the past, there is a script to dump a text file.
It is also in the PreProcessing dir.

## Machine Learning

The above setup is image-based so fully connected networks and CNN are quite natural.
It would also be interesting to try RNNs on the entire list of particles.
This would require some tweaking of the output, but not much.
This feature is being worked on to be added in near future.
