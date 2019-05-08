# Example
This is an example of using the pipeline, with even generator.
It's also a program used in paper.

## Requirements
* [Pythia8](http://home.thep.lu.se/~torbjorn/Pythia.html) - Version 8.2 assumed here.
  If you are using different version you might need to modify --tunePP and --colorMode flags.
* [Boost](https://www.boost.org/)
* [FastJet](http://fastjet.fr/)
* [FastJet Contrib](https://fastjet.hepforge.org/contrib/)

An example of working setup is:
* Pythia 8.235
* Boost 1.68.0
* FastJet 3.3.2
* FastJet Contrib 1.039

## Setup
If you haven't modified setup.sh yet, make sure it's pointing to:
* PYTHIA8LOCATION - points to root folder of pythia8 installation
* path in setup_ROOT points to ROOT executable (on linux it's *thisroot.sh* file)
* FASTJETLOCATION - points to root folder of fastjet installation
* BOOSTINCDIR - points to include folder of boost installation
* BOOSTLIBLOCATION - points to lib folder of boost installation

## Docs
Setting for *myexample.exe* are:
* --help                     : produce help message
* --NEvents: Number of Events. Default is 10000.
* --NPixels: Number of Pixels. Default is 25.
* --Debug: Debug flag. Default is 0. Unless you are modifying code do not change.
* --OutFile: output file name (.root will be added automatically). Default is test.
* --InFile: input file name (.root will be added automatically). Default is test_input.
* --pTMin: Upper bound of cut on pT of leading jet. Default is 300.
* --pTMax: Lower bound of cut on pT of leading jet. Default is 600.
* --etaMax: Upper bound of cut on absolute value of eta of leading jet. Default is 2.
* --massMin: Lower bound of cut on mass GeV of leading jet. Default is 100.
* --massMax: Upper bound of cut on mass GeV of leading jet. Default is 150.
* --colorMode: Type of model used for colour reconnection by pythia8. Default is 0.
* --tunePP: Choice of tune to pp/ppbar data used by pythia8. Default is 14.

## Example usage
To run always start with:

```
source setup.sh
make
./myexample.exe --NEvents 600000 --NPixels 65 --colorMode 2 --InFile h_gg/h_gg --OutFile h_gg_col_2.root
```

## Modification/Advanced Usage
In order to modify behaviour of shower/modify settings change *myexample.C*  
To make changes to the analysis and the output file, change *myexampleAnalysis.cc*  
Allowed modification are based on version of pythia you are using