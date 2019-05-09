# Preprocess
This directory, as the name suggests, is intended to preprocess images,
which then can be fed into the neural network.
The main file if jetconverter_withpull, which will output a root TTrees, numpy arrays or hdf5 files,
which then can be fed into respective framework.

## Documentation
* --verbose: Verbose output
* --signal: String to search for in filenames to indicate a signal file
* --dump: ROOT file *prefix* to dump all this into (writes to TTree images). The *.root* extension will be added.
* --hdf5: Filename *prefix* to write out the data. The .h5 extension will be added to the end
* --save: Filename *prefix* to write out the data. The .npy ext will be added to the end
* --plot: File prefix that will be part of plotting file names. If plot not chosen, no plots will be generated.
* --pixelSize: Size of one side of the image, in pixels,
* --chunk: Number of files to chunk together.
* --pTmin: Lower bound of cut on pT of leading jet
* --pTmax: Upper bound of cut on pT of leading jet
* --etaMax: Upper bound of cut on eta of leading jet
* --massMin: Lower bound of cut on mass of leading jet
* --massMax: Upper bound of cut on mass of leading jet
* --applyCuts: Whether to apply all the cuts or not
* --rotate: Whether to rotate the image or not
* --norm: Whether to normalize the image or not

## Example
```
python jetconverter_withpull.py --hdf5 h_gg --rotate True --norm True --applyCuts False --pixelSize 65 inputs/h_gg/
```