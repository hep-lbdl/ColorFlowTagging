#!/usr/bin/env python
'''
file: jetconverter.py
author: Luke de Oliveira, Aug 2015 

This file takes files (*.root) produced by the event-gen
portion of the jet-simulations codebase and converts them into 
a more usable format.

This also dumps everything (optionally) into a ROOT file.

'''

from argparse import ArgumentParser
import logging
import sys
sys.path.append("../jettools")

import numpy as np
import h5py

from jettools import plot_mean_jet
from processing import buffer_to_jet, is_signal

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def perfectsquare(n):
    '''
    I hope this is self explanatory...
    '''
    return n % n ** 0.5 == 0


if __name__ == '__main__':
    tree = None
    parser = ArgumentParser()
    parser.add_argument('--verbose',
                        action='store_true',
                        help='Verbose output')

    parser.add_argument('--signal',
                        default='wprime',
                        help='String to search for in\
                         filenames to indicate a signal file')

    parser.add_argument('--dump',
                        default=None,
                        help='ROOT file *prefix* to dump all this into (writes to TTree `images`). the .root '
                             'extension will be added.')
    parser.add_argument('--hdf5',
			default=None,
			help='Filename *prefix* to write out the data. (a .h5 ext will be added to the end)')
    parser.add_argument('--save',
                        default=None,
                        help='Filename *prefix* to write out the data. (a .npy ext will be added to the end)')
    parser.add_argument('--plot',
                        help='File prefix that\
                         will be part of plotting filenames.')
    parser.add_argument('--pixelSize', default=25, type=int, help='size of one side of the image, in pixels')

    parser.add_argument('--chunk', default=10, type=int, help='number of files to chunk together')

    parser.add_argument('--pTmin', default=300.0, type=float, help='Lower bound of cut on pT of leading jet')

    parser.add_argument('--pTmax', default=600.0, type=float, help='Upper bound of cut on pT of leading jet')

    parser.add_argument('--etaMax', default=2.0, type=float, help='Upper bound of cut on eta of leading jet')

    parser.add_argument('--massMin', default=100.0, type=float, help='Lower bound of cut on mass of leading jet')

    parser.add_argument('--massMax', default=150.0, type=float, help='Upper bound of cut on mass of leading jet')

    parser.add_argument('--applyCuts', default=True, type=str, help='Whether to apply all the cuts or not')

    parser.add_argument('files', nargs='*', help='Files to pass in')

    
    args = parser.parse_args()

    # -- check for logic errors
    if len(args.files) < 1:
        logger.error('Must pass at least one file in -- terminating with error.')
        exit(1)

    if (args.save is None) and (args.dump is None) and (args.hdf5 is None):
        logger.error('Must write to NPY, HDF5 and/or ROOT file.')
        exit(1)

    signal_match = args.signal
    files = args.files
    savefile = args.save
    hdf5 = args.hdf5
    plt_prefix = ''

    ptj_min = args.pTmin
    ptj_max = args.pTmax
    eta_max = args.etaMax
    mass_min = args.massMin
    mass_max = args.massMax
    if args.applyCuts in ['1', 'True', 'true', 'y', 'Y', 'T', 't']:
        apply_cuts = True
    else:
        apply_cuts = False
    
    if args.plot:
        plt_prefix = args.plot

    try:  # -- see if this rootpy business works
        from rootpy.io import root_open
        from rootpy.tree import Tree, TreeModel, FloatCol, FloatArrayCol
    except ImportError:
        raise ImportError('rootpy (www.rootpy.org) not installed\
         -- install, then try again!')


    # -- create buffer for the tree
    class JetImage(TreeModel):
        '''
        Buffer for Jet Image
        '''
        # -- START BUFFER

        # -- raveled image
        image = FloatArrayCol(args.pixelSize ** 2)

        # -- 1 if signal, 0 otherwise
        signal = FloatCol()

        # -- kinematics
        jet_pt = FloatCol()
        jet_eta = FloatCol()
        jet_phi = FloatCol()
        jet_m = FloatCol()
        jet_delta_R = FloatCol()

        # -- NSJ
        tau_32 = FloatCol()
        tau_21 = FloatCol()
        tau_1 = FloatCol()
        tau_2 = FloatCol()
        tau_3 = FloatCol()
        pull1 = FloatCol()
        pull2 = FloatCol()
        # -- END BUFFER


    pix_per_side = -999
    entries = []
    img_entries = []

    CHUNK_MAX = int(args.chunk)
    print('chunk max is {}'.format(CHUNK_MAX))
    N_CHUNKED = 0
    CURRENT_CHUNK = 0

    import glob

    # -- hack
    # files = args.files # --
    # files = glob.glob(args.files[0] + '/*.root')
    files = glob.glob(args.files[0] + '/*.root')

    print('length: {}'.format(len(files)), args.files[0] + '/*.root')

    ROOTfile = None

    for i, fname in enumerate(files):
        logger.info('ON: CHUNK #{}'.format(CURRENT_CHUNK))
        try:
            if args.dump and ROOTfile is None:
                logger.info('Making ROOT file: {}'.format(args.dump + '-chnk{}.root'.format(CURRENT_CHUNK)))
                ROOTfile = root_open(args.dump + '-chnk{}.root'.format(CURRENT_CHUNK), "recreate")
                tree = Tree('images', model=JetImage)
        except Exception:
            continue

        logger.info('({} of {}) working on file: {}'.format(i, len(files), fname))
        try:
            with root_open(fname) as f:
                df = f.EventTree.to_array()

                n_entries = df.shape[0]

                pix = df[0]['Intensity'].shape[0]

                if not perfectsquare(pix):
                    raise ValueError('shape of image array must be square.')

                if (pix_per_side > 1) and (int(np.sqrt(pix)) != pix_per_side):
                    raise ValueError('all files must have same sized images.')

                pix_per_side = int(np.sqrt(pix))
                logger.info('Pixel Size: {}'.format(pix_per_side))

                tag = is_signal(fname, signal_match)
                logger.info('Logging as {}'.format(tag))
                print(not apply_cuts)
                for jet_nb, jet in enumerate(df):
                    if jet_nb % 1000 == 0:
                        logger.info('processing jet {} of {} for file {}'.format(
                            jet_nb, n_entries, fname
                        )
                        )
                    if ((np.abs(jet['LeadingEta']) < eta_max) & (jet['LeadingPt'] > ptj_min) & (
                           jet['LeadingPt'] < ptj_max) & (jet['LeadingM'] < mass_max) & (
                           jet['LeadingM'] > mass_min)) or not apply_cuts:
			if not (np.abs(jet['LeadingEta']) < eta_max) & (jet['LeadingPt'] > ptj_min) & (jet['LeadingPt'] < ptj_max) & (jet['LeadingM'] < mass_max) & (jet['LeadingM'] > mass_min):

                        buf = buffer_to_jet(jet, args.pixelSize, tag, max_entry=100000)
                        if args.dump:
                            tree.image = buf[0].ravel()  # .astype('float32')
                            tree.signal = buf[1]
                            tree.jet_pt = buf[2]
                            tree.jet_eta = buf[3]
                            tree.jet_phi = buf[4]
                            tree.jet_m = buf[5]
                            tree.jet_delta_R = buf[6]
                            tree.tau_32 = buf[7]
                            tree.tau_21 = buf[8]
                            tree.tau_1 = buf[9]
                            tree.tau_2 = buf[10]
                            tree.tau_3 = buf[11]
                            tree.pull1 = buf[12]
                            tree.pull2 = buf[13]
                        if savefile is not None:
                            entries.append(buf)
                        if hdf5 is not None:
                            img_entries.append(buf[0])
                            entries.append(buf[1:])
                        if args.dump:
                            tree.fill()

            # -- Check for chunking
            N_CHUNKED += 1
            # -- we've reached the max chunk size
            if N_CHUNKED >= CHUNK_MAX:
                logger.info('{} files chunked, max is {}'.format(N_CHUNKED, CHUNK_MAX))
                N_CHUNKED = 0
                # -- clear the env, and reset
                if args.dump:
                    tree.write()
                    ROOTfile.close()
                    ROOTfile = None
                    tree = None
                CURRENT_CHUNK += 1
        except KeyboardInterrupt:
            logger.info('Skipping file {}'.format(fname))
        except AttributeError:
            logger.info('Skipping file {} for compatibility reasons'.format(fname))
    if args.dump and tree is not None:
        tree.write()
        ROOTfile.close()

    if hdf5 is not None:

        with h5py.File(hdf5 + '.h5', mode='w') as hf:
            entries = np.array(entries)
            t = hf.create_group('meta_variables')
            _buf_names = ['signal', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_mass', 'jet_delta_R', 'tau_32', 'tau_21', 
                          'tau_1', 'tau_2', 'tau_3', 'pull1', 'pull2']
            for idx, meta in enumerate(_buf_names):
                t.create_dataset(meta, data=entries[:, idx])
            del entries
            hf.create_dataset('images', data=img_entries)
            del img_entries

    if savefile is not None:
        # -- datatypes for outputted file.
        _bufdtype = [('image', 'float32', (pix_per_side, pix_per_side)),
                     ('signal', 'float32'),
                     ('jet_pt', 'float32'),
                     ('jet_eta', 'float32'),
                     ('jet_phi', 'float32'),
                     ('jet_mass', 'float32'),
                     ('jet_delta_R', 'float32'),
                     ('tau_32', 'float32'),
                     ('tau_21', 'float32'),
                     ('tau_1', 'float32'),
                     ('tau_2', 'float32'),
                     ('tau_3', 'float32'),
                     ('pull1', 'float32'),
                     ('pull2', 'float32'), ]
        # print(entries)
        # print(_bufdtype)
        df = np.array(entries, dtype=_bufdtype)
        logger.info('saving to file: {}'.format(savefile))
        np.save(savefile, df)

        if plt_prefix != '':
            logger.info('plotting...')
            plot_mean_jet(df[df['signal'] == 0], title="Average Jet Image, Background").savefig(plt_prefix + '_bkg.pdf')
            plot_mean_jet(df[df['signal'] == 1], title="Average Jet Image, Signal").savefig(plt_prefix + '_signal.pdf')
