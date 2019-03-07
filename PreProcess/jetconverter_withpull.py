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

def change_angles(x):
    return (x > np.pi) * (2 * np.pi - x) + x * (x <= np.pi)


if __name__ == '__main__':
    tree_charged = None
    tree_standard = None
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
                        help='File prefix that will be part of plotting file names.')
    parser.add_argument('--pixelSize', default=65, type=int, help='size of one side of the image, in pixels')

    parser.add_argument('--chunk', default=10, type=int, help='number of files to chunk together')

    parser.add_argument('--pTmin', default=300.0, type=float, help='Lower bound of cut on pT of leading jet')

    parser.add_argument('--pTmax', default=600.0, type=float, help='Upper bound of cut on pT of leading jet')

    parser.add_argument('--etaMax', default=2.0, type=float, help='Upper bound of cut on eta of leading jet')

    parser.add_argument('--massMin', default=100.0, type=float, help='Lower bound of cut on mass of leading jet')

    parser.add_argument('--massMax', default=150.0, type=float, help='Upper bound of cut on mass of leading jet')

    parser.add_argument('--applyCuts', default=False, type=str, help='Whether to apply all the cuts or not')

    parser.add_argument('--rotate', default=True, type=str, help='Whether to rotate the image or not')

    parser.add_argument('--norm', default=True, type=str, help='Whether to normalize the image or not')

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

    truth = ['1', 'True', 'true', 'y', 'Y', 'T', 't']

    if args.applyCuts in truth:
        apply_cuts = True
    else:
        apply_cuts = False
    
    if args.rotate in truth:
        rotate = True
    else:
        rotate = False

    if args.norm in truth:
        normalize = True
    else:
        normalize = False

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
        """
        Buffer for Jet Image
        """
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

        # standard vs charged
        s_vs_c_dR = FloatCol()
        # -- END BUFFER


    pix_per_side = -999
    entries_standard = []
    entries_charged = []
    img_entries_standard = []
    img_entries_charged = []

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
                tree_standard = Tree('images', model=JetImage)
                tree_charged = Tree('images', model=JetImage)
        except Exception:
            continue

        logger.info('({} of {}) working on file: {}'.format(i, len(files), fname))
        try:
            with root_open(fname) as f:
                df_s = f.StandardEventTree.to_array()
                df_c = f.ChargedEventTree.to_array()

            assert df_s.shape[0] == df_c.shape[0]

            n_entries = df_s.shape[0]

            pix_s = df_s[0]['Intensity'].shape[0]
            pix_c = df_c[0]['Intensity'].shape[0]

            assert pix_s == pix_c

            if not perfectsquare(pix_s):
                raise ValueError('shape of image array must be square.')

            if (pix_per_side > 1) and (int(np.sqrt(pix_s)) != pix_per_side):
                raise ValueError('all files must have same sized images.')

            pix_per_side = int(np.sqrt(pix_s))
            logger.info('Pixel Size: {}'.format(pix_per_side))

            tag = is_signal(fname, signal_match)
            logger.info('Logging as {}'.format(tag))

            for jet_nb, jet_s in enumerate(df_s):
                jet_c = df_c[jet_nb]
                if jet_nb % 1000 == 0:
                    logger.info('processing jet {} of {} for file {}'.format(
                        jet_nb, n_entries, fname
                    )
                    )
                if (((np.abs(jet_s['LeadingEta']) < eta_max) & (jet_s['LeadingPt'] > ptj_min) & (
                       jet_s['LeadingPt'] < ptj_max) & (jet_s['LeadingM'] < mass_max) & (
                       jet_s['LeadingM'] > mass_min)) and
                    ((np.abs(jet_c['LeadingEta']) < eta_max) & (jet_c['LeadingPt'] > ptj_min) & (
                       jet_c['LeadingPt'] < ptj_max) & (jet_c['LeadingM'] < mass_max) & (
                       jet_c['LeadingM'] > mass_min))) or not apply_cuts:

                    buf_c = list(buffer_to_jet(jet_c, args.pixelSize, tag, max_entry=100000, rotate=rotate, normalize=normalize))
                    buf_s = list(buffer_to_jet(jet_s, args.pixelSize, tag, max_entry=100000, rotate=rotate, normalize=normalize))
                    if buf_s is None or buf_c is None:
                        continue
                    buf_s[4] = change_angles(buf_s[4])
                    buf_c[4] = change_angles(buf_c[4])
                    dR = np.hypot(buf_s[3] - buf_c[3], buf_s[4] - buf_s[4])
                    if args.dump:
                        tree_standard.image = buf_s[0].ravel()  # .astype('float32')
                        tree_standard.signal = buf_s[1]
                        tree_standard.jet_pt = buf_s[2]
                        tree_standard.jet_eta = buf_s[3]
                        tree_standard.jet_phi = buf_s[4]
                        tree_standard.jet_m = buf_s[5]
                        tree_standard.jet_delta_R = buf_s[6]
                        tree_standard.tau_32 = buf_s[7]
                        tree_standard.tau_21 = buf_s[8]
                        tree_standard.tau_1 = buf_s[9]
                        tree_standard.tau_2 = buf_s[10]
                        tree_standard.tau_3 = buf_s[11]
                        tree_standard.pull1 = buf_s[12]
                        tree_standard.pull2 = buf_s[13]
                        tree_standard.s_vs_c_dR = dR

                        tree_charged.image = buf_c[0].ravel()  # .astype('float32')
                        tree_charged.signal = buf_c[1]
                        tree_charged.jet_pt = buf_c[2]
                        tree_charged.jet_eta = buf_c[3]
                        tree_charged.jet_phi = buf_c[4]
                        tree_charged.jet_m = buf_c[5]
                        tree_charged.jet_delta_R = buf_c[6]
                        tree_charged.tau_32 = buf_c[7]
                        tree_charged.tau_21 = buf_c[8]
                        tree_charged.tau_1 = buf_c[9]
                        tree_charged.tau_2 = buf_c[10]
                        tree_charged.tau_3 = buf_c[11]
                        tree_charged.pull1 = buf_c[12]
                        tree_charged.pull2 = buf_c[13]
                        tree_charged.s_vs_c_dR = dR

                    if savefile is not None:
                        entries_standard.append(buf_s + [dR])
                        entries_charged.append(buf_c + [dR])
                    if hdf5 is not None:
                        img_entries_standard.append(buf_s[0])
                        entries_standard.append(buf_s[1:] + [dR])
                        img_entries_charged.append(buf_c[0])
                        entries_charged.append(buf_c[1:] + [dR])
                    if args.dump:
                        tree_standard.fill()
                        tree_charged.fill()

            # -- Check for chunking
            N_CHUNKED += 1
            # -- we've reached the max chunk size
            if N_CHUNKED >= CHUNK_MAX:
                logger.info('{} files chunked, max is {}'.format(N_CHUNKED, CHUNK_MAX))
                N_CHUNKED = 0
                # -- clear the env, and reset
                if args.dump:
                    tree_standard.write()
                    tree_charged.write()
                    ROOTfile.close()
                    ROOTfile = None
                    tree_standard = None
                    tree_charged = None
	    CURRENT_CHUNK += 1
        except KeyboardInterrupt:
            logger.info('Skipping file {}'.format(fname))
        except AttributeError:
            logger.info('Skipping file {} for compatibility reasons'.format(fname))
    if args.dump and tree_standard is not None and tree_charged is not None:
        tree_charged.write()
        tree_standard.write()
        ROOTfile.close()

    if hdf5 is not None:
        with h5py.File(hdf5 + '_standard.h5', mode='w') as hf:
            entries = np.array(entries_standard)
            del entries_standard

            t = hf.create_group('meta_variables')
            _buf_names = ['signal', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_mass', 'jet_delta_R', 'tau_32', 'tau_21', 
                          'tau_1', 'tau_2', 'tau_3', 'pull1', 'pull2', 'dR_standard_vs_charged']
            for idx, meta in enumerate(_buf_names):
                t.create_dataset(meta, data=entries[:, idx])
            del entries
            hf.create_dataset('images', data=img_entries_standard)
            del img_entries_standard

        with h5py.File(hdf5 + '_charged.h5', mode='w') as hf:
            entries = np.array(entries_charged)
            t = hf.create_group('meta_variables')
            _buf_names = ['signal', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_mass', 'jet_delta_R', 'tau_32', 'tau_21',
                          'tau_1', 'tau_2', 'tau_3', 'pull1', 'pull2', 'dR_standard_vs_charged']
            for idx, meta in enumerate(_buf_names):
                t.create_dataset(meta, data=entries[:, idx])
            del entries
            del entries_charged
            hf.create_dataset('images', data=img_entries_charged)
            del img_entries_charged

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
                     ('pull2', 'float32'),
                     ('dR_standard_vs_charged', 'float32')]

        df_s = np.array(entries_standard, dtype=_bufdtype)
        df_c = np.array(entries_standard, dtype=_bufdtype)
        logger.info('saving to file: {}'.format(savefile))
        np.save(savefile + '_standard', df_s)
        np.save(savefile + '_charged', df_s)

        if plt_prefix != '':
            logger.info('plotting...')
            plot_mean_jet(df_s[df_s['signal'] == 0], title="Average Standard Jet Image, Background").\
                savefig(plt_prefix + '_bkg.pdf')
            plot_mean_jet(df_c[df_c['signal'] == 1], title="Average Charged Jet Image, Signal").\
                savefig(plt_prefix + '_signal.pdf')
