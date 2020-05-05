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
import glob

import numpy as np
import h5py

from tqdm import tqdm

try:  # -- see if this rootpy business works
    from rootpy.io import root_open
    from rootpy.tree import Tree, TreeModel, FloatCol, FloatArrayCol
except ImportError:
    raise ImportError('rootpy (www.rootpy.org) not installed\
        -- install, then try again!')

sys.path.append("../jettools")
from jettools import plot_mean_jet
from processing import buffer_to_jet, is_signal


def perfectsquare(n):
    return n % n ** 0.5 == 0


def change_angles(x):
    return (x > np.pi) * (2 * np.pi - x) + x * (x <= np.pi)


def parse_from_command_line():
    """
    Loads arguments from command line.
    """
    parser = ArgumentParser()
    parser.add_argument('--verbose',
                        action='store_true',
                        help='Verbose output')

    parser.add_argument('--hdf5',
                        default=None,
                        help='Filename *prefix* to write out the data. (a .h5 ext will be added to the end)')

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

    parser.add_argument('--pixelated', default=True, type=str, help='Whether to use pixaleted variables or not')

    parser.add_argument('--numSamples', default=None, type=int, help='Number of samples to generate')

    parser.add_argument('files', nargs='*', help='Files to pass in')

    args = parser.parse_args()

    # -- check for logic errors
    if len(args.files) < 1:
        raise ValueError('Must pass at least one file in -- terminating with error.')

    files = args.files
    hdf5 = args.hdf5

    ptj_min = args.pTmin
    ptj_max = args.pTmax
    eta_max = args.etaMax
    mass_min = args.massMin
    mass_max = args.massMax
    pixelated = args.pixelated

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

    if hdf5 is None:
        raise ValueError('--hdf5 has to be filled with file name which will be filled with name to which save the images')

    return files, hdf5, ptj_min, ptj_max, eta_max, mass_min, mass_max, pixelated, apply_cuts, rotate, normalize, args.chunk, args.pixelSize, args.numSamples


def process_root_file(file_path, key, logger, chunk, apply_cuts, ptj_min, ptj_max, eta_max, mass_min, mass_max, rotate, normalize, pixelated, pixelSize, num_samples):
    pix_per_side = -1
    with root_open(file_path, 'r') as f:
        tag = is_signal(file_path)

        tree = getattr(f, key)

        all_images = []
        all_meta_vars = []

        for event_idx, event in enumerate(tqdm(tree, leave=False)):

            if num_samples is not None and len(all_images) >= num_samples:
                # This ensures, that only so many samples will be generated
                break

            if ((np.abs(event['LeadingEta']) < eta_max) & (event['LeadingPt'] > ptj_min) & (
                    event['LeadingPt'] < ptj_max) & (event['LeadingM'] < mass_max) & (
                    event['LeadingM'] > mass_min)) or not apply_cuts:

                processed_jet = list(buffer_to_jet(event, pixelSize, tag, max_entry=100000, rotate=rotate, normalize=normalize, pixelated=pixelated))
                processed_jet[4] = change_angles(processed_jet[4])

                # Extract processed jet
                image = processed_jet[0]
                meta_vars = map(float, processed_jet[2:])

                all_images.append(image)
                all_meta_vars.append(meta_vars)

    return np.array(all_images), np.array(all_meta_vars)


def main(files, hdf5, ptj_min, ptj_max, eta_max, mass_min, mass_max, pixelated, apply_cuts, rotate, normalize, chunk, pixelSize, num_samples):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    tree_charged = None
    tree_standard = None

    entries_standard = []
    entries_charged = []
    img_entries_standard = []
    img_entries_charged = []

    CHUNK_MAX = int(chunk)
    print('chunk max is {}'.format(CHUNK_MAX))
    N_CHUNKED = 0
    CURRENT_CHUNK = 0

    files = glob.glob(files[0] + '/*.root')

    print('length: {}'.format(len(files)), files[0] + '/*.root')

    ROOTfile = None
    _buf_names = [
        'jet_pt',
        'jet_eta', 'jet_phi',
        'jet_mass',
        'jet_delta_R',
        'tau_32', 'tau_21',
        'tau_1', 'tau_2', 'tau_3',
        'pull1', 'pull2',
        'ec_1', 'ec_2', 'ec_3'
    ]

    use_idx = len(files) > 1

    for i, fname in enumerate(files):
        logger.info('({} of {}) working on file: {}'.format(i, len(files), fname))

        try:
            standard_images, standard_meta_vars = process_root_file(fname, 'StandardEventTree', logger, chunk, apply_cuts, ptj_min, ptj_max, eta_max, mass_min, mass_max, rotate, normalize, pixelated, pixelSize, num_samples)

            if use_idx:
                standard_save_name = '{}_{}_standard.h5'.format(hdf5, i + 1)
                charged_save_name = '{}_{}_charged.h5'.format(hdf5, i + 1)
            else:
                standard_save_name = '{}_standard.h5'.format(hdf5)
                charged_save_name = '{}_charged.h5'.format(hdf5)

            print('Saving Standard File')
            with h5py.File(standard_save_name, mode='w') as hf:
                t = hf.create_group('meta_variables')
                for idx, meta in enumerate(_buf_names):
                    t.create_dataset(meta, data=standard_meta_vars[:, idx])
                hf.create_dataset('images', data=standard_images)
            print('Finished Saving Standard File')

            del standard_images, standard_meta_vars

            charged_images, charged_meta_vars = process_root_file(fname, 'ChargedEventTree', logger, chunk, apply_cuts, ptj_min, ptj_max, eta_max, mass_min, mass_max, rotate, normalize, pixelated, pixelSize, num_samples)

            print('Saving Charged File')
            with h5py.File(charged_save_name, mode='w') as hf:
                t = hf.create_group('meta_variables')
                for idx, meta in enumerate(_buf_names):
                    t.create_dataset(meta, data=charged_meta_vars[:, idx])
                hf.create_dataset('images', data=charged_images)
            print('Finished Saving Charged File')

            del charged_images, charged_meta_vars

        except KeyboardInterrupt:
            logger.info('Skipping file {}'.format(fname))
        except AttributeError:
            logger.info('Skipping file {} for compatibility reasons'.format(fname))

if __name__ == "__main__":
    args = parse_from_command_line()
    main(*args)
