'''
processing.py
author: Luke de Oliveira, July 2015 

Simple utilities for processing the junk that comes out of the ntuple event generation.
'''

import numpy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
from jettools import rotate_jet, flip_jet, plot_mean_jet
from ROOT import*

def sum_image(image):
    out = 0.
    for i in range(len(image)):
        for j in range(len(image)):
            out+=image[i][j]
            pass
        pass
    return out

def ben_rotate(image,angle,pixels):
    arange = 1.25;
    hold = TH2F("", "", pixels, -arange, arange, pixels, -arange, arange);
    out_image=[]
    vecs=[]
    for i2 in range(1,hold.GetNbinsX()+1):
        for j2 in range(1,hold.GetNbinsY()+1):
            etaval = hold.GetXaxis().GetBinCenter(i2);
            phival = hold.GetYaxis().GetBinCenter(j2);
            ImageMass_hold = TLorentzVector();
            ImageMass_hold.SetPtEtaPhiM(image[i2-1][j2-1],etaval*np.cos(angle)-phival*np.sin(angle),phival*np.cos(angle)+etaval*np.sin(angle),0.)
            hold.Fill(ImageMass_hold.Eta(),ImageMass_hold.Phi(),ImageMass_hold.Pt())
            pass
        pass
    for i2 in range(1,hold.GetNbinsX()+1):
        out_image+=[[]]
        for j2 in range(1,hold.GetNbinsY()+1):
            out_image[i2-1]+=[hold.GetBinContent(i2,j2)]
            pass
        pass
    return out_image

def image_mass(image,pixels):
    arange = 1.25;
    hold = TH2F("", "", pixels, -arange, arange, pixels, -arange, arange);
    counter = 0;
    ImageMass = TLorentzVector(0.,0.,0.,0.); 
    for i2 in range(1,hold.GetNbinsX()+1):
        for j2 in range(1,hold.GetNbinsY()+1):
            etaval = hold.GetXaxis().GetBinCenter(i2);
            phival = hold.GetYaxis().GetBinCenter(j2);
            ImageMass_hold = TLorentzVector();
            ImageMass_hold.SetPtEtaPhiM(image[i2-1][j2-1],etaval,phival,0.);   
            ImageMass+=ImageMass_hold;
            counter+=1
            pass
        pass
    return ImageMass.M()*4000

def angle_from_vec(v1, v2):
    # cosang = np.dot(v1, v2)
    cosang = v1
    sinang = la.norm(np.cross(np.array([1.0, 0.0]), np.array([v1, v2])))
    return np.arctan2(sinang, cosang)


def buffer_to_jet(entry, pix, tag=0, side='r', max_entry=None, rotate=True, normalize=True, pixelated=True):
    """
    Takes an *single* entry from an structured ndarray, i.e., X[i],
    and a tag = {0, 1} indicating if its a signal entry or not.
    The parameter 'side' indicates which side of the final
    jet image we want the highest energy.

    The `entry` must have the following fields (as produced by event-gen)
        * Intensity
        * PCEta, PCPhi
        * LeadingPt
        * LeadingEta
        * LeadingPhi
        * SubLeadingEta
        * SubLeadingPhi
        * LeadingM
        * DeltaR
        * Tau32
        * Tau21
        * Tau{n} for n = 1, 2, 3
        * pull1
        * pull2
        * ec_{n} for n = 1, 2, 3
    """

    accessed = False
    if (entry['SubLeadingEta'] < -10) | (entry['SubLeadingPhi'] < -10):
        accessed = True
        e, p = (entry['PCEta'], entry['PCPhi'])
    else:
        e, p = (entry['SubLeadingEta'], entry['SubLeadingPhi'])


    if rotate:

        angle = np.arctan2(p, e) + 2.0 * np.arctan(1.0)
        if (-np.sin(angle) * e + np.cos(angle) * p) > 0:
            angle += -4.0 * np.arctan(1.0)

        image = flip_jet(rotate_jet(np.array(entry['Intensity']), -angle, pix, normalizer=4000.0), side) # change between (-angle <-> -4*np.arctan(1.0)) if testing
    else:
        image = flip_jet(rotate_jet(np.array(entry['Intensity']), 0.0, pix, normalizer=4000.0), side)

    if normalize:
        e_norm = np.linalg.norm(image)
    	img = (image / e_norm).astype('float32')
    else:
        img = (image).astype('float32')

    if pixelated:
        return (
            img,
            np.float32(tag),
            np.float32(entry['LeadingPt']), np.float32(entry['LeadingEta']),
            np.float32(entry['LeadingPhi']), np.float32(entry['LeadingM']),
            np.float32(entry['DeltaR']),
            np.float32(entry['Tau32']), np.float32(entry['Tau21']),
            np.float32(entry['Tau1']), np.float32(entry['Tau2']), np.float32(entry['Tau3']),
            np.float32(entry['pull1']), np.float32(entry['pull2']),
            np.float32(entry['ec_1']), np.float32(entry['ec_2']), np.float32(entry['ec_3'])
        )
    else:
        return (
            img,
            np.float32(tag),
            np.float32(entry['LeadingPt_nopix']), np.float32(entry['LeadingEta_nopix']),
            np.float32(entry['LeadingPhi_nopix']), np.float32(entry['LeadingM_nopix']),
            np.float32(entry['DeltaR_nopix']),
            np.float32(entry['Tau32_nopix']), np.float32(entry['Tau21_nopix']),
            np.float32(entry['Tau1_nopix']), np.float32(entry['Tau2_nopix']), np.float32(entry['Tau3_nopix']),
            np.float32(entry['pull1_nopix']), np.float32(entry['pull2_nopix']),
            np.float32(entry['ec_1_nopix']), np.float32(entry['ec_2_nopix']), np.float32(entry['ec_3_nopix'])
        )

def is_signal(f, matcher = 'wprime'):
    """
    Takes as input a filename and a string to match. If the 
    'matcher' string is found in the filename, the file is 
    taken to be a signal file.
    """
    key = matcher.lower().replace(' ', '').replace('-', '')
    if key in f.lower().replace(' ', '').replace('-', ''):
        return 1.0
    return 0.0
