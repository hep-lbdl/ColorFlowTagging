'''
processing.py
author: Luke de Oliveira, July 2015 

Simple utilities for processing the junk that comes out of the ntuple event generation.
'''

import numpy.linalg as la
import numpy as np
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

#def reset_image(image):
#    arange = 1.25;
#    pixels      = 25;
#    hold = TH2F("", "", pixels, -arange, arange, pixels, -arange, arange);
#    counter = 0;
#    eta0 = hold.GetYaxis().FindBin(0)
#    phi0 = hold.GetXaxis().FindBin(1)
#    phi1 = hold.GetXaxis().FindBin(-1)
#    for i2 in range(1,hold.GetNbinsX()+1):
#        for j2 in range(1,hold.GetNbinsY()+1):
#            image[counter]=0
#            if (j2==eta0 and i2==phi0):
#                image[counter]=1
#                pass
#            if (j2==eta0 and i2==phi1):
#                image[counter]=1
#                pass
#            counter+=1
#            pass
#        pass
#    return image

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


def buffer_to_jet(entry, pix, tag = 0, side = 'r', max_entry = None):
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
    """


    if (entry['SubLeadingEta'] < -10) | (entry['SubLeadingPhi'] < -10):
        e, p = (entry['PCEta'], entry['PCPhi'])
    else:
        e, p = (entry['SubLeadingEta'], entry['SubLeadingPhi'])
    
    angle = np.arctan2(p, e) + 2.0 * np.arctan(1.0)

    if (-np.sin(angle) * e + np.cos(angle) * p) > 0:
        angle += -4.0 * np.arctan(1.0)
    
    '''    
    print "angleeee:,",angle
    for kk in range(-1000,1000):
        angle = 4*np.arctan(1.)*kk/1000
        #image_map = reset_image(np.array(entry['Intensity']))
        image_map = np.array(entry['Intensity'])
        image = flip_jet(rotate_jet(image_map, 0., pix, normalizer=4000.0), side)
        before = image_mass(image, pix)
        #print sum_image(image)
        image = flip_jet(rotate_jet(image_map, -angle, pix, normalizer=4000.0), side)
        #print sum_image(image)
        after = image_mass(image,pix)
        image = flip_jet(rotate_jet(image_map, 0., pix, normalizer=4000.0), side)
        image = ben_rotate(image,angle,pix)
        #print sum_image(image)
        after2 = image_mass(image,pix)
        print "angle, before, after, after(Ben), diff, diff",angle,before,after,after2,after-before,after-after2
        pass

    exit(1)
    '''

    image = flip_jet(rotate_jet(np.array(entry['Intensity']), -4*np.arctan(1.0), pix, normalizer=4000.0), side) # change between (-angle <-> -4*np.arctan(1.0)) if testing
    e_norm = np.linalg.norm(image)
    #e_norm = 1./4000.0
    return ((image / e_norm).astype('float32'), np.float32(tag), 
        np.float32(entry['LeadingPt']), np.float32(entry['LeadingEta']), 
        np.float32(entry['LeadingPhi']), np.float32(entry['LeadingM']), np.float32(entry['DeltaR']),
        np.float32(entry['Tau32']), np.float32(entry['Tau21']), np.float32(entry['Tau1']), np.float32(entry['Tau2']), np.float32(entry['Tau3']), np.float32(entry['pull1']), np.float32(entry['pull2']))


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
