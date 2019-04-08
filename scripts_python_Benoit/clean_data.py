#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Benoit Epinat

modified by Wilfried Mercier 5th April 2019
"""

import os, glob
import numpy as np
import astropy.io.ascii as ascii
import astropy.io.fits as fits
#import pyfits as fits
#import asciitable as ascii
#import copy
import logging
#import ipdb
from astropy import constants as ct

from sys import exit

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('Clean')

'''
This module is to be used to clean kinematics map created using camel.
'''


def create_mask(image, thrl=None, thru=None):
    '''
    This function creates a mask from one image using a lower and an upper threshold (True everywhere thrl<=image<=thru)
    
    Parameters
    ----------
    image: numpy array
        image used to create the mask
    thrl: float
        lower threshold
    thru: float
        upper threshold
        
    Returns a boolean mask
    '''
    if thrl is None:
        thrl = np.nanmin(image)
    if thru is None:
        thru = np.nanmax(image)
    logger.info('create_mask: lower threshold % s' % (str(thrl)) )
    logger.info('create_mask: upper threshold % s' % (str(thru)) )
    return ((image <= thru) & (image >= thrl))


def apply_mask(mask, image):
    '''
    This function applies a mask to an image and puts nan in the masked area
    
    Parameters
    ----------
    mask: numpy array
        array containing the index of the pixels to be masked
    image: numpy array
        image to be masked
        
    Returns the masked image
    '''
    image[mask] = np.nan
    return image


def clean_galaxy(path, outputpath, name, lsfw, fraction, data_mask='snr', thrl=None, thru=None, option='', line='', clean=None):
    '''
    This function cleans the maps created by camel for a given galaxy
    
    Parameters
    ----------
    path: string
        path where the input data are stored
    ouputpath: string
        path where the ouput data will be stored
    name: string
        name of the galaxy
    lsfw: float
        spectral resolution in km/s (sigma)
    fraction: float
        fraction for a lower threshold on the velocity dispersion map
    data_mask: string
        basename of the map used for threshold (e.g. snr for signal to noise ratio map)
    thrl: float
        lower threshold for cleaning
    thru: float
        upper threshold for cleaning
    option: string
        option of camel to find the files to clean (e.g. '_ssmooth')
    line: string
        line used (suffixe, e.g. '_Ha')
    clean: string
        name of the manually cleaned map
    
    XXX utiliser des listes si on veut ajouter des fichiers pour masquer?
    '''
    
    #Checking path exists
    if not(os.path.isdir(path)):
        logger.info('clean_galaxy: path % s does not exist' % (str(path)) )
        return
 
    #Checking output path exists
    if not(os.path.isdir(outputpath)):
        logger.info('clean_galaxy: path % s does not exist' % (str(outputpath)) )
        return
    
    
    smin = lsfw * fraction
    logger.info('clean_galaxy: dispersion threshold % s' % (str(smin)) )
    
    files = glob.glob(path + name + option + '*.fits')
    #fim0 = glob.glob(path + name + '/' + name + option + '_disp_*[pn]' + line + '.fits')
    fim0 = glob.glob(path + name + option + '_disp_*[pn].fits')
    fim1 = glob.glob(path + name + option + '_' + data_mask + '_*[pn].fits')
    
    try:
        hdul0 = fits.open(fim0[0])
        im0 = hdul0[0].data
        logger.info('clean_galaxy: using % s' % (str(fim0[0])) )
    except:
        logger.info('clean_galaxy: % s not found' % (str(path + name + option + '_disp_*[pn].fits')) )
        return
    try:
        hdul1 = fits.open(fim1[0])
        im1 = hdul1[0].data
        logger.info('clean_galaxy: using % s' % (str(fim1[0])) )
    except:
        logger.info('clean_galaxy: % s not found' % (str(path + name + option + '_' + data_mask + '_*[pn].fits')) )
        return
    
    #True where im0>=smin
    mask0 = create_mask(im0, thrl=smin)
    
    #True where thrl<=im1<=thru
    mask1 = create_mask(im1, thrl=thrl, thru=thru)
    
    #Keep positions where the values are out of bounds
    mask = np.where((np.logical_not(mask0)) | (np.logical_not(mask1)))
    
    #ipdb.set_trace()
    if clean is not None:
        fcl = glob.glob(path + clean)
        try:
            hducl = fits.open(fcl[0])
            imcl = hducl[0].data
            logger.info('clean_galaxy: using % s' % (str(fcl[0])) )    
            maskcl = create_mask(imcl, thrl=None, thru=None)
            mask2 = np.where((np.logical_not(mask0)) | (np.logical_not(mask1)) | (np.logical_not(maskcl)))
        except:
            logger.info('clean_galaxy: % s not found' % (str(path + clean)) )
            clean = None
    
    for fim in files:
        if 'clean' in fim:
            continue
        hdul = fits.open(fim)
        im = hdul[0].data
        if im.ndim == 3:
            continue
        
        thr = 0
        if thrl is not None:
            thr = thrl
        if thru is not None:
            thr = thru
            
        hdul[0].data = apply_mask(mask, im)
        fimcl = fim.split('.fits')[0].split('/')[-1] + '_clean%3.1f.fits' %thr
        hdul.writeto(outputpath + fimcl, overwrite=True)

        if clean is not None:
            hdul[0].data = apply_mask(mask2, im)
            fimcl = fim.split('.fits')[0].split('/')[-1] + '_mclean%3.1f.fits' %thr
            hdul.writeto(outputpath + fimcl, overwrite=True)
        logger.info('output written in %s' %(outputpath+fimcl))



def clean_setofgalaxies(path='/home/wilfried/ST2/', outfilename='list_output_folders', filename='list_gal', 
                        fraction=1., data_mask='snr', thrl=None, thru=None, option='_ssmooth', line='', 
                        clean='clean.fits'):
    '''
    This function cleans the maps created by camel for a list of galaxies
    
    Parameters
    ----------
    path: string
        path where the data are stored
    filename: string
        name of the file containing the list of galaxies and the associated spectral resolution in km/s (sigma)
    outfilename: string
        name of of the file containing the list of the output folders names
    fraction: float
        fraction for a lower threshold on the velocity dispersion map
    data_mask: string
        basename of the map used for threshold (e.g. snr for signal to noise ratio map)
    thrl: float
        lower threshold for cleaning
    thru: float
        upper threshold for cleaning
    option: string
        option of camel to find the files to clean (e.g. '_ssmooth')
    line: string
        line used (suffixe, e.g. '_Ha')
    clean: string
        name of the manually cleaned map
    '''
    
    cat = ascii.read(filename)
    
    for ligne in cat:
        name        = ligne[0].split('/')[-1].split('.config')[0]
        path        = ligne[0].split(name+'.config')[0]
        outpath     = '../outputs/MUSE/' + path.split('../data/')[1]
        clean_galaxy(path, outpath, name, ligne[1], fraction, data_mask=data_mask, thru=thru, thrl=thrl, line=line, option=option, clean=clean)

def compute_velres(z, lbda0, a2=5.835e-8, a1=-9.080e-4, a0=5.983):
    '''
    This function computes the resolution in terms of velocity sigma from the line restframe wavelength, the redshift of the source and from MUSE LSF model: FWHM(lbda) = a2 * lbda ** 2 + a1 * lbda + a0
    
    Parameters
    ----------
    z: float
        redshift of the galaxy
    lbda0: float
        rest frame wavelength of the line used to infer kinematics (in Angstroms)
    a2: float
        lambda ** 2 coefficient of the variation of LSF FWHM with respect to lambda
    a1: float
    '''
    '''
        lambda ** 1 coefficient of the variation of LSF FWHM with respect to lambda
    a0: float
        lambda ** 0 coefficient of the variation of LSF FWHM with respect to lambda
    
    Returns
    -------
    Returns the observed wavelength, the LSF FWHM in Angstroms and the LSF dispersion in km/s, assuming a Gaussian shape for the LSF profile.
    '''
    lbda = lbda0 * (1 + z)
    fwhm = a2 * lbda ** 2 + a1 * lbda + a0
    velsig = fwhm / (lbda * 2 * np.sqrt(2 * np.log(2))) * ct.c.value * 1e-3
    return lbda, fwhm, velsig
    
    
def velres_setofgalaxies(inname, outname, lbda0, a2=5.835e-8, a1=-9.080e-4, a0=5.983):
    '''
    This function computes the resolution in velocity for a list of galaxies
    
    Parameters
    ----------
    inname: string
        input file name containing the list of galaxies. In this version, the redshift is in the name itself, as well as the line used.
    outname: string
        output file name containing the list of galaxies and the spectral resolution
    lbda0: float
        rest frame wavelength of the line used to infer kinematics (in Angstroms)
    a2: float
        lambda ** 2 coefficient of the variation of LSF FWHM with respect to lambda
    a1: float
        lambda ** 1 coefficient of the variation of LSF FWHM with respect to lambda
    a0: float
        lambda ** 0 coefficient of the variation of LSF FWHM with respect to lambda
    '''
    
    cat = ascii.read(inname)
    f = open(outname, 'w')
    
    for ligne in cat:
        z    = float(ligne[1])
        lbda, fwhm, velsig = compute_velres(z, lbda0, a2=a2, a1=a1, a0=a0)
        line = '{0:100} {1:5.1f} \n'.format(ligne[0], velsig)
        f.write(line)
    f.close()
   
def main():
    #spectral: FWHM = 2.5 A => sigma = FWHM/(2.*np.sqrt(2.*np.log(2))) = 1.06165 A
    # sigma = 2.5/(2.*np.sqrt(2.*np.log(2.))) / (lbda * (z+1)) * ct.c*1e-3 (km/s)
    
    ##Groups
    #path0 = '/media/bepinat/WD2To/data/Instruments/MUSE/groups/kinematics/'
    #groups = ['gr28_shallow', 'gr32', 'gr83', 'gr116', 'gr28_best_seeing', 'gr28_deep']
    #groups = ['gr28_shallow']
    #for group in groups:
        #path = path0 + group + '/o2/'
        #clean_setofgalaxies(path, thru=None, thrl=5, fraction=1., filename='../clean_o2.txt', clean='clean.fits', option='_ssmooth', line='_OII3729', data_mask='snr')
        #path = path0 + group + '/o3hb/'
        #clean_setofgalaxies(path, thru=None, thrl=5, fraction=1., filename='../clean_o3hb.txt', clean='clean.fits', option='_ssmooth', line='_OIII5007', data_mask='snr')
        
    #UDF
    path        = '/home/wilfried/ST2/'
    scripts     = 'scripts_python_Benoit/'
    inname      = path + scripts + 'list_gal'
    inname      = path + scripts + 'test'
    outname     = path + scripts + 'clean_o2'
    lbda0       = 3729.  # OII wavelength at restframe in Angstroms
    velres_setofgalaxies(inname, outname, lbda0)
    
    fraction = 0.8
    clean_setofgalaxies(path=path, thru=None, thrl=5, fraction=fraction, filename=outname, clean='clean.fits', option='_ssmooth', line='_OII3729', data_mask='snr')
    

if __name__ == "__main__":
    main()
