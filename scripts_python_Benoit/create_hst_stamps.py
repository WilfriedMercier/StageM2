#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import astropy.wcs as wcs
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import numpy as np
import astropy.coordinates as coord
import astropy.units as u
import copy
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import pylab as plb
import ipdb


def extract_psf(image_file, rastr, decstr, size=2., outname='psf.fits'):
    '''This function enables to extract an image centered on given coordinates from a large image
    
    Parameters
    ----------
    image_file: string
        name of the input image (fits format)
    rastr: string
        right ascention in sexagesimal (hh:mm:ss.ss)
    decstr: string
        declination in sexagesimal (dd:mm:ss.s)
    size: float
        half size of the small images in arcsec
    outname: string
        output fits file name
    '''
    hdul = fits.open(image_file)
    im = hdul[0].data
    hdr = hdul[0].header
    ra = coord.Angle(rastr, unit=u.hour)
    dec = coord.Angle(decstr, unit=u.degree)
    radeg = np.array(ra.degree)
    decdeg = np.array(dec.degree)
    hst_pix = np.sqrt(hdr['CD1_1'] ** 2 + hdr['CD1_2'] ** 2) * 3600
    npix = np.round(size / hst_pix)
    w = wcs.WCS(hdr, hdul)
    xc, yc = w.wcs_world2pix(radeg, decdeg, 0)
    hdr['CRPIX1'] = hdr['CRPIX1'] - xc + npix
    hdr['CRPIX2'] = hdr['CRPIX2'] - yc + npix
    im1 = im[yc - npix:yc + npix, xc - npix:xc + npix]
    hdu = fits.PrimaryHDU(im1, hdr)
    hdulist = fits.HDUList(hdu)
    hdulist.writeto(outname, overwrite=True)


def extract_stamps_udf(image_file, gal_list, size=2., factor=1., pathout='./', suffixe='f814_comb'):
    '''This function enables to extract several images centered on galaxies from a large image
    
    Parameters
    ----------
    image_file: string
        name of the input image (fits format)
    gal_list: string
        name of the galaxy list (columns: ID, z, Flag, RA, DEC, I_AB)
    size: float
        half size of the small images in arcsec
    factor: float
        factor by which multiplying the data for using GALFIT
    pathout: string
        path where to write the extracted images
    '''
    
    imname = str.split(image_file, '/')[-1]
    hdul = fits.open(image_file)
    im = hdul[0].data * factor
    hdr = hdul[0].header
    hdr['UZERO'] = -2.5 * np.log10(1. / factor)
    cat = ascii.read(gal_list)
    
    radeg = cat['ra']
    decdeg = cat['dec']
    
    hst_pix = np.sqrt(hdr['CD1_1'] ** 2 + hdr['CD1_2'] ** 2) * 3600
    npix = np.round(size / hst_pix)
    
    w = wcs.WCS(hdr, hdul)
    xc, yc = w.wcs_world2pix(radeg, decdeg, 0)
    
    for i in range(np.size(cat['ID'])):
        hdrc = copy.deepcopy(hdr)
        hdrc['CRPIX1'] = hdr['CRPIX1'] - xc[i] + npix
        hdrc['CRPIX2'] = hdr['CRPIX2'] - yc[i] + npix
        im1 = im[np.int(np.round(yc[i] - npix)):np.int(np.round(yc[i] + npix)), np.int(np.round(xc[i] - npix)):np.int(np.round(xc[i] + npix))]
        hdu = fits.PrimaryHDU(im1, hdrc)
        hdulist = fits.HDUList(hdu)
        hdulist.writeto(pathout + str(cat['ID'][i]) + '_' + suffixe + '.fits', overwrite=True)



def extract_stamps_groups(image_file, gal_list, size=2., factor=1., pathout='./'):
    '''This function enables to extract several images centered on galaxies from a large image
    
    Parameters
    ----------
    image_file: string
        name of the input image (fits format)
    gal_list: string
        name of the galaxy list (columns: ID, z, Flag, RA, DEC, I_AB)
    size: float
        half size of the small images in arcsec
    factor: float
        factor by which multiplying the data for using GALFIT
    pathout: string
        path where to write the extracted images
    '''
    
    imname = str.split(image_file, '/')[-1]
    hdul = fits.open(image_file)
    im = hdul[0].data * factor
    hdr = hdul[0].header
    hdr['UZERO'] = -2.5 * np.log10(1. / factor)
    cat = ascii.read(gal_list)
    
    radeg = cat['ra']
    decdeg = cat['dec']
    
    hst_pix = np.sqrt(hdr['CD1_1'] ** 2 + hdr['CD1_2'] ** 2) * 3600
    npix = np.round(size / hst_pix)
    
    w = wcs.WCS(hdr, hdul)
    xc, yc = w.wcs_world2pix(radeg, decdeg, 0)
    
    for i in range(np.size(cat['ID'])):
        hdrc = copy.deepcopy(hdr)
        hdrc['CRPIX1'] = hdr['CRPIX1'] - xc[i] + npix
        hdrc['CRPIX2'] = hdr['CRPIX2'] - yc[i] + npix
        im1 = im[np.int(np.round(yc[i] - npix)):np.int(np.round(yc[i] + npix)), np.int(np.round(xc[i] - npix)):np.int(np.round(xc[i] + npix))]
        hdu = fits.PrimaryHDU(im1, hdrc)
        hdulist = fits.HDUList(hdu)
        nm=image_file.split('_')[1].split('.fits')[0]
        hdulist.writeto(pathout + str(cat['ID'][i]) + '_' + nm + '.fits', overwrite=True)

def extract_stamps(image_file, gal_list, size=2., factor=1., pathout='./'):
    '''This function enables to extract several images centered on galaxies from a large image
    
    Parameters
    ----------
    image_file: string
        name of the input image (fits format)
    gal_list: string
        name of the galaxy list (columns: ID, ID_Jarle, RA_Jarle, DEC_Jarle, z_jarle, log(M*)_BC03, log(SFR)_CL01
    size: float
        half size of the small images in arcsec
    factor: float
        factor by which multiplying the data for using GALFIT
    pathout: string
        path where to write the extracted images
    '''
    
    imname = str.split(image_file, '/')[-1]
    hdul = fits.open(image_file)
    im = hdul[0].data * factor
    hdr = hdul[0].header
    hdr['UZERO'] = -2.5 * np.log10(1. / factor)
    c = np.genfromtxt(gal_list, names=True, dtype=('O', 'O', 'O', 'O', 'f8', 'f8', 'f8'))
    
    ra = coord.Angle(list(c['RA_Jarle']), unit=u.hour)  # create an Angle object
    dec = coord.Angle(list(c['DEC_Jarle']), unit=u.degree)  # create an Angle object
    radeg = np.array(ra.degree)
    decdeg = np.array(dec.degree)
    
    hst_pix = np.sqrt(hdr['CD1_1'] ** 2 + hdr['CD1_2'] ** 2) * 3600
    npix = np.round(size / hst_pix)
    
    w = wcs.WCS(hdr, hdul)
    xc, yc = w.wcs_world2pix(radeg, decdeg, 0)
    
    for i in range(np.size(c['ID'])):
        hdrc = copy.deepcopy(hdr)
        hdrc['CRPIX1'] = hdr['CRPIX1'] - xc[i] + npix
        hdrc['CRPIX2'] = hdr['CRPIX2'] - yc[i] + npix
        im1 = im[yc[i] - npix:yc[i] + npix, xc[i] - npix:xc[i] + npix]
        hdu = fits.PrimaryHDU(im1, hdrc)
        hdulist = fits.HDUList(hdu)
        hdulist.writeto(pathout + c['ID_Jarle'][i] + '_' + imname, overwrite=True)
        #hdulist.writeto('PSF' + c['ID'][i] + '_f814_comb.fits', overwrite=True)


def create_modeled_psf(file1, file2, fileout='psf.fits'):
    '''This function enables to write a fits image containing the psf modeled by galfit.
    
    Parameters
    ----------
    file1: string
        name of the stamp image of a real psf (to recover the heaer)
    file2: string
        name of GALFIT output (to recover the model)
    fileout: string
        name of the output file
    '''
    
    hdul1 = fits.open(file1)
    hdr = hdul1[0].header
    hdul2 = fits.open(file2)
    data = hdul2[2].data
    hdu = fits.PrimaryHDU(data, hdr)
    hdulist = fits.HDUList(hdu)
    hdulist.writeto(fileout, overwrite=True)


def display_hst_models(file1, fileout='test.pdf', title='title'):
    '''This function enables to display an image, the associated GALFIT model and residuals.
    
    Parameters
    ----------
    file1: string
        name of GALFIT file that contains the model
    fileout: string
        name of the output file
    title: string
        title of the output image
    '''
    hdul = fits.open(file1)
    data = hdul[1].data
    model = hdul[2].data
    res = hdul[3].data
    fig = plt.figure(figsize=(12, 3))
    fig.suptitle(title)
    ax = fig.add_subplot(131)
    #norm = mpl.colors.LogNorm(vmin=np.min([data, model]), vmax=np.max([data, model]))
    imgplot = plt.imshow(data, origin='lower', cmap='spectral', interpolation='nearest', vmin=np.min([data, model]), vmax=np.max([data, model]))
    plt.colorbar(fraction=0.05, shrink=1.)
    ax = fig.add_subplot(132)
    imgplot = plt.imshow(model, origin='lower', cmap='spectral', interpolation='nearest', vmin=np.min([data, model]), vmax=np.max([data, model]))
    plt.colorbar(fraction=0.05, shrink=1.)
    ax = fig.add_subplot(133)
    imgplot = plt.imshow(res, origin='lower', cmap='spectral', interpolation='nearest')
    plt.colorbar(fraction=0.05, shrink=1.)
    plt.savefig(fileout)
    #fig.clf()
    plt.close()


def display_all_hst_models(path, pathout, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt'):
    '''This function enables to display an image, the associated GALFIT model and residuals for all models in path with the basename
    
    Parameters
    ----------
    path: string
        path where to find the GALFIT models
    pathout: string
        path where to write the images
    basename: string
        suffixe attached to the GALFIT models
    filecorr: string
        name of the file containing the corresponding name
    '''
    cat = np.genfromtxt(path + filecorr, dtype=('i', 'i'), names=True)
    files = glob.glob(path + '*' + basename)
    for file1 in files:
        filen = file1.split('/')[-1]
        name = filen.split(basename)[0]
        namejb = str(cat['JB'][np.where(cat['TC'] == int(name))][0])
        namecorr = 'JB_' + namejb + ' / TC_' + name
        fileout = pathout + name + '_f814_comb.pdf'
        display_hst_models(file1, fileout=fileout, title=namecorr)


def display_all_hst_models_groups(path, pathout, basename='_f814_comb_out.fits'):
    '''This function enables to display an image, the associated GALFIT model and residuals for all models in path with the basename
    
    Parameters
    ----------
    path: string
        path where to find the GALFIT models
    pathout: string
        path where to write the images
    basename: string
        suffixe attached to the GALFIT models
    '''
    files = glob.glob(path + '*' + basename)
    for file1 in files:
        filen = file1.split('/')[-1]
        gr = file1.split('/')[-2].split('_')[-1]
        name = filen.split(basename)[0]
        fileout = pathout + name + '_f814_comb.pdf'
        display_hst_models(file1, fileout=fileout, title=gr.upper() + ' #ID' + name)


def main():
    '''
    '''
    
    #------#
    # UDF  #
    #------#
#    factor = 1e4
#    path = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/camel/'
#    gal_list = path + 'udf_mos_c042_e030_spectrocat.txt'
#    path_hst = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/hst/'
#    pathout = path_hst + 'stamps/'
#    bands = ['f105w', 'f125w', 'f140w', 'f160w', 'f435w', 'f606w', 'f775w', 'f814w', 'f850lp']
    #band = 'f160w'
#    sz = 2.
#    for band in bands:
        #image_file = glob.glob(path_hst + '*_hst_*_' + band + '_*.fits')[0]
#        image_file = glob.glob(path_hst + 'hlsp_xdf_hst_*_' + band + '_*.fits')[0]
#        print(image_file)
#        extract_stamps_udf(image_file, gal_list, size=sz, factor=factor, pathout=pathout, suffixe=band)
    #gal_list = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho_largerstamps.txt'
    #sz = 3.
    #extract_stamps(image_file, gal_list, size=sz, factor=factor, pathout=pathout)
    
    
    #------#
    # HDFS #
    #------#
    #factor = 1e4
    #image_file = '/media/bepinat/WD2To/data/Instruments/MUSE/data/Commissioning/hdf/f814_comb.fits'
    #pathout = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_stamps/'
    #gal_list = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/interz_els_HDFS_v031c_vfinal_forbenoit.txt'
    #sz = 2.
    #extract_stamps(image_file, gal_list, size=sz, factor=factor, pathout=pathout)
    #gal_list = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho_largerstamps.txt'
    #sz = 3.
    #extract_stamps(image_file, gal_list, size=sz, factor=factor, pathout=pathout)
    #gal_list = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/HDFS_stars_extract.txt'
    #sz = 1.
    #extract_stamps(image_file, gal_list, size=sz, factor=factor)
    
    #for imname in ('f300_comb.fits', 'f450_comb.fits', 'f606_comb.fits', 'f814_comb.fits'):
        #factor = 1e4
        #pathim = '/media/bepinat/WD2To/data/Instruments/MUSE/data/Commissioning/hdf/'
        #image_file = pathim + imname
        #pathout = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_davor/'
        #gal_list = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_davor/davor_list.txt'
        #sz = 4.
        #extract_stamps(image_file, gal_list, size=sz, factor=factor, pathout=pathout)
    
    ##file1 = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psf1/10_f814_comb_out.fits'
    ##display_hst_models(file1, title=file1)
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffat/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffat_bulgefree/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfstack/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfstack_bulgefree/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psftinytim/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psftinytim_bulgefree/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffattrujillo_bulgefree/'
    #display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    
    ##path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psf_tinytim/'
    ##display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')
    ##path = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfJR/'
    ##display_all_hst_models(path, path, basename='_f814_comb_out.fits', filecorr='names_correspondances.txt')

    #--------#
    # groups #
    #--------#
    
    path = '/home/wilfried/ST2/data/hst/'
    groups = ['CGr114', 'CGr23', 'CGr26', 'CGr28', 'CGr30',
              'CGr32', 'CGr34', 'CGr51', 'CGr61', 'CGr79',
              'CGr84', 'CGr84-N']
    groups = [groups[7]]
    factor = 1e4
    szg = 3.0
#    szs = 1.
    for group in groups:
        pathgr = path + group + '/'
#        pathpsf = path + 'acs_mosaic_2.0_' + group + '/psf/'
        image_file = pathgr + 'HST_' + group + '.fits'
        gal_list = pathgr + 'catalog_field_galaxies_' + group + '.txt'
        #star_list = path + group + '_starcatalog.txt'
        extract_stamps_groups(image_file, gal_list, pathout=pathgr, size=szg, factor=factor)
        #extract_stamps_groups(image_file, star_list, pathout=pathpsf, size=szs, factor=factor)
        #display_all_hst_models_groups(pathgr, pathgr, basename='_f814_comb_out.fits')


if __name__ == "__main__":
    main()
