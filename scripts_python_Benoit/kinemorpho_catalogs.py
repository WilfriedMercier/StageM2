#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import numpy as np
#import asciitable
#import pyfits as pf
from astropy import constants as ct
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import astropy.wcs as wcs
import astropy.coordinates as coord
import astropy.units as u
import matplotlib as mpl
from matplotlib import pyplot as plt
import copy
import ipdb

def convert_udfcat_kincat(udfcat, outcat, filedir, lbda0=3729, prefix='udf_mos_c042_e030_'):
    hdu = fits.open(udfcat)
    cat = hdu[1].data
    f = open(outcat, 'w')
    result = []
    line = '%35s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %7s \n'%('gal', 'x', 'y', 'pa', 'i', 'vs', 'vm', 'd', 'sig', 'psfx', 'psfz', 'smooth')   
    f.write(line)
    for obj in cat:
        try:
            dossier = glob.glob(filedir + prefix + str(obj['ID']) + '_*')[0]
        except:
            print(filedir + prefix + str(obj['ID']) + '_*')
            continue
        print(dossier)
        name = dossier.split('/')[-1]
        kinim = dossier + '/' + name + '_ssmooth_vel_common.fits'
        hdu = fits.open(kinim)
        hdr = hdu[0].header
        hdr.remove('WCSAXES')
        wk = wcs.WCS(hdr)
        ra = obj['RA']
        dec = obj['DEC']
        print(ra, dec)
        #ipdb.set_trace()
        xc, yc = wk.wcs_world2pix(ra, dec, 0)
        # Je prends pour l'instant la bande F160W pour toutes les galaxies
        pa = obj['gfit_pa_h']
        inc = np.degrees(np.arccos(obj['gfit_q_h']))  # dans certains = nan ou -999, dans certains cas =1, 0.5 à la décimale près...
        pixsize = np.sqrt(hdr['CD1_1'] ** 2 + hdr['CD1_2'] ** 2) * 3600  # taille du pixel MUSE en arcsec
        psfx = obj['GF2D_camel_oii_PSF_Moffat_fwhm'] / pixsize  #attention, FWHM PSF ajustée par une Moffat, habituellement, Gaussienne utilisée
        z = float(name.split(sep='Z_')[1])
        lbda, fwhm, velsig = compute_velres(z, lbda0)
        psfz = velsig
        result.append([name, xc, yc, pa, inc, psfx, psfz])
    for obj in result:
        line = '%35s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.1f \n'%(obj[0], obj[1], obj[2], obj[3], obj[4], 0, 80, 2., 0, obj[5], obj[6], 2)
        f.write(line)
    f.close()
    return result
 #'gfit_f_h'
 #'gfit_mag_h'
 #'gfit_emag_h'
 #'gfit_sma_h'
 #'gfit_esma_h'
 #'gfit_n_h'
 #'gfit_en_h'
 #'gfit_q_h'
 #'gfit_eq_h'
 #'gfit_pa_h'
 #'gfit_epa_h'
 #'gfit_snr_h'

def compute_velres(z, lbda0, a2=5.835e-8, a1=-9.080e-4, a0=5.983):
    '''
    This function computes the resolution in terms of velocity sigma from the line restframe wavelength, the redshift of the source and from MUSE LSF model: FWHM(lbda) = a2 * lbda ** 2 + a1 * lbda + a0. Default is UDF mosaic parameters
    
    Parameters
    ----------
    z: float
        redshift of the galaxy
    lbda0: float
        rest frame wavelength of the line used to infer kinematics (in Angstroms)
    a2: float
        lambda ** 2 coefficient of the variation of LSF FWHM with respect to lambda
    a1: float
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


def convert_morphocat_radec(pathm, pathc):
    cat = ascii.read(pathm + 'recap_morpho_params.txt')
    catc = ascii.read(pathc + 'names_correspondances_tc_jb.txt')
    index=np.argsort(cat['ID'])
    result2 = []
    f2 = open(pathm + 'morph2kin2.txt', 'w')
    line = '%10s %10s %15s %15s %10s %10s %10s %10s \n'%('ID', 'ID2', 'RA', 'DEC', 'PA', 'B/A', 'Rd', 'Rb')
    f2.write(line)
    for ind in index:
        morphim = pathm + str(cat['ID'][ind]) + '.fits'  #'_f814_comb.fits'
        print(morphim)
        wm = wcs.WCS(morphim)
        idjb = catc['JB'][catc['TC'] == cat['ID'][ind]][0]
        ra, dec = wm.wcs_pix2world(cat['X'][ind] - 1, cat['Y'][ind] - 1, 0)  # -1 car premier pixel est 0 (et non 1)
        result2.append([idjb, cat['ID'][ind], ra, dec, cat['pa_d'][ind], cat['b/a_d'][ind], cat['R_d'][ind], cat['R_b'][ind]])
    for obj in result2:
        line = '%10i %10i %15.8f %15.8f %10.2f %10.2f %10.2f %10.2f \n'%(obj[0], obj[1], obj[2], obj[3], obj[4], obj[5], obj[6], obj[7])
        #print('%30s %8.2f %8.2f %8.2f %8.2f'%(obj[0], obj[1], obj[2], obj[3], obj[4]))
        f2.write(line)
    f2.close()
    return result2

def convert_morphocat_kincat(pathm, pathk):
    cat = ascii.read(pathm + 'recap_morpho_params.txt')
    print(cat)
    print(cat['ID'])
    index=np.argsort(cat['ID'])
    result = []
    f = open(pathk + 'morph2kin.txt', 'w')
    line = '%30s %8s %8s %8s %8s \n'%('ID', 'X', 'Y', 'PA', 'INC')
    f.write(line)
    for ind in index:
        morphim = pathm + str(cat['ID'][ind]) + '.fits' #_f814_comb.fits #'.fits'
        print(morphim)
#        dossiers = glob.glob(pathk + '*[obj, _]' + str(cat['ID'][ind]) + '_*')
#        dossiers = glob.glob(pathk + str(cat['ID'][ind]) + '_*') #use this for one galaxy and if the morphology files are named starting with the group instead of the ID number
        dossiers = glob.glob(pathk + '*[obj, _]' + cat['ID'][ind].split('_')[0] + '_*') #change to this line for groups 84A and B and 51A and B
        print(dossiers)
        #ipdb.set_trace()
        if np.size(dossiers) == 0: continue
        #if cat['ID'][ind] != 70: continue
        wm = wcs.WCS(morphim)
        #ra, dec = wm.wcs_pix2world(cat['X'][ind] - 2, cat['Y'][ind] - 2, 0)  # -1 pour la PSF (hdfs), -1 car premier pixel est 0 (et non 1)
        ra, dec = wm.wcs_pix2world(cat['X'][ind] - 1, cat['Y'][ind] - 1, 0)  # -1 car premier pixel est 0 (et non 1)
        for dossier in dossiers:
            print(dossier)
            name = dossier.split('/')[-1]
            kinim = dossier + '/' + name + '_ssmooth_vel_common.fits'
            wk = wcs.WCS(kinim)
            xc, yc = wk.wcs_world2pix(ra, dec, 0)
            result.append([name, xc, yc, cat['pa_d'][ind], np.degrees(np.arccos(cat['b/a_d'][ind]))])
    for obj in result:
        line = '%30s %8.2f %8.2f %8.2f %8.2f \n'%(obj[0], obj[1], obj[2], obj[3], obj[4])
        f.write(line)
    f.close()
    return result


def recap_kinem_params(path, basename='hdfs_v*', model='slp', suff1='_xyi', suff2='_mclean5.0'):
    '''
    Creates an output file 'recap_kinematics_parameters.txt' containing the output kinematics parameters from the models.
    
    Parameters
    ----------
    path: string
        path of the subdirectories
    basename: string
        suffixe
    model: string
        string corresponding to the model used (exp, iso, slp, ata)
    suff1: string
        string corresponding to the suffixe for parameters fixed
    suff2: string
        string corresponding to the suffixe for the cleaning used
    '''
    gals = glob.glob(path + basename)
    gals.sort()
    f = open(path + 'recap_kinematics_parameters_2_' + model + suff1 + suff2 + '.txt', 'w')
    line = '%38s %6s %5s %6s %5s %6s %5s %6s %5s %6s %5s %6s %5s %6s %5s %6s %5s %6s %5s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n'%('ID', 'X', 'DX', 'Y', 'DY', 'VS', 'DVS', 'PA', 'DPA','INC','DINC', 'VC', 'DVC', 'RC', 'DRC', 'SIG', 'DSIG', 'CHI2', 'DOF', 'MEANRESV', 'MEDRESV', 'STDRESV', 'MINRESV', 'MAXRESV', 'MEANRESD', 'MEDRESD', 'STDRESD', 'MINRESD', 'MAXRESD', 'RLAST', 'RLASTSYM')
    
    f.write(line)
    for gal in gals:
        filen = gal + '/' + gal.split('/')[-1] + '_parameters_red_' + model + suff1 + suff2 + '.txt'
        filen2 = gal + '/' + gal.split('/')[-1] + '_parameters_residual_' + model + suff1 + suff2 + '.txt'
        filen3 = gal + '/' + gal.split('/')[-1] + '_vmax_map_rlast' + suff2 + '.txt'
        cube = gal + '/' + gal.split('/')[-1] + '_ssmooth_cube.fits'
#        print("gal", gal, "\n", gal.split('/')[-1], "\n", filen, "\n", filen2, "\n", filen3, "\n", cube)
        line = '%38s '%(gal.split('/')[-1])
        try:
            data = np.genfromtxt(filen, names=True, usecols=(0,1,2,3,4,5,6,7))
            data1 = np.genfromtxt(filen, names=True, usecols=(8,9), skip_footer=True)
            data2 = np.genfromtxt(filen2, names=True)
            data3 = np.genfromtxt(filen3, usecols=(6,7))
            hdul = fits.open(cube)
            pixsize = np.sqrt(hdul[0].header['CD1_1'] ** 2 + hdul[0].header['CD1_2'] ** 2) * 3600
            #pixsize = np.max([np.abs(hdul[0].header['CDELT1']), np.sqrt(hdul[0].header['CD1_1'] ** 2 + hdul[0].header['CD1_2'] ** 2)]) * 3600
            
            for i in range(8):
                line += '%6.1f %5.1f '%(data[0][i], data[1][i])
            line += '%6.2f %5i '%(data1['chi2'], data1['dof'])
            for i in data2.item():
                line += '%8.1f '%(i)
            for i in list(data3):
                line += '%8.1f '%(i / pixsize)
            line += ' \n'
            f.write(line)
        except:
            print('1')
            #for i in range(8):
                #line += '%6.1f %5.1f '%(-99.9, -99.9)
            #line += '%6.2f %5i '%(-99.9, -99.9)
            #for i in range(12):
                #line += '%8.1f '%(-99.9)
            #line += ' \n'
            #f.write(line)
    f.close()

def compare_pa(morphof, kinemf):
    '''
    '''
    #Objectif: faire un plot qui compare les pa morpho (galfit) vs kinematics en indiquant le nom de objets
    #We need to open morpho parameters file + kinem parameters file + correspondance name?
    
    morpho = np.genfromtxt(morphof, names=True)
    dt = morpho.dtype.descr
    dt[0] = (dt[0][0], 'O')
    dt = np.dtype(dt)
    morpho = np.genfromtxt(morphof, names=True, dtype=dt)

    kin = np.genfromtxt(kinemf, names=True)
    dt = kin.dtype.descr
    dt[0] = (dt[0][0], 'O')
    dt = np.dtype(dt)
    kin = np.genfromtxt(kinemf, names=True, dtype=dt)
    indm = list()
    indk = list()
    for i in range(kin.shape[0]):
        idsplit = kin['ID'][i].split(b'_')
        name = idsplit[2].split(b'obj')[1]
        line = idsplit[3]
        red = idsplit[1]
        cc = np.where(morpho['ID'] == name)[0]
        if cc.size != 0:
            indm.append(cc[0])
            indk.append(i)
        else:
            print(name)
    dpa = np.mod((morpho['pa_d'][np.array(indm)] - kin['PA'][np.array(indk)]), 180)
    inc = kin['INC'][np.array(indk)]
    names = morpho['ID'][np.array(indm)]
    cc = np.where(dpa > 90)
    dpa[cc] -= 180
    plt.plot(inc, dpa, 'bo', [0, 90], [30, 30], 'r', [0, 90], [-30, -30], 'r', [0, 90], [0, 0], 'k--')
    for i, txt in enumerate(names):
        plt.text(inc[i], dpa[i], txt.decode('utf8'))
    plt.axis([0.,90.,-90.0,90.0])
    plt.yticks(np.arange(-90, 91, 30))
    plt.ylabel(r'$\Delta PA$ (Deg)')
    plt.xlabel(r'$i$ (Deg)')
    plt.savefig('deltaPA_inc.pdf')
    plt.close()
    #plt.show()
    
def inclination_histo(filen):
    '''
    Plots inclinations from galfit as a function of the inclinations from Trujillo et al.
    Plots also histogram of inclinations from Galfit and Trujillo
    
    Parameters
    ----------
    filen: string
        filename containing both Trujillo and Galfit inclinations
    '''
    inc = np.genfromtxt(filen, names=True)
    plt.plot(inc['inct'], inc['incg'], 'x', [0, 90], [0, 90])
    plt.show()
    cc=np.where(inc['inct'] != -1)
    plt.hist(inc['inct'][cc], bins=9, range=[0, 90])
    plt.show()
    plt.hist(inc['incg'], bins=9, range=[0, 90])
    plt.show()

    #Objective: plot an histogram of inclinations + compare inclinations from Galfit and Trujillo + galpak?
    
def read_galfit_output(galfile):
    '''
    Reads an output from Galfit and store the result inside a structure
    
    Parameters
    ----------
    galfile: string
        name of the galfit output
    '''
    lines = open(galfile, 'r')
    ncomp = -1
    dtype = {'names':('type', 'X', 'Y', 'Mag', 'R_e', 'n', 'b/a', 'pa'), 'formats':('O','f8','f8','f8','f8','f8','f8','f8')}
    result = np.ndarray(1, dtype=dtype)
    for line in lines:
        lsplit = line.split(')')
        #print(lsplit[0].strip().isdigit())
        if lsplit[0].strip().isdigit():
            num = int(lsplit[0])
            params = lsplit[1].split('#')[0]
            if num == 0:
                typec = params.strip()
                if typec == 'sersic':
                    ncomp += 1
            if ncomp == -1: continue
            if typec != 'sersic': continue
            if num == 0:
                if ncomp > 0:
                    result = np.append(result, np.ndarray(1, dtype=dtype))
                result[ncomp]['type'] = typec
                continue
            p = np.array([np.float(param) for param in params.split()])
            if num == 1:
                result[ncomp]['X'] = p[0]
                result[ncomp]['Y'] = p[1]
            if num == 3:
                result[ncomp]['Mag'] = p[0]
            if num == 4:
                result[ncomp]['R_e'] = p[0]
            if num == 5:
                result[ncomp]['n'] = p[0]
            if num == 9:
                result[ncomp]['b/a'] = p[0]
            if num == 10:
                result[ncomp]['pa'] = p[0]
    return result
    
def recap_morpho_params0(path, basename='_f814_comb_galfit.txt'):
    '''
    Reads all outputs from Galfit in a given directory and write all morphological parameters (bulge + disk) in a file
    
    Parameters
    ----------
    path: string
        path where to find the GALFIT models
    basename: string
        suffixe attached to the GALFIT models
    '''
    f = open(path + 'recap_morpho_params.txt', 'w')
    line = '%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n'%('ID', 'X', 'Y', 'Mag_d', 'R_d', 'n_d', 'b/a_d', 'pa_d', 'Mag_b', 'R_b', 'n_b', 'b/a_b', 'pa_b')
    f.write(line)
    
    files = glob.glob(path + '*' + basename)
    for file1 in files:
        filen = file1.split('/')[-1]
        name = filen.split(basename)[0]
        result = read_galfit_output(file1)
        line = '%10s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n'%(name, result[0]['X'], result[0]['Y'], result[0]['Mag'], result[0]['R_e'], result[0]['n'], result[0]['b/a'], result[0]['pa'], result[1]['Mag'], result[1]['R_e'], result[1]['n'], result[1]['b/a'], result[1]['pa'])
        f.write(line)
    f.close()

def read_galfit_log(galfile):
    '''
    Reads an output from Galfit and store the result inside a structure
    
    Parameters
    ----------
    galfile: string
        name of the galfit output log file
    '''
    lines = open(galfile, 'r')
    ncomp = -1
    dtype = {'names':('type', 'X', 'dX', 'Y', 'dY', 'Mag', 'dMag', 'R_e', 'dR_e', 'n', 'dn', 'b/a', 'db/a', 'pa', 'dpa'), 'formats':('O','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')}
    result = np.ndarray(1, dtype=dtype)
    prev = ''
    for line in lines:
        if ('sersic' in prev) | ('sersic' in line):
            line = line.replace(',', '').replace(':', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace('{', '').replace('}', '').replace('*', '')
            lsplit = line.split()
            prev = lsplit[0]
            if prev == 'sersic':
                ncomp += 1
                if ncomp > 0:
                    result = np.append(result, np.ndarray(1, dtype=dtype))
                p = np.array([np.float(param) for param in lsplit[1:]])
                d = ''
                result[ncomp]['type'] = prev
            else:
                p = np.array([np.float(param) for param in lsplit])
                d = 'd'
            result[ncomp][d + 'X'] = p[0]
            result[ncomp][d + 'Y'] = p[1]
            result[ncomp][d + 'Mag'] = p[2]
            result[ncomp][d + 'R_e'] = p[3]
            result[ncomp][d + 'n'] = p[4]
            result[ncomp][d + 'b/a'] = p[5]
            result[ncomp][d + 'pa'] = p[6]
        if ('Chi^2 =' in line):
            line = line.replace(',', '')
            lsplit = line.split()
            chi2 = float(lsplit[2])
            ndof = int(lsplit[5])
        if ('Chi^2/nu =' in line):
            lsplit = line.split()
            chi2r = float(lsplit[2])
    return result, chi2, ndof, chi2r

def recap_morpho_params(path, basename='_comb.log', ipath='*', outname=''): #basename='_f814_comb_galfit.log'
    '''
    Reads all outputs from Galfit in a given directory and write all morphological parameters (bulge + disk) in a file
    
    Parameters
    ----------
    path: string
        path where to find the GALFIT models
    basename: string
        suffixe attached to the GALFIT models
    ipath: string
        intermediate path name, including shell-style wildcards
    outname: string
        suffixe attached to the output file
    '''
    f = open(path + 'recap_morpho_params' + outname + '.txt', 'w')
    line = '%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n'%('ID', 'X', 'dX', 'Y', 'dY', 'Mag_d', 'dMag_d', 'R_d', 'dR_d', 'n_d', 'dn_d', 'b/a_d', 'db/a_d', 'pa_d', 'dpa_d', 'Mag_b', 'dMag_b', 'R_b', 'dR_b', 'n_b', 'dn_b', 'b/a_b', 'db/a_b', 'pa_b', 'dpa_b', 'chi2', 'ndof', 'chi2r')
    f.write(line)
    
    files = glob.glob(path + ipath + basename)
    for file1 in files:
        filen = file1.split('/')[-1]
        name = filen.split(basename)[0]
        result, chi2, ndof, chi2r = read_galfit_log(file1)
        line = '%10s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10i %10.4f \n'%(name, result[0]['X'], result[0]['dX'], result[0]['Y'], result[0]['dY'], result[0]['Mag'], result[0]['dMag'], result[0]['R_e'], result[0]['dR_e'], result[0]['n'], result[0]['dn'], result[0]['b/a'], result[0]['db/a'], result[0]['pa'], result[0]['dpa'], result[1]['Mag'], result[1]['dMag'], result[1]['R_e'], result[1]['dR_e'], result[1]['n'], result[1]['dn'], result[1]['b/a'], result[1]['db/a'], result[1]['pa'], result[1]['dpa'], chi2, ndof, chi2r)
        f.write(line)
    f.close()

def read_galfitpsf_log(galfile):
    '''
    Reads an output from Galfit and store the result inside a structure
    
    Parameters
    ----------
    galfile: string
        name of the galfit output log file
    '''
    lines = open(galfile, 'r')
    ncomp = -1
    dtype = {'names':('type', 'X', 'dX', 'Y', 'dY', 'Mag', 'dMag', 'R_e', 'dR_e', 'n', 'dn', 'b/a', 'db/a', 'pa', 'dpa'), 'formats':('O','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')}
    result = np.ndarray(1, dtype=dtype)
    prev = ''
    prev2 = ''
    mag_s = -1
    dmag_s = -1
    for line in lines:
        if ('moffat' in prev) | ('moffat' in line) | ('gaussian' in prev) | ('gaussian' in line):
            prev2 = prev
            line = line.replace(',', '').replace(':', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace('{', '').replace('}', '').replace('*', '')
            lsplit = line.split()
            prev = lsplit[0]
            if ('moffat' in prev) | ('gaussian' in prev):
                ncomp += 1
                if ncomp > 0:
                    result = np.append(result, np.ndarray(1, dtype=dtype))
                p = np.array([np.float(param) for param in lsplit[1:]])
                d = ''
                result[ncomp]['type'] = prev
            else:
                p = np.array([np.float(param) for param in lsplit])
                d = 'd'
            if ('moffat' in prev) | ('moffat' in prev2):
                result[ncomp][d + 'X'] = p[0]
                result[ncomp][d + 'Y'] = p[1]
                result[ncomp][d + 'Mag'] = p[2]
                result[ncomp][d + 'R_e'] = p[3]
                result[ncomp][d + 'n'] = p[4]
                result[ncomp][d + 'b/a'] = p[5]
                result[ncomp][d + 'pa'] = p[6]
            if ('gaussian' in prev) | ('gaussian' in prev2):
                result[ncomp][d + 'X'] = p[0]
                result[ncomp][d + 'Y'] = p[1]
                result[ncomp][d + 'Mag'] = p[2]
                result[ncomp][d + 'R_e'] = p[3]
                result[ncomp][d + 'n'] = -1
                result[ncomp][d + 'b/a'] = p[4]
                result[ncomp][d + 'pa'] = p[5]
        if ('sky' in prev) | ('sky ' in line):
            line = line.replace(',', '').replace(':', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace('{', '').replace('}', '').replace('*', '')
            lsplit = line.split()
            prev = lsplit[0]
            if prev == 'sky':
                mag_s = float(lsplit[3])
            else:
                dmag_s = float(lsplit[0])
        if ('Chi^2 =' in line):
            line = line.replace(',', '')
            lsplit = line.split()
            chi2 = float(lsplit[2])
            ndof = int(lsplit[5])
        if ('Chi^2/nu =' in line):
            lsplit = line.split()
            chi2r = float(lsplit[2])
    return result, mag_s, dmag_s, chi2, ndof, chi2r

def recap_morphopsf_params(path, basename='_comb.log', ipath='*', outname=''): #basename='_f814_comb.log'
    '''
    Reads all outputs from Galfit in a given directory and write all morphological parameters (bulge + disk) in a file
    
    Parameters
    ----------
    path: string
        path where to find the GALFIT models
    basename: string
        suffixe attached to the GALFIT models
    ipath: string
        intermediate path name, including shell-style wildcards
    outname: string
        suffixe attached to the output file
    '''
    f = open(path + 'recap_morphopsf_params' + outname + '.txt', 'w')
    if ('_mg' in outname):
        line = '%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n'%('ID', 'X1', 'dX1', 'Y1', 'dY1', 'Mag1', 'dMag1', 'FWHM1', 'dFWHM1', 'n1', 'dn1', 'b/a1', 'db/a1', 'pa1', 'dpa1', 'X2', 'dX2', 'Y2', 'dY2', 'Mag2', 'dMag2', 'FWHM2', 'dFWHM2', 'n2', 'dn2', 'b/a2', 'db/a2', 'pa2', 'dpa2', 'Mag_s', 'dMag_s', 'chi2', 'ndof', 'chi2r')
    else:
        line = '%20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n'%('ID', 'X', 'dX', 'Y', 'dY', 'Mag', 'dMag', 'FWHM', 'dFWHM', 'n', 'dn', 'b/a', 'db/a', 'pa', 'dpa', 'Mag_s', 'dMag_s', 'chi2', 'ndof', 'chi2r')
    f.write(line)
    
    files = glob.glob(path + ipath + basename)
    for file1 in files:
        filen = file1.split('/')[-1]
        name = filen.split(basename)[0]
        result, mag_s, dmag_s, chi2, ndof, chi2r = read_galfitpsf_log(file1)
        if ('_mg' in outname):
            line = '%20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10i %10.4f \n'%(name, result[0]['X'], result[0]['dX'], result[0]['Y'], result[0]['dY'], result[0]['Mag'], result[0]['dMag'], result[0]['R_e'], result[0]['dR_e'], result[0]['n'], result[0]['dn'], result[0]['b/a'], result[0]['db/a'], result[0]['pa'], result[0]['dpa'], result[1]['X'], result[1]['dX'], result[1]['Y'], result[1]['dY'], result[1]['Mag'], result[1]['dMag'], result[1]['R_e'], result[1]['dR_e'], result[1]['n'], result[1]['dn'], result[1]['b/a'], result[1]['db/a'], result[1]['pa'], result[1]['dpa'], mag_s, dmag_s, chi2, ndof, chi2r)
        else:
            line = '%20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10i %10.4f \n'%(name, result[0]['X'], result[0]['dX'], result[0]['Y'], result[0]['dY'], result[0]['Mag'], result[0]['dMag'], result[0]['R_e'], result[0]['dR_e'], result[0]['n'], result[0]['dn'], result[0]['b/a'], result[0]['db/a'], result[0]['pa'], result[0]['dpa'], mag_s, dmag_s, chi2, ndof, chi2r)
        f.write(line)
    f.close()

def extract_lines(cat_zurich, cat_group, outname, radius=1.5):
    '''
    Extract in a large catalog the line corresponding to the closest object for each object in a small catalog
    
    Parameters
    ----------
    cat_zurich: ndarray
        large catalog (from COSMOS)
    cat_group: ndarray
        small catalog (from MUSE)
    outname: string
        name of the output catalog
    radius: float
        radius around which to search (arcsec)
    '''
    f = open(outname, 'w')
    title = '%18s%18s'%('COSMOS_ID', 'distance')
    for name in cat_zurich.dtype.names[:-1]:
        title += '%18s'%name
    title += ' \n'
    f.write(title)
    #print(title)
    for obj in cat_group:
        dist = 3600 * (((cat_zurich['RA'] - obj['RA']) / obj['DEC']) ** 2 + (cat_zurich['DEC'] - obj['DEC']) ** 2) ** 0.5
        cond = np.where(dist < radius)
        if np.shape(cond)[1] < 1:
            line = '%18i%18s'%(obj['COSMOS_ID'], 'no_correspondance')
            #print('%18i: no correspondance'%obj['COSMOS_ID'])
        elif np.shape(cond)[1] > 1:
            line = '%18i%18s'%(obj['COSMOS_ID'], str(np.shape(cond)[1]) + '_correspondances')
            for i in range(np.shape(cond)[1]):
                line += '%18.7f%18.7f'%(cat_zurich['SequentialID'][cond][i], dist[cond][i])
            line += ' \n'
            #print('%18i: %i correspondances'%(obj['COSMOS_ID'], np.shape(cond)[1]))
            f.write(line)
            line = '%18i%18.7f'%(obj['COSMOS_ID'], dist[cond][0])
            for name in cat_zurich.dtype.names[:-1]:
                line += '%18.7f'%cat_zurich[name][cond][0]
        else:
            line = '%18i%18.7f'%(obj['COSMOS_ID'], dist[cond][0])
            for name in cat_zurich.dtype.names[:-1]:
                line += '%18.7f'%cat_zurich[name][cond][0]
        line += '\n'
            #print(line)
        f.write(line)
    f.close()

def main():
    '''
    '''
    
    #-------#
    #  UDF  #
    #-------#
    
    #udfcat = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/mosaic_c042_e030_withz_o2_forkineanalysis_withMERGED.fits'
    #outcat = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/o2/input_fit_o2.txt'
    #filedir = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/o2/'
    
    #convert_udfcat_kincat(udfcat, outcat, filedir)
    
    #model = 'slp'
    #paths = ['/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/o2/']
    #for path in paths:
        #recap_kinem_params(path, basename='udf_mos_c042_e030_*', model=model, suff1='_xyi', suff2='_mclean5.0')
    
    #---------#
    #  UDF10  #
    #---------#
    
    #udfcat = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/mosaic_c042_e030_withz_o2_forkineanalysis_withMERGED.fits'
    #outcat = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/udf10/o2/input_fit_o2.txt'
    #filedir = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/udf10/o2/'
    
    #convert_udfcat_kincat(udfcat, outcat, filedir, prefix='udf10_c042_e031_')
    
#    model = 'slp'
#    paths = ['/media/bepinat/WD2To/data/Instruments/MUSE/analyse/UDF/data/camel/udf10/o2/']
#    for path in paths:
#        recap_kinem_params(path, basename='udf10_c042_e031_*', model=model, suff1='_xyi', suff2='_clean5.0')

    #--------#
    #  HDFS  #
    #--------#
    
    #galfile = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psf1/31_f814_comb_galfit.txt'
    #result = read_galfit_output(galfile)
    
    #galfile = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psf1/31_f814_comb_galfit.log'
    #result = read_galfit_log(galfile)
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffat/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffat_bulgefree/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfstack/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfstack_bulgefree/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psftinytim/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psftinytim_bulgefree/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #path = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffattrujillo_bulgefree/'
    #recap_morpho_params(path, basename='_f814_comb.log')
    
    #filen = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/tables/inclinations_comparison.txt'
    #inclination_histo(filen)
    
    #pathm = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/morpho/morpho_f814_psfmoffat/'
    #pathk = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/v031c/han2/'
    #res = convert_morphocat_kincat(pathm, pathk)
    #pathk = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/v031c/o3hb/'
    #res = convert_morphocat_kincat(pathm, pathk)
    #pathk = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/v031c/o2/'
    #res = convert_morphocat_kincat(pathm, pathk)
    
    #model = 'slp'
    #paths = ['/home/bepinat/Instruments/MUSE/analyse/Commissioning/v031c/o2/', '/home/bepinat/Instruments/MUSE/analyse/Commissioning/v031c/o3hb/', '/home/bepinat/Instruments/MUSE/analyse/Commissioning/v031c/han2/']
    #for path in paths:
        #recap_kinem_params(path, basename='hdfs_v*', model=model)
    
    #morphof = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/tables/recap_morpho_params.txt'
    #kinemf = '/home/bepinat/Instruments/MUSE/analyse/Commissioning/tables/recap_kinematics_parameters_best.txt'
    #compare_pa(morphof, kinemf)
    
    #---------------------------------------------------------#
    # Groups Snapsots and Deep fields with right nomenclature #
    #---------------------------------------------------------#
##    path = '/home/vabril/MUSE/Groups_Snapshots/'
##    pathmorpho = '/home/vabril/MUSE/Groups_Snapshots/morphology/'
#    
##    new paths for deep fields with right nomenclature and additional objects included
#    path = '/home/vabril/MUSE/Groups_Deep/'
#    pathmorpho = '/home/vabril/MUSE/Groups_Deep/morphology/'
##    pathmorpho = '/home/vabril/MUSE/Emission_region_CGr32-M1/galaxy_473/Morpho_modeling/'
#    
#    #groups = ['84-North']
#    #groups = ['23','26', '28', '61', '79','32-M1', '32-M2', '32-M3']
#    groups = ['32-M1']
#    #groups = ['32-M3'] #'79','32-M1', '32-M2', '32-M3'
#    #groups = ['30','34','84','114' ]
#    #groups = ['84',] 
#    
#    #section for one galaxy (not listed)
#    
##    pathm = pathmorpho
##    pathk = '/home/vabril/MUSE/Emission_region_CGr32-M1/galaxy_473/o2/'
##    pathgr = pathmorpho
#    
##    recap_morpho_params(pathgr, basename='_comb.log', outname='')
##    res = convert_morphocat_kincat(pathm, pathk)
##    print(res)
#    
#    for group in groups:
#        pathgr = pathmorpho + 'HST_CGr'+ group + '/'
#        #pathgr = pathmorpho + 'HST_CGr'+ group + '/CGr84_z0.68/'
#        #pathgr = pathmorpho + 'HST_CGr'+ group +'/CGr'+group+'A/'
#        #pathgr = pathmorpho + 'HST_CGr'+ group +'/CGr'+group+'B/'
#        
#        #recap_morpho_params(pathgr, basename='_comb.log', outname='CGr'+group+'A') #outname='CGr'+group+'B',  #outname='CGr'+group , outname='CGr'+group+'A'
#        recap_morpho_params(pathgr, basename='_comb.log', outname='')
#
#        gr = group.split('_')[0]
#        #pathm = pathmorpho + 'HST_CGr' + gr  +'/CGr'+group+'A/'  # A or B for group 84-North
#        
#        pathm = pathmorpho + 'HST_CGr' + gr + '/'
#        #pathm = pathmorpho + 'HST_CGr' + gr + '/CGr84_z0.68/'
##        pathk = path + 'CGr' + group + '/o2/'
#        #pathk = path + 'CGr' + group + '_best_seeing/o2/'
#        pathk = path + 'CGr' + group + '_deep/o2/'
#        #pathk = path + 'CGr' + group + '_mid_deep/CGr84_mdeep_A/o2/' #A/ or B/ for the two structures in field CGr84
#        #pathk = path + 'CGr' + group + '_snshot/o2/'
##        pathk = path + 'CGr' + group + 'B' + '/o2/' # A or B for group 84-North
#        
#        
#        res = convert_morphocat_kincat(pathm, pathk)
#        
##        pathk = path + 'CGr' + group + '/o3hb/'
#        #pathk = path + 'CGr' + group + '_best_seeing/o3hb/'
#        pathk = path + 'CGr' + group + '_deep/o3hb/'
##        pathk = path + 'CGr' + group + '_mid_deep/CGr84_mdeep_A/o3hb/'
#        #pathk = path + 'CGr' + group + '_snshot/o3hb/'
##        pathk = path + 'CGr' + group + 'B' + '/o3hb/' # A or B for group 84-North
#        
#        res = convert_morphocat_kincat(pathm, pathk)
##        pathk = path + 'CGr' + group + '/n2ha/'
##        res = convert_morphocat_kincat(pathm, pathk)
##        pathk = path + 'CGr' + group + '/s2/'
##        res = convert_morphocat_kincat(pathm, pathk)
#
#        print(pathm, pathk)
    
    ##--------#
    ## Groups #
    ##--------#
##    path = '/media/bepinat/WD2To/data/Instruments/MUSE/groups/morpho/acs_mosaic_2.0_'
##    pathkin = '/media/bepinat/WD2To/data/Instruments/MUSE/groups/kinematics/'
#    
#    path = '/home/vabril/MUSE/groups/morpho/acs_mosaic_2.0_'
#    pathkin = '/home/vabril/MUSE/groups/kinematics/'
#    #groups = ['gr28_shallow', 'gr32', 'gr83', 'gr116', 'gr28_best_seeing', 'gr28_deep']
#    groups = ['gr32', 'gr32_bs', 'gr83', 'gr116', 'gr28_best_seeing','gr28_shallow']
#    #groups = ['gr32_middeep', 'gr83_middeep', 'gr28_deep']
#    #groups = ['gr28_deep']
#    #groups = ['gr28_best_seeing', 'gr28_deep']
#    for group in groups:
#        #pathgr = path + group + '/'
#        
#        #recap_morpho_params(pathgr, basename='_f814_comb.log', outname='')
#        
#        gr = group.split('_')[0]
#        
#        pathgr = path + gr + '/'
#        recap_morpho_params(pathgr, basename='_f814_comb2.log', outname='')
#        pathm = path + gr + '/'
#        pathk = pathkin + group + '/o2/'
#        res = convert_morphocat_kincat(pathm, pathk)
#        pathk = pathkin + group + '/o3hb/'
#        res = convert_morphocat_kincat(pathm, pathk)
        
        ###################
        # PSF pour le GTO #
        ###################
        
#        path_psf = pathkin + group + '/psf_gto/'
#        gr = group.split('_')[0]
#        ipath = gr + '*5000/*'
#        recap_morphopsf_params(path_psf, basename='_sky_g.log', ipath=ipath, outname='_5000_g')
#        #recap_morphopsf_params(path_psf, basename='_sky_mg.log', ipath=ipath, outname='_5000_mg')
#        recap_morphopsf_params(path_psf, basename='_sky_m.log', ipath=ipath, outname='_5000_m')
#        ipath = gr + '*7000/*'
#        recap_morphopsf_params(path_psf, basename='_sky_g.log', ipath=ipath, outname='_7000_g')
#        #recap_morphopsf_params(path_psf, basename='_sky_mg.log', ipath=ipath, outname='_7000_mg')
#        recap_morphopsf_params(path_psf, basename='_sky_m.log', ipath=ipath, outname='_7000_m')
#        ipath = gr + '*9000/*'
#        recap_morphopsf_params(path_psf, basename='_sky_g.log', ipath=ipath, outname='_9000_g')
#        #recap_morphopsf_params(path_psf, basename='_sky_mg.log', ipath=ipath, outname='_9000_mg')
#        recap_morphopsf_params(path_psf, basename='_sky_m.log', ipath=ipath, outname='_9000_m')
#        
        ###################################
        ## PSF pour l'analyse cinématique #
        ##################################
#    pathkin = '/home/vabril/MUSE/Groups_Snapshots/'
#    pathkin = '/home/vabril/MUSE/Emission_region_CGr32-M1/galaxy_473/o2/'
#    pathkin = '/home/vabril/MUSE/groups/kinematics/'
    pathkin = '/home/wilfried/ST2/outputs/MUSE/'
     
    #groups = ['23','26', '28', '51', '61', '79','32-M1', '32-M2', '32-M3','84-NorthA','84-NorthB']
    #groups = ['23', '61']   #Activate the options for calculating the psf for n2ha and o2 bands
    #groups = ['26', '28',  '79','32-M1', '32-M2', '32-M3','84-NorthA'] 
    groups = ['30'] #'32-M1', '32-M2', '32-M3', '79'
    #groups = ['gr32_middeep','gr83_middeep', 'gr28_deep']
    #groups = ['gr32', 'gr83', 'gr116', 'gr28_best_seeing','gr28_shallow']
#    groups = ['51'] 
    
    for group in groups:
#        #path_star = pathkin + group + '/stars/'
#        path_star = pathkin + 'CGr' + group + '/stars/' # '/stars_A/' or '/stars_B/' for group 51
#        #path_star = pathkin + 'CGr' + group + '_deep/stars/'
#        #path_star = pathkin + 'CGr' + group + 'A/stars/' #Activate only for group 84-North
#        #path_star = pathkin + 'CGr' + group +'/'+'CGr84_mdeep_A/stars/'
#        
#        gr = group.split('_')[0]
#        print(gr)
#        ipath = 'CGr'+ gr + '*o2/*'
#        #ipath = gr + '*o2/*'
#        recap_morphopsf_params(path_star, basename='_comb_sky_g.log', ipath=ipath, outname='_o2_g')
#        recap_morphopsf_params(path_star, basename='_comb_sky_mg.log', ipath=ipath, outname='_o2_mg')
#        recap_morphopsf_params(path_star, basename='_comb_sky_m.log', ipath=ipath, outname='_o2_m')
#        ipath = 'CGr' +gr + '*o3hb/*'
#        #ipath = gr + '*o3hb/*'
#        recap_morphopsf_params(path_star, basename='_comb_sky_g.log', ipath=ipath, outname='_o3hb_g')
#        recap_morphopsf_params(path_star, basename='_comb_sky_mg.log', ipath=ipath, outname='_o3hb_mg')
#        recap_morphopsf_params(path_star, basename='_comb_sky_m.log', ipath=ipath, outname='_o3hb_m')
#        
#        ipath = 'CGr'+ gr + '*s2/*'
#        recap_morphopsf_params(path_star, basename='_comb_sky_g.log', ipath=ipath, outname='_s2_g')
#        recap_morphopsf_params(path_star, basename='_comb_sky_mg.log', ipath=ipath, outname='_s2_mg')
#        recap_morphopsf_params(path_star, basename='_comb_sky_m.log', ipath=ipath, outname='_s2_m')
#        ipath = 'CGr' +gr + '*n2ha/*'
#        recap_morphopsf_params(path_star, basename='_comb_sky_g.log', ipath=ipath, outname='_n2ha_g')
#        recap_morphopsf_params(path_star, basename='_comb_sky_mg.log', ipath=ipath, outname='_n2ha_mg')
#        recap_morphopsf_params(path_star, basename='_comb_sky_m.log', ipath=ipath, outname='_n2ha_m')

        gr = group.split('_')[0]
        print(group)
#        pathk = pathkin + 'CGr'+ group + '/o2/'
        pathk = pathkin
        pathk = pathkin + 'CGr'+ group + '_d/o2/'
        #pathk = pathkin + 'CGr'+ group + '_mid_deep/CGr84_mdeep_A/o2/'
        basename = 'CGr'+ gr + '_deep_67_o2'
#        basename = 'CGr'+ gr + '*_*' #Activate for Group 51 Snapshots to include 51B and 51 A in the same output
        print(basename)
        recap_kinem_params(pathk, basename=basename)
        
        
#        pathk = pathkin +'CGr'+ group + '/o3hb/'
#        pathk = pathkin +'CGr'+ group + '_deep/o3hb/'
        #pathk = pathkin +'CGr'+ group + '_mid_deep/CGr84_mdeep_A/o3hb/'
#        recap_kinem_params(pathk, basename=basename)
#        
#        pathk = pathkin + 'CGr'+ group + '/s2/'
#        basename = 'CGr'+ gr + '_*'
#        recap_kinem_params(pathk, basename=basename)
#        pathk = pathkin +'CGr'+ group + '/n2ha/'
#        recap_kinem_params(pathk, basename=basename)
#        
#        gr = group.split('_')[0]
#        pathk = pathkin + group + '/o2/'
#        basename = gr + '_*'
#        recap_kinem_params(pathk, basename=basename)
#        pathk = pathkin +group+ '/o3hb/'
#        recap_kinem_params(pathk, basename=basename)
    
    #path = '/media/bepinat/WD2To/data/Instruments/MUSE/groups/'
    #morpho_file = 'catalogues_irsa/cosmos_morph_zurich_1.0_lisible.tbl'
    #cat_zurich = ascii.read(path + morpho_file)
    #radius = 1.5  # arcsec
    #groups = ['gr28', 'gr32', 'gr83', 'gr116']
    #for group in groups:
        #group_file = 'morpho/C' + group.upper() + '_zmeasurement_reconciled.txt'
        #cat_group = ascii.read(path + group_file)
        #outname = path + 'catalogues_irsa/cosmos_morph_zurich_' + group + '.tbl'
        #extract_lines(cat_zurich, cat_group, outname, radius=radius)




if __name__ == "__main__":
    main()
    
