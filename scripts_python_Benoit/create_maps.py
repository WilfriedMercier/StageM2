#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import astropy.wcs as wcs
from astropy.io import fits
from astropy.io import ascii
import numpy as np
import astropy.coordinates as coord
import astropy.units as u
import copy

from scipy.ndimage.filters import gaussian_filter

import matplotlib as mpl
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
from matplotlib import cm
from matplotlib import rc
from matplotlib import gridspec
mpl.style.use('classic')
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
rc('figure', figsize=(6.0, 4.5))

from reproject import reproject_interp

import ipdb


def hst_interp_wcs(hsthdu, hsthdr, fluxhdr):
    hsthdr2 = copy.deepcopy(fluxhdr)
    fac = np.sqrt((fluxhdr['CD1_1'] ** 2 + fluxhdr['CD1_2'] ** 2) / (hsthdr['CD1_1'] ** 2 + hsthdr['CD1_2'] ** 2))
    #print(fac)
    hsthdr2['CD1_1'] = fluxhdr['CD1_1'] / fac
    hsthdr2['CD1_2'] = fluxhdr['CD1_2'] / fac
    hsthdr2['CD2_1'] = fluxhdr['CD2_1'] / fac
    hsthdr2['CD2_2'] = fluxhdr['CD2_2'] / fac
    hsthdr2['CRVAL1'] = fluxhdr['CRVAL1']
    hsthdr2['CRVAL2'] = fluxhdr['CRVAL2']
    hsthdr2['CRPIX1'] = fluxhdr['CRPIX1'] * fac
    hsthdr2['CRPIX2'] = fluxhdr['CRPIX2'] * fac
    hsthdr2['NAXIS1'] = fluxhdr['NAXIS1'] * fac
    hsthdr2['NAXIS2'] = fluxhdr['NAXIS2'] * fac
    hstmap, footprint = reproject_interp(hsthdu[0], hsthdr2)
    hstwcs = wcs.WCS(hsthdr2)
    #print(hstmap.shape)
    return hstmap, hsthdr2, hstwcs, fac


def plot_hst(gs, hstwcs, hstmap, fac, xmin, xmax, ymin, ymax, vmin, vmax, ra0, dec0):
    
    xh, yh = hstwcs.wcs_world2pix(ra0, dec0, 0)
    
    xch = xh[0]
    ych = yh[0]
    
    xmah1 = xh[1]
    xmah2 = xh[2]
    ymah1 = yh[1]
    ymah2 = yh[2]
    
    axhstmap = plt.subplot(gs, projection=hstwcs)
    #axhstmap = fig.add_subplot(3, 5, 1, projection=fluxwcs)
    
    #rav, decv = fluxwcs.wcs_pix2world([xmin, xmax], [ymin, ymax], 0)
    #xhst, yhst = hstwcs.wcs_world2pix(rac, decv, 0)
    
    #axhstmap.set_xlim(xhst[0], xhst[1])
    #axhstmap.set_ylim(yhst[0], yhst[1])
    
    delta = np.round(fac/2)
    xminh = -0.5+delta
    xmaxh = hstmap.shape[1]-0.5+delta
    yminh = -0.5+delta
    ymaxh = hstmap.shape[0]-0.5+delta
    
    axhstmap.set_xlim(xminh, xmaxh)
    axhstmap.set_ylim(yminh, ymaxh)
    
    hstmap[hstmap < vmin] = vmin
    
    imfluxmap = axhstmap.imshow(hstmap, norm=LogNorm(vmin=vmin, vmax=vmax), cmap=plt.cm.gray_r, origin='lower', interpolation='nearest')  # Logarithmic scale
    ra = axhstmap.coords['ra']
    dec = axhstmap.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    #cb = plt.colorbar(imfluxmap, orientation='horizontal', fraction=0.04, pad=0.)
    #cbarC.set_label(---)
    
    plt.plot([xch], [ych], 'w+', mew=1)
    plt.plot([xmah1, xmah2], [ymah1, ymah2], 'w')
    
    
    #axhstmap = fig.add_subplot(3, 5, 1, projection=hstwcs)
    #axhstmod = fig.add_subplot(3, 5, 6, projection=hstwcs)
    #axhstres = fig.add_subplot(3, 5, 11, projection=hstwcs)


def paper_map(hst, flux, snr, vf, vfm, vfr, sig, sigm, sigr, name, z, xc, yc, vsys, pa, rc, lsf, psf, pathout, deltapa=-42):
    
    # File reading
    hst435hdu = fits.open(hst['f435w'])
    hst814hdu = fits.open(hst['f814w'])
    hst160hdu = fits.open(hst['f160w'])
    
    hst814map = hst814hdu[0].data
    hst814hdr = hst814hdu[0].header
    
    hst435map = hst435hdu[0].data
    hst435hdr = hst435hdu[0].header
    
    hst160map = hst160hdu[0].data
    hst160hdr = hst160hdu[0].header
    
    #sigma = 3. # this depends on how noisy your data is, play with it!
    #data = gaussian_filter(hstmap, sigma)
    
    #hstmap = hsthdu[1].data
    #hstmod = hsthdu[2].data
    #hstres = hsthdu[3].data
    #hsthdr = hsthdu[1].header
    
    fluxhdu = fits.open(flux[0])
    fluxmap = fluxhdu[0].data * 0
    fluxhdr = fluxhdu[0].header
    for file in flux:
        temphdu = fits.open(file)
        fluxmap += temphdu[0].data
    
    snrhdu = fits.open(snr)
    snrmap = snrhdu[0].data
    snrhdr = snrhdu[0].header
    
    vfhdu = fits.open(vf)
    vfmap = vfhdu[0].data
    vfhdr = vfhdu[0].header
    
    vfmhdu = fits.open(vfm)
    vfmmap = vfmhdu[0].data
    vfmhdr = vfmhdu[0].header
    
    vfrhdu = fits.open(vfr)
    vfrmap = vfrhdu[0].data
    vfrhdr = vfrhdu[0].header
    
    sighdu = fits.open(sig)
    sigmap = sighdu[0].data
    sighdr = sighdu[0].header
    
    sigmhdu = fits.open(sigm)
    sigmmap = sigmhdu[0].data
    sigmhdr = sigmhdu[0].header
    
    sigrhdu = fits.open(sigr)
    sigrmap = sigrhdu[0].data
    sigrhdr = sigrhdu[0].header
    
    # Defining wcs
    fluxhdr['WCSAXES'] = 2
    vfhdr['WCSAXES'] = 2
    sighdr['WCSAXES'] = 2
    snrhdr['WCSAXES'] = 2
    hst435hdr['WCSAXES'] = 2
    hst814hdr['WCSAXES'] = 2
    hst160hdr['WCSAXES'] = 2
    fluxwcs = wcs.WCS(fluxhdr)
    vfwcs = wcs.WCS(vfhdr)
    sigwcs = wcs.WCS(sighdr)
    snrwcs = wcs.WCS(snrhdr)
    
    hst435map, hst435hdr2, hst435wcs, fac435 = hst_interp_wcs(hst435hdu, hst435hdr, fluxhdr)
    hst814map, hst814hdr2, hst814wcs, fac814 = hst_interp_wcs(hst814hdu, hst814hdr, fluxhdr)
    hst160map, hst160hdr2, hst160wcs, fac160 = hst_interp_wcs(hst160hdu, hst160hdr, fluxhdr)
    #hstwcs = wcs.WCS(hsthdr)
    
    sigma = 3.  # this depends on how noisy your data is, play with it!
    data = gaussian_filter(hst160map, sigma)
    
    # positions of center and major axis
    
    xma1 = xc + rc * np.sin(np.radians(pa+deltapa))
    yma1 = yc - rc * np.cos(np.radians(pa+deltapa))
    xma2 = xc - rc * np.sin(np.radians(pa+deltapa))
    yma2 = yc + rc * np.cos(np.radians(pa+deltapa))
    
    #print(xma1, xma2, yma1, yma2, pa)
    
    ra0, dec0 = fluxwcs.wcs_pix2world([xc, xma1, xma2], [yc, yma1, yma2], 0)
    
    #rac, decc = fluxwcs.wcs_pix2world(xc, yc, 0)
    #xch, ych= hstwcs.wcs_world2pix(rac, decc, 0)
    
    xmin = -0.5
    xmax = vfmap.shape[1]-0.5
    ymin = -0.5
    ymax = vfmap.shape[0]-0.5
    
    # Initializing figure
    fig = plt.figure(figsize=(10.08, 6.27))  # sums of width in gridspec + intervals in subplots_adjust, eventually multiply by some factor
    
    plt.figtext(0.5, 0.25, 'ID ' + name.split('_o2')[0].split('_')[-1])
    plt.figtext(0.5, 0.2, r'$z={:.5f}$'.format(z))
    
    gs = gridspec.GridSpec(4, 5, width_ratios=[1,1,1,1,1], height_ratios=[1,1,1,0.075])
    
    # Valeurs F160W
    vmin160 = 2e-4
    vmax160 = 2e0
    
    # Valeurs F814W
    vmin814 = 1e-4
    vmax814 = 1e-1
    
    # Valeurs F435W
    vmin435 = 1e-4
    vmax435 = 1e-1
    
    plot_hst(gs[0], hst160wcs, hst160map, fac160, xmin, xmax, ymin, ymax, vmin160, vmax160, ra0, dec0)
    plot_hst(gs[5], hst814wcs, hst814map, fac814, xmin, xmax, ymin, ymax, vmin814, vmax814, ra0, dec0)
    plot_hst(gs[10], hst435wcs, hst435map, fac435, xmin, xmax, ymin, ymax, vmin435, vmax435, ra0, dec0)
    
    # Flux map

    axflux = plt.subplot(gs[1], projection=fluxwcs)
    #axflux = fig.add_subplot(3, 5, 2, projection=fluxwcs)
    
    fluxmap[fluxmap == 0] = np.nan
    vmin = np.nanmin(fluxmap)
    vmax = np.nanmax(fluxmap)
    
    axflux.set_xlim(xmin, xmax)
    axflux.set_ylim(ymin, ymax)
    
    imfluxmap = axflux.imshow(fluxmap, norm=LogNorm(vmin=vmin, vmax=vmax), cmap=plt.cm.gray_r, origin='lower', interpolation='nearest')  # Logarithmic scale
    ra = axflux.coords['ra']
    dec = axflux.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    #cb = plt.colorbar(imfluxmap, orientation='horizontal', fraction=0.04, pad=0.)
    #cbarC.set_label(---)
    
    # SNR map
    
    axsnr = plt.subplot(gs[2], projection=snrwcs)
    #axsnr = fig.add_subplot(3, 5, 3, projection=snrwcs)
    
    vmin = np.nanmin(snrmap)
    vmax = np.nanmax(snrmap)
    
    axsnr.set_xlim(xmin, xmax)
    axsnr.set_ylim(ymin, ymax)
    
    imsnrmap = axsnr.imshow(snrmap, norm=LogNorm(vmin=vmin, vmax=vmax), cmap=plt.cm.gray_r, origin='lower', interpolation='nearest')  # Logarithmic scale
    ra = axsnr.coords['ra']
    dec = axsnr.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    #cb = plt.colorbar(imsnrmap, orientation='horizontal', fraction=0.04, pad=0.)
    #cbarC.set_label(---)

    
    # VF map

    axvf = plt.subplot(gs[3], projection=vfwcs)
    #axvf = fig.add_subplot(3, 5, 4, projection=vfwcs)
    val = np.max([np.abs(np.nanmax(vfmmap - vsys)), np.abs(np.nanmin(vfmmap - vsys))])
    
    vmin = -val
    vmax = val
    
    axvf.set_xlim(xmin, xmax)
    axvf.set_ylim(ymin, ymax)
    
    imvf = axvf.imshow(vfmap - vsys, vmin=vmin, vmax=vmax, cmap=plt.cm.seismic, origin='lower', interpolation='nearest')
    ra = axvf.coords['ra']
    dec = axvf.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    trans = axvf.get_transform(hst160wcs)
    axvf.contour(data, levels=np.logspace(-3, 0., 7), colors='k', transform=trans)
    
    plt.plot([xc], [yc], 'k+', mew=1)
    plt.plot([xma1, xma2], [yma1, yma2], 'k')
    
    circle = Circle(xy=(psf*1.25, psf*1.25), radius=psf, edgecolor='grey', fc='0.8', lw=0.5)
    axvf.add_patch(circle)
    #Ellipse(xy=(xreg[i], yreg[i]), width=2 * rereg[i] * 2.2, height=2 * rereg[i] * 2.2 * np.cos(np.radians(cat2['INC'][i])), angle=90+cat2['PA'][i], edgecolor=col, fc='None', lw=lw, ls=ls, zorder=zorder)
    #plt.plot([psf*1.5], [psf*1.5], 'ok', ms=psf)
    
    # VF model
    
    axvfm = plt.subplot(gs[8], projection=vfwcs)
    #axvfm = fig.add_subplot(3, 5, 9, projection=vfwcs)
    
    axvfm.set_xlim(xmin, xmax)
    axvfm.set_ylim(ymin, ymax)
    
    imvfm = axvfm.imshow(vfmmap - vsys, vmin=vmin, vmax=vmax, cmap=plt.cm.seismic, origin='lower', interpolation='nearest')
    ra = axvfm.coords['ra']
    dec = axvfm.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    
    # VF residuals
    
    axvfr = plt.subplot(gs[13], projection=vfwcs)
    #axvfr = fig.add_subplot(3, 5, 14, projection=vfwcs)
    
    axvfr.set_xlim(xmin, xmax)
    axvfr.set_ylim(ymin, ymax)
    
    imvfr = axvfr.imshow(vfrmap, vmin=vmin, vmax=vmax, cmap=plt.cm.seismic, origin='lower', interpolation='nearest')
    ra = axvfr.coords['ra']
    dec = axvfr.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    axcbvf = plt.subplot(gs[18])
    
    pos_axvfr = axvfr.get_position()
    pos_axcbvf = axcbvf.get_position()
    axcbvf.set_position([pos_axvfr.bounds[0], pos_axcbvf.bounds[1], pos_axvfr.width, pos_axcbvf.height])
    
    stepvf = np.ceil(vmax/3/10) * 10
    ticksvf = np.arange(-4 * stepvf,vmax + 1, stepvf)
    
    cb = plt.colorbar(imvfm, orientation='horizontal', cax=axcbvf, ticks=ticksvf)
    cb.set_label(r'km s$^{-1}$')


    # Velocity dispersion

    #lsff = np.sqrt(sigmap**2-sigmmap**2-sigrmap**2)
    #cond = np.logical_not(np.isnan(lsff))
    #print(np.mean(lsff[cond]), np.std(lsff[cond]))

    axsig = plt.subplot(gs[4], projection=sigwcs)
    #axsig = fig.add_subplot(3, 5, 5, projection=sigwcs)
    
    vmin = 0
    cond = np.logical_not(np.isnan(np.sqrt(sigmap**2 - lsf**2)))
    ss = np.sqrt(sigmap**2 - lsf**2)[cond]
    ss.sort()
    vmax = ss[-3]
    #print(vmax)
    
    axsig.set_xlim(xmin, xmax)
    axsig.set_ylim(ymin, ymax)
    
    imsig = axsig.imshow(np.sqrt(sigmap**2 - lsf**2), vmin=vmin, vmax=vmax, cmap=plt.cm.CMRmap, origin='lower', interpolation='nearest')
    ra = axsig.coords['ra']
    dec = axsig.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    

    # Velocity dispersion model
    
    axsigm = plt.subplot(gs[9], projection=sigwcs)
    #axsigm = fig.add_subplot(3, 5, 10, projection=sigwcs)
    
    axsigm.set_xlim(xmin, xmax)
    axsigm.set_ylim(ymin, ymax)
    
    imsigm = axsigm.imshow(np.sqrt(sigmmap**2), vmin=vmin, vmax=vmax, cmap=plt.cm.CMRmap, origin='lower', interpolation='nearest')
    ra = axsigm.coords['ra']
    dec = axsigm.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    
    # Velocity dispersion residuals
    
    axsigr = plt.subplot(gs[14], projection=sigwcs)
    #axsigr = fig.add_subplot(3, 5, 15, projection=sigwcs)
    
    axsigr.set_xlim(xmin, xmax)
    axsigr.set_ylim(ymin, ymax)
    
    imsigr = axsigr.imshow(sigrmap, vmin=vmin, vmax=vmax, cmap=plt.cm.CMRmap, origin='lower', interpolation='nearest')
    ra = axsigr.coords['ra']
    dec = axsigr.coords['dec']
    #dec.set_ticklabel(rotation=90)
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    #dec.set_axislabel('Declination (J2000)')
    #ra.set_axislabel('Right ascention (J2000)', minpad=0.4)
    
    axcbsig = plt.subplot(gs[19])
    
    pos_axsigr = axsigr.get_position()
    pos_axcbsig = axcbsig.get_position()
    axcbsig.set_position([pos_axsigr.bounds[0], pos_axcbsig.bounds[1], pos_axsigr.width, pos_axcbsig.height])
    
    stepsig = np.round(vmax/4/5) * 5
    tickssig = np.arange(0,vmax + 1, stepsig)
    cb = plt.colorbar(imsigm, orientation='horizontal', cax=axcbsig, ticks=tickssig)
    cb.set_label(r'km s$^{-1}$')
    
    #cb = plt.colorbar(imsigm, orientation='horizontal', pad=0.)
    
    fig.subplots_adjust(hspace=0.05, wspace=0.05, left=0.18, right=0.98, top=0.9, bottom=0.1)  # spacing left+right vs top+bottom must be identical to keep good aspect, otherwize chanfe figure size
    
    # Additional map?
    
    #ax32 = fig.add_subplot(3, 5, 6, projection=fluxwcs)
    fig.savefig(pathout + name.split('_Z_')[0] + '.pdf', bbox_inches='tight')
    plt.close()


def main():
    '''
    '''
    pathhst = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/hst/'
    #hst = pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f814w_v1_sci.fits'
    #hst = pathhst + 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits'
    
    hst = {'f105w0':pathhst + 'hlsp_hudf12_hst_wfc3ir_udfmain_f105w_v1.0_drz.fits',
           'f125w0':pathhst + 'hlsp_hudf12_hst_wfc3ir_udfmain_f125w_v1.0_drz.fits',
           'f140w0':pathhst + 'hlsp_hudf12_hst_wfc3ir_udfmain_f140w_v1.0_drz.fits',
           'f160w0':pathhst + 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits',
           'f105w':pathhst + 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f105w_v1_sci.fits',
           'f125w':pathhst + 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f125w_v1_sci.fits',
           'f140w':pathhst + 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits',
           'f160w':pathhst + 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f160w_v1_sci.fits',
           'f435w':pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f435w_v1_sci.fits',
           'f606w':pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f606w_v1_sci.fits',
           'f775w':pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f775w_v1_sci.fits',
           'f814w':pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f814w_v1_sci.fits',
           'f850lp':pathhst + 'hlsp_xdf_hst_acswfc-30mas_hudf_f850lp_v1_sci.fits'}
    
    line = 'o2'
    #name = 'udf_mos_c042_e030_1_o2_Z_0.621992'
    
    mosaic = True
    #mosaic = False
    
    if mosaic:
        path = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/camel/'
        #file_recap = path + line + '/recap_kinematics_parameters_2_slp_xyi_mclean5.0_work.txt'
        file_recap = path + line + '/recap_kinematics_parameters_2_slp_xyi_mclean5.0.txt'
        file_input = path + line + '/input_fit_o2_v3.txt'
        file_input = path + line + '/input_fit_o2_v4.txt'
        pathout = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/pdf/udfmosaic/'
        deltapa = -42
    else:
        path = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/camel/udf10/'
        file_recap = path + line + '/recap_kinematics_parameters_2_slp_xyi_mclean5.0.txt'
        file_input = path + line + '/input_fit_o2_v2.txt'
        file_input = path + line + '/input_fit_o2_v3.txt'
        pathout = '/home/bepinat/Instruments/MUSE/analyse/UDF/data/pdf/udf10/'
        deltapa = 0
        
    cat_recap = ascii.read(file_recap)
    cat_input = ascii.read(file_input)
    
    for obj in cat_recap:
        name = obj['ID']
        print(name)
        z = float(name.split('_Z_')[-1])
        ind, = np.where(cat_input['gal'] == name)
        psf = cat_input['psfx'][ind[0]]  # pixels MUSE (FWHM)
        smooth = cat_input['smooth'][ind[0]]  # pixels MUSE (FWHM)
        psff = np.sqrt(psf**2 + smooth**2)
        #print(psff)
        lsf = cat_input['psfz'][ind[0]]  # km/s  (sigma)
        xc = obj['X']
        yc = obj['Y']
        pa = obj['PA']  # par rapport au nord dans l'image MUSE
        vsys = obj['VS']
        
        rc = obj['RLAST']  # pixels
        #rc = obj[]  # a prendre dans la morpho
        
        #flux = glob.glob(path + line + '/' + name + '/' + name + '_ssmooth_flux_common_mclean5.0.fits')
        flux = list(set(glob.glob(path + line + '/' + name + '/' + name + '_ssmooth_flux_common_*.fits')) - set(glob.glob(path + line + '/' + name + '/' + name + '_ssmooth_flux_common_*clean5.0.fits')))
    
        #snr = path + line + '/' + name + '/' + name + '_ssmooth_snr_common_mclean5.0.fits'
        snr = path + line + '/' + name + '/' + name + '_ssmooth_snr_common.fits'
    
        vf = path + line + '/' + name + '/' + name + '_ssmooth_vel_common_mclean5.0.fits'
        vfm = path + line + '/' + name + '/' + name + '_modv_slp_xyi_mclean5.0.fits'
        vfr = path + line + '/' + name + '/' + name + '_resv_slp_xyi_mclean5.0.fits'
    
        sig = path + line + '/' + name + '/' + name + '_ssmooth_disp_common_mclean5.0.fits'
        sigm = path + line + '/' + name + '/' + name + '_modd_slp_xyi_mclean5.0.fits'
        sigr = path + line + '/' + name + '/' + name + '_resd_slp_xyi_mclean5.0.fits'
    
        paper_map(hst, flux, snr, vf, vfm, vfr, sig, sigm, sigr, name, z, xc, yc, vsys, pa, rc, lsf, psff, pathout, deltapa=deltapa)
        #return

if __name__ == "__main__":
    main()
