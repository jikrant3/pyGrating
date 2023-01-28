import os
import numpy as np

from scipy import interpolate
from scipy.signal import find_peaks

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

from astropy import units as u
from astropy import constants as const
from astropy.convolution import convolve, Gaussian2DKernel, Gaussian1DKernel
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.nddata import StdDevUncertainty
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, QTable

from specutils.spectra import Spectrum1D

######################################################################
######################################################################

# Calibration parameters of FUV Grating1 -1 order # Dewangan+2022 priv. communication
dict_FUV_G1_M1 ={'DETECTOR':'FUV', 'FILTERID':'Grating1', 'ORDER':-1,
                 'PIX_MIN':-322*u.pix, 'PIX_MAX':-212*u.pix, 
                 'C0':-18.0335*u.Angstrom,'C1':-5.8332*u.Angstrom/u.pix,     # -5.5690588440719395
                 'EFFECTIVE_AREA_FILE':'data/fuv_grating1m1_effarea_2021.fits'}

# Calibration parameters of FUV Grating1 -2 order # Dewangan+2022 priv. communication
dict_FUV_G1_M2 ={'DETECTOR':'FUV', 'FILTERID':'Grating1', 'ORDER':-2,
                 'PIX_MIN':-629*u.pix, 'PIX_MAX':-413*u.pix, 
                 'C0':45.52709688716913*u.Angstrom,'C1':-2.789558814560311*u.Angstrom/u.pix,
                 'EFFECTIVE_AREA_FILE':'data/fuv_grating1m2_effarea_9nov22_2.fits'}

# Calibration parameters of FUV Grating2 -2 order # Dewangan+2022 priv. communication
dict_FUV_G2_M2 ={'DETECTOR':'FUV', 'FILTERID':'Grating2', 'ORDER':-2,
                 'PIX_MIN':-626*u.pix, 'PIX_MAX':-430*u.pix, 
                 'C0':31.2*u.Angstrom,'C1':-2.812*u.Angstrom/u.pix,
                 'EFFECTIVE_AREA_FILE':'data/fuv_grating2m2_effarea_12nov22_2.fits'}

# Calibration parameters of NUV Grating2 -1 order  # Tandon+2020
dict_NUV_G_M1  ={'DETECTOR':'NUV', 'FILTERID':'Grating', 'ORDER':-1,
                 'PIX_MIN':-545*u.pix, 'PIX_MAX':-336*u.pix, 
                 'C0':18.1*u.Angstrom,'C1':-5.5868*u.Angstrom/u.pix,
                 'EFFECTIVE_AREA_FILE':'data/nuv_gratingm1_effarea_2020.fits'}

######################################################################
######################################################################

class GratingImage:
    '''
    Class to visualise CCDLAB reduced UVIT Grating image and identifying sources
    '''
    def __init__(self, file_name, crop=True, verbose=True):
        '''
        Initialises the image.
        
        Parameters
        ----------
        file_name : str
            Name of the Grating file. The spectral dispersion direction should should be exactly horizantal.
        crop : bool
            Cropping in to central circular region of the image.
        '''
        self.file_name    = file_name
        hdu               = fits.open(self.file_name)
        header            = hdu[0].header
        self.data         = hdu[0].data
        self.DETECTOR     = header['DETECTOR']
        self.FILTERID     = header['FILTERID']
        self.EXPTIME      = header['RDCDTIME']*u.s
        self.OBJECT       = header['OBJECT']
        self.OBS_ID       = header['OBS_ID']

        if verbose: print('%s\n\tOBS_ID   = %s\n\tDETECTOR = %s\n\tFILTERID = %s\n\tRDCDTIME = %s'%(self.OBJECT, self.OBS_ID, self.DETECTOR, self.FILTERID, self.EXPTIME))
        
        # Smoothing data for plotting
        self.data_smooth     = convolve(self.data, Gaussian2DKernel(3))
        if crop: self.crop_to_circle()

    def crop_to_circle(self, rmax=2000):
        '''
        Crop to central circular region.
        UVIT pixel scale: 1 pixel = 0.41684142 arcsec --> 2000 pix ~ 13.9 arcmin
        
        Parameters
        ----------
        rmax : int
            Radius of the usable region in pixels

        Returns
        -------
        data: numpy.ndarray
            Cropped image data 
        data_smooth: numpy.ndarray
            Cropped smooth image data 
        '''
        lx, ly                 = self.data.shape
        X, Y                   = np.ogrid[0:lx, 0:ly]
        mask                   = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > rmax**2
        self.data[mask]        = np.nan
        self.data_smooth[mask] = np.nan

    def plot(self, ax=None):
        '''
        Plot the image. Option to show detected sources if available.
        '''
        if ax==None: fig, ax = plt.subplots(figsize=(6,6))
        median               = np.median(self.data_smooth[self.data_smooth>0])
        ax.imshow(self.data_smooth, norm=PowerNorm(0.1, vmin=median/2, vmax=median*20), origin='lower')
        ax.set(xlabel='x [pix]', ylabel='y [pix]', title='%s %s %s'%(self.OBJECT, self.DETECTOR, self.FILTERID))
        if hasattr(self, 'ypix_list'):
            ax.scatter(np.zeros(len(self.ypix_list)), self.ypix_list, marker=1, c='r')
        
    def identify_sources(self, prominence=3., plot=True):
        '''
        Detects sources by averaging across dispersion axis and identifies peaks based on prominance. 
        
        Parameters
        ----------
        prominence : float
            Multiplication factor for measuring the prominance for the peak detection. [prominence * standard_deviation]
        
        Returns
        -------
        ypix_list: list
            y pixels of all detected sources.
        '''
        row_wise_counts   = np.log10(np.nanmean(self.data_smooth, axis=1))
        mean, median, std = sigma_clipped_stats(data=row_wise_counts, sigma=3.0,maxiters=20)
        row_wise_counts[row_wise_counts<mean-5*std] = mean-5*std        # Replaces extremly low counts 
        self.ypix_list, _ = find_peaks(row_wise_counts, height=mean-3*std, prominence=prominence*std)
        if plot:
            fig, ax = plt.subplots(figsize=(20,5), ncols=2)
            [axi.set_axis_off() for axi in ax.ravel()]
            ax[0] = fig.add_axes([0.05, 0.06, 0.22, 0.88])
            ax[1] = fig.add_axes([0.32, 0.06, 0.65, 0.88])
            
            self.plot(ax=ax[0])
            
            ax[1].plot(row_wise_counts)
            ax[1].axhline(mean, c='r', ls='--')
            ax[1].axhline(mean+std, c='red', ls=':')
            ax[1].axhline(mean-std, c='red', ls=':')
            ax[1].plot(self.ypix_list, row_wise_counts[self.ypix_list], "rx", ms=10)
            ax[1].set(xlabel='y [pix]', ylabel='Average log(counts)')
            for i in range(len(self.ypix_list)):
                ax[1].annotate(self.ypix_list[i], (self.ypix_list[i], row_wise_counts[self.ypix_list][i] + std), 
                            rotation=90, ha='center', va='bottom')
            ax[1].grid()

######################################################################
######################################################################
            
class GratingSpectrum:
    '''
    Class to get spectrum from UVIT Grating image
    '''
    def __init__(self, image, YPIX, YPIX_bkg, xpix=None, cross_disp_buffer=25):
        '''
        Initialise the spectrum object based on y value
        
        Parameters
        ----------
        image : GratingImage
            Grating image
        YPIX : int
            y value of the source[pix]
        YPIX_bkg : int
            y value of the background region[pix]
        xpix : int
            x value of the [pix]
        cross_disp_buffer : int
            Buffer to get 2d spectrum. The 2d spetrum is made using region of "YPIX-cross_disp_buffer" to "YPIX+cross_disp_buffer" [pix]
        '''
        self.YPIX              = YPIX
        self.YPIX_bkg          = YPIX_bkg
        self.xpix              = xpix
        self.cross_disp_buffer = cross_disp_buffer
        # Inherit metadata from GratingImage
        self.EXPTIME           = image.EXPTIME
        self.DETECTOR          = image.DETECTOR
        self.FILTERID          = image.FILTERID
        self.OBJECT            = image.OBJECT
        self.OBS_ID            = image.OBS_ID
        
        # Get 2d spectrum
        self.spec_2d_bkg = image.data[YPIX_bkg-cross_disp_buffer:YPIX_bkg+cross_disp_buffer+1,:]*u.ct
        self.spec_2d     = image.data[YPIX-cross_disp_buffer:YPIX+cross_disp_buffer+1,:]*u.ct
        # Get 1d spectrum and errors
        self.spec_1d     = np.sum(self.spec_2d, axis=0)/u.pix
        self.fractional_error_Poisson = self.spec_1d.value**-0.5
        self.spec_1d_bkg = np.sum(self.spec_2d_bkg, axis=0)/u.pix
        self.spec_1d_net = self.spec_1d - self.spec_1d_bkg
        self.calculate_xpix()
        # pixels_from_zeroth
        pixels_obs       = np.arange(0.,len(self.spec_1d_net),1.)
        self.pixels      = (pixels_obs-self.xpix)*u.pix
        
        # get orders and wavelength/flux calibration
        if (self.DETECTOR=='FUV') & (self.FILTERID=='Grating1'):
            self.spectrum_m1, self.spectrum_cps_m1, self.spectrum_ct_m1 = self.get_calibrated_spectrum(dict_FUV_G1_M1)
            self.spectrum_m2, self.spectrum_cps_m2, self.spectrum_ct_m2 = self.get_calibrated_spectrum(dict_FUV_G1_M2)
        if (self.DETECTOR=='FUV') & (self.FILTERID=='Grating2'):
            self.spectrum_m2, self.spectrum_cps_m2, self.spectrum_ct_m2 = self.get_calibrated_spectrum(dict_FUV_G2_M2)
        if (self.DETECTOR=='NUV') & (self.FILTERID=='Grating'):
            self.spectrum_m1, self.spectrum_cps_m1, self.spectrum_ct_m1 = self.get_calibrated_spectrum(dict_NUV_G_M1)
            
    def calculate_xpix(self, plot=False):
        '''
        Calculate the x pixel of the 0th order
        
        Returns
        -------
        xpix : float
            Moffat fitted x value of the 0th order [pix]
        '''
        if self.xpix!=None:
            '''
            Finds the maximum value near the initial guess (Default buffer of 100 pixel around the guess).
            '''
            buffer               = 100 # half width of region around the initial guess
            order_0_neighborhood = self.spec_1d_net[self.xpix-buffer:self.xpix+buffer]
            self.pix             = self.xpix-buffer+order_0_neighborhood.argmax()

        if self.xpix==None:
            '''
            If the initial guess is not given: 
            Finds a peak with high prominance (>1/2th of the max count) and small width (<20 pix).
            '''
            peaks, _        = find_peaks(self.spec_1d_net, prominence=np.nanmax(self.spec_1d_net.value)/2.,
                                         height=0, width=(2,20))
            if len(peaks)==1: self.xpix = peaks[0]
            else:
                '''
                If more than 1 peak is found, requires a initial guess
                '''
                fig, ax = plt.subplots(figsize =(20,6), ncols=1, nrows=3)
                [axi.set_axis_off() for axi in ax.ravel()]
                ax[0] = fig.add_axes([0.05, 0.85, 0.90, 0.1])
                ax[1] = fig.add_axes([0.05, 0.75, 0.90, 0.1])
                ax[2] = fig.add_axes([0.05, 0.05, 0.90, 0.7])

                ax[0].imshow(self.spec_2d_bkg, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
                ax[1].imshow(self.spec_2d, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
                
                ax[2].plot(self.spec_1d_net, c='b')
                ax[2].plot(peaks, self.spec_1d_net[peaks], "rx", ms=10)
                ax[2].grid()
                ax[2].set_xlim(ax[0].get_xlim())
                raise AttributeError('Multiple peaks look similar to the zeroth order.\n' +
                                     '                Provide initial guess of the x pixel of the zeroth order.')
        buffer      = 30 # half width of region around the initial guess
        order_0_neighborhood = self.spec_1d_net[self.xpix-buffer:self.xpix+buffer]
        xrange      = np.arange(self.xpix-buffer,self.xpix+buffer)
        index_notna = ~(np.isnan(order_0_neighborhood) | np.isnan(xrange))  # To ignore nan's if present
        
        fit_init    = models.Moffat1D(amplitude=np.nanmax(order_0_neighborhood), x_0=self.xpix-buffer+order_0_neighborhood.argmax(), gamma=4.)
        fitter      = fitting.LevMarLSQFitter()
        fit_init.amplitude.min = np.nanmax(order_0_neighborhood)/2.
        # print(fit_init)
        fit_final   = fitter(fit_init, xrange[index_notna], order_0_neighborhood[index_notna])
        self.xpix   = fit_final.x_0.value
        # print(fit_final)
        if plot:
            fig, ax = plt.subplots(figsize =(20,6), ncols=1, nrows=3)
            [axi.set_axis_off() for axi in ax.ravel()]
            ax[0] = fig.add_axes([0.05, 0.85, 0.90, 0.1])
            ax[1] = fig.add_axes([0.05, 0.75, 0.90, 0.1])
            ax[2] = fig.add_axes([0.05, 0.05, 0.90, 0.7])
            
            ax[0].imshow(self.spec_2d_bkg, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
            ax[1].imshow(self.spec_2d, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
            
            ax[2].plot(xrange, order_0_neighborhood, label='Observed spetrum')
            ax[2].plot(xrange, fit_final(xrange), label='Moffat fit')
            ax[2].legend()
            ax[2].grid()
            ax[2].set_xlim(ax[0].get_xlim())
            
    def get_calibrated_spectrum(self, dict_cal, plot=False):
        '''
        Flux and waveength calibration
        
        Parameters
        ----------
        dict_cal : dict
            Calibration parameters of the grating and the order
            
        Returns
        -------
        spectrum : Spectrum1D
            Flux calibrated spectrum in erg/s/cm2/A. Wavelength is in Angstrom
        spectrum_cps : Spectrum1D
            Spectrum in counts/s/A. Wavelength is in Angstrom
        spectrum_ct : Spectrum1D
            Spectrum in counts/pix. Wavelength is in pixels
        '''
        # wavelength calibration
        pixel_mask    = (self.pixels>dict_cal['PIX_MIN']) & (self.pixels<dict_cal['PIX_MAX'])
        wave_pix      = self.pixels[pixel_mask]
        wave_angstrom = (dict_cal['C0'] + dict_cal['C1']*wave_pix)
        
        # interpolation of the effective area curve
        table_ea       = Table.read(dict_cal['EFFECTIVE_AREA_FILE'])
        f              = interpolate.interp1d(table_ea['spectral_axis'], table_ea['flux'])
        effective_area = f(wave_angstrom) * table_ea['flux'].unit
        
        # Flux conversion: 
        # counts/pix  --> counts/s/pix  --> counts/s/A --> erg/s/cm2/A
        spec_ct_pix         = self.spec_1d_net[pixel_mask]
        spec_ct_s_pix       = spec_ct_pix/self.EXPTIME
        spec_ct_s_angstrom  = spec_ct_s_pix/abs(dict_cal['C1'])
        hc_by_lambda        = (const.h * const.c / wave_angstrom).cgs /u.ct
        spec_erg_s_angstrom = hc_by_lambda * spec_ct_s_angstrom / effective_area

        # Poisson error
        fractional_error_Poisson = self.fractional_error_Poisson[pixel_mask]
        e_spec_ct_pix            = fractional_error_Poisson * spec_ct_pix
        e_spec_ct_s_angstrom     = fractional_error_Poisson * spec_ct_s_angstrom
        e_spec_erg_s_angstrom    = fractional_error_Poisson * spec_erg_s_angstrom

        meta_data = {'DETECTOR':dict_cal['DETECTOR'], 'FILTERID':dict_cal['FILTERID'], 'ORDER':dict_cal['ORDER'],
                     'OBJECT':self.OBJECT, 'OBS_ID':self.OBS_ID, 'YPIX':self.YPIX}
        spectrum_ct  = Spectrum1D(spectral_axis = wave_pix+self.xpix*u.pix, 
                                  flux          = spec_ct_pix, 
                                  uncertainty   = StdDevUncertainty(e_spec_ct_pix),
                                  meta          = meta_data.copy())
        spectrum_cps = Spectrum1D(spectral_axis = wave_angstrom, 
                                  flux          = spec_ct_s_angstrom, 
                                  uncertainty   = StdDevUncertainty(e_spec_ct_s_angstrom),
                                  meta          = meta_data.copy())
        spectrum     = Spectrum1D(spectral_axis = wave_angstrom, 
                                  flux          = spec_erg_s_angstrom, 
                                  uncertainty   = StdDevUncertainty(e_spec_erg_s_angstrom),
                                  meta          = meta_data.copy())
        return spectrum, spectrum_cps, spectrum_ct

    def plot(self):
        '''
        Plot flux calibrated spectra (order -1 and -2)
        '''
        fig, ax = plt.subplots(figsize=(10,6))
        ax.set_title(f'{self.OBJECT}  y={self.YPIX}')
        if hasattr(self, 'spectrum_m1'): plot_spectrum(self.spectrum_m1, ax=ax, color='b')
        if hasattr(self, 'spectrum_m2'): plot_spectrum(self.spectrum_m2, ax=ax, color='r')
        
    def plot_all(self, save=False, folder_path=''):
        '''
        Plot the 2d spectrum, 1d spectrum, wavelength calibrated count spectra and flux calibrated spectra.
        
        Returns
        -------
        fig : Figure
            matplotlib.figure.Figure
        ax : ndarray
            Array of matplotlib.axes._subplots.AxesSubplot
        '''
        fig, ax = plt.subplots(figsize =(20,10), ncols=1, nrows=5)
        [axi.set_axis_off() for axi in ax.ravel()]
        ax[0] = fig.add_axes([0.05, 0.94, 0.90, 0.05])
        ax[1] = fig.add_axes([0.05, 0.89, 0.90, 0.05])
        ax[2] = fig.add_axes([0.05, 0.77, 0.90, 0.12])
        ax[3] = fig.add_axes([0.05, 0.52, 0.90, 0.22])
        ax[4] = fig.add_axes([0.05, 0.10, 0.90, 0.42])
        plt.setp(ax[0].get_xticklabels(),visible=False)
        plt.setp(ax[1].get_xticklabels(),visible=False)
        plt.setp(ax[3].get_xticklabels(),visible=False)
        
        ax[0].set_title(f'{self.OBJECT}  y={self.YPIX}')

        ax[0].text(0.01,0.99, 'Background 2d spectrum',transform=ax[0].transAxes, va='top', ha='left')
        ax[0].imshow(self.spec_2d_bkg, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
        
        ax[1].text(0.01,0.99, 'Source 2d spectrum',transform=ax[1].transAxes, va='top', ha='left')
        ax[1].imshow(self.spec_2d, origin='lower', aspect='auto', norm=PowerNorm(0.1), interpolation='gaussian')
        
        ax[2].text(0.01,0.99, 'Source 1d spectrum',transform=ax[2].transAxes, va='top', ha='left')
        ax[2].plot(self.spec_1d_net, c='0.5')
        ax[2].axvline(self.xpix, c="k", ls=':', zorder=0)
        ax[2].set_xlim(ax[0].get_xlim())
        if hasattr(self, 'spectrum_ct_m1'): plot_spectrum(self.spectrum_ct_m1, ax=ax[2], color='b')
        if hasattr(self, 'spectrum_ct_m2'): plot_spectrum(self.spectrum_ct_m2, ax=ax[2], color='r')
        
        if (self.DETECTOR=='FUV') & (self.FILTERID=='Grating1'):
            ax[2].axvline(dict_FUV_G1_M1['PIX_MIN'].value+self.xpix, c='b', ls=':', zorder=0)
            ax[2].axvline(dict_FUV_G1_M1['PIX_MAX'].value+self.xpix, c='b', ls=':', zorder=0)
            ax[2].axvline(dict_FUV_G1_M2['PIX_MIN'].value+self.xpix, c='r', ls=':', zorder=0)
            ax[2].axvline(dict_FUV_G1_M2['PIX_MAX'].value+self.xpix, c='r', ls=':', zorder=0)
        if (self.DETECTOR=='FUV') & (self.FILTERID=='Grating2'):
            ax[2].axvline(dict_FUV_G2_M2['PIX_MIN'].value+self.xpix, c='r', ls=':', zorder=0)
            ax[2].axvline(dict_FUV_G2_M2['PIX_MAX'].value+self.xpix, c='r', ls=':', zorder=0)
        if (self.DETECTOR=='NUV') & (self.FILTERID=='Grating'):
            ax[2].axvline(dict_NUV_G_M1['PIX_MIN'].value+self.xpix, c='b', ls=':', zorder=0)
            ax[2].axvline(dict_NUV_G_M1['PIX_MAX'].value+self.xpix, c='b', ls=':', zorder=0)
        
            
        ax[3].text(0.01,0.99, 'Source CPS spectrum',transform=ax[3].transAxes, va='top', ha='left')
        if hasattr(self, 'spectrum_cps_m1'): plot_spectrum(self.spectrum_cps_m1, ax=ax[3], color='b')
        if hasattr(self, 'spectrum_cps_m2'): plot_spectrum(self.spectrum_cps_m2, ax=ax[3], color='r')
        
        ax[4].text(0.01,0.99, 'Flux calibrated spectrum',transform=ax[4].transAxes, va='top', ha='left')
        if hasattr(self, 'spectrum_m1'): plot_spectrum(self.spectrum_m1, ax=ax[4], color='b')
        if hasattr(self, 'spectrum_m2'): plot_spectrum(self.spectrum_m2, ax=ax[4], color='r')
        
        if save:
            if not os.path.exists(folder_path): os.makedirs(folder_path)
            file_name = folder_path + '%s_%s_%s_%d.jpg'%(self.OBJECT, self.DETECTOR, self.FILTERID, self.YPIX)
            plt.savefig(file_name, dpi=300, bbox_inches='tight')
            plt.close()
        return fig, ax
        
    def save_all(self, folder_path='',overwrite=False):
        '''
        Save spectra if available.
        
        Parameters
        ----------
        folder_path : str
            Folder in which to save the files. Default is working directory.
        '''
        self.plot_all(save=True, folder_path=folder_path)

        for spectrum_name in ['spectrum_m1', 'spectrum_m2']:
            if hasattr(self, spectrum_name):
                spectrum = getattr(self, spectrum_name)
                file_name = folder_path + '%s_%s_%s_%d_m%d.fits'%(spectrum.meta['OBJECT'], spectrum.meta['DETECTOR'], 
                                                                  spectrum.meta['FILTERID'], self.YPIX, abs(spectrum.meta['ORDER']))
                save_spectrum(spectrum, file_name, overwrite=overwrite)

######################################################################
######################################################################

def plot_spectrum(spec, ax=None, color='k', label=None):
    '''
    Plot specutils spec1D.

    Parameters
    ----------
    spec : Spectrum1D
        Spectum1D object with errors and UVIT metadata.
    '''
    if ax   ==None: fig, ax = plt.subplots(figsize=(20,6))
    if label==None: label='%s %s order=%d'%(spec.meta['DETECTOR'],spec.meta['FILTERID'],spec.meta['ORDER'])

    ax.plot(spec.spectral_axis, spec.flux, drawstyle='steps-mid', c=color, label=label)
    ax.fill_between(spec.spectral_axis.value, spec.flux.value-spec.uncertainty.array, spec.flux.value+spec.uncertainty.array, 
                    color=color, alpha=0.2, step='mid')
        
    ax.set(ylabel=spec.flux.unit, xlabel=spec.spectral_axis.unit)
    ax.legend()
    return ax

def save_spectrum(spec, file_name, overwrite=False):
    '''
    Save Spectrum1D object
    
    Parameters
    ----------
    spec : Spectrum1D
        Spectrum1D object
    file_name : str
        File name
        
    Returns
    -------
    qt_spec : QTable
        Astropy QTable .fits file including units and metadata
    '''
    qt_spec = QTable([spec.spectral_axis.quantity, spec.flux, spec.uncertainty.array*spec.uncertainty.unit],
                         names=('spectral_axis', 'flux', 'uncertainty'),
                         meta =spec.meta)
    
    qt_spec.write(file_name, overwrite=overwrite)
    print('Saved %s'%file_name)
    
def load_spectrum(file_name):
    '''
    Load Spectrum1D object from a QTable .fits file
    
    Parameters
    ---------
    file_name : str
        File name with path
    '''
    qt_spec = Table.read(file_name)
    spec    = Spectrum1D(spectral_axis     = qt_spec['spectral_axis'].quantity, 
                              flux         = qt_spec['flux'].quantity, 
                              uncertainty  = StdDevUncertainty(qt_spec['uncertainty'].quantity),
                              meta         = qt_spec.meta)
    return spec