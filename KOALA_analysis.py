#!/usr/bin/python
# -*- coding: utf-8 -*-
# # KOALA analysis pipeline
# Brief description of the KOALA class
# by Ángel López and Yago

from astropy.io import fits
from astropy.wcs import WCS

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sys
import copy

from scipy import interpolate, signal, optimize

# COLOR

# rgb
fuego_color_map = colors.LinearSegmentedColormap.from_list("fuego", ((0.25, 0, 0),  (0.5,0,0),    (1, 0, 0), (1, 0.5, 0), (1, 0.75, 0), (1, 1, 0), (1, 1, 1)), N=256, gamma=1.0)
fuego_color_map.set_bad('lightgray')
plt.register_cmap(cmap=fuego_color_map)

projo = [0.25, 0.5, 1, 1.0, 1.00, 1, 1]
pverde = [0.00, 0.0, 0, 0.5, 0.75, 1, 1]
pazul = [0.00, 0.0, 0, 0.0, 0.00, 0, 1]


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / N


#def weighted_mean_and_variance(vector, weight):
#    total_weight = np.sum(weight)
#    mean = np.sum(vector*weight) / total_weight
#    v_squared = np.sum(vector**2*weight) / total_weight
#    return mean, v_squared-mean**2


def cumulaive_Moffat(r2, L_star, alpha2, beta):
    return L_star*(1 - np.power(1+(r2/alpha2), -beta))


def fit_Moffat(r2_growth_curve, F_growth_curve,
               F_guess, r2_half_light, r_out, plot=False):
    """
    Fits a Moffat profile to a flux growth curve
    as a function of radius squared,
    cutting at to r_out (in units of the half-light radius),
    provided an initial guess of the total flux and half-light radius squared.
    """
    index_cut = np.searchsorted(r2_growth_curve, r2_half_light*r_out**2)
    fit, cov = optimize.curve_fit(cumulaive_Moffat,
                                  r2_growth_curve[:index_cut], F_growth_curve[:index_cut],
                                  p0=(F_guess, r2_half_light, 1)
                                  )
    if plot:
        print "Best-fit: L_star =", fit[0]
        print "          alpha =", np.sqrt(fit[1])
        print "          beta =", fit[2]
        r_norm = np.sqrt(np.array(r2_growth_curve) / r2_half_light)
        plt.plot(r_norm, cumulaive_Moffat(np.array(r2_growth_curve),
                                          fit[0], fit[1], fit[2])/fit[0], ':')

    return fit


# -----------------------------------------------------------------------------
class Interpolated_cube(object):
    """
    Constructs a map by accumulating RSS with given offsets.
    """

    # -------------------------------------------------------------------------
    def __init__(self, RSS, pixel_size_arcsec, kernel_size_arcsec,
                 centre_deg=[], size_arcsec=[]):
        self.RSS = RSS
        self.n_wave = RSS.n_wave
        self.pixel_size_arcsec = pixel_size_arcsec
        self.kernel_size_arcsec = kernel_size_arcsec
        self.kernel_size_pixels = kernel_size_arcsec/pixel_size_arcsec  # must be a float number!

        if len(centre_deg) == 2:
            self.RA_centre_deg = centre_deg[0]
            self.DEC_centre_deg = centre_deg[1]
        else:
            RA_min, RA_max, DEC_min, DEC_max = coord_range([RSS])
#            print RA_min, RA_max, DEC_min, DEC_max
            self.RA_centre_deg = (RA_min + RA_max)/2.
            self.DEC_centre_deg = (DEC_min + DEC_max)/2.
        xoffset_centre_arcsec = (self.RA_centre_deg-RSS.RA_centre_deg)*3600.
        yoffset_centre_arcsec = (self.DEC_centre_deg-RSS.DEC_centre_deg)*3600.

        if len(size_arcsec) == 2:
            self.n_cols = np.int(size_arcsec[0]/pixel_size_arcsec)
            self.n_rows = np.int(size_arcsec[1]/pixel_size_arcsec)
        else:
            self.n_cols = 2 * \
                (np.int(np.nanmax(np.abs(RSS.offset_RA_arcsec-xoffset_centre_arcsec))/pixel_size_arcsec) +
                 np.int(self.kernel_size_pixels + 2))
            self.n_rows = 2 * \
                (np.int(np.nanmax(np.abs(RSS.offset_DEC_arcsec-yoffset_centre_arcsec))/pixel_size_arcsec) +
                 np.int(self.kernel_size_pixels + 2))


#        print x_min, x_centre_arcsec, x_max, y_min, y_centre_arcsec, y_max,
#        print self.n_rows, self.n_cols, self.RA_centre_deg, self.DEC_centre_deg

##        self.n_cols = int((np.max(RSS.offset_RA_arcsec) - np.min(RSS.offset_RA_arcsec) + 3*kernel_size_arcsec) / pixel_size_arcsec + 1)
##        self.n_rows = int((np.max(RSS.offset_DEC_arcsec) - np.min(RSS.offset_DEC_arcsec) + 3*kernel_size_arcsec) / pixel_size_arcsec + 1)
##        self.RA_centre_deg = (np.max(RSS.offset_RA_arcsec) + np.min(RSS.offset_RA_arcsec)) / 2.
##        self.DEC_centre_deg = (np.max(RSS.offset_DEC_arcsec) + np.min(RSS.offset_DEC_arcsec)) / 2.
#        self.RA_centre_deg = 0.
#        self.DEC_centre_deg = 0.
#        self.n_cols = 2*int((np.max(np.abs(RSS.offset_RA_arcsec)) + 2*kernel_size_arcsec) / pixel_size_arcsec + 1)
#        self.n_rows = 2*int((np.max(np.abs(RSS.offset_DEC_arcsec)) + 2*kernel_size_arcsec) / pixel_size_arcsec + 1)

        self._weighted_I = np.zeros((self.n_wave, self.n_rows, self.n_cols))
        self._weight = np.zeros_like(self._weighted_I)

        print "\n> Smooth cube, (RA, DEC)_centre = ({}, {}) degree" \
            .format(self.RA_centre_deg, self.DEC_centre_deg)
        print "  size = {} columns (RA) x {} rows (DEC); {:.2f} x {:.2f} arcsec" \
            .format(self.n_cols, self.n_rows, self.n_cols*pixel_size_arcsec, self.n_rows*pixel_size_arcsec)
        sys.stdout.write("  Adding {} spectra...       ".format(RSS.n_spectra))
        sys.stdout.flush()
        output_every_few = np.sqrt(RSS.n_spectra)+1
        next_output = -1
        for i in range(RSS.n_spectra):
            if i > next_output:
                sys.stdout.write("\b"*6)
                sys.stdout.write("{:5.2f}%".format(i*100./RSS.n_spectra))
                sys.stdout.flush()
                next_output = i + output_every_few
            offset_rows = (RSS.offset_DEC_arcsec[i]-yoffset_centre_arcsec) / pixel_size_arcsec
            offset_cols = (RSS.offset_RA_arcsec[i]-xoffset_centre_arcsec) / pixel_size_arcsec

            corrected_intensity = (RSS.intensity[i] - RSS.sky_emission) \
                * RSS.extinction_correction / RSS.relative_throughput[i]

            self.add_spectrum(corrected_intensity, offset_rows, offset_cols)
        self.data = self._weighted_I / self._weight
        sys.stdout.write("\b"*6)
        sys.stdout.write(" DONE!\n")

        self.trace_peak()

    # -------------------------------------------------------------------------
    def add_spectrum(self, intensity, offset_rows, offset_cols):
        """
        Add one single spectrum to the datacube

        Parameters
        ----------
        intensity: np.array(float)
          Spectrum.
        offset_rows, offset_cols: float
          Offset with respect to the image centre, in pixels.
        kernel_FWHM_pixels: float
          FWHM of the interpolating kernel, in pixels
        """
#        print "---"
        kernel_centre_x = .5*self.n_cols + offset_cols
        x_min = int(kernel_centre_x - self.kernel_size_pixels)
        x_max = int(kernel_centre_x + self.kernel_size_pixels) + 1
        n_points_x = x_max-x_min
        x = np.linspace(x_min-kernel_centre_x, x_max-kernel_centre_x, n_points_x) / self.kernel_size_pixels
        x[0] = -1.
        x[-1] = 1.
        weight_x = np.diff((3.*x - x**3 + 2.) / 4)
#        print x_min, kernel_centre, x_max

        kernel_centre_y = .5*self.n_rows + offset_rows
        y_min = int(kernel_centre_y - self.kernel_size_pixels)
        y_max = int(kernel_centre_y + self.kernel_size_pixels) + 1
        n_points_y = y_max-y_min
        y = np.linspace(y_min-kernel_centre_y, y_max-kernel_centre_y, n_points_y) / self.kernel_size_pixels
        y[0] = -1.
        y[-1] = 1.
        weight_y = np.diff((3.*y - y**3 + 2.) / 4)
#        print y_min, kernel_centre, y_max

        if x_min < 0 or x_max >= self.n_cols or y_min < 0 or y_max >= self.n_rows:
            print "\n***************************************"
            print "WARNING: Spectra outside field of view:"
            print x_min, kernel_centre_x, x_max
            print y_min, kernel_centre_y, y_max
            print "***************************************"
        else:
            bad_wavelengths = np.argwhere(np.isnan(intensity))
            intensity[bad_wavelengths] = 0.
            ones = np.ones_like(intensity)
            ones[bad_wavelengths] = 0.
            self._weighted_I[:, y_min:y_max-1, x_min:x_max-1] += intensity[:, np.newaxis, np.newaxis] * weight_y[np.newaxis, :, np.newaxis] * weight_x[np.newaxis, np.newaxis, :]
            self._weight[:, y_min:y_max-1, x_min:x_max-1] += ones[:, np.newaxis, np.newaxis] * weight_y[np.newaxis, :, np.newaxis] * weight_x[np.newaxis, np.newaxis, :]

#        print self._weight[y_min:y_max-1, x_min:x_max-1]

#        Old stuff:
#        dx = np.arange(self.n_cols) - image_centre - offset_cols
#        image_centre = .5*(self.n_rows-1)
#        dy = np.arange(self.n_rows) - image_centre - offset_rows
#        deltas = np.meshgrid(dx, dy)
#        normalised_distance_from_fibre_squared = (deltas[0]**2 + deltas[1]**2) / kernel_FWHM_pixels**2
#        fibre_weights = self._kernel(normalised_distance_from_fibre_squared)
#        
#        self._weight += fibre_weights
#        self._weighted_I = self._weighted_I + np.outer(intensity, fibre_weights).reshape((self.n_wave, self.n_rows, self.n_cols))

    # -------------------------------------------------------------------------
    def plot_wavelength(self, wavelength,
                        norm=colors.PowerNorm(gamma=1./4.),
                        save_file=""):
        interpolated_map = self.data[np.searchsorted(self.RSS.wavelength, wavelength)]

        plt.figure(figsize=(12, 12))
        plt.imshow(interpolated_map, origin='lower', interpolation='none',
                   norm=norm, cmap=fuego_color_map,
                   extent=(-.5*self.n_cols*self.pixel_size_arcsec,
                           0.5*self.n_cols*self.pixel_size_arcsec,
                           -.5*self.n_rows*self.pixel_size_arcsec,
                           0.5*self.n_rows*self.pixel_size_arcsec)
                   )
        plt.colorbar()
        plt.contour(interpolated_map,
                    extent=(-.5*self.n_cols*self.pixel_size_arcsec,
                            0.5*self.n_cols*self.pixel_size_arcsec,
                            -.5*self.n_rows*self.pixel_size_arcsec,
                            0.5*self.n_rows*self.pixel_size_arcsec)
                    )
#        plt.gray()
        plt.gca().invert_xaxis()
        plt.title("{} - {} $\AA$".format(self.RSS.description, wavelength))
        plt.xlabel("$\Delta$ RA [arcsec]")
        plt.ylabel("$\Delta$ DEC [arcsec]")
        if save_file == "":
            plt.show()
        else:
            plt.savefig(save_file)
        plt.close()

    # -------------------------------------------------------------------------
    def plot_weight(self):
        interpolated_map = np.mean(self._weight, axis=0)

        plt.figure(figsize=(12, 12))
        plt.imshow(interpolated_map, origin='lower', interpolation='none',
                   extent=(-.5*self.n_cols*self.pixel_size_arcsec,
                           0.5*self.n_cols*self.pixel_size_arcsec,
                           -.5*self.n_rows*self.pixel_size_arcsec,
                           0.5*self.n_rows*self.pixel_size_arcsec)
                   )
#        plt.gray()
        plt.colorbar()
        plt.gca().invert_xaxis()
        plt.show()

#    # -------------------------------------------------------------------------
#    def find_peak_at_wavelength(self, index, tolerance_in_pixels):
#        x = np.arange(self.n_cols)
#        y = np.arange(self.n_rows)
#    
#        weight = self.data[index]**2
#        x_peak, var_x = weighted_mean_and_variance(x, np.nansum(weight, axis=0))
#        y_peak, var_y = weighted_mean_and_variance(y, np.nansum(weight, axis=1))
#    
#        delta = tolerance_in_pixels**2 + 1
#        while delta > tolerance_in_pixels**2:
#            new_weight = weight * np.exp(-.5*((x-x_peak)**2/(var_x+4))[np.newaxis, :] - .5*((y-y_peak)**2/(var_y+4))[:, np.newaxis])
#            new_x, new_var_x = weighted_mean_and_variance(x, np.nansum(new_weight, axis=0))
#            new_y, new_var_y = weighted_mean_and_variance(y, np.nansum(new_weight, axis=1))
#            delta = (new_x-x_peak)**2 + (new_y-y_peak)**2
#    #        print x_peak, y_peak, var_x, var_y, np.sqrt(delta)
#    #        print new_x, new_y, new_var_x, new_var_y, np.sqrt(delta)
#            x_peak = new_x
#            y_peak = new_y
#    #        var_x = new_var_x
#    #        var_y = new_var_y
#            weight = new_weight
#
#        return x_peak, y_peak
#
#    # -------------------------------------------------------------------------
#    def trace_peak(self, tolerance_in_pixels=0.01, plot=True):
#        print "\n> Tracing intensity peak over all wavelengths with {:.3f} pixels tolerance...".format(tolerance_in_pixels)
##        peak = np.array([self.find_peak_at_wavelength(wavelength, tolerance_in_pixels)
##                        for wavelength in range(self.n_wave)])
#        peak = []
#        sys.stdout.write("  Adding {} wavelengths...       ".format(self.n_wave))
#        sys.stdout.flush()
#        output_every_few = np.sqrt(self.n_wave)+1
#        next_output = -1
#        for wavelength in range(self.n_wave):
#            if wavelength > next_output:
#                sys.stdout.write("\b"*6)
#                sys.stdout.write("{:5.2f}%".format(wavelength*100./self.n_wave))
#                sys.stdout.flush()
#                next_output = wavelength + output_every_few
#            peak.append(self.find_peak_at_wavelength(wavelength,
#                                                     tolerance_in_pixels))
#        sys.stdout.write("\b"*6)
#        sys.stdout.write(" DONE!\n")
#        peak = np.array(peak)
#        self.x_peak = peak[:, 0]
#        self.y_peak = peak[:, 1]
#
#        x = (self.x_peak-.5*self.n_cols)*self.pixel_size_arcsec
#        y = (self.y_peak-.5*self.n_rows)*self.pixel_size_arcsec
#
#        print "  Peak offset in arcsec from image centre:"
#        print "  (RA, DEC) = ({:.2f}+-{:.2f}, {:.2f}+-{:.2f})"\
#              .format(np.nanmedian(x),
#                      (np.nanpercentile(x, 84)-np.nanpercentile(x, 16))/2,
#                      np.nanmedian(y),
#                      (np.nanpercentile(y, 84)-np.nanpercentile(y, 16))/2
#                      )
#        print "  DONE!"
#
#        if plot:
#            plt.figure(figsize=(10, 5))
#            wl = self.RSS.wavelength
#            x = self.x_peak - np.nanmedian(self.x_peak)
#            y = self.y_peak - np.nanmedian(self.y_peak)
#            plt.plot(wl, x, 'k.', alpha=0.2)
#            plt.plot(wl, y, 'r.', alpha=0.2)
#            plt.plot(wl, signal.medfilt(x, 5), 'k-')
#            plt.plot(wl, signal.medfilt(y, 5), 'r-')
#            hi = np.max([np.nanpercentile(x, 90), np.nanpercentile(y, 90)])
#            lo = np.min([np.nanpercentile(x, 10), np.nanpercentile(y, 10)])
#            plt.ylim(lo, hi)
#            plt.show()
#            plt.close()

    # -------------------------------------------------------------------------
    def trace_peak(self, plot=False):
        print "\n> Tracing intensity peak over all wavelengths..."
        x = np.arange(self.n_cols)
        y = np.arange(self.n_rows)
        weight = np.nan_to_num(self.data)
        mean_iamge = np.nanmean(weight, axis=0)
        mean_iamge /= np.nanmean(mean_iamge)
        weight *= mean_iamge[np.newaxis, :, :]
        xw = x[np.newaxis, np.newaxis, :] * weight
        yw = y[np.newaxis, :, np.newaxis] * weight
        w = np.nansum(weight, axis=(1, 2))
        self.x_peak = np.nansum(xw, axis=(1, 2)) / w
        self.y_peak = np.nansum(yw, axis=(1, 2)) / w
        if plot:
            plt.figure(figsize=(10, 5))
            wl = self.RSS.wavelength
            x = (self.x_peak-np.nanmedian(self.x_peak))*self.pixel_size_arcsec
            y = (self.y_peak-np.nanmedian(self.y_peak))*self.pixel_size_arcsec
            plt.plot(wl, x, 'k.', alpha=0.2)
            plt.plot(wl, y, 'r.', alpha=0.2)
            odd_number = 2 * int(np.sqrt(cube.n_wave) / 2) + 1
            plt.plot(wl, signal.medfilt(x, odd_number), 'k-')
            plt.plot(wl, signal.medfilt(y, odd_number), 'r-')
            hi = np.max([np.nanpercentile(x, 90), np.nanpercentile(y, 90)])
            lo = np.min([np.nanpercentile(x, 10), np.nanpercentile(y, 10)])
            plt.ylim(lo, hi)
            plt.ylabel("offset [arcsec]")
            plt.xlabel("wavelength [$\AA$]")
            plt.title(self.RSS.description)
            plt.show()
            plt.close()
        print "  Peak coordinates, in pixel: ({:.2f}+-{:.2f}, {:.2f}+-{:.2f})"\
              .format(np.nanmedian(self.x_peak), np.std(self.x_peak),
                      np.nanmedian(self.y_peak), np.std(self.y_peak))
        print "  DONE!"

    # -------------------------------------------------------------------------
    def growth_curve_between(self, lambda_min, lambda_max, plot=False):
        index_min = np.searchsorted(self.RSS.wavelength, lambda_min)
        index_max = np.searchsorted(self.RSS.wavelength, lambda_max)
        intensity = np.nanmean(self.data[index_min:index_max, :, :], axis=0)
        x_peak = np.median(self.x_peak[index_min:index_max])
        y_peak = np.median(self.y_peak[index_min:index_max])
        x = np.arange(self.n_cols) - x_peak
        y = np.arange(self.n_rows) - y_peak
        r2 = np.sum(np.meshgrid(x**2, y**2), axis=0)
        sorted_by_distance = np.argsort(r2, axis=None)

        F_growth_curve = []
        r2_growth_curve = []
        total_flux = 0.
        for spaxel in sorted_by_distance:
            index = np.unravel_index(spaxel, (self.n_rows, self.n_cols))
            I = intensity[index]
    #        print spaxel, r2[index], L, total_flux, np.isnan(L)
    #        if np.isnan(L) == False and L > 0:
            if np.isnan(I) == False:
                total_flux += I  # TODO: Properly account for solid angle...
                F_growth_curve.append(total_flux)
                r2_growth_curve.append(r2[index])

        F_guess = np.max(F_growth_curve)
        r2_half_light = np.interp(.5*F_guess, F_growth_curve, r2_growth_curve)

        if plot:
            r_norm = np.sqrt(np.array(r2_growth_curve) / r2_half_light)
            F_norm = np.array(F_growth_curve) / F_guess
            print "Flux guess =", F_guess, np.nansum(intensity), np.nansum(intensity)/F_guess
            print "Half-light radius:", np.sqrt(r2_half_light)*self.pixel_size_arcsec
            print "Light within 2, 3, 4 half-lght radii:", np.interp([2, 3, 4], r_norm, F_norm)
            plt.plot(r_norm, F_norm, '-')

        return r2_growth_curve, F_growth_curve, F_guess, r2_half_light

    # -------------------------------------------------------------------------
    def half_light_spectrum(self, r_max=1, plot=False):
        r2_growth_curve, F_growth_curve, flux, r2_half_light = self.growth_curve_between(0, 1e30)
        print "\n> Computing growth-curve spectrum..."
        intensity = []
        smooth_x = signal.medfilt(self.x_peak, 11)
        smooth_y = signal.medfilt(self.y_peak, 11)
        for l in range(self.n_wave):
            wavelength = self.RSS.wavelength[l]
            if l % (self.n_wave/10+1) == 0:
                print "  {:.2f} Angstroms (wavelength {}/{})..." \
                      .format(wavelength, l+1, self.n_wave)
            x = np.arange(self.n_cols) - smooth_x[l]
            y = np.arange(self.n_rows) - smooth_y[l]
            r2 = np.sum(np.meshgrid(x**2, y**2), axis=0)
            spaxels = np.where(r2 < r2_half_light*r_max**2)
            intensity.append(np.nansum(self.data[l][spaxels]))
        if plot:
            plt.plot(self.RSS.wavelength, intensity, 'b.', alpha=0.2)
            plt.plot(self.RSS.wavelength, signal.medfilt(intensity,11), 'r-')
            plt.ylim(np.nanpercentile(intensity, 1), np.nanpercentile(intensity, 99))
    #        plt.xlim(6200, 6400)
            plt.xlim(6200, 7400)
    #        plt.xlim(7200, 7400)
        return np.array(intensity)

    # -------------------------------------------------------------------------
    def do_response_curve(self, filename, plot=False):
        lambda_cal, flux_cal, delta_lambda = np.loadtxt(filename, usecols=(0,1,3), unpack=True)
        lambda_min = lambda_cal - delta_lambda/2
        lambda_max = lambda_cal + delta_lambda/2
        measured_counts = np.array([self.fit_Moffat_between(lambda_min[i],
                                                            lambda_max[i])[0]
                                    if lambda_cal[i] > 6300 and
                                    lambda_cal[i] < 7400
                                    else np.NaN
                                    for i in range(len(lambda_cal))])

        self.response_wavelength = lambda_cal
        self.response_curve = measured_counts / flux_cal
        scale = np.nanmedian(self.response_curve)
        if plot:
            plt.figure(figsize=(12, 8))
            self.half_light_spectrum(5, True)
            plt.plot(lambda_cal, measured_counts, 'k+', ms=30, mew=3)
            plt.plot(lambda_cal, flux_cal*scale, 'k*-')
            plt.plot(lambda_cal, flux_cal*self.response_curve, 'k:')
            plt.show()

            plt.figure(figsize=(12, 8))
            plt.plot(lambda_cal, self.response_curve, 'k-')
            plt.show()

    # -------------------------------------------------------------------------
    def fit_Moffat_between(self, lambda_min=0, lambda_max=1e30, r_fit=5, plot=False):
        r2_growth_curve, F_growth_curve, flux, r2_half_light = self.growth_curve_between(lambda_min, lambda_max, plot)
        flux, alpha, beta = fit_Moffat(r2_growth_curve, F_growth_curve,
                                       flux, r2_half_light, r_fit, plot)
        r2_half_light = alpha * (np.power(2., 1./beta) - 1)
        if plot:
            print "Moffat fit: Flux = {:.3e},".format(flux), \
                "HWHM = {:.3f},".format(np.sqrt(r2_half_light)*self.pixel_size_arcsec), \
                "beta = {:.3f}".format(beta)

        return flux, np.sqrt(r2_half_light)*self.pixel_size_arcsec, beta

    # -------------------------------------------------------------------------
    def total_spectrum(self):
        x = np.arange(self.n_cols)
        y = np.arange(self.n_rows)
        
        self.data


# -----------------------------------------------------------------------------
class RSS(object):
    """
    Collection of row-stacked spectra (RSS).

    Attributes
    ----------
    wavelength: np.array(float)
      Wavelength, in Angstrom.
    intensity: np.array(float)
      Intensity :math:`I_\lambda` per unit wavelength.
    variance: np.array(float)
      Variance :math:`\sigma^2_\lambda` per unit wavelength
      (note the square in the definition of the variance).
    """

    # -------------------------------------------------------------------------
    def __init__(self):
        self.description = "Undefined row-stacked spectra (RSS)"
        
        self.n_spectra = 0
        self.n_wave = 0

        self.wavelength = np.zeros((0))
        self.intensity = np.zeros((0, 0))
        self.variance = np.zeros_like(self.intensity)

        self.RA_centre_deg = 0.
        self.DEC_centre_deg = 0.
        self.offset_RA_arcsec = np.zeros((0))
        self.offset_DEC_arcsec = np.zeros_like(self.offset_RA_arcsec)

    # -------------------------------------------------------------------------
    def set_data(self, wavelength, intensity, variance, offset_RA_arcsec, offset_DEC_arcsec):
        self.wavelength = wavelength
        self.n_wave = len(wavelength)

        if variance.shape != intensity.shape:
            print "\n* ERROR: * the intensity and variance matrices are", \
                  intensity.shape, "and", variance.shape, "respectively\n"
            raise ValueError
        n_dim = len(intensity.shape)
        if n_dim == 2:
            self.intensity = intensity
            self.variance = variance
        elif n_dim == 1:
            self.intensity = intensity.reshape((1, self.n_wave))
            self.variance = variance.reshape((1, self.n_wave))
        else:
            print "\n* ERROR: * the intensity matrix supplied has", \
                  n_dim, "dimensions\n"
            raise ValueError

        self.n_spectra = self.intensity.shape[0]
        self.n_wave = len(self.wavelength)
        print "  {} spectra with {} wavelengths" \
              .format(self.n_spectra, self.n_wave), \
              "between {:.2f} and {:.2f} Angstrom" \
              .format(self.wavelength[0], self.wavelength[-1])
        if self.intensity.shape[1] != self.n_wave:
            print "\n* ERROR: * spectra have", self.intensity.shape[1], \
                  "wavelengths rather than", self.n_wave
            raise ValueError
        if len(offset_RA_arcsec) != self.n_spectra or \
           len(offset_DEC_arcsec) != self.n_spectra:
            print "\n* ERROR: * offsets (RA, DEC) = ({},{})" \
                  .format(len(self.offset_RA_arcsec),
                          len(self.offset_DEC_arcsec)), \
                  "rather than", self.n_spectra
            raise ValueError
        else:
            self.offset_RA_arcsec = offset_RA_arcsec
            self.offset_DEC_arcsec = offset_DEC_arcsec

        self.sky_emission = np.zeros(self.n_wave)
#        self.find_sky_emission()
        self.relative_throughput = np.ones(self.n_spectra)
        self.airmass = 0
        self.extinction_correction = np.ones(self.n_wave)

    # -------------------------------------------------------------------------
    def plot_spectrum(self, spectrum_number):
        plt.plot(self.wavelength, self.intensity[spectrum_number])
        error = 3*np.sqrt(self.variance[spectrum_number])
        plt.fill_between(self.wavelength,
                         self.intensity[spectrum_number]-error,
                         self.intensity[spectrum_number]+error, alpha=.1)

    # -------------------------------------------------------------------------
    def plot_spectra(self, list_spectra='all', wavelength_range=[0]):
        if list_spectra == 'all':
            list_spectra = range(self.n_spectra)
        if len(wavelength_range) == 2:
            plt.xlim(wavelength_range[0], wavelength_range[1])

        for i in list_spectra:
            self.plot_spectrum(i)

    # -------------------------------------------------------------------------
    def flux_between(self, lambda_min, lambda_max, list_spectra=[]):
        index_min = np.searchsorted(self.wavelength, lambda_min)
        index_max = np.searchsorted(self.wavelength, lambda_max)+1
        if len(list_spectra) == 0:
            list_spectra = range(self.n_spectra)

        n_spectra = len(list_spectra)
        fluxes = np.empty(n_spectra)
        variance = np.empty(n_spectra)
        for i in range(n_spectra):
            fluxes[i] = np.nanmean(self.intensity[list_spectra[i],
                                                  index_min:index_max])
            variance[i] = np.nanmean(self.variance[list_spectra[i],
                                                   index_min:index_max])

        return fluxes*(lambda_max-lambda_min), variance*(lambda_max-lambda_min)
#        WARNING: Are we overestimating errors?

    # -------------------------------------------------------------------------
    def median_between(self, lambda_min, lambda_max, list_spectra=[]):
        index_min = np.searchsorted(self.wavelength, lambda_min)
        index_max = np.searchsorted(self.wavelength, lambda_max)+1
        if len(list_spectra) == 0:
            list_spectra = range(self.n_spectra)

        n_spectra = len(list_spectra)
        medians = np.empty(n_spectra)
        for i in range(n_spectra):
            medians[i] = np.nanmedian(self.intensity[list_spectra[i],
                                                     index_min:index_max])
        return medians

    # -------------------------------------------------------------------------
    def line_flux(self,
                  left_min, left_max,
                  line_min, line_max,
                  right_min, right_max,
                  list_spectra=[]):
        if len(list_spectra) == 0:
            list_spectra = range(self.n_spectra)

        line, var_line = self.flux_between(line_min, line_max, list_spectra)
        left, var_left = self.flux_between(left_min, left_max,
                                           list_spectra)/(left_max-left_min)
        right, var_right = self.flux_between(right_min, right_max,
                                             list_spectra)/(left_max-left_min)
        wavelength_left = (left_min+left_max)/2
        wavelength_line = (line_min+line_max)/2
        wavelength_right = (right_min+right_max)/2
        continuum = left + \
            (right-left)*(wavelength_line-wavelength_left) \
            / (wavelength_right-wavelength_left)
        var_continuum = (var_left+var_right)/2

        return line - continuum*(line_max-line_min), \
            var_line+var_continuum*(line_max-line_min)
#           WARNING: Are we overestimating errors?

#    # -------------------------------------------------------------------------
#    def find_sky_emission(self, reference_min, reference_max,
#                          filename='sky_lines.txt'):
#        print "\n> Identifying sky spaxels based on wavelenth interval:", \
#              reference_min, reference_max
#        print "  Sky lines in file '{}':".format(filename)
#        self.sky_line_flux = []
#        self.sky_line_variance = []
#        sky_lines_list = np.loadtxt(filename)
#        for line in sky_lines_list:
#            print "   ", line[2], "-", line[3], "Angstrom"
#            flux_sky, var_sky = self.line_flux(line[0], line[1],
#                                               line[2], line[3],
#                                               line[4], line[5])
##            flux_sky[flux_sky < 0] = np.NaN
#            self.sky_line_flux.append(flux_sky)
#            self.sky_line_variance.append(var_sky)
#
#        flux_ratio = self.flux_between(reference_min, reference_max)[0] \
#            / np.abs(np.nansum(self.sky_line_flux, axis=0))
#        sorted_by_flux = np.argsort(flux_ratio)

    # -------------------------------------------------------------------------
    def find_sky_emission(self, plot=False):
        print "\n> Identifying sky spaxels..."
        I = np.nanmedian(self.intensity, axis=0)
        I2 = np.nanmean(self.intensity**2, axis=0)
        var = np.sqrt(I2-I**2)
        weight_sky = I/np.nanmedian(I)
        weight_obj = var/np.nanmedian(var)
        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(self.wavelength, weight_sky, 'c-', label='sky')
            plt.plot(self.wavelength, weight_obj, 'k-', label='object')
            plt.yscale('log')
            plt.ylabel("weight")
            plt.xlabel("wavelength [$\AA$]")
            plt.title(rss.description)
            plt.legend(frameon=False)
            plt.show()
            plt.close()

        flux_sky = np.nanmean(self.intensity*weight_sky[np.newaxis, :], axis=1)
        flux_object = np.nanmean(self.intensity*weight_obj[np.newaxis, :], axis=1)
        flux_ratio = flux_object / flux_sky
        sorted_by_flux = np.argsort(flux_ratio)
        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(flux_ratio[sorted_by_flux], 'r-')
            plt.plot(flux_sky[sorted_by_flux], 'c-')
            plt.plot(flux_object[sorted_by_flux], 'k-')
            plt.yscale('log')
            plt.show()

        n = 10
        minimum_delta = 1e30
        optimal_n = n
        while n <= self.n_spectra:
            sky = np.nanmedian(self.intensity[sorted_by_flux[:n]], axis=0)
            upper = np.nanpercentile(self.intensity[sorted_by_flux[:n]], 84, axis=0)
            lower = np.nanpercentile(self.intensity[sorted_by_flux[:n]], 16, axis=0)
            med_error = np.sqrt(np.nanmedian(self.variance[sorted_by_flux[:n]], axis=0))
            delta = np.nanmax((upper-lower)/med_error) \
                * np.nanmax(((upper-sky)-(sky-lower))/med_error) \
                * np.nanmax(upper-sky) / n
            if delta < minimum_delta:
                minimum_delta = delta
                optimal_n = n
#            print n, delta
#            plt.semilogy(n, delta, 'k*')
            n += n/5 + 1
        print " ", optimal_n, "spaxels identified as sky"
        self.sky_spaxels = sorted_by_flux[:optimal_n]
        self.sky_emission = np.nanmedian(self.intensity[sorted_by_flux[:optimal_n]], axis=0)
        if plot:
            self.RSS_map(flux_ratio, None, self.sky_spaxels)

    # -------------------------------------------------------------------------
    def find_relative_throughput(self):
        """
        Determine the relative transmission of each spectrum
        from the median flux in sky lines (assumed constant).
        """
        line_throughput = []
        line_throughput_var = []
        for line_flux, line_var in zip(self.sky_line_flux, self.sky_line_variance):
            norm = np.nanmedian(line_flux)
            line_throughput.append(line_flux/norm)
            line_throughput_var.append(line_var/norm)
        self.relative_throughput = np.nanmedian(line_throughput, axis=0)
        relative_throughput_var = np.nanmean(line_throughput_var, axis=0)
        global_mean_throughput = np.nanmean(self.relative_throughput)
        global_var_throughput = np.nanvar(self.relative_throughput)
        print "\n> Relative throughput: {:.3f} +- {:.3f} using {} skypl lines" \
              .format(np.nanmean(self.relative_throughput),
                      np.nanvar(self.relative_throughput),
                      len(self.sky_line_flux))
        print "                       (min, max) = ({:.3f}, {:.3f}), {} NaNs" \
              .format(np.nanmin(self.relative_throughput),
                      np.nanmax(self.relative_throughput),
                      np.isnan(self.relative_throughput).sum())
        print "  overall flux correction: {:.3f}" \
              .format(np.nansum(self.intensity*self.relative_throughput[:, np.newaxis]) /
                      np.nansum(self.intensity))

        global_var_throughput *= 3.
        new_inverse_variance = 1/relative_throughput_var + 1/global_var_throughput
        self.relative_throughput /= relative_throughput_var
        self.relative_throughput += global_mean_throughput/global_var_throughput
        self.relative_throughput /= new_inverse_variance
        flux_correction = np.nansum(self.intensity*self.relative_throughput[:, np.newaxis]) / np.nansum(self.intensity)
        self.relative_throughput /= flux_correction
        print "  Relative throughput: {:.3f} +- {:.3f} using {} skypl lines" \
              .format(np.nanmean(self.relative_throughput),
                      np.nanvar(self.relative_throughput),
                      len(self.sky_line_flux))
        print "                       (min, max) = ({:.3f}, {:.3f}), {} NaNs" \
              .format(np.nanmin(self.relative_throughput),
                      np.nanmax(self.relative_throughput),
                      np.isnan(self.relative_throughput).sum())
        print "  overall flux correction: {:.3f}" \
              .format(np.nansum(self.intensity*self.relative_throughput[:, np.newaxis]) /
                      np.nansum(self.intensity))

#    # -------------------------------------------------------------------------
#    def RSS_map(self, variable, norm=colors.LogNorm(), list_spectra='all'):
#        if list_spectra == 'all':
#            list_spectra = range(self.n_spectra)
#
#        plt.figure(figsize=(10, 10))
#        plt.scatter(self.offset_RA_arcsec[list_spectra],
#                    self.offset_DEC_arcsec[list_spectra],
#                    c=variable[list_spectra], norm=norm, s=50)
#        plt.gray()
#        plt.colorbar()
#        plt.gca().invert_xaxis()
#        plt.show()
    # -------------------------------------------------------------------------
    def RSS_map(self, variable, list_spectra=[],
                new_figure=True, norm=colors.LogNorm()):
        """
        Plot map showing the offsets, coloured by variable.
        """
        if len(list_spectra) == 0:
            list_spectra = range(self.n_spectra)

        if(new_figure):
            plt.figure(figsize=(10, 10))
        plt.scatter(self.offset_RA_arcsec[list_spectra],
                    self.offset_DEC_arcsec[list_spectra],
                    c=variable[list_spectra], cmap=fuego_color_map, norm=norm,
                    s=75, marker="h")
        plt.title(self.description+" - RSS map")
        plt.xlabel("$\Delta$ RA [arcsec]")
        plt.ylabel("$\Delta$ DEC [arcsec]")
        plt.colorbar()
        plt.gca().invert_xaxis()
        if(new_figure):
            plt.show()
            plt.close()

    # -------------------------------------------------------------------------
    def do_extinction_curve(self, observatory_file='ssoextinct.dat', plot=False):
        data_observatory = np.loadtxt(observatory_file, unpack=True)
        extinction_curve_wavelenghts = data_observatory[0]
        extinction_curve = data_observatory[1]
        extinction_corrected_airmass = 10**(0.4*self.airmass*extinction_curve)

#        print data_observatory[1]
#        print " "
#        for i in range(len(extinction_curve)):
#            print extinction_curve_wavelenghts[i],extinction_curve[i],extinction_corrected_airmass[i]

        tck = interpolate.splrep(extinction_curve_wavelenghts, extinction_corrected_airmass, s=0)
        self.extinction_correction = interpolate.splev(self.wavelength, tck, der=0)
#        self.intensity_corrected_extinction = self.intensity * self.extinction_correction

        # Plot  
        if plot == True:
            plt.figure(figsize=(10, 5))
            plt.plot(extinction_curve_wavelenghts, extinction_corrected_airmass, '+')
            plt.xlim(np.min(self.wavelength),np.max(self.wavelength))
            cinco_por_ciento = 0.05 * (np.max(self.extinction)- np.min(self.extinction))
            plt.ylim(np.min(self.extinction)-cinco_por_ciento,np.max(self.extinction)+cinco_por_ciento)
            plt.plot(self.wavelength,self.extinction, "g")
            plt.minorticks_on() 
            plt.title('Correction for extinction using airmass = '+str(self.airmass))

        print "  Airmass = ", self.airmass
        print "  Observatory file with extinction curve :", observatory_file
        print "  Correction for extinction using airmass obtained!"


# -----------------------------------------------------------------------------
class KOALA_RSS(RSS):
    """
    This class reads the FITS files returned by
    `2dfdr
    <https://www.aao.gov.au/science/software/2dfdr>`_
    and performs basic analysis tasks (see description under each method).

    Parameters
    ----------
    filename : string
      FITS file returned by 2dfdr, containing the Raw Stacked Spectra.
      The code makes sure that it contains 1000 spectra
      with 2048 wavelengths each.

    Example
    -------
    >>> pointing1 = KOALA_RSS('data/16jan20058red.fits')
    > Reading file "data/16jan20058red.fits" ...
      2048 wavelength points between 6271.33984375 and 7435.43408203
      1000 spaxels
      These numbers are the right ones for KOALA!
      DONE!
    """

    # -------------------------------------------------------------------------
    def __init__(self, filename):

        # Create RSS object
        super(KOALA_RSS, self).__init__()

        print "\n> Reading file", '"'+filename+'"', "..."
        RSS_fits_file = fits.open(filename)  # Open file

#        General info:
        self.object = RSS_fits_file[0].header['OBJECT']
        self.description = self.object + ' - ' + filename
        self.RA_centre_deg = RSS_fits_file[2].header['CENRA'] * 180/np.pi
        self.DEC_centre_deg = RSS_fits_file[2].header['CENDEC'] * 180/np.pi
#        WARNING: Something is probably wrong/inaccurate here!
#                 Nominal offests between pointings are totally wrong!

        # Read spaxel positions on sky (converting radians to deg)
        # and good/bad ('P'/'N') flag
        all_spaxels = range(len(RSS_fits_file[2].data))
        quality_flag = [RSS_fits_file[2].data[i][11] for i in all_spaxels]
        good_spaxels = [i for i in all_spaxels if quality_flag[i] == 'P']
#        print quality_flag
#        print np.where(quality_flag!='P')

        # Create wavelength, intensity, and variance arrays
        wcsKOALA = WCS(RSS_fits_file[0].header)
        index_wave = np.arange(RSS_fits_file[0].header['NAXIS1'])
        wavelength = wcsKOALA.dropaxis(1).wcs_pix2world(index_wave, 0)[0]
        intensity = RSS_fits_file[0].data[good_spaxels]
        variance = RSS_fits_file[1].data[good_spaxels]
        offset_RA_arcsec = np.array([RSS_fits_file[2].data[i][5]
                                     for i in good_spaxels])
        offset_DEC_arcsec = np.array([RSS_fits_file[2].data[i][6]
                                      for i in good_spaxels])

        self.set_data(wavelength, intensity, variance,
                      offset_RA_arcsec, offset_DEC_arcsec)
        # Check that dimensions match KOALA numbers
        if self.n_wave != 2048 and len(all_spaxels) != 1000:
            print "\n *** WARNING *** : These numbers are NOT the standard ones for KOALA"

        # KOALA-specific stuff
        self.ID = np.array([RSS_fits_file[2].data[i][0] for i in good_spaxels])
        self.PA = RSS_fits_file[0].header['TEL_PA']
        self.grating = RSS_fits_file[0].header['GRATID']
        # Check RED / BLUE arm for AAOmega
        if (RSS_fits_file[0].header['SPECTID'] == "RD"):
            AAOmega_Arm = "RED"
        if (RSS_fits_file[0].header['SPECTID'] == "BL"):
            AAOmega_Arm = "BLUE"

        # Check if NARROW (spaxel_size = 0.7 arcsec)
        # or WIDE (spaxel_size=1.25) field of view
        # (if offset_max - offset_min > 31 arcsec in both directions)
        if np.max(self.offset_RA_arcsec)-np.min(self.offset_RA_arcsec) > 31 and \
           np.max(self.offset_DEC_arcsec)-np.min(self.offset_DEC_arcsec) > 31:
            self.spaxel_size = 1.25
            field = "WIDE"
        else:
            self.spaxel_size = 0.7
            field = "NARROW"

        # Print information from header
        print "  This is a KOALA '{}' file,".format(AAOmega_Arm), \
              "using grating '{}' in AAOmega".format(self.grating)
        print "  Object:", self.object
        print "  Field of view:", field, \
              "(spaxel size =", self.spaxel_size, "arcsec)"
        print "  Center position: (RA, DEC) = ({:.3f}, {:.3f}) degree" \
              .format(self.RA_centre_deg, self.DEC_centre_deg)
        print "  Position angle (PA) = {:.2f}".format(self.PA)

        # Get airmass and do extinction
        ZD = (RSS_fits_file[0].header['ZDSTART'] +
              RSS_fits_file[0].header['ZDEND']) / 2
        self.airmass = 1 / np.cos(np.radians(ZD))
        self.do_extinction_curve('ssoextinct.dat')

        print "  DONE!"

#   --------------------------------------------------------------------

def plot_star(pointing):
    star = np.argmax(pointing.flux_between(6635, 6645))
    plt.figure(figsize=(10, 5))
    plt.plot(pointing.wavelength, pointing.intensity[star], 'y-')
    plt.plot(pointing.wavelength, pointing.intensity[star]-pointing.sky_emission, 'r-')
    plt.plot(pointing.wavelength, pointing.sky_emission, 'c-')
    plt.xlim(6200, 6400)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot(pointing.wavelength, pointing.intensity[star], 'y-')
    plt.plot(pointing.wavelength, pointing.intensity[star]-pointing.sky_emission, 'r-')
    plt.plot(pointing.wavelength, pointing.sky_emission, 'c-')
    plt.xlim(6800, 7100)
    plt.show()


def plot_object(cube, index):
    plt.figure(figsize=(12, 12))

#        interpolated_map = np.nansum(cube.data, axis=0)
    interpolated_map = cube.data[index]
    plt.imshow(interpolated_map, origin='lower', interpolation='none', norm=colors.LogNorm())

    I_max = np.unravel_index(np.nanargmax(interpolated_map), interpolated_map.shape)
    print 'maximum_I:', I_max[1], I_max[0], interpolated_map.shape
    plt.plot(I_max[1], I_max[0], 'w+', ms=50, mew=3)

    weight = interpolated_map**2
    total_weight = np.nansum(weight)
    x_cm = np.sum(np.arange(interpolated_map.shape[1])*np.nansum(weight, axis=0)) / total_weight
    y_cm = np.sum(np.arange(interpolated_map.shape[0])*np.nansum(weight, axis=1)) / total_weight
    print 'initial cm:', x_cm, y_cm
    plt.plot(x_cm, y_cm, 'bx', ms=50, mew=3)
    
    x_cm, y_cm = find_cm(cube, index)
    print 'final cm:', x_cm, y_cm
    plt.plot(x_cm, y_cm, 'wx', ms=50, mew=3)

    plt.colorbar()
    plt.gca().invert_xaxis()
    plt.show()
    plt.close()


#    Find maximum intensity:
#    I_max = np.array([np.unravel_index(np.nanargmax(cube.data[wavelength]), (cube.n_rows, cube.n_cols)) for wavelength in range(cube.n_wave)])
#    x_max = I_max[:, 1]
#    y_max = I_max[:, 0]


# def growth_curve_at_wavelength(cube, wavelength):
#    median_filter_half_width = 75
#    lambda_min = max(wavelength-median_filter_half_width, 0)
#    lambda_max = min(wavelength+median_filter_half_width+1, cube.n_wave)
#    x_peak = np.median(cube.x_peak[lambda_min:lambda_max])
#    y_peak = np.median(cube.y_peak[lambda_min:lambda_max])
#
#    x = np.arange(cube.n_cols) - x_peak
#    y = np.arange(cube.n_rows) - y_peak
#    r2 = np.sum(np.meshgrid(x**2, y**2), axis=0)
#
#    sorted_by_distance = np.argsort(r2, axis=None)
#    print r2.shape, sorted_by_distance.shape
#
#    L_growth_curve = []
#    r_growth_curve = []
#    total_luminosity = 0.
#    for spaxel in sorted_by_distance:
#        index = np.unravel_index(spaxel, (cube.n_rows, cube.n_cols))
#        L = cube.data[wavelength][index]
##        print spaxel, r2[index], L, total_luminosity, np.isnan(L)
#        if np.isnan(L) == False:
#            total_luminosity += L
#            L_growth_curve.append(total_luminosity)
#            r_growth_curve.append(np.sqrt(r2[index]))
#
#    plt.plot(r_growth_curve, L_growth_curve)


#def half_light_spectrum(cube, plot=False):
#    r2_growth_curve, F_growth_curve, flux, r2_half_light = cube.growth_curve_between(0, 1e30)
#    print "\n> Computing growth-curve spectrum..."
#    intensity = []
#    smooth_x = signal.medfilt(cube.x_peak, 11)
#    smooth_y = signal.medfilt(cube.y_peak, 11)
#    for l in range(cube.n_wave):
#        wavelength = cube.RSS.wavelength[l]
#        if l % (cube.n_wave/10+1) == 0:
#            print "  {:.2f} Angstroms (wavelength {}/{})..." \
#                  .format(wavelength, l+1, cube.n_wave)
#        x = np.arange(cube.n_cols) - smooth_x[l]
#        y = np.arange(cube.n_rows) - smooth_y[l]
#        r2 = np.sum(np.meshgrid(x**2, y**2), axis=0)
#        spaxels = np.where(r2 < r2_half_light*25)
#        intensity.append(np.nansum(cube.data[l][spaxels]))
#    if plot:
#        plt.plot(cube.RSS.wavelength, intensity, 'b.', alpha=0.2)
#        plt.plot(cube.RSS.wavelength, signal.medfilt(intensity,11), 'r-')
##        plt.xlim(6200, 6400)
#        plt.xlim(6200, 7400)
##        plt.xlim(7200, 7400)
#        plt.ylim(np.nanpercentile(intensity, 1), np.nanpercentile(intensity, 99))
#        lambda_cal, flux_cal, delta_lambda = np.loadtxt('data/FLUX_CAL/feg21.dat', usecols=(0,1,3), unpack=True)
#        points = np.array([cube.fit_Moffat_between(wl-25, wl+25)[0] if wl > 6300 and wl < 7400 else np.NaN for wl in lambda_cal])
#        scale = np.nanmedian(points/flux_cal)
#        plt.plot(lambda_cal, flux_cal*scale, 'k*-')
#        plt.plot(lambda_cal, points, 'k+', ms=30, mew=3)
#    return np.array(intensity)


def KOALA_offsets(s, pa):
    print "\n> Offsets towards North and East between pointings," \
        "according to KOALA manual, for pa =", pa, "degrees"
    pa *= np.pi/180
    print "  a -> b :", s*np.sin(pa), -s*np.cos(pa)
    print "  a -> c :", -s*np.sin(60-pa), -s*np.cos(60-pa)
    print "  b -> d :", -np.sqrt(3)*s*np.cos(pa), -np.sqrt(3)*s*np.sin(pa)


def offset_between_cubes(cube1, cube2, plot=True):
    x = (cube2.x_peak - cube2.n_cols/2. + cube2.RA_centre_deg*3600./cube2.pixel_size_arcsec) \
        - (cube1.x_peak - cube1.n_cols/2. + cube1.RA_centre_deg*3600./cube1.pixel_size_arcsec)
    y = (cube2.y_peak - cube2.n_rows/2. + cube2.DEC_centre_deg*3600./cube2.pixel_size_arcsec) \
        - (cube1.y_peak - cube1.n_rows/2. + cube1.DEC_centre_deg*3600./cube1.pixel_size_arcsec)
    delta_RA_pix = np.nanmedian(x)
    delta_DEC_pix = np.nanmedian(y)
#    weight = np.nansum(cube1.data+cube2.data, axis=(1, 2))
#    total_weight = np.nansum(weight)
#    print "--- lambda=", np.nansum(cube1.RSS.wavelength*weight) / total_weight
#    delta_RA_pix = np.nansum(x*weight) / total_weight
#    delta_DEC_pix = np.nansum(y*weight) / total_weight
    delta_RA_arcsec = delta_RA_pix * cube1.pixel_size_arcsec
    delta_DEC_arcsec = delta_DEC_pix * cube1.pixel_size_arcsec
    print '(delta_RA, delta_DEC) = ({:.3f}, {:.3f}) arcsec' \
        .format(delta_RA_arcsec, delta_DEC_arcsec)
#    delta_RA_headers = (cube2.RSS.RA_centre_deg - cube1.RSS.RA_centre_deg) * 3600
#    delta_DEC_headers = (cube2.RSS.DEC_centre_deg - cube1.RSS.DEC_centre_deg) * 3600
#    print '                        ({:.3f}, {:.3f}) arcsec according to headers!!!???' \
#        .format(delta_RA_headers, delta_DEC_headers)
#    print 'difference:             ({:.3f}, {:.3f}) arcsec' \
#        .format(delta_RA-delta_RA_headers, delta_DEC-delta_DEC_headers)

    if plot:
        x -= delta_RA_pix
        y -= delta_DEC_pix
        smooth_x = signal.medfilt(x, 151)
        smooth_y = signal.medfilt(y, 151)

        plt.figure(figsize=(10, 5))
        wl = cube1.RSS.wavelength
        plt.plot(wl, x, 'k.', alpha=0.1)
        plt.plot(wl, y, 'r.', alpha=0.1)
        plt.plot(wl, smooth_x, 'k-')
        plt.plot(wl, smooth_y, 'r-')
    #    plt.plot(wl, x_max-np.nanmedian(x_max), 'g-')
    #    plt.plot(wl, y_max-np.nanmedian(y_max), 'y-')
        plt.ylim(-1.6, 1.6)
        plt.show()
        plt.close()

    return delta_RA_arcsec, delta_DEC_arcsec


def compare_cubes(cube1, cube2, line):
    l = np.searchsorted(cube1.RSS.wavelength, line)
    map1 = cube1.data[l]
    map2 = cube2.data[l]
    scale = np.nanmedian(map1+map2)*3

    plt.figure(figsize=(12, 8))
    plt.imshow(map1-map2, vmin=-scale, vmax=scale, cmap=plt.cm.get_cmap('RdBu'))
    plt.colorbar()
    plt.contour(map1, colors='w', linewidths=2, norm=colors.LogNorm())
    plt.contour(map2, colors='k', linewidths=1, norm=colors.LogNorm())
    plt.title("{:.2f} AA".format(line))
    plt.show()
    plt.close()

#    plt.figure(figsize=(12, 8))
#    plt.imshow(map1)
#    plt.colorbar()
#    plt.contour(map2)
#    plt.plot(cube1.x_peak[l], cube1.y_peak[l], 'w+', ms=150)
#    plt.plot(cube2.x_peak[l], cube2.y_peak[l], 'wx', ms=150)
#    plt.show()
#    plt.close()


def plot_response(cube, calibration_star_cubes):
    plt.figure(figsize=(12, 8))
    wavelength = cube.RSS.wavelength
    mean_curve = np.zeros_like(wavelength)
    for star in calibration_star_cubes:
        good = np.where(~np.isnan(star.response_curve))
        wl = star.response_wavelength[good]
        R = star.response_curve[good]
        mean_curve += np.interp(wavelength, wl, R)
        plt.plot(star.response_wavelength, star.response_curve,
                 label=star.RSS.description)
        print np.nanmean(star.response_curve)
    mean_curve /= len(calibration_star_cubes)
    plt.plot(wavelength, mean_curve, 'k.', label='mean response curve')
    plt.legend(frameon=False)
    plt.show()


def coord_range(rss_list):
    RA = [rss.RA_centre_deg+rss.offset_RA_arcsec/3600. for rss in rss_list]
    RA_min = np.nanmin(RA)
    RA_max = np.nanmax(RA)
    DEC = [rss.DEC_centre_deg+rss.offset_DEC_arcsec/3600. for rss in rss_list]
    DEC_min = np.nanmin(DEC)
    DEC_max = np.nanmax(DEC)
    return RA_min, RA_max, DEC_min, DEC_max


#def create_cube(rss):
#    rss.find_sky_emission(6635, 6645)
#    rss.find_relative_throughput()
#    cube = Interpolated_cube(rss, .3, 1.5)
#    return(cube)


if __name__ == "__main__":

    print "\nTesting KOALA RSS class..."
#    KOALA_offsets(1.25, 120)

#   --------------------------------------------------------------------
#    Calibration stars
#   --------------------------------------------------------------------

    stars_RSS = []
    stars_RSS.append(KOALA_RSS('data/NO_THROUGHPUT/16jan20046red.fits'))
    stars_RSS.append(KOALA_RSS('data/NO_THROUGHPUT/16jan20052red.fits'))
    stars_RSS.append(KOALA_RSS('data/NO_THROUGHPUT/16jan20064red.fits'))
    stars_cubes = [Interpolated_cube(rss, .3, 1.5) for rss in stars_RSS]
    stars_cubes[0].do_response_curve('data/FLUX_CAL/feg21.dat')
    stars_cubes[1].do_response_curve('data/FLUX_CAL/fhilt600.dat')
    stars_cubes[2].do_response_curve('data/FLUX_CAL/ffeige56.dat')

    plot_response(stars_cubes[0], stars_cubes)


#   --------------------------------------------------------------------
#    Pointings
#   --------------------------------------------------------------------

#    pointings_RSS = []
#    pointings_RSS.append(KOALA_RSS('data/16jan20058red.fits'))
#    pointings_RSS.append(KOALA_RSS('data/16jan20059red.fits'))
#    pointings_RSS.append(KOALA_RSS('data/16jan20060red.fits'))
#    pointings_cubes = []
#    for rss in pointings_RSS:
#        rss.find_sky_emission(6635, 6645)
#        rss.find_relative_throughput()
#        pointings_cubes.append(Interpolated_cube(rss, .3, 1.5))

#   --------------------------------------------------------------------

#    rss1 = KOALA_RSS('data/16jan20058red.fits')
#    rss1.find_sky_emission(6635, 6645)
#    cube1 = Interpolated_cube(rss1, .3, 1.5)
#    rss2 = KOALA_RSS('data/16jan20059red.fits')
#    rss2.find_sky_emission(6635, 6645)
#    cube2 = Interpolated_cube(rss2, .3, 1.5)


#    dx, dy = offset_between_cubes(cube1, cube2)
#    rss1.RA_centre_deg += dx/3600.
#    rss1.DEC_centre_deg += dy/3600.
#    cube1 = \
#        Interpolated_cube(rss1,
#                          cube2.pixel_size_arcsec,
#                          cube2.kernel_size_arcsec,
#                          [cube2.RA_centre_deg, cube2.DEC_centre_deg],
#                          [cube2.n_cols*cube2.pixel_size_arcsec,
#                           cube2.n_rows*cube2.pixel_size_arcsec])
#

#   --------------------------------------------------------------------
#    Other tests with real data
#   --------------------------------------------------------------------

#    plt.figure(figsize=(12,8))
#    intensity = half_light_spectrum(cube1s, True)
#    plt.show()


#    plt.figure(figsize=(12,8))
#    cube1s.fit_Moffat_between(6560, 6600, plot=True)
#    growth_curve(cube2s, 6560, 6600)
#    growth_curve(cube3s, 6560, 6600)
#    growth_curve(cube1s)
#    growth_curve(cube2s)
#    growth_curve(cube3s)
#    plt.show()

#    print cube1s.growth_curve_between(6750, 6850)[2], \
#        cube1s.growth_curve_between(6850, 6950)[2], \
#        cube1s.growth_curve_between(6950, 7050)[2]
#    print cube2s.growth_curve_between(6750, 6850)[2], \
#        cube2s.growth_curve_between(6850, 6950)[2], \
#        cube2s.growth_curve_between(6950, 7050)[2]
#    print cube3s.growth_curve_between(6750, 6850)[2], \
#        cube3s.growth_curve_between(6850, 6950)[2], \
#        cube3s.growth_curve_between(6950, 7050)[2]

#    plot_object(cube1s)
#    plot_object(cube2s)
#    plot_object(cube3s)
#    trace_object(cube1s)
#    trace_object(cube2s)
#    trace_object(cube3s)

#    pointing1 = KOALA_RSS('data/16jan20058red.fits')
#    pointing2 = KOALA_RSS('data/16jan20059red.fits')
#    pointing3 = KOALA_RSS('data/16jan20060red.fits')
#
#    pointing1.find_sky_emission(6635, 6645)
#    pointing2.find_sky_emission(6635, 6645)
#    pointing3.find_sky_emission(6635, 6645)
#
#    cube1 = Interpolated_cube(pointing1, .3, 1.5)
#    cube2 = Interpolated_cube(pointing2, .3, 1.5)
#    cube3 = Interpolated_cube(pointing3, .3, 1.5)

#    cube1.trace_peak(plot=True)
#    cube2.trace_peak(plot=True)
#    cube3.trace_peak(plot=True)

#    trace_object(cube2)
#    trace_object(cube3)
#    offset_between_cubes(cube1, cube2)
#    offset_between_cubes(cube2, cube3)
#    offset_between_cubes(cube3, cube1)

#   --------------------------------------------------------------------

#    plot_object(cube3, 621)

#    pointing1.RSS_map(pointing1.flux_between(6635, 6645))

#    plt.figure(figsize=(10, 5))
#    index_Ha = np.searchsorted(pointing.wavelength, 6639.5)
#    print index
#    cube.plot_wavelength(index)
#    plt.gca().invert_xaxis()
#    plt.colorbar()

#    med_throughput = np.mean([pointing1.throughput, pointing2.throughput, pointing3.throughput], axis=0)
#    plt.plot(pointing1.ID, pointing1.throughput/med_throughput, 'r.')
#    plt.plot(pointing2.ID, pointing2.throughput/med_throughput, 'g.')
#    plt.plot(pointing3.ID, pointing3.throughput/med_throughput, 'b.')

#    pointing1.plot_spectra([160],[6600,6700])
#    sky_spaxels = pointing1.find_sky(6635, 6645)
#    pointing1.RSS_map(pointing1.flux_between(6635, 6645) / np.nansum(pointing1.sky, axis=0), None, sky_spaxels)
#    sky1 = np.nanmedian(pointing1.intensity[sky_spaxels], axis=0)
#    sky2 = np.nanmedian(pointing1.throughput_corrected[sky_spaxels], axis=0)
#    plt.figure(figsize=(10, 5))
#    plt.plot(pointing1.wavelength, 100*(sky1/sky2-1), 'c.')
#    plt.plot(pointing1.wavelength, sky1-sky2, 'r.')
#    plt.show()
#
#    plt.figure(figsize=(10, 5))
#    plt.plot(pointing1.wavelength, sky1, 'c-')
#    plt.xlim(6600, 7400)
#    plt.show()
#
#    plt.figure(figsize=(10, 5))
#    plt.imshow(pointing1.intensity-sky1[np.newaxis, :], vmin=0, vmax=10)
#    plt.colorbar()
#    plt.show()
#
#    plt.figure(figsize=(10, 5))
#    plt.imshow(pointing1.throughput_corrected-sky1[np.newaxis, :], vmin=0, vmax=10)
#    plt.colorbar()
#    plt.show()


#    pointing2.find_sky(6635, 6645)
#    pointing3.find_sky(6635, 6645)


#    for i in range(pointing1.n_spectra):
#        print i+1, pointing1.offset_RA_arcsec[i], pointing1.offset_DEC_arcsec[i]
#        
#    plt.plot(pointing1.offset_RA_arcsec, pointing1.offset_DEC_arcsec, '+')
#    
#    plt.figure(figsize=(20,5))
#    pointing1.plot_spectra([156])
#    plt.xlim([6500,6800])
##    plt.ylim([10,1e4])
#    plt.yscale('log')
#    plt.show()

#    cube = Interpolated_cube(60, 60, 2048)
#    cube.add_spectra(pointing1, pointing1.offset_DEC_arcsec, pointing1.offset_RA_arcsec, 1)
#    cube.add_spectra(pointing2, pointing2.offset_DEC_arcsec+(pointing2.DEC_centre_deg-pointing1.DEC_centre_deg), pointing2.offset_RA_arcsec+(pointing2.RA_centre_deg-pointing1.RA_centre_deg), 1)
#    cube.add_spectra(pointing3, pointing3.offset_DEC_arcsec+(pointing3.DEC_centre_deg-pointing1.DEC_centre_deg), pointing3.offset_RA_arcsec+(pointing3.RA_centre_deg-pointing1.RA_centre_deg), 1)
#
#    plt.figure(figsize=(10, 5))
#    index_Ha = np.searchsorted(pointing1.wavelength, 6639.5)
#    print index
#    cube.plot_wavelength(index)
#    plt.gca().invert_xaxis()
#    plt.colorbar()

#    pointing1.RSS_map(pointing1.throughput, None)
#    pointing1.RSS_map(pointing1.sky[0], None)
#    pointing1.RSS_map(pointing1.sky[1], None)
#    pointing1.RSS_map(pointing1.sky[2], None)

#    star1 = KOALA_RSS('data/NO_THROUGHPUT/16jan20046red.fits')
#    star2 = KOALA_RSS('data/NO_THROUGHPUT/16jan20052red.fits')
#    star3 = KOALA_RSS('data/NO_THROUGHPUT/16jan20064red.fits')
##    star1 = KOALA_RSS('data/16jan20046red.fits')
##    star2 = KOALA_RSS('data/16jan20052red.fits')
##    star3 = KOALA_RSS('data/16jan20064red.fits')

#    star1.find_sky_emission(6635, 6645)
#    star2.find_sky_emission(6635, 6645)
#    star3.find_sky_emission(6635, 6645)
#
#    star1.find_relative_throughput()
#    star2.find_relative_throughput()
#    star3.find_relative_throughput()
#
#    cube1s = Interpolated_cube(star1, .3, 1.5)
#    cube2s = Interpolated_cube(star2, .3, 1.5)
#    cube3s = Interpolated_cube(star3, .3, 1.5)
#
#    cube1s.trace_peak(plot=True)
#    cube2s.trace_peak(plot=True)
#    cube3s.trace_peak(plot=True)
#
#    cube1s.do_response_curve('data/FLUX_CAL/feg21.dat', True)
#    cube2s.do_response_curve('data/FLUX_CAL/fhilt600.dat', True)
#    cube3s.do_response_curve('data/FLUX_CAL/ffeige56.dat', True)

#   --------------------------------------------------------------------
#    Tests with mock data
#   --------------------------------------------------------------------

#    wl = np.linspace(0, 2*np.pi, 100)
#    spec = np.sin(wl)**2
##    mock_spectrum = RSS(wl, spec)
##    mock_spectrum.plot_spectra()
#
#    smooth_map = Interpolated_cube(50, 100, len(wl))
#    smooth_map.add_spectrum(spec, -25, -20, 5)
#    smooth_map.add_spectrum(2*spec, 0, 0, 5)
#    smooth_map.add_spectrum(-spec, -25, 20, 5)

#   --------------------------------------------------------------------

#    wl = np.linspace(0, 2*np.pi, 20)
#    spec = np.sin(wl)**2
#    spec = np.outer([1, 2, -1], spec)
#    
#    mock1 = RSS("mock1", wl, spec, .05*spec)
#    mock1.offset_RA_arcsec = np.array([-25., 0., -25.])
#    mock1.offset_DEC_arcsec = np.array([5., 0., 10.])
#    mock1.RA_centre_deg = 12./3600
#    mock1.DEC_centre_deg = 2./3600
#    cube1 = Interpolated_cube(mock1, .2, 5)
#    
#    mock2 = RSS("mock1", wl, spec, .05*spec)
#    mock2.offset_RA_arcsec = mock1.offset_RA_arcsec
#    mock2.offset_DEC_arcsec = mock1.offset_DEC_arcsec
#    cube2 = Interpolated_cube(mock2, .2, 5)
#
##    mock1.plot_spectra()
##    plt.show()
##    plt.close()
##    cube1.plot_wavelength(1., None)
#
#    dx, dy = offset_between_cubes(cube1, cube2)
#
#    mock1.RA_centre_deg += dx/3600.
#    mock1.DEC_centre_deg += dy/3600.
#    
#    cube1 = Interpolated_cube(mock1,
#                              cube2.pixel_size_arcsec, cube2.kernel_size_arcsec,
#                              [cube2.RA_centre_deg, cube2.DEC_centre_deg],
#                              [cube2.n_cols*cube2.pixel_size_arcsec,
#                               cube2.n_rows*cube2.pixel_size_arcsec])
#
#    dx, dy = offset_between_cubes(cube1, cube2)

#   --------------------------------------------------------------------
#    New peak finder
#   --------------------------------------------------------------------

#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20046red.fits')
#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20052red.fits')
#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20064red.fits')
#    rss = KOALA_RSS('data/16jan20058red.fits')
#    rss = KOALA_RSS('data/16jan20059red.fits')
#    rss = KOALA_RSS('data/16jan20060red.fits')
#    cube = create_cube(rss)

#    x = np.arange(cube.n_cols)
#    y = np.arange(cube.n_rows)
#    weight = np.nan_to_num(cube.data)
#    mean_iamge = np.nanmean(weight, axis=0)
#    mean_iamge /= np.nanmean(mean_iamge)
#    weight *= mean_iamge[np.newaxis, :, :]
#    xw = x[np.newaxis, np.newaxis, :] * weight
#    yw = y[np.newaxis, :, np.newaxis] * weight
#    w = np.nansum(weight, axis=(1, 2))
#    x_peak = np.nansum(xw, axis=(1, 2)) / w
#    y_peak = np.nansum(yw, axis=(1, 2)) / w
#
#    plt.figure(figsize=(10, 5))
#    wl = cube.RSS.wavelength
#    x = (x_peak-np.nanmedian(x_peak))*cube.pixel_size_arcsec
#    y = (y_peak-np.nanmedian(y_peak))*cube.pixel_size_arcsec
#    plt.plot(wl, x, 'k.', alpha=0.2)
#    plt.plot(wl, y, 'r.', alpha=0.2)
#    odd_number = 2 * int(np.sqrt(cube.n_wave) / 2) + 1
#    plt.plot(wl, signal.medfilt(x, odd_number), 'k-')
#    plt.plot(wl, signal.medfilt(y, odd_number), 'r-')
##    plt.plot(wl[7:-7], running_mean(x, 15), 'k--')
##    plt.plot(wl[7:-7], running_mean(y, 15), 'r--')
#    hi = np.max([np.nanpercentile(x, 90), np.nanpercentile(y, 90)])
#    lo = np.min([np.nanpercentile(x, 10), np.nanpercentile(y, 10)])
##    plt.ylim(lo, hi)
#    plt.ylim(-.2, .2)
#    plt.ylabel("offset [arcsec]")
#    plt.xlabel("wavelength [$\AA$]")
#    plt.title(cube.RSS.description)
#    plt.show()
#    plt.close()
#    print "  Peak coordinates, in pixel: ({:.2f}+-{:.2f}, {:.2f}+-{:.2f})"\
#          .format(np.nanmedian(x_peak), np.std(x_peak),
#                  np.nanmedian(y_peak), np.std(y_peak))

#   --------------------------------------------------------------------
#    New sky subtraction
#   --------------------------------------------------------------------

#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20046red.fits')
#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20052red.fits')
#    rss = KOALA_RSS('data/NO_THROUGHPUT/16jan20064red.fits')
#    rss = KOALA_RSS('data/16jan20058red.fits')
#    rss = KOALA_RSS('data/16jan20059red.fits')
#    rss = KOALA_RSS('data/16jan20060red.fits')


#    plt.plot(flux_sky[sort_by_flux]-flux_object[sort_by_flux], 'r-')
#    plt.show()

#   --------------------------------------------------------------------
#    Registration tests with real data
#   --------------------------------------------------------------------

#    pointings_RSS = []
#    pointings_RSS.append(KOALA_RSS('data/16jan20058red.fits'))
#    pointings_RSS.append(KOALA_RSS('data/16jan20059red.fits'))
#    pointings_RSS.append(KOALA_RSS('data/16jan20060red.fits'))
#    pointings_cubes = [Interpolated_cube(rss, .3, 1.5) for rss in pointings_RSS]

#   --------------------------------------------------------------------

#    x01, y01 = offset_between_cubes(pointings_cubes[0], pointings_cubes[1])
#    x12, y12 = offset_between_cubes(pointings_cubes[1], pointings_cubes[2])
#    x20, y20 = offset_between_cubes(pointings_cubes[2], pointings_cubes[0])
#    print x01, x12, x20, x01+x12+x20
#    print y01, y12, y20, y01+y12+y20
#
#    pointings_RSS[0].RA_centre_deg += (x01-x20)/3/3600.
#    pointings_RSS[0].DEC_centre_deg += (y01-y20)/3/3600.
#    pointings_RSS[1].RA_centre_deg += (x12-x01)/3/3600.
#    pointings_RSS[1].DEC_centre_deg += (y12-y01)/3/3600.
#    pointings_RSS[2].RA_centre_deg += (x20-x12)/3/3600.
#    pointings_RSS[2].DEC_centre_deg += (y20-y12)/3/3600.
#
#    RA_min, RA_max, DEC_min, DEC_max = coord_range(pointings_RSS)
#    RA_centre_deg = (RA_min + RA_max)/2.
#    DEC_centre_deg = (DEC_min + DEC_max)/2.
#    RA_size_arcsec = 1.1*(RA_max - RA_min)*3600.
#    DEC_size_arcsec = 1.1*(DEC_max - DEC_min)*3600.
#    pointings_cubes = [Interpolated_cube(rss, .3, 1.5,
#                                         centre_deg=[RA_centre_deg, DEC_centre_deg],
#                                         size_arcsec=[RA_size_arcsec, DEC_size_arcsec])
#                       for rss in pointings_RSS]

#   --------------------------------------------------------------------

#    pointings_cubes[0].plot_wavelength(6640, save_file="Ha-0.pdf")
#    pointings_cubes[1].plot_wavelength(6640, save_file="Ha-1.pdf")
#    pointings_cubes[2].plot_wavelength(6640, save_file="Ha-2.pdf")

#    compare_cubes(pointings_cubes[0], pointings_cubes[1], 6640)
#    compare_cubes(pointings_cubes[1], pointings_cubes[2], 6640)
#    compare_cubes(pointings_cubes[2], pointings_cubes[0], 6640)

#   --------------------------------------------------------------------

#    plt.figure(figsize=(12, 12))
#    plt.imshow(np.nanmedian(pointings_cubes[0].data/pointings_cubes[1].data, axis=0))
#    plt.colorbar()
#    plt.show()
#    plt.close()
#
#    plt.figure(figsize=(12, 12))
#    plt.imshow(np.nanmedian(pointings_cubes[2].data/pointings_cubes[1].data, axis=0))
#    plt.colorbar()
#    plt.show()
#    plt.close()

#   --------------------------------------------------------------------

#    norm0 = np.nanmean(pointings_cubes[0].data, axis=(1, 2))
#    norm1 = np.nanmean(pointings_cubes[1].data, axis=(1, 2))
#    norm2 = np.nanmean(pointings_cubes[2].data, axis=(1, 2))
#    wl = pointings_cubes[0].RSS.wavelength
#    plt.figure(figsize=(10, 5))
#    plt.plot(wl, norm0, 'k:', label=pointings_cubes[0].RSS.description)
#    plt.plot(wl, norm1, 'r--', label=pointings_cubes[1].RSS.description)
#    plt.plot(wl, norm2, 'b-', label=pointings_cubes[2].RSS.description)
#    plt.yscale('log')
#    plt.ylim(3, 300)
#    plt.legend(frameon=False)
#    plt.show()
#    plt.close()
#
#    norm0 = np.nansum(pointings_cubes[0].data)
#    norm1 = np.nansum(pointings_cubes[1].data)
#    norm2 = np.nansum(pointings_cubes[2].data)
#    print norm0, norm1, norm2
#    norm = np.mean([norm0, norm1, norm2])
#    pointings_cubes[0].data *= norm/norm0
#    pointings_cubes[1].data *= norm/norm1
#    pointings_cubes[2].data *= norm/norm2
#    compare_cubes(pointings_cubes[0], pointings_cubes[1], 6640)
#    compare_cubes(pointings_cubes[1], pointings_cubes[2], 6640)
#    compare_cubes(pointings_cubes[2], pointings_cubes[0], 6640)
    
#   --------------------------------------------------------------------

#    combined_data = np.array([cube.data for cube in pointings_cubes])
#    median_cube = np.nanmedian(combined_data, axis=0)
#    wl = pointings_cubes[0].RSS.wavelength
#    SED = np.nanmean(median_cube, axis=(1, 2))
#    plt.figure(figsize=(10, 5))
#    plt.plot(wl, SED, 'k-')
##    plt.yscale('log')
#    plt.ylim(-50, 250)
##    plt.xlim(6250, 6500)
#    plt.show()
#    plt.close()
#
#    combined_cube = copy.deepcopy(pointings_cubes[0])
#    combined_cube.data = median_cube
#    combined_cube.plot_wavelength(6640, save_file="Ha-median.pdf")
#    combined_cube.data = np.nanmean(combined_data, axis=0)
#    combined_cube.plot_wavelength(6640, save_file="Ha-mean.pdf")

#   --------------------------------------------------------------------

#    plt.figure(figsize=(12, 12))
#    plt.imshow(np.nanmedian(pointings_cubes[0].data/pointings_cubes[1].data, axis=0))
#    plt.colorbar()
#    plt.show()
#    plt.close()
#
#    plt.figure(figsize=(12, 12))
#    plt.imshow(np.nanmedian(pointings_cubes[2].data/pointings_cubes[1].data, axis=0))
#    plt.colorbar()
#    plt.show()
#    plt.close()

# -----------------------------------------------------------------------------
#                                                    ... Paranoy@ Rulz! ;^D
# -----------------------------------------------------------------------------
