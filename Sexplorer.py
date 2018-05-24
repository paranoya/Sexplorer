#!/usr/bin/python
# -*- coding: utf-8 -*-
# First attempt to show a KOALA RSS file, based on the
# script to estimate AAT focus by Angel R. Lopez-Sanchez (AAO/MQU)
# Yago Ascasibar (AAT, 19/05/2018)

from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import KOALA_analysis as KOALA
import copy

# Tkinter is for python 2; tkinter is for python 3
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFont
    import tkFileDialog
else:
    import tkinter as tk
    from tkinter import filedialog as tkFileDialog

matplotlib.use("TkAgg")
color_normalisation = {'linear': matplotlib.colors.Normalize,
                       'log': matplotlib.colors.LogNorm,
                       }


class MainApp(tk.Frame):
    # -------------------------------------------------------------------------

    def __init__(self, parent):

        # --- set up the window

        tk.Frame.__init__(self, parent)
        self.parent_window = parent
        self.parent_window.title('SELGIFS explorer')
        self.parent_window.bind('<Escape>', lambda f: root.destroy())

        # --- buttons

        self.controls = tk.Frame(self)
        tk.Button(self.controls,
                  text='New RSS file',
                  command=self.open_file).grid(row=0, sticky='ew')
        tk.Button(self.controls,
                  text='Sky subtraction',
                  command=self.sky_subtraction).grid(row=1, sticky='ew')
        self.controls.grid(row=0, column=0)

        # --- plots

        self.plotting_area = tk.Frame(self)
        self.figure = {}
        self.canvas = {}

        self.add_panel('RA - DEC', (10, 8), 0, 0)
        self.add_panel('Spectra', (20, 8), 0, 1)
        self.plotting_area.grid(row=0, column=1)

        self.active_fibres = []
        self.canvas['RA - DEC'].mpl_connect('button_press_event', self.pick_fibre)

        # --- variables and parameters

        self.parameters = {}
        self.active_fibres = []

#        self.open_file('')
        self.open_file('22may20047red.fits')

        self.update_parameters()

    def add_panel(self, plot_name, size, row, col):
        fig = plt.figure(num=plot_name, figsize=size, dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=self.plotting_area)
        canvas.get_tk_widget().grid(row=row, column=col)
        self.figure[plot_name] = fig
        self.canvas[plot_name] = canvas

    # -------------------------------------------------------------------------

    def open_file(self, filename=False):
        if(filename is False):
            filename = \
                tkFileDialog.askopenfilename(initialdir=".",
                                             title="Select file",
                                             filetypes=(("FITS files", "*.fits"),
                                                        ("all files", "*.*")))
        self.rss = KOALA.KOALA_RSS(filename)
        self.rss.intensity -= self.rss.sky_emission

        self.plot_range_original = {}
        l_min = self.rss.wavelength[0]
        l_max = self.rss.wavelength[-1]
        value = self.rss.flux_between(l_min, l_max)[0]
        self.plot_range_original['RA - DEC'] = {
            'wavelength': (l_min, l_max),
            'offset RA': (np.nanmin(self.rss.offset_RA_arcsec),
                          np.nanmax(self.rss.offset_RA_arcsec)),
            'offset DEC': (np.nanmin(self.rss.offset_DEC_arcsec),
                           np.nanmax(self.rss.offset_DEC_arcsec)),
            'value': (np.nanmin(value), np.nanmax(value)),
            }

        self.normalised_spectrum = np.nansum(self.rss.intensity, axis=0)
        np.savetxt(filename+'_wavelength.txt', self.rss.wavelength)
        np.savetxt(filename+'_spectrum.txt', self.normalised_spectrum)
        np.savetxt(filename+'_sky.txt', self.rss.sky_emission)
        self.normalised_spectrum *= self.rss.n_wave/np.nansum(self.normalised_spectrum)
        
        self.plot_range_original['Spectra'] = {
            'wavelength': (l_min, l_max),
            'value': (np.nanpercentile(self.normalised_spectrum, 1),
                      np.nanpercentile(self.normalised_spectrum, 99)),
            }

        self.plot_range_current = copy.deepcopy(self.plot_range_original)
#        self.update_plots()

    # -------------------------------------------------------------------------

    def sky_subtraction(self):
        self.rss.find_sky_emission()
        self.rss.intensity -= self.rss.sky_emission
        self.update_parameters()

    # -------------------------------------------------------------------------

    def plot_clicked(self, plot):
        print("PLOT CLICK --------------------", plot)
        self.dialog_window = tk.Toplevel(self)
        self.dialog_window.title('Configure plot:'+plot)

        tk.Label(self.dialog_window, text="plot").grid(row=0, column=1, columnspan=2)
        tk.Label(self.dialog_window, text="data").grid(row=0, column=3, columnspan=2)
        tk.Label(self.dialog_window, text="min").grid(row=1, column=1)
        tk.Label(self.dialog_window, text="max").grid(row=1, column=2)
        tk.Label(self.dialog_window, text="min").grid(row=1, column=3)
        tk.Label(self.dialog_window, text="max").grid(row=1, column=4)
        row = 2
        self.min_values = {}
        self.max_values = {}
        for thing in self.plot_range_original[plot]:
            tk.Label(self.dialog_window, text=thing).grid(row=row, column=0)
            self.min_values[thing] = tk.Entry(self.dialog_window)
            self.min_values[thing].grid(row=row, column=1)
            self.min_values[thing].delete(0, tk.END)
            self.min_values[thing].insert(0, self.plot_range_current[plot][thing][0])
            self.max_values[thing] = tk.Entry(self.dialog_window)
            self.max_values[thing].grid(row=row, column=2)
            self.max_values[thing].delete(0, tk.END)
            self.max_values[thing].insert(0, self.plot_range_current[plot][thing][1])
            tk.Label(self.dialog_window,
                     text=self.plot_range_original[plot][thing][0]).grid(row=row, column=3)
            tk.Label(self.dialog_window,
                     text=self.plot_range_original[plot][thing][1]).grid(row=row, column=4)
            row = row+1

        tk.Button(self.dialog_window, text="OK",
                  command=lambda x=plot: self.plot_changed(x)).grid(row=row, column=1)

    def plot_changed(self, plot):
        for thing in self.plot_range_original[plot]:
            self.plot_range_current[plot][thing] = \
                (np.float(self.min_values[thing].get()),
                 np.float(self.max_values[thing].get()))
        self.dialog_window.destroy()
        self.update_plots()

    def update_parameters(self, new_values={}):
        print('WWWW')
        self.set_parameter(new_values, 'RA_min',
                           np.nanmin(self.rss.offset_RA_arcsec))
        self.set_parameter(new_values, 'RA_max',
                           np.nanmax(self.rss.offset_RA_arcsec))
        self.set_parameter(new_values, 'DEC_min',
                           np.nanmin(self.rss.offset_DEC_arcsec))
        self.set_parameter(new_values, 'DEC_max',
                           np.nanmax(self.rss.offset_DEC_arcsec))

        l_min = self.rss.wavelength[0]
        l_max = self.rss.wavelength[-1]
        self.set_parameter(new_values, 'wl_min', l_min)
        self.set_parameter(new_values, 'wl_max', l_max)
        self.set_parameter(new_values, 'band_min', l_min)
        self.set_parameter(new_values, 'band_max', l_max)

        map_values = self.rss.flux_between(self.parameters['band_min'],
                                           self.parameters['band_max'])[0]
        map_min = np.nanpercentile(map_values, 1)
        map_med = np.nanpercentile(map_values, 50)
        map_max = np.nanpercentile(map_values, 99)
        self.set_parameter(new_values, 'map_values', map_values)
        self.set_parameter(new_values, 'map_min', map_min)
        self.set_parameter(new_values, 'map_max', map_max)
        self.set_parameter(new_values, 'exponent',
                           np.log(0.5)/np.log((map_med-map_min)/(map_max-map_min)))

        spec = np.nansum(self.rss.intensity, axis=0)/np.nansum(map_values)
        self.parameters['spec_min'] = np.nanpercentile(spec, 1)
        self.parameters['spec_max'] = np.nanpercentile(spec, 99)

        self.update_plots()

    def set_parameter(self, new_values, key, auto_value):
        """
        Update `self.parameters` to `new_values`.
        If the `key` is not defined in `new_values` (or it is blank),
        set it to the `auto_value` provided
        """
#        print('>>>>', key, new_values.get(key, 'NOPE'))
        value = new_values.get(key, '')
        if(value == ''):
            self.parameters[key] = auto_value
#            print('    a=', auto_value)
        else:
            self.parameters[key] = value
#            print('    v=', value)
#        print('    ', self.parameters.get(key, 'NOPE'))

    # -------------------------------------------------------------------------

    def update_plots(self):
        self.update_RA_DEC_plot()
        self.update_Spectra_plot()

    def update_Spectra_plot(self):
        print('sss')
        wl_min = self.parameters['wl_min']
        wl_max = self.parameters['wl_max']
        spec_min = self.parameters['spec_min']
        spec_max = self.parameters['spec_max']
        band_min = self.parameters['band_min']
        band_max = self.parameters['band_max']
        map_values = self.parameters['map_values']

        plt.figure('Spectra')
        plt.clf()
        plot_axis = plt.axes((.05, .1, .94, .85))  # left,bottom,width,height
        plot_axis.set_xlim(wl_min, wl_max)
        plot_axis.set_ylim(spec_min, spec_max)
        plot_axis.set_xlabel("$\lambda [\AA]$")
        plot_axis.set_ylabel("normalised intensity")

        plot_axis.axvspan(band_min, band_max, alpha=0.1, color='c')

        total_spec = np.nansum(self.rss.intensity, axis=0)
        norm_total = np.nansum(map_values)
        total_spec /= norm_total
        plot_axis.plot(self.rss.wavelength, total_spec, 'b:',
                       label='total ({:.3f})'.format(norm_total))

        if(len(self.active_fibres) > 0):
            active_spec = np.nansum(self.rss.intensity[self.active_fibres], axis=0)
            norm_active = np.nansum(map_values[self.active_fibres])
            active_spec /= norm_active
            if(len(self.active_fibres) > 3):
                lbl = 'fibres ({:.3f})'.format(norm_active)
            else:
                lbl = 'fibres {} ({:.3f})'.format(self.active_fibres, norm_active)
            plot_axis.plot(self.rss.wavelength, active_spec, 'k-', label=lbl)

        self.canvas['Spectra'].draw()

    def update_RA_DEC_plot(self):
        print('rrr')
        RA_min = self.parameters['RA_min']
        RA_max = self.parameters['RA_max']
        DEC_min = self.parameters['DEC_min']
        DEC_max = self.parameters['DEC_max']

        map_values = self.parameters['map_values']
        map_min = self.parameters['map_min']
        map_max = self.parameters['map_max']
        exponent = self.parameters['exponent']

        fig = plt.figure('RA - DEC')
        plt.clf()
        plot_axis = plt.axes((.1, .08, .7, .7/.8))
        cbar_axis = plt.axes((.85, .08, .05, .7/.8))
        plot_axis.cla()

        plot_axis.set_xlim(RA_min, RA_max)
        plot_axis.set_ylim(DEC_min, DEC_max)
        plot_axis.set_title(self.rss.description+" - RSS map")
        plot_axis.set_xlabel("$\Delta$ RA [arcsec]")
        plot_axis.set_ylabel("$\Delta$ DEC [arcsec]")

        plot_axis.set_aspect('equal')
        plot_axis.invert_xaxis()

        bbox = plot_axis.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        plot_area_pt2 = bbox.width*bbox.height*fig.dpi**2
        plot_area_arcsec2 = (RA_max-RA_min)*(DEC_max-DEC_min)
        arcsec2_to_pt2 = plot_area_pt2/plot_area_arcsec2

        colour_map = plot_axis.scatter(self.rss.offset_RA_arcsec,
                                       self.rss.offset_DEC_arcsec,
                                       c=map_values,
                                       norm=matplotlib.colors.PowerNorm(exponent),
                                       vmin=map_min, vmax=map_max,
                                       # cmap=fuego_color_map, norm=norm,
                                       s=.6*arcsec2_to_pt2,
                                       # ad-hoc scaling for KOALA
                                       # (should be fibre_area;
                                       # maybe colorbar?)
                                       marker="h")
        plot_axis.plot(self.rss.offset_RA_arcsec[self.active_fibres],
                       self.rss.offset_DEC_arcsec[self.active_fibres],
                       linestyle='none',
                       markerfacecolor='none',
                       markeredgecolor='black',
                       markeredgewidth=2,
                       markersize=np.sqrt(arcsec2_to_pt2),
                       marker='h')

        fig.colorbar(colour_map, cax=cbar_axis)

        self.canvas['RA - DEC'].draw()

    def pick_fibre(self, event):
        try:
            fibre = np.argmin(
                (self.rss.offset_RA_arcsec-event.xdata)**2 +
                (self.rss.offset_DEC_arcsec-event.ydata)**2)
            print('Fibre', fibre,
                  'ID =', self.rss.ID[fibre],
                  'is the closest to', event.xdata, event.ydata)
            if(fibre in self.active_fibres):
                self.active_fibres.remove(fibre)
            else:
                self.active_fibres.append(fibre)
        except:
            print('Clicked outside plot')
        self.update_plots()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    root = tk.Tk()
#    root.geometry("1000x730+10+10")
    root.resizable(0, 0)
    for font_name in tkFont.names():
        font = tkFont.nametofont(font_name)
        font.configure(size=12)
    Sexplorer_window = MainApp(root).pack(side=tk.TOP)
    root.mainloop()
    plt.close()

# -----------------------------------------------------------------------------
#                                                    ... Paranoy@ Rulz! ;^D
# -----------------------------------------------------------------------------
