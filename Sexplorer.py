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
    import ttk
else:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as tkFileDialog

matplotlib.use("TkAgg")
color_normalisation = {'linear': matplotlib.colors.Normalize,
                       'log': matplotlib.colors.LogNorm,
                       }


class MainApp(ttk.Frame):
    # -------------------------------------------------------------------------

    def __init__(self, parent):

        # --- set up the window

        ttk.Frame.__init__(self, parent)
        self.parent_window = parent
        self.parent_window.title('SELGIFS explorer')
        self.parent_window.bind('<Escape>', lambda f: root.destroy())

        # --- buttons

        self.controls = ttk.Frame(self)
        ttk.Button(self.controls,
                   text='New RSS file',
                   command=self.open_file).pack(fill=tk.X)
#        ttk.Button(self.controls,
#                   text='Fibre throughput',
#                   command=self.fibre_throughput).pack(fill=tk.X)
        ttk.Button(self.controls,
                   text='Sky subtraction',
                   command=self.sky_subtraction).pack(fill=tk.X)
        ttk.Button(self.controls,
                   text='Configure plots',
                   command=self.configure_plots).pack(fill=tk.X)

        self.active_fibres = []
        self.pick_many_fibres = tk.IntVar()
        ttk.Checkbutton(self.controls,
                        text="Pick many fibre",
                        variable=self.pick_many_fibres).pack(fill=tk.X)

        self.controls.pack(side=tk.LEFT)

        # --- plots

        self.plotting_area = ttk.Frame(self)
        self.figure = {}
        self.canvas = {}

        self.add_panel('RA - DEC', (9, 8), 0, 0)
        self.add_panel('Spectra', (19, 8), 0, 1)
        self.plotting_area.pack(side=tk.LEFT)

        self.active_fibres = []
        self.canvas['RA - DEC'].mpl_connect('button_press_event', self.pick_fibre)
#        self.canvas['RA - DEC'].mpl_connect('key_press_event', self.RA_DEC_key_press)
#        self.canvas['RA - DEC'].mpl_connect('key_press_release', self.RA_DEC_key_release)
        self.canvas['Spectra'].mpl_connect('button_press_event', self.pick_band)

        # --- variables and parameters

        self.parameters = {}

#        self.open_file('')
        self.open_file('22may20047red.fits')

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
        self.update_parameters()

    # -------------------------------------------------------------------------

    def fibre_throughput(self):
        self.rss.find_relative_throughput()
        self.rss.intensity /= self.rss.relative_throughput
        self.update_parameters()

    def sky_subtraction(self):
        self.rss.find_sky_emission()
        self.rss.intensity -= self.rss.sky_emission
        self.update_parameters()

    # -------------------------------------------------------------------------

    def configure_plots(self):
        self.dialog_window = tk.Toplevel(self)
        self.dialog_window.title('Plot parameters')

        row = 0
        ttk.Label(self.dialog_window,
                  text='Enter your preferred values for plot ranges').grid(
                  row=row, column=0, columnspan=4)
        row = row+1
        ttk.Label(self.dialog_window,
                  text='Leave it blank for automatic guess').grid(
                  row=row, column=0, columnspan=4)
        row = row+1
        ttk.Separator(self.dialog_window, orient=tk.HORIZONTAL).grid(
                      row=row, column=0, columnspan=4, sticky='nsew')
        row = row+1

        entry = {}
        ttk.Label(self.dialog_window, text='RA').grid(row=row, column=0)
        entry['RA_min'] = ttk.Entry(self.dialog_window)
        entry['RA_min'].grid(row=row, column=1)
        entry['RA_min'].insert(0, self.parameters['RA_min'])
        entry['RA_max'] = ttk.Entry(self.dialog_window)
        entry['RA_max'].grid(row=row, column=2)
        entry['RA_max'].insert(0, self.parameters['RA_max'])
        ttk.Label(self.dialog_window, text='Coordinate range').grid(row=row, column=3)
        row = row+1
        ttk.Label(self.dialog_window, text='DEC').grid(row=row, column=0)
        entry['DEC_min'] = ttk.Entry(self.dialog_window)
        entry['DEC_min'].grid(row=row, column=1)
        entry['DEC_min'].insert(0, self.parameters['DEC_min'])
        entry['DEC_max'] = ttk.Entry(self.dialog_window)
        entry['DEC_max'].grid(row=row, column=2)
        entry['DEC_max'].insert(0, self.parameters['DEC_max'])
        row = row+1
        ttk.Label(self.dialog_window, text='wavelength').grid(row=row, column=0)
        entry['wl_min'] = ttk.Entry(self.dialog_window)
        entry['wl_min'].grid(row=row, column=1)
        entry['wl_min'].insert(0, self.parameters['wl_min'])
        entry['wl_max'] = ttk.Entry(self.dialog_window)
        entry['wl_max'].grid(row=row, column=2)
        entry['wl_max'].insert(0, self.parameters['wl_max'])
        ttk.Label(self.dialog_window, text='Range of spectrum plot').grid(row=row, column=3)
        row = row+1
        ttk.Label(self.dialog_window, text='Continuum band').grid(row=row, column=0)
        entry['band_min'] = ttk.Entry(self.dialog_window)
        entry['band_min'].grid(row=row, column=1)
        entry['band_min'].insert(0, self.parameters['band_min'])
        entry['band_max'] = ttk.Entry(self.dialog_window)
        entry['band_max'].grid(row=row, column=2)
        entry['band_max'].insert(0, self.parameters['band_max'])
        ttk.Label(self.dialog_window, text='To compute mean intensities').grid(row=row, column=3)
        row = row+1
        ttk.Label(self.dialog_window, text='Colour map').grid(row=row, column=0)
        entry['map_min'] = ttk.Entry(self.dialog_window)
        entry['map_min'].grid(row=row, column=1)
        entry['map_min'].insert(0, self.parameters['map_min'])
        entry['map_max'] = ttk.Entry(self.dialog_window)
        entry['map_max'].grid(row=row, column=2)
        entry['map_max'].insert(0, self.parameters['map_max'])
        ttk.Label(self.dialog_window, text='Intensity range in colour map').grid(row=row, column=3)
        row = row+1
        ttk.Label(self.dialog_window, text='exponent').grid(row=row, column=0)
        entry['exponent'] = ttk.Entry(self.dialog_window)
        entry['exponent'].grid(row=row, column=1)
        entry['exponent'].insert(0, self.parameters['exponent'])
        ttk.Label(self.dialog_window, text='for colour scale').grid(row=row, column=3)
        row = row+1
        ttk.Label(self.dialog_window, text='Spectra').grid(row=row, column=0)
        entry['spec_min'] = ttk.Entry(self.dialog_window)
        entry['spec_min'].grid(row=row, column=1)
        entry['spec_min'].insert(0, self.parameters['spec_min'])
        entry['spec_max'] = ttk.Entry(self.dialog_window)
        entry['spec_max'].grid(row=row, column=2)
        entry['spec_max'].insert(0, self.parameters['spec_max'])
        ttk.Label(self.dialog_window, text='Intensity range in spectrum plot').grid(row=row, column=3)
        row = row+1
#        for thing in self.parameters:
#            ttk.Label(self.dialog_window, text=thing).grid(row=row, column=0)
#            entry[thing] = ttk.Entry(self.dialog_window)
#            entry[thing].grid(row=row, column=1)
##            entry[thing].delete(0, ttk.END)
#            entry[thing].insert(0, self.parameters[thing])
#            row = row+1

        ttk.Button(self.dialog_window, text="OK",
                  command=lambda x=entry: self.config_OK(entry=x)).grid(row=row, column=0)
        self.dialog_window.bind('<Return>', lambda event, x=entry: self.config_OK(event, x))

        ttk.Button(self.dialog_window, text="Cancel",
                  command=lambda x=entry: self.config_OK(x)).grid(row=row, column=1)
        self.dialog_window.bind('<Escape>', lambda f: self.dialog_window.destroy())

    def config_OK(self, event=None, entry={}):
        values = {}
        for thing in entry:
            values[thing] = entry[thing].get()
        self.dialog_window.destroy()
        self.update_parameters(values)

    def update_parameters(self, new_values={}):
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
        
        self.reset_map(new_values)

        spec = np.nansum(self.rss.intensity, axis=0)/np.nansum(self.map_values)
        self.set_parameter(new_values, 'spec_min', np.nanpercentile(spec, 1))
        self.set_parameter(new_values, 'spec_max', np.nanpercentile(spec, 99))

        self.update_plots()

    def reset_map(self, new_values={}):
        self.map_values = \
            self.rss.flux_between(self.parameters['band_min'],
                                  self.parameters['band_max'])[0] / \
            (self.parameters['band_max']-self.parameters['band_min'])

        map_min = np.nanpercentile(self.map_values, 1)
        map_med = np.nanpercentile(self.map_values, 50)
        map_max = np.nanpercentile(self.map_values, 99)
        self.set_parameter(new_values, 'map_min', map_min)
        self.set_parameter(new_values, 'map_max', map_max)
        self.set_parameter(new_values, 'exponent',
                           np.log(0.5)/np.log((map_med-map_min)/(map_max-map_min)))

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
            self.parameters[key] = float(value)
#            print('    v=', value)
#        print('    ', self.parameters.get(key, 'NOPE'))

    # -------------------------------------------------------------------------

    def update_plots(self):
        self.update_RA_DEC_plot()
        self.update_Spectra_plot()

    def update_Spectra_plot(self):
        wl_min = self.parameters['wl_min']
        wl_max = self.parameters['wl_max']
        spec_min = self.parameters['spec_min']
        spec_max = self.parameters['spec_max']
        band_min = self.parameters['band_min']
        band_max = self.parameters['band_max']
#        map_values = self.map_values

        plt.figure('Spectra')
        plt.clf()
        plot_axis = plt.axes((.05, .1, .94, .85))  # left,bottom,width,height
        plot_axis.set_xlim(wl_min, wl_max)
        plot_axis.set_ylim(spec_min, spec_max)
        plot_axis.set_xlabel("$\lambda [\AA]$")
        plot_axis.set_ylabel("normalised intensity")

        plot_axis.axvspan(band_min, band_max, alpha=0.1, color='c')

        total_spec = np.nansum(self.rss.intensity, axis=0)
        norm_total = np.nansum(self.map_values)
        total_spec /= norm_total
        plot_axis.plot(self.rss.wavelength, total_spec, 'b:',
                       label='total ({:.3e})'.format(norm_total))

        if(len(self.active_fibres) > 0):
            active_spec = np.nansum(self.rss.intensity[self.active_fibres], axis=0)
            norm_active = np.nansum(self.map_values[self.active_fibres])
            active_spec /= norm_active
            if(len(self.active_fibres) > 10):
                lbl = 'fibres ({:.3e})'.format(norm_active)
            else:
                lbl = 'fibres {} ({:.3e})'.format(self.active_fibres, norm_active)
            plot_axis.plot(self.rss.wavelength, active_spec, 'k-', label=lbl)

        plot_axis.legend()
        self.canvas['Spectra'].draw()

    def pick_band(self, event):
        try:
#            print(event.xdata)
#            if(event.xdata < .5*(self.parameters['band_min']+self.parameters['band_max'])):
#                print(.5*(self.parameters['band_min']+self.parameters['band_max']))
#                self.parameters['band_min'] = event.xdata
#            else:
#                print('BAND')
#                self.parameters['band_max'] = event.xdata
            if(event.xdata < self.parameters['band_min']):
                self.parameters['band_min'] = event.xdata
            elif(event.xdata > self.parameters['band_max']):
                self.parameters['band_max'] = event.xdata
            elif(event.button == 1):
                self.parameters['band_min'] = event.xdata
            else:
                self.parameters['band_max'] = event.xdata
        except:
            print('Clicked outside plot')
        self.reset_map()
        self.update_plots()

    def update_RA_DEC_plot(self):
        RA_min = self.parameters['RA_min']
        RA_max = self.parameters['RA_max']
        DEC_min = self.parameters['DEC_min']
        DEC_max = self.parameters['DEC_max']

#        map_values = self.parameters['map_values']
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
                                       c=self.map_values,
                                       norm=matplotlib.colors.PowerNorm(exponent, clip=False),
                                       vmin=map_min, vmax=map_max,
                                       cmap=matplotlib.cm.terrain,
                                       s=.5*arcsec2_to_pt2,
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
                       markersize=.6*np.sqrt(arcsec2_to_pt2),
                       marker='h')

        fig.colorbar(colour_map, cax=cbar_axis)

        self.canvas['RA - DEC'].draw()

    def pick_fibre(self, event):
        try:
            fibre = np.argmin(
                (self.rss.offset_RA_arcsec-event.xdata)**2 +
                (self.rss.offset_DEC_arcsec-event.ydata)**2)
#            print('Fibre', fibre,
#                  'ID =', self.rss.ID[fibre],
#                  'is the closest to', event.xdata, event.ydata)
            print(self.pick_many_fibres.get())
            if(self.pick_many_fibres.get() == 1):
                if(fibre in self.active_fibres):
                    self.active_fibres.remove(fibre)
                else:
                    self.active_fibres.append(fibre)
            else:
                self.active_fibres = [fibre]
        except:
            print('Clicked outside plot')
        self.update_plots()

#    def RA_DEC_key_press(self, event):
#        print('press', self.pick_many_fibres)
#        self.pick_many_fibres = (event.key == 'shift' or
#                                 event.key == 'control')
#        print(event.key, self.pick_many_fibres)
#
#    def RA_DEC_key_release(self, event):
#        print('release')
#        self.pick_many_fibres = False

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
    Sexplorer_window = MainApp(root).pack()
    root.mainloop()
    plt.close()

# -----------------------------------------------------------------------------
#                                                    ... Paranoy@ Rulz! ;^D
# -----------------------------------------------------------------------------
