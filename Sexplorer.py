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

    def __init__(self, parent):

        # --- set up the window

        tk.Frame.__init__(self, parent)
        self.parent_window = parent
        self.parent_window.title('SELGIFS explorer')
        self.parent_window.bind('<Escape>', lambda f: root.destroy())

        self.file_button = tk.Button(self, text='New RSS file',
                                     command=self.file_clicked)
        self.file_button.grid(row=0, column=0)

        # --- plotting area

        self.figure = {}
        self.canvas = {}
        
        self.figure['RA - DEC'] = plt.figure(figsize=(10, 8))  # to account for ~25% colorbar
        self.canvas['RA - DEC'] = \
            FigureCanvasTkAgg(self.figure['RA - DEC'], master=self)
        self.canvas['RA - DEC'].get_tk_widget().grid(row=0, column=1, rowspan=8)
        self.figure['RA - DEC'].canvas.mpl_connect('button_press_event', self.onclick)
        self.active_fibre = False

        self.figure['Spectra'] = plt.figure(figsize=(20, 8))
        self.canvas['Spectra'] = \
            FigureCanvasTkAgg(self.figure['Spectra'], master=self)
        self.canvas['Spectra'].get_tk_widget().grid(row=0, column=2, rowspan=8)

        # --- control buttons

        tk.Label(self, text="Configure plots:").grid(row=1, column=0)
        row_controls = 2
        for plot in self.figure:
            tk.Button(self, text=plot,  # state=tk.DISABLED,
                      command=lambda x=plot: self.plot_clicked(x)).grid(
                      row=row_controls, column=0)
            row_controls += 1

        self.file_clicked('22may20047red.fits')

    def onclick(self, event):
        try:
            self.active_fibre = np.argmin(
                (self.rss.offset_RA_arcsec-event.xdata)**2 +
                (self.rss.offset_DEC_arcsec-event.ydata)**2)
            print('Fiber', self.active_fibre,
                  'ID =', self.rss.ID[self.active_fibre],
                  'is the closest to', event.xdata, event.ydata)
        except:
            self.active_fibre = False
        self.update_plots()

    def update_plots(self):
        self.update_RA_DEC_plot()
        self.update_Spectra_plot()

    def update_Spectra_plot(self):
        plot = 'Spectra'
        l_min = self.plot_range_current[plot]['wavelength'][0]
        l_max = self.plot_range_current[plot]['wavelength'][1]
        v_min = self.plot_range_current[plot]['value'][0]
        v_max = self.plot_range_current[plot]['value'][1]
        band_min = self.plot_range_current['RA - DEC']['wavelength'][0]
        band_max = self.plot_range_current['RA - DEC']['wavelength'][1]

        fig = self.figure[plot]
        plt.figure(fig.number)
        plt.clf()
        plot_axis = plt.axes((.05, .1, .94, .85))  # left,bottom,width,height

        plot_axis.plot(self.rss.wavelength, self.normalised_spectrum, 'b:')
        plot_axis.axvspan(band_min, band_max, alpha=0.1, color='c')
#        plot_axis.plot(self.rss.wavelength, self.rss.sky_emission/np.nanmean(self.rss.sky_emission), 'r-')

        if(self.active_fibre is not False):
            active_spectrum = np.array(self.rss.intensity[self.active_fibre])
            active_spectrum *= self.rss.n_wave/np.nansum(active_spectrum)
            plot_axis.plot(self.rss.wavelength, active_spectrum, 'k-')

        plot_axis.set_xlim(l_min, l_max)
        plot_axis.set_ylim(v_min, v_max)
        plot_axis.set_xlabel("$\lambda [\AA]$")
        plot_axis.set_ylabel("NORMALISED intensity")

        self.canvas[plot].draw()

    def update_RA_DEC_plot(self):
        plot = 'RA - DEC'
        RA_min = self.plot_range_current[plot]['offset RA'][0]
        RA_max = self.plot_range_current[plot]['offset RA'][1]
        DEC_min = self.plot_range_current[plot]['offset DEC'][0]
        DEC_max = self.plot_range_current[plot]['offset DEC'][1]
        l_min = self.plot_range_current[plot]['wavelength'][0]
        l_max = self.plot_range_current[plot]['wavelength'][1]
        v_min = self.plot_range_current[plot]['value'][0]
        v_max = self.plot_range_current[plot]['value'][1]

        fig = self.figure[plot]
        plt.figure(fig.number)
        plt.clf()
        plot_axis = plt.axes((.1, .08, .7, .7/.8))  # left, bottom, width, height
        cbar_axis = plt.axes((.85, .08, .05, .7/.8))  # left, bottom, width, height
        plot_axis.cla()

        plot_axis.set_aspect('equal')
        bbox = plot_axis.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        plot_area_pt2 = bbox.width*bbox.height*fig.dpi**2
        plot_area_arcsec2 = (RA_max-RA_min)*(DEC_max-DEC_min)
        arcsec2_to_pt2 = plot_area_pt2/plot_area_arcsec2
        value = self.rss.flux_between(l_min, l_max)[0]

        absolute_value = np.absolute(value)
        v_min = np.nanpercentile(absolute_value, 1)
        v_med = np.nanpercentile(absolute_value, 50)
        v_max = np.nanpercentile(absolute_value, 99)
        exponent = np.log(0.5)/np.log((v_med-v_min)/(v_max-v_min))

        if(self.active_fibre is not False):
            plot_axis.plot(self.rss.offset_RA_arcsec[self.active_fibre],
                           self.rss.offset_DEC_arcsec[self.active_fibre],
                           markerfacecolor='none',
                           markeredgecolor='black',
                           markeredgewidth=2,
                           markersize=np.sqrt(arcsec2_to_pt2),
                           marker='h')
        colour_map = plot_axis.scatter(self.rss.offset_RA_arcsec,
                                       self.rss.offset_DEC_arcsec,
                                       c=value,
                                       norm=matplotlib.colors.PowerNorm(exponent),
                                       vmin=v_min, vmax=v_max,
                                       # cmap=fuego_color_map, norm=norm,
                                       s=.6*arcsec2_to_pt2,
                                       # ad-hoc scaling for KOALA
                                       # (should be fibre_area;
                                       # maybe colorbar?)
                                       marker="h")
        plot_axis.invert_xaxis()

        plot_axis.set_xlim(RA_min, RA_max)
        plot_axis.set_ylim(DEC_min, DEC_max)
        plot_axis.set_title(self.rss.description+" - RSS map")
        plot_axis.set_xlabel("$\Delta$ RA [arcsec]")
        plot_axis.set_ylabel("$\Delta$ DEC [arcsec]")
        fig.colorbar(colour_map, cax=cbar_axis)

        self.canvas[plot].draw()

    def file_clicked(self, filename=False):
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
        self.update_plots()

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
