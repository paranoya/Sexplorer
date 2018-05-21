#!/usr/bin/python
# -*- coding: utf-8 -*-
# First attempt to show a KOALA RSS file, based on the
# script to estimate AAT focus by Angel R. Lopez-Sanchez (AAO/MQU)
# Yago Ascasibar (AAT, 19/05/2018)

import numpy as np
from astropy.io import fits
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
# Tkinter is for python 2; tkinter is for python 3
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFont
#    import tkMessageBox
    import tkFileDialog
    import tkSimpleDialog
else:
    import tkinter as tk
    from tkinter import messagebox as tkMessageBox
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

        # --- define the file access widgets

        self.open_file_button = tk.Button(self, text='File:',
                                          command=self.open_file_clicked)
        self.open_file_button.grid(row=0, column=0, sticky="sw")
        self.filename = tk.Label(self, text="None")
        self.filename.grid(row=0, column=1, columnspan=3, sticky="sw")

        tk.Label(self, text="Object:").grid(row=1, column=1, sticky="ne")
        self.object = tk.Label(self, text="None")
        self.object.grid(row=1, column=2, sticky="nw")
        self.t_exp_label = tk.Label(self, text="t_exp = 0 s")
        self.t_exp_label.grid(row=1, column=3, sticky="nw")

        # --- define the plotting area

        self.figure = plt.figure(figsize=(18, 18))
        self.colorbar = False
        self.figure_canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.figure_canvas.get_tk_widget().grid(column=4, row=0,
                                                rowspan=10, columnspan=20,
                                                sticky="nesw")

        # --- define the widgets to control what's plotted

        tk.Label(self, text="Image range:").grid(row=12, column=3, sticky="w")

        self.xmin_button = tk.Button(self, text='x min', command=self.set_xmin)
        self.xmin = tk.Label(self, text='0.')
        self.xmin_file = tk.Label(self, text='0.')
        self.xmin_button.grid(row=11, column=4, sticky="e")
        self.xmin.grid(row=11, column=5, sticky="w")
        self.xmin_file.grid(row=12, column=5, sticky="w")

        self.xmax_button = tk.Button(self, text='x max', command=self.set_xmax)
        self.xmax = tk.Label(self, text='1.')
        self.xmax_file = tk.Label(self, text='1.')
        self.xmax.grid(row=11, column=6, sticky="e")
        self.xmax_button.grid(row=11, column=7, sticky="w")
        self.xmax_file.grid(row=12, column=6, sticky="e")

        self.ymin_button = tk.Button(self, text='y min', command=self.set_ymin)
        self.ymin = tk.Label(self, text='0.')
        self.ymin_file = tk.Label(self, text='0.')
        self.ymin_button.grid(row=11, column=8, sticky="e")
        self.ymin.grid(row=11, column=9, sticky="w")
        self.ymin_file.grid(row=12, column=9, sticky="w")

        self.ymax_button = tk.Button(self, text='y max', command=self.set_ymax)
        self.ymax = tk.Label(self, text='1.')
        self.ymax_file = tk.Label(self, text='1.')
        self.ymax.grid(row=11, column=10, sticky="e")
        self.ymax_button.grid(row=11, column=11, sticky="w")
        self.ymax_file.grid(row=12, column=10, sticky="e")

        self.vmin_button = tk.Button(self, text='v min', command=self.set_vmin)
        self.vmin = tk.Label(self, text='0.')
        self.vmin_file = tk.Label(self, text='0.')
        self.vmin_button.grid(row=11, column=12, sticky="e")
        self.vmin.grid(row=11, column=13, sticky="w")
        self.vmin_file.grid(row=12, column=13, sticky="w")

        self.vmax_button = tk.Button(self, text='v max', command=self.set_vmax)
        self.vmax = tk.Label(self, text='1.')
        self.vmax_file = tk.Label(self, text='1.')
        self.vmax.grid(row=11, column=14, sticky="e")
        self.vmax_button.grid(row=11, column=15, sticky="w")
        self.vmax_file.grid(row=12, column=14, sticky="e")

        self.plot_scale = tk.StringVar(self)
        self.plot_scale.set(color_normalisation.keys()[0])
        tk.OptionMenu(self, self.plot_scale, 'linear', 'log',
                      command=self.update_plot).grid(row=11, column=16)

    def update_plot(self, scale=None): # Need to specify that because of OptionMenu
        plt.clf()
        plt.title(self.object['text'])
        plt.xlim(float(self.xmin['text']), float(self.xmax['text']))
        plt.ylim(float(self.ymin['text']), float(self.ymax['text']))
        plt.imshow(self.data, origin='lower', interpolation='none',
                   vmin=float(self.vmin['text']),
                   vmax=float(self.vmax['text']),
                   norm=color_normalisation[self.plot_scale.get()](),
                   )
        plt.colorbar()
        self.figure_canvas.draw()

    def open_file_clicked(self):
        filename = \
            tkFileDialog.askopenfilename(initialdir=".",
                                         title="Select file",
                                         filetypes=(("FITS files", "*.fits"),
                                                    ("all files", "*.*")))
#        print('Opening file "'+filename+'"')
        print(filename)

        self.hdu = fits.open(filename)
        self.filename.configure(text=filename)
        hdr = self.hdu[0].header
#        self.hdu.info()
        self.object.configure(text=hdr['object'])
        print(hdr['object'])
        
        try:
            self.t_exp = hdr['exptime']
        except:
            try:
                self.t_exp = hdr['totalexp']
            except:
                self.t_exp = 'unknown'
        self.t_exp_label.configure(text='t_exp = '+str(self.t_exp)+' s')

        self.data = self.hdu[0].data

        self.xmin.configure(text='0')
        self.xmin_file.configure(text='0')
        self.xmax.configure(text=hdr['naxis1'])
        self.xmax_file.configure(text=hdr['naxis1'])
        self.ymin.configure(text='0')
        self.ymin_file.configure(text='0')
        self.ymax.configure(text=hdr['naxis2'])
        self.ymax_file.configure(text=hdr['naxis2'])
        self.vmin.configure(text=np.nanmin(self.data))
        self.vmin_file.configure(text=np.nanmin(self.data))
        self.vmax.configure(text=np.nanmax(self.data))
        self.vmax_file.configure(text=np.nanmax(self.data))

        self.update_plot()

    def set_xmin(self):
        xmin = tkSimpleDialog.askfloat('Set xmin',
                                       'Enter minimum value for x axis',
                                       initialvalue=self.xmin['text'])
        self.xmin.configure(text=xmin)
        self.update_plot()

    def set_xmax(self):
        xmax = tkSimpleDialog.askfloat('Set xmax',
                                       'Enter maximum value for x axis',
                                       initialvalue=self.xmax['text'])
        self.xmax.configure(text=xmax)
        self.update_plot()

    def set_ymin(self):
        ymin = tkSimpleDialog.askfloat('Set ymin',
                                       'Enter minimum value for y axis',
                                       initialvalue=self.ymin['text'])
        self.ymin.configure(text=ymin)
        self.update_plot()

    def set_ymax(self):
        ymax = tkSimpleDialog.askfloat('Set ymax',
                                       'Enter maximum value for y axis',
                                       initialvalue=self.ymax['text'])
        self.ymax.configure(text=ymax)
        self.update_plot()

    def set_vmin(self):
        vmin = tkSimpleDialog.askfloat('Set vmin',
                                       'Enter minimum value for color map',
                                       initialvalue=self.vmin['text'])
        self.vmin.configure(text=vmin)
        self.update_plot()

    def set_vmax(self):
        vmax = tkSimpleDialog.askfloat('Set vmax',
                                       'Enter maximum value for color map',
                                       initialvalue=self.vmax['text'])
        self.vmax.configure(text=vmax)
        self.update_plot()

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
