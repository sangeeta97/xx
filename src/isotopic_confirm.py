import os
import pymzml
import json
import pandas as pd
import statistics
from scipy import interpolate
import numpy as np
from numpy.linalg import norm
import copy
from collections import defaultdict
import operator
from collections import defaultdict as ddict
import math
from threading import Thread
from collections import defaultdict
from collections import namedtuple
import plotly.express as px
import tempfile
from functools import wraps
import time
from .formula import isotopic_table
import pandas as pd
from bokeh.io import output_notebook
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, NumeralTickFormatter, Label
from bokeh.palettes import Category10
from molmass import Formula
from bokeh.plotting import figure, output_file, show
from bokeh.resources import CDN
from bokeh.embed import file_html



def add_axis_labels(p):
    p.xaxis.axis_label = 'Fragment m/z'
    p.xaxis.axis_label_text_font_size = '10pt'
    p.xaxis.major_label_text_font_size = '9pt'

    p.yaxis.axis_label = 'Intensity'
    p.yaxis.axis_label_text_font_size = '10pt'
    p.yaxis.major_label_text_font_size = '9pt'
    p.yaxis.formatter = NumeralTickFormatter(format='0.')


def create_p(width=800, height=300,
            main_title='IM EIC plot'):
    tooltips = [
        ('m/z','@mz{0.0000}'),
        ('Int','@intensity')
        ]
    p = figure(
        plot_width=width, plot_height=height,
        title = main_title,
        tools = 'xwheel_zoom,xpan,box_zoom,undo,reset',
        tooltips=tooltips
        )
    return p

'''
ion_dict : dictionary of ion and the corresponding monoisotopic mass value for the ion
selected_ions : dictionary of ion and the corresponding MF to isotopic confirmation
found_molecular_formula : if the molecular formula is found based on the accurate mass of primary ions, keys are the molecular formula and values are the retention time for those MF.
dt_present : dataframe of combination of rt and ion_type with selected peak in drift time space.
rt_found_ion: key is rt and value is list of ion types found at that rt
ions_data[ion][i]: first key is ion_type and second key is rt and value is a raw dataframe.
dt_peaks[ion][i]: first key is ion_type and second key is rt and value is picked peak.
'''


class Xic_Eic:

    def __init__(self, data):
        self.data = data



    def theoretical_function(self, z):
        mf = self.data.selected_ions[z]
        df = isotopic_table("spectrum", mf)
        xx = df['relative_mass'].values
        xx = xx.astype(float)
        yy = df['intensity'].values
        yy = yy.astype(float)
        return xx, yy




    def plotdt_IC(self):
        for b in self.data.rt_found_ion:
            print(b)
            zz = self.data.rt_found_ion[b]
            df_all = []
            for z in zz:
                df = self.data.ions_data[z][b]
                x= df['drift_time'].values
                y= df['intensity'].values
                f = interpolate.interp1d(x, y, kind = "cubic")
                xnew = np.arange(df['drift_time'].min(), df['drift_time'].max(), 0.05)
                ynew = f(xnew)
                data = {"drift_time": xnew, "intensity": ynew}
                df1 = pd.DataFrame(data = data)
                df1['ion_type'] = z
                df_all.append(df1)
            df = pd.concat(df_all)
            b = str(b)
            fig = px.line(df, x="drift_time", y="intensity", color='ion_type', title=f'overlay_IM plot of {self.data.molecular_formula}_{b}')
            # fig.write_html(f"{self.data.molecular_formula}_{b}_overlay.html")
            fig.write_html(os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{b}_IMoverlay.html"))




    def isotope_distribution_overlay(self, ion, rt):
        df = self.data.dt_peaks[ion][rt]
        run = pymzml.run.Reader(self.data.mzml,  MS1_Precision=8e-6)
        spec_list = np.array([spec for spec in run])
        spec_index = df["spec_number"].values
        spec_index = spec_index[0]
        target =  df['theoretical_mass'].values[0]
        spec = spec_list[spec_index]
        ff= self.data.extract_decorator(spec)
        test_mz = np.array([target, target + 1.00335, target + 1.00335 + 1.00335])
        mz_index = np.searchsorted(ff.mz, test_mz)
        result_mz = ff.mz[mz_index]
        result_intensity = ff.intensity[mz_index]
        result_intensity = result_intensity.astype(float)
        combine = result_intensity.max()
        result_intensity = np.divide(result_intensity, combine)
        result_intensity = np.multiply(result_intensity, 100)
        result_mz = result_mz.astype(float)
        df_exp = pd.DataFrame({"mz": result_mz, "intensity": result_intensity})
        _, theoretical_intensity = self.theoretical_function(ion)
        theoretical_intensity = theoretical_intensity[0:3]
        df_thr = pd.DataFrame({"mz": test_mz, "intensity": theoretical_intensity})
        sources = [ColumnDataSource(df_exp), ColumnDataSource(df_thr)]
        p = create_p()
        p.vbar(x = "mz", width = 0.1, bottom = 0, top = "intensity", color="navy", fill_alpha = 0, source = sources[-1])
        p.vbar(x = "mz", width = 0.001, bottom = 0, top = "intensity", color = '#324ea8', fill_alpha = 1, source = sources[0])
        add_axis_labels(p)
        html = file_html(p, CDN, "my plot")
        f = open(os.path.join(self.data.temp_spectrum, f"{self.data.molecular_formula}_{ion}_{rt}.html"), "w")
        f.write(html)
        f.close()

        # export_png(p, filename =




    def ticplot(self, ion, rt):
        df = self.data.ions_data[ion][rt]
        df = df.sort_values(by=['drift_time'])
        x= df['drift_time'].values
        y= df['intensity'].values
        f = interpolate.interp1d(x, y, kind = "cubic")
        xnew = np.arange(df['drift_time'].min(), df['drift_time'].max(), 0.05)
        ynew = f(xnew)
        fig = px.line(x=xnew, y=ynew, labels={'x':'drift_time', 'y':'intensity'}, title=f'IM plot of {self.data.molecular_formula}_{rt}_{ion}')
        # fig.write_html(f"{self.data.molecular_formula}_{rt}_{ion}.html")
        fig.write_html(os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{ion}_{rt}_IM.html"))
