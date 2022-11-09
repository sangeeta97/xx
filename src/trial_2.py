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
from ms_peak_picker import pick_peaks
from threading import Thread
from collections import defaultdict
from collections import namedtuple
import plotly.express as px
from peaks import *

example_file = "20220808_Mutanobactin.mzML"

run = pymzml.run.Reader(example_file,  MS1_Precision=8e-6)

mass_all = defaultdict(lambda: defaultdict(list))

dict1 = {721.4317: 6.330, 735.4473: 6.905, 721.4858: 6.621}




def peak_dt(mm, name):
    ll = mm.sort_values(by=['drift_time'])
    xnew = ll['drift_time'].values
    ynew = ll['intensity'].values
    mz = ll['mz'].values
    noise = estimate_noise(ynew)
    baseline = estimate_baseline(ynew, noise)
    start, peaks, end = detect_peaks(ynew, noise, baseline)
    peak_intensities  = np.array([ynew[x] for x in peaks])
    max_intensity = max(peak_intensities)
    peaks_index = peak_intensities > max_intensity/10
    peaks = peaks[peaks_index]
    start = start[peaks_index]
    end = end[peaks_index]
    peak_mid = xnew[peaks]
    mz_mid = mz[peaks]
    peak_start = xnew[start]
    peak_end = xnew[end]
    dict_all = {"compound": [name]*mz_mid.size, "mz_top": mz_mid, "dt_start": peak_start, "dt_mid": peak_mid, "dt_end": peak_end}
    df = pd.DataFrame(data = dict_all)
    df.to_csv(f"peaks{name}.csv")








#adduct mass for all
find_mass = {287.0555: [721.4317, 735.4473, 721.4858]}



def extract_decorator(spectrum, target):
    mz = spectrum.mz[spectrum.i > 20.0]
    intensity = spectrum.i[spectrum.i > 20.0]
    whole = namedtuple('Whole', ['mz', 'intensity', 'dt', 'rt'])
    gg = round(spectrum.get("MS:1002476"), 2)
    if mz.size > 0:
        li = [mz, intensity, gg, round(spectrum.scan_time_in_minutes(), 3)]
        mm = whole._make(li)
        return mm

def all_rt(mm, name):
    ll = []
    pp= mm[['rt', 'intensity']].groupby('rt')['intensity'].apply(max)
    rt= np.array(pp.index.values)
    intensity= pp.values
    f = interpolate.interp1d(rt, intensity)
    xnew = np.arange(rt.min(), rt.max(), 0.005)
    ynew = f(xnew)
    fig = px.line(x=xnew, y=ynew, labels={'x':'rt', 'y':'intensity'})
    fig.write_html(f"{name}_rt.html")
    noise = estimate_noise(intensity)
    baseline = estimate_baseline(intensity, noise)
    start, peaks, end = detect_peaks(intensity, noise, baseline)
    peak_intensities  = np.array([intensity[x] for x in peaks])
    max_intensity = max(peak_intensities)
    peaks_index = peak_intensities > max_intensity/10
    peaks = peaks[peaks_index]
    start = start[peaks_index]
    end = end[peaks_index]
    return rt[peaks]

    # peak_list = pick_peaks(xnew, ynew, fit_type="quadratic")
    # if len(peak_list) > 0:
    #     for pp in peak_list:
    #         if pp.full_width_at_half_max < 0.5:
    #             less = pp.mz - 0.5
    #             more = pp.mz + 0.5
    #             df0 = mm.loc[(mm['rt'] >= less) & (mm['rt'] <= more), ]
    #             ll.append(df0)
    # return ll



def ticplot(df1, name):
    df1 = df1.sort_values(by=['drift_time'])
    x= df1['drift_time'].values
    y= df1['intensity'].values
    f = interpolate.interp1d(x, y, kind = "cubic")
    xnew = np.arange(df1['drift_time'].min(), df1['drift_time'].max(), 0.05)
    ynew = f(xnew)
    fig = px.line(x=xnew, y=ynew, labels={'x':'drift_time', 'y':'intensity'})
    fig.write_html(f"{name}.html")


def extract_mzml(spectrum):

    df1 = pd.DataFrame({"mz": spectrum.mz, "intensity": spectrum.i, "drift_time": np.array(spectrum.get("MS:1002476") * spectrum.mz.size)})
    return df1


from functools import wraps
import time


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper

@timeit
def final():
    spec_one = (spec for spec in run if spec.ms_level == 1)
    tt= []
    gg = defaultdict(list)
    target = np.array(find_mass[287.0555])
    for i, spec in enumerate(spec_one):
        ff= extract_decorator(spec, target)
        if bool(ff):
            for p in target:
                condition = np.absolute(ff.mz- p) < 0.003
                result_mz = np.extract(condition, ff.mz)
                result_intensity = np.extract(condition, ff.intensity)
                if result_mz.any():
                    ll = result_mz.size
                    mass_all[p]['spec_index'].append(i)
                    mass_all[p]['mz'].append(result_mz[0])
                    mass_all[p]['intensity'].append(result_intensity[0])
                    mass_all[p]['drift_time'].append(ff.dt)
                    mass_all[p]['rt'].append(ff.rt)



    for z in mass_all.keys():
        print(z)

        zz = mass_all[z]
        mm = pd.DataFrame(data = zz)
        mm.columns = ['index', 'mz', 'intensity', 'drift_time', 'rt']
        dd = []
        pp = all_rt(mm, z)
        for y in pp:
            print(y)
            rt_range_min = y - 0.2
            print(rt_range_min)
            rt_range_max = y + 0.2
            print(rt_range_max)
            ll= mm.loc[(mm['rt'].values >= rt_range_min) & (mm['rt'].values <= rt_range_max), ]
            print(ll)
            dd.append(ll)

        # dd = dd.append(df)
        ll = pd.concat(dd)
        lt = ll[['drift_time', 'intensity', 'mz']].sort_values(by = ['intensity']).drop_duplicates(subset = ['drift_time'], keep = 'last')
        print("yes")
        print(lt['mz'].min())
        print(lt['mz'].max())
        lt.to_csv(f'final{z}.csv')

        ticplot(lt, z)
        peak_dt(lt, z)
        # ll = lt.sort_values(by=['drift_time'])
        # xnew = ll['drift_time'].values
        # ynew = ll['intensity'].values
        # noise = estimate_noise(ynew)
        # baseline = estimate_baseline(ynew, noise)
        # start, peaks, end = detect_peaks(ynew, noise, baseline)
        # peak_intensities  = np.array([ynew[x] for x in peaks])
        # max_intensity = max(peak_intensities)
        # peaks_index = peak_intensities > max_intensity/10
        # peaks = peaks[peaks_index]
        # start = start[peaks_index]
        # end = end[peaks_index]
        # print(xnew[peaks])
        # print(xnew[start])
        # print(xnew[end])


        # rt_range = dict1[z]
        # rt_range_min = rt_range - 0.500
        # rt_range_max = rt_range + 0.500
        # ll= mm.loc[(mm['rt'] >= rt_range_min) & (mm['rt'] <= rt_range_max), ]
        # lt = ll[['drift_time', 'intensity']].sort_values(by = ['intensity']).drop_duplicates(subset = ['drift_time'], keep = 'last')
        # lt.to_csv(f'final{z}.csv')
        #
        # ticplot(lt, z)
        # ll = lt.sort_values(by=['drift_time'])
        # xnew = ll['drift_time'].values
        # ynew = ll['intensity'].values
        # noise = estimate_noise(ynew)
        # baseline = estimate_baseline(ynew, noise)
        # start, peaks, end = detect_peaks(ynew, noise, baseline)
        # peak_intensities  = np.array([ynew[x] for x in peaks])
        # max_intensity = max(peak_intensities)
        # peaks_index = peak_intensities > max_intensity/5
        # peaks = peaks[peaks_index]
        # start = start[peaks_index]
        # end = end[peaks_index]
        #
        #
        # print(xnew[peaks])
        # print(xnew[start])
        # print(xnew[end])



        # pp= mm[['rt', 'intensity']].groupby('rt')['intensity'].apply(max)
        # rt= np.array(pp.index)
        # intensity= pp.values
        # f = interpolate.interp1d(rt, intensity)
        # xnew = np.arange(rt.min(), rt.max(), 0.1)
        # ynew = f(xnew)
        # peak_list = pick_peaks(xnew, ynew, fit_type="quadratic")
        # if len(peak_list) > 0:
        #     for pp in peak_list:
        #         if pp.full_width_at_half_max < 0.5:
        #             gg[z].append(round(pp.mz, 1))






if __name__ == "__main__":
    xx = final()
