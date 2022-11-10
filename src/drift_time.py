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
from .peak import *
import tempfile
from functools import wraps
import time
import matplotlib.pyplot as plt
from scipy.signal import peak_widths
from .isotopic_confirm import Xic_Eic
import plotly.figure_factory as ff
import numpy as np





#https://towardsdatascience.com/interactive-mass-spectra-with-bokeh-3b9163881b12

#https://spectrum-utils.readthedocs.io/en/latest/plotting.html




full_result = defaultdict(lambda: defaultdict(lambda: None))

#

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




def find_fwhm(half_max, start, end, xnew, ynew):
    fwhm_all = []
    for i, j, k in zip(start, end, half_max):
        peak_range = np.arange(i, j, 1)
        peak_range_y = ynew[peak_range]
        peak_range_x = xnew[peak_range]
        peak_range_y = np.ravel(peak_range_y)
        bool_array = peak_range_y > k
        print(bool_array)
        peak_range_x = peak_range_x[bool_array]
        print(peak_range_x)
        first = peak_range_x.min()
        last = peak_range_x.max()
        fwhm = last - first
        print(fwhm)
        fwhm_all.append(fwhm)
    fwhm_all = np.array(fwhm_all)
    fwhm_all = np.ravel(fwhm_all)
    return fwhm_all





class Parse_MF:
    '''
    molecular formula parsing and giving table of drift time and rt information

    '''
    def __init__(self, data):
        self.data = data
        self.monoPPM = self.data.ppm_values
        self.mass_df = self.data.ions
        self.mzml = self.data.primary_data["mzml"]
        self.primary_ion = self.data.primary_data['primary_ion']
        self.isotope = self.data.isotope_ratio
        self.temp_spectrum = self.data.temp_spectrum
        self.temp_drift = self.data.temp_drift
        self.run_all()




    def run_all(self):
        for x, y in self.mass_df.items():
            drift_cal = Drift_Time(self.mzml, x, y.values, list(y.index), self.monoPPM, self.isotope, self.primary_ion, self.temp_spectrum, self.temp_drift, self.data.message, self.data.isotope_match, self.data.spec_list)
            drift_cal.launch()



    def dataframe_all(self):
        kk = []
        k = []
        for x, y in full_result.items():
            dd = []
            td = []
            for c, t in y.items():
                dd.append(t)
                td.append(c)
                df1 = pd.concat(dd, keys = td)
            df1 = df1.reset_index()
            df1 = df1.drop(['level_1'], axis = 1)
            df1 = df1.rename({"level_0": "ion_type"}, axis = 1)
            kk.append(df1)
            k.append(x)
        df_final = pd.concat(kk, keys = k)
        df_final = df_final.reset_index()
        df_final = df_final.drop(['level_1'], axis = 1)
        df_final = df_final.rename({"level_0": "molecular_formula"}, axis = 1)
        return df_final



'''
This class parse the mzml file, extract the exact rt value for each molecular formula, if compound is not found it flags that as well.
it confirms the presence of molecular formula using its primary ion by spectral isotope pattern matching.
primary ion gives the rt window in which all the ions search get limited.
Also, it generates ion dataframe for each molecular formula if ion is being found, again confirmed by spectral isotope matching.

The spectral confirmation give spectral ids which give idea of possible dt values (putting dt values accuracy upto two decimal points)
which later being confirmed by peak picking on the dt/intensity axis.
'''





class Drift_Time:
    '''
    parse the mzml to produce ion dataframe for the ion+molecular formula studied.The ion dataframe for each ion species corresponding to a molecular formula have similar rt like primary ions
    and its own dt, intensity and mz values. mzml is path to the file, mass_values is exact mass of the ions, ions is the string of ions
    types, ppm is mass accuracy defined by the user, primary_ion is the primary ion type for the MF, and temp is path to the TemporaryDirectory.

    '''

    def __init__(self, mzml, molecular_formula, mass_values, ions, PPM, isotope, primary_ion, temp_spectrum, temp_drift, message, isotope_match, spec_list):
        self.run = pymzml.run.Reader(mzml,  MS1_Precision=8e-6)
        self.molecular_formula = molecular_formula
        self.mass_values = mass_values
        self.ions = ions
        self.PPM1 = PPM[2]
        self.PPM2 = PPM[1]
        self.isotope = isotope
        self.mzml = mzml
        self.primary_ion = primary_ion
        self.temp_spectrum = temp_spectrum
        self.temp_drift = temp_drift
        self.ion_dict = {x:y for x, y in zip(ions, mass_values)}
        self.selected_ions = defaultdict(lambda: None)
        self.rt_found_ion = defaultdict(list)
        self.ions_data = defaultdict(lambda: defaultdict(lambda: None))
        self.found_molecular_formula = defaultdict(lambda: None)
        self.dt_peaks = defaultdict(lambda: defaultdict(lambda: None))
        self.mass_all = defaultdict(lambda: defaultdict(list))
        self.ion_present = list()
        self.status = False
        self.peaks = False
        self.message = message
        self.isotope_match = isotope_match
        self.spec_list = spec_list
        self.formula_define()
        '''
        Attributes:
        MF: the studies molecular formula
        mass_values : The array of monoisotopic mass of ions of the molecular formula to be studied.
        ions : The array of string ions to be studied
        ppm : mass accuracy selected for monoisotopic mass match
        isotope : table of abundance ratio for different isotopes
        mzml : path to the mzml file
        primary_ion : primary ion for the MF used to confirm the rt space for each MF
        temp_spectrum : tempory folder to store spectrum except drift time related spectra
        temp_drift : temporary folder to store drift time related spectra
        ion_dict : dictionary of ion and the corresponding monoisotopic mass value for the ion
        selected_ions : dictionary of ion and the corresponding MF to isotopic confirmation
        found_molecular_formula : if the molecular formula is found based on the accurate mass of primary ions, keys are the molecular formula and values are the retention time for those MF.
        dt_present : dataframe of combination of rt and ion_type with selected peak in drift time space.
        rt_found_ion: key is rt and value is list of ion types found at that rt
        ions_data[ion][i]: first key is ion_type and second key is rt and value is a raw dataframe.
        dt_peaks[ion][i]: first key is ion_type and second key is rt and value is picked peak.

        '''


    def formula_define(self):
        for y in self.ions:
            ion_not = ["[M+H]+", "[M-H]-", "[M-]", "[M+]", "+electron", "-electron"]
            x = '' if y in ion_not else y
            x = x.strip("+")
            x = x.strip("-")
            ions = self.molecular_formula + x.strip()
            self.selected_ions[y] = ions


    def launch(self):
        self.ion_all_extract()
        if len(self.mass_all[self.primary_ion]['mz']) < 10:
            self.message['warning'].append("The molecular formula is not present in the sample")
            return None
        self.rt_window_formula()
        ss = self.found_molecular_formula[self.molecular_formula]
        if ss.size > 0:
            self.ion_dataframe()
            self.plotrt_EIC()
            plot_object = Xic_Eic(self)
            plot_object.plotdt_IC()
            for z in self.rt_found_ion:
                b = self.rt_found_ion[z]
                for g in b:
                    plot_object.isotope_distribution_overlay(g, z)
                    plot_object.ticplot(g, z)
        else:
            self.message['warning'].append("The molecular formula is not present in the sample")
            return None




    def extract_decorator(self, spectrum):
        mz = spectrum.mz[spectrum.i > 20.0]
        intensity = spectrum.i[spectrum.i > 20.0]
        whole = namedtuple('Whole', ['mz', 'intensity', 'dt', 'rt'])
        gg = round(spectrum.get("MS:1002476"), 2)
        if mz.size > 0:
            li = [mz, intensity, gg, round(spectrum.scan_time_in_minutes(), 3)]
            mm = whole._make(li)
            return mm

    @timeit
    def ion_all_extract(self):
        spec_one = [spec for spec in self.run if spec.ms_level == 1]
        self.spec_list['spec_list'] = spec_one
        tt= []
        gg = defaultdict(list)
        target = self.mass_values
        for i, spec in enumerate(spec_one):
            ff= self.extract_decorator(spec)
            if bool(ff):
                for t, p in enumerate(target):
                    condition = np.absolute(ff.mz- p) < self.PPM1/1000000 * p
                    result_mz = np.extract(condition, ff.mz)
                    result_intensity = np.extract(condition, ff.intensity)
                    k = self.ions[t]
                    if result_mz.any():
                        ll = result_mz.size
                        self.mass_all[k]['spec_index'].append(i)
                        self.mass_all[k]['mz'].append(result_mz[0])
                        self.mass_all[k]['intensity'].append(result_intensity[0])
                        self.mass_all[k]['drift_time'].append(ff.dt)
                        self.mass_all[k]['rt'].append(ff.rt)
                        self.mass_all[k]['theoretical_mass'].append(p)



    @timeit
    def detect_rt(self):
        zz = self.mass_all[self.primary_ion]
        mm = pd.DataFrame(data = zz)
        mm.columns = ['index', 'mz', 'intensity', 'drift_time', 'rt', 'theoretical_mass']
        confirm_p = mm['index'].values
        pp= mm[['rt', 'intensity']].groupby('rt')['intensity'].apply(max)
        rt= np.array(pp.index.values)
        intensity= pp.values
        noise = estimate_noise(intensity)
        baseline = estimate_baseline(intensity, noise)
        start, peaks, end = detect_peaks(intensity, noise, baseline)
        peak_intensities  = np.array([intensity[x] for x in peaks])
        max_intensity = max(peak_intensities)
        peaks_index = peak_intensities > max_intensity/3
        peaks = peaks[peaks_index]
        start = start[peaks_index]
        end = end[peaks_index]
        rt_peaks = rt[peaks]
        spec_index = mm[mm['rt'].isin(rt_peaks)]
        print(spec_index)
        spec_index = spec_index.sort_values(by = ['intensity'], ascending=False)
        spec_index = spec_index.drop_duplicates(subset=['rt'])
        spec_index = spec_index['index'].values
        return rt[peaks], spec_index

    @timeit
    def plotrt_EIC(self):
        ions = []
        for x in self.rt_found_ion:
            yy = self.rt_found_ion[x]
            ions.extend(yy)
        ions = set(ions)
        for x in self.rt_found_ion:
            df_all = []
            for ion in ions:
                zz = self.mass_all[ion]
                mm = pd.DataFrame(data = zz)
                mm.columns = ['index', 'mz', 'intensity', 'drift_time', 'rt', 'theoretical_mass']
                df = mm.sort_values(by=['rt'])
                g= df['rt'].values
                y= df['intensity'].values
                f = interpolate.interp1d(g, y, kind = "linear")
                xnew = np.arange(df['rt'].min(), df['rt'].max(), 0.001)
                ynew = f(xnew)
                data = {"rt": xnew, "intensity": ynew}
                df1 = pd.DataFrame(data = data)
                df1['ion_type'] = ion
                df_all.append(df1)
            df = pd.concat(df_all)
            fig = px.line(df, x="rt", y="intensity", color='ion_type', title=f'rt EIC plot of {self.molecular_formula}')
            # fig.write_html(f"{self.molecular_formula}_rt_overlay.html")
            fig.write_html(os.path.join(self.temp_spectrum, f"{self.molecular_formula}_rt_overlay.html"))




    @timeit
    def plot3d_EIC(self):
        ions = []
        for x in self.rt_found_ion:
            yy = self.rt_found_ion[x]
            ions.extend(yy)
        ions = set(ions)
        for x in self.rt_found_ion:
            for ion in ions:
                zz = self.mass_all[ion]
                mm = pd.DataFrame(data = zz)
                mm.columns = ['index', 'mz', 'intensity', 'drift_time', 'rt', 'theoretical_mass']
                df = mm.sort_values(by=['rt'])
                x= df['rt'].values
                z= df['intensity'].values
                y = df['drift_time'].values
                fig = px.line(df, x="rt", y="intensity", color='ion_type', title=f'rt EIC plot of {self.molecular_formula}')
                # fig.write_html(f"{self.molecular_formula}_rt_overlay.html")
                fig.write_html(os.path.join(self.temp_spectrum, f"{self.molecular_formula}_rt_overlay.html"))




    @timeit
    def rt_window_formula(self):
        pp, spec_index = self.detect_rt()
        final_index, rr = self.spectrum_confirm(spec_index, self.primary_ion)
        if np.any(final_index):
            self.status = True
            pp = pp[final_index]
            pp = np.array(pp)
            self.found_molecular_formula[self.molecular_formula] = pp
        else:
            self.message.append(f'{self.molecular_formula} not found')
            self.found_molecular_formula[self.molecular_formula] = np.array([])
            return np.array([])


    def spectrum_confirm(self, spec_index, ion_type):
        spec_list = np.array(self.spec_list['spec_list'])
        n = spec_index.size
        spec_peaks = spec_list[spec_index]
        spec_peaks = np.array([spec_peaks])
        spec_peaks = np.ravel(spec_peaks)
        target1 = self.ion_dict[ion_type]
        target2 = target1 + 1.00335
        mz1 = []
        mz2 = []
        intensity_ratio = []
        for spec in spec_peaks:
            ff= self.extract_decorator(spec)
            condition = np.absolute(ff.mz- target1) < self.PPM1/1000000 * target1
            result_mz1 = np.extract(condition, ff.mz)
            result_intensity1 = np.extract(condition, ff.intensity)
            condition = np.absolute(ff.mz- target2) < self.PPM2/1000000 * target2
            result_mz2 = np.extract(condition, ff.mz)
            result_intensity2 = np.extract(condition, ff.intensity)
            intensity_r = result_intensity2/result_intensity1 * 100
            mz1.append(result_mz1)
            mz2.append(result_mz2)
            intensity_ratio.append(intensity_r)
        test_case = self.isotope[self.molecular_formula]
        test_case = test_case['intensity'].values[1]
        test_case = float(test_case)
        intensity_ratio = np.array(intensity_ratio)
        kd = np.absolute(intensity_ratio - test_case)
        kd = np.ravel(kd)
        dk = kd < test_case/5
        dk = np.ravel(dk)
        ratio_bool = dk
        intensity_ratio = np.ravel(intensity_ratio)
        print(intensity_ratio)
        print(ratio_bool)
        final_ratio = intensity_ratio[ratio_bool]
        final_ratio = np.ravel(final_ratio)
        spec_peaks = np.array(spec_peaks)
        if np.any(final_ratio):
            return ratio_bool, spec_peaks[ratio_bool]
        else:
            return np.array([]), np.array([])



    def spectrum_confirm_ion(self, spec_index, ion_type, rt):
        spec_list = np.array(self.spec_list['spec_list'])
        n = spec_index.size
        spec_peaks = spec_list[spec_index]
        spec_peak = np.array([spec_peaks])
        spec = np.ravel(spec_peaks)
        target1 = self.ion_dict[ion_type]
        target2 = target1 + 1.00335
        ff= self.extract_decorator(spec)
        condition = np.absolute(ff.mz- target1) < self.PPM1/1000000 * target1
        result_mz1 = np.extract(condition, ff.mz)
        result_intensity1 = np.extract(condition, ff.intensity)
        condition = np.absolute(ff.mz- target2) < self.PPM2/1000000 * target2
        result_mz2 = np.extract(condition, ff.mz)
        result_intensity2 = np.extract(condition, ff.intensity)
        intensity_r = result_intensity2/result_intensity1 * 100
        rt = ff.rt if not rt else rt
        intensity_ratio = np.array(intensity_r)
        test_case = self.isotope[self.molecular_formula]
        test_case = test_case['intensity'].values[1]
        test_case = float(test_case)
        kd = np.absolute(intensity_ratio - test_case)
        kd = np.ravel(kd)
        dk = kd < test_case/2
        dk = np.ravel(dk)
        ratio_bool = dk
        intensity_ratio = np.ravel(intensity_ratio)
        print(intensity_ratio)
        print(ratio_bool)
        final_ratio = np.extract(ratio_bool, intensity_ratio)
        if np.any(final_ratio):
            label = f"{self.data.molecular_formula}_{ion_type}_{rt}"
            self.isotope_match[label]['mz'] = np.ravel(np.array(result_mz1[0], result_mz2[0]))
            self.isotope_match[label]['intensity'] = np.ravel(np.array(result_intensity1[0], result_intensity2[0]))
            return ratio_bool, spec[ratio_bool]
        else:
            return np.array([]), np.array([])



    def ion_dataframe(self):

        for ion in self.mass_all.keys():
            zz = self.mass_all[ion]
            mm = pd.DataFrame(data = zz)
            mm.columns = ['index', 'mz', 'intensity', 'drift_time', 'rt', 'theoretical_mass']
            pp = self.found_molecular_formula[self.molecular_formula]
            self.rt = pp
            rt_range_min = pp - 0.1
            rt_range_max = pp + 0.1
            lc = []
            for i in range(len(pp)):
                ll= mm.loc[(mm['rt'].values >= rt_range_min[i]) & (mm['rt'].values <= rt_range_max[i]), ]
                if ll.index.size > 5:
                    print(ion)
                    lt = ll[['index' ,'drift_time', 'intensity', 'mz', 'rt', 'theoretical_mass']].sort_values(by = ['intensity'], ascending = False).drop_duplicates(subset = ['drift_time'], keep = "first")
                    spec_index = lt['index'].values[0]
                    num = 5 if ion == self.primary_ion else 2
                    ratio_bool, _ = self.spectrum_confirm_ion(spec_index, ion, rt = i)
                    if ratio_bool.size > 0:
                        rt = pp[i]
                        self.ions_data[ion][rt] = lt
                        self.rt_found_ion[rt].append(ion)
                        self.ion_present.append(ion)
                        df = self.peak_dt(lt, ion, rt)
                        if self.peaks:
                            lc.append(df)
                        if len(lc) > 0:
                            dd = pd.DataFrame()
                            lc.append(dd)
                        dt_peaks = pd.concat(lc)
                        full_result[self.molecular_formula][ion] = dt_peaks




    def peak_dt(self, df, z, i):
        ll = df.sort_values(by=['drift_time'])
        xnew = ll['drift_time'].values
        ynew = ll['intensity'].values
        mz = ll['mz'].values
        spec = ll['index'].values
        rt = ll['rt'].values
        theoretical_mass = ll['theoretical_mass'].values
        noise = estimate_noise(ynew)
        baseline = estimate_baseline(ynew, noise)
        start, peaks, end = detect_peaks(ynew, noise, baseline)
        peak_intensities  = np.array([ynew[x] for x in peaks])
        max_intensity = max(peak_intensities)
        peaks_index = peak_intensities > max_intensity/20
        peaks = peaks[peaks_index]
        if peaks.size > 0:
            self.peaks = True
            start = start[peaks_index]
            end = end[peaks_index]
            end = end - 1
            peak_mid = xnew[peaks]
            spec_number = spec[peaks]
            mz_mid = mz[peaks]
            rt_mid = rt[peaks]
            intensity_mid = ynew[peaks]
            half_max = intensity_mid * 0.5
            theoretical_mass = theoretical_mass[peaks]
            peak_start = xnew[start]
            peak_end = xnew[end]
            fwhm = find_fwhm(half_max, start, end, xnew, ynew)
            dict_all = {"mz_top": mz_mid, "peak_start": peak_start, "peak_end": peak_end, "dt_mid": peak_mid, "fwhm": fwhm, "spec_number": spec_number, "rt_mid": rt_mid, "theoretical_mass": theoretical_mass}
            df = pd.DataFrame(data = dict_all)
            df['rt'] = i
            self.dt_peaks[z][i] = df
            df['number of conformers'] = df['dt_mid'].size
            df['Error(PPM)'] = np.absolute(df['mz_top'].values - df['theoretical_mass'].values)
            df['Error(PPM)'] = np.divide(df['Error(PPM)'].values, df['theoretical_mass'].values)
            df['Error(PPM)'] = df['Error(PPM)'] * 1000000
            print(np.divide(peak_mid, fwhm))
            df['resolving_power'] = np.divide(peak_mid, fwhm)
            df['resolving_power'] = df['resolving_power'].map(lambda x: round(x, 2))

            final_result = df[["mz_top", "Error(PPM)", "number of conformers", "dt_mid", "fwhm", "rt_mid", "resolving_power", "rt"]]
            final_result.columns = ["mz_measured", "Error(PPM)", "#conformer", "drift_time", "fwhm", "retention_time", "resolving_power", "rt"]
            return final_result
        else:
            self.peaks = False
            return None
