
# from . import ccs
from . import formula
# from . import isotopic_confirm
from . import drift_time
from . import ccs
import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict

parameters = {"abundance_combobox": 17, "c13_combobox": 8, "mono_combobox": 5}

#making c12 as mono 100% and then calculate the c13 as pertage should be idearlly, and then subtract and then calculate the percentage of that difference.


#The primary data dict contains keys for all the primary analysis, keys are
#primary_ion
#drift_gas
#formula
#mzml
#calibration
#beta
#tfix
#buffer_text


#The secondary data dict contains keys for all the advanced settings.

#checked_ions
#positive_mass
#negative_mass
#mono_combobox
#c13_combobox
#abundance_combobox

class Final:

    def __init__(self, primary_data, secondary_data, optional_data, mass_values, function_values, temp_drift, temp_spectrum, message):

        self.primary_data =  primary_data
        self.secondary_data = secondary_data
        self.optional_data = optional_data
        self.ion_species = False
        self.mass_values = mass_values
        self.function_values = function_values
        self.temp_drift = temp_drift
        self.temp_spectrum = temp_spectrum
        self.message = message
        self.ppm_values = list(parameters.values())

        print(self.secondary_data)



    def run(self):
        self.unpack_primary()
        print(self.message)
        if not self.message['warning']:
            self.unpack_secondary()
            self.formula_df = formula.Molecular_Formula(self)
            self.ions, self.isotope_ratio = self.formula_df.run()
            self.drift_time_df = drift_time.Parse_MF(self)
            df9 = self.drift_time_df.dataframe_all()
            self.ccs_df = ccs.CCS(self, df9)
            df22 = self.ccs_df.finish_analysis()
            return df22






    def unpack_primary(self):
        xx = list(self.primary_data.values())
        tt = [bool(x) for x in xx]
        bool_value = np.all(tt)
        print(bool_value)
        if not bool_value:
            bb = self.optional_data['formula']
            x = isinstance(bb, list)
            print(x)
            if x:
                self.primary_data['buffer_text'] = self.optional_data['formula']
            bb = self.optional_data['calibration']
            x = isinstance(bb, dict)
            print(x)
            if x:
                print("doing")
                self.primary_data['tfix'] = bb['TFix']
                self.primary_data['beta'] = bb['Beta']
            xx = list(self.primary_data.values())
            print(xx)
            tt = [bool(x) for x in xx]
            bool_value = np.all(tt)
            if not bool_value:
                self.message['warning'].append("The compulsory fields should be filled, the mzml file,  tfix, beta value, (or .xml cal file should be uploaded) MF (or formula file should be upload), drift gas must be provided")






    def unpack_secondary(self):

        bool_value = bool(list(self.secondary_data.values()))
        if bool_value:
            self.secondary_data = OrderedDict(sorted(self.secondary_data.items()))
            print(self.secondary_data)

            self.ppm_values = [self.secondary_data.get(x, 0.0) if self.secondary_data.get(x, 0.0) != 0.0 else parameters.get(x, 0.0) for x in parameters.keys()]
            print(self.ppm_values)
            self.ppm_values = [float(x) for x in self.ppm_values]
            if bool(self.secondary_data.get("checked_ions", 0)):
                self.ion_species = True
                keys_to_search = set(self.secondary_data["checked_ions"])
                print(keys_to_search)
                mass_set = set(list(self.mass_values.keys()))
                print(mass_set)
                function_set = set(list(self.function_values.keys()))
                print(function_set)
                mass_set = mass_set.intersection(keys_to_search)
                print(mass_set)
                self.mass_values= {k: self.mass_values[k] for k in mass_set}
                function_set = function_set.intersection(keys_to_search)
                self.function_values= {k: self.function_values[k] for k in function_set}
