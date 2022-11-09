from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys
from PyQt5.uic import loadUiType
import datetime
from xlrd import *
from xlsxwriter import *
import sqlite3
from collections import defaultdict
import numpy as np
import pandas as pd
import io
import src
import tempfile
import pathlib
import asyncio
from PyQt5.QtCore import QThread
import contextlib
import time
import matplotlib.pyplot as plt
from src.formula import isotopic_table
from file_parser import *


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)



ui,_ = loadUiType(resource_path('ccs.ui'))

#primary_ion
#drift_gas
#formula
#mzml
#calibration
#beta
#tfix
#buffer_text

p_mass = {"[M+H]+": 1.007276, "[M+]": 0.000548579909, "[M-H]-": 1.007276, "[M-]": 0.000548579909}


theoretical_isotope = {"+electron": '', "-electron": '', "+Cl":  "Cl",  "+Br": "Br", "+CH3COO": "CH3COO", "+HCOO": "HCOO", "+H": "H", "-H": '', "+Na": 'Na', "+K": 'K', "+NH4": 'NH4'}

mass_values = {"+electron": 0.000548579909, "-electron": 0.000548579909, "+Cl":  34.969402,  "+Br": 78.918885, "+CH3COO": 59.013851, "+HCOO": 44.998201, "+H": 1.007276, "-H": 1.007276, "+Na": 22.989218, "+K": 38.963158, "+NH4": 18.033823}

function_values = {"+electron": np.subtract, "-electron": np.add, "+Cl": np.add,  "+Br": np.add, "+CH3COO": np.add, "+HCOO": np.add, "+H": np.add, "-H": np.subtract, "+Na": np.add, "+K": np.add, "+NH4": np.add}



# class PlainTextEdit(QPlainTextEdit):
#     def __init__(self, parent=None):
#         QPlainTextEdit.__init__(self, parent)
#         regexp = QRegExp("[A-Z][a-z]?\d*|\(.*?\)\d+")
#         self.validator= QRegExpValidator(regexp)
#
#     def keyPressEvent(self, event):
#         state = self.validator.validate(event.text(), 0)
#         if state[0] == QValidator.Acceptable:
#             QPlainTextEdit.keyPressEvent(self, event)




class MyAbstract(QThread):
    """Base export thread"""
    done = pyqtSignal(object)
    fail = pyqtSignal(object)
    loop = asyncio.get_event_loop()

    def __init__(self, func, parent=None):
        super().__init__(parent)
        self.func= func

    def run(self):
        try:
            result= self.loop.run_until_complete(self.func())

        except Exception as exception:
            print(exception)
            self.fail.emit(exception)
        else:
            self.done.emit(result)


class MainApp(QMainWindow , ui):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.primary_ion.currentTextChanged.connect(self.primary_ion_text)
        self.drift_gas.currentTextChanged.connect(self.drift_gas_text)
        self.mono_combobox.currentTextChanged.connect(self.mono_combobox_text)
        self.c13_combobox.currentTextChanged.connect(self.c13_combobox_text)
        self.abundance_combobox.currentTextChanged.connect(self.abundance_combobox_text)
        self.upload_mzml.clicked.connect(self.open_file_mzml)
        self.upload_formula.clicked.connect(self.open_file_formula)
        self.upload_calibration.clicked.connect(self.open_file_calibration)
        # self.findccs.clicked.connect(self.download3)
        self.findccs.clicked.connect(self.find_ccs)
        self.tfix.textEdited.connect(self.onTextEdited)
        self.beta.textEdited.connect(self.onTextEdited)
        self.formula.textChanged.connect(self.plaintextdata)
        self.positive_listwidget.itemClicked.connect(self.remove1)
        self.negative_listwidget.itemClicked.connect(self.remove2)

        self.progressbar.setStyle(QStyleFactory.create("windows"))
        self.progressbar.setRange(0, 1)
        self.progress = MyAbstract(self.find_ccs)
        self.delete_row.clicked.connect(self.table_deleterow)
        self.overlay_im.clicked.connect(self.table_xic)
        self.overlay_eic.clicked.connect(self.table_EIC)
        self.im.clicked.connect(self.table_im)
        self.distribution.clicked.connect(self.table_spectrum)
        self.download.clicked.connect(self.download_result)



        self.positive_electron.stateChanged.connect(lambda:self.btnstate(self.positive_electron))
        self.positive_hydrogen.stateChanged.connect(lambda:self.btnstate(self.positive_hydrogen))
        self.positive_sodium.stateChanged.connect(lambda:self.btnstate(self.positive_sodium))
        self.positive_potassium.stateChanged.connect(lambda:self.btnstate(self.positive_potassium))
        self.positive_ammonium.stateChanged.connect(lambda:self.btnstate(self.positive_ammonium))

        self.negative_electron.stateChanged.connect(lambda:self.btnstate(self.negative_electron))
        self.negative_hydrogen.stateChanged.connect(lambda:self.btnstate(self.negative_hydrogen))
        self.negative_chlorine.stateChanged.connect(lambda:self.btnstate(self.negative_chlorine))
        self.negative_bromine.stateChanged.connect(lambda:self.btnstate(self.negative_bromine))
        self.negative_formate.stateChanged.connect(lambda:self.btnstate(self.negative_formate))
        self.negative_acetate.stateChanged.connect(lambda:self.btnstate(self.negative_acetate))

        self.positive_pushbutton.clicked.connect(lambda:self.mass_added("positive_pushbutton"))
        self.negative_pushbutton.clicked.connect(lambda:self.mass_added("negative_pushbutton"))
        #
        reg_ex = QRegExp("[A-Z][a-z]?\d*|\(.*?\)\d+")
        input_validator = QRegExpValidator(reg_ex, self.formula)
        self.primary_data = {"primary_ion": "", "drift_gas": "Nitrogen", "mzml": "", "beta": "", "tfix": "", "buffer_text": ""}
        self.optional_data = defaultdict(lambda: None)
        self.secondary_data = defaultdict(list)
        self.message = defaultdict(list)






    @staticmethod
    def write_df_to_qtable(df, table):
        headers = list(df)
        table.setRowCount(df.shape[0])
        table.setColumnCount(df.shape[1])
        table.setHorizontalHeaderLabels(headers)

        # getting data from df is computationally costly so convert it to array first
        df_array = df.values
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                table.setItem(row, col, QTableWidgetItem(str(df_array[row,col])))



    def btnstate(self, b):
        if b.isChecked() == True:
            self.secondary_data["checked_ions"].append(b.text())
        elif b.isChecked() == False:
            index = self.secondary_data["checked_ions"].index(b.text())
            self.secondary_data["checked_ions"].pop(index)



    @pyqtSlot()
    def download3(self):

        def fail(exception):
            if not bool(self.n):
                self.n += 1
                print(str(self.message['warning']))
                QMessageBox.warning(self, "Warning", f"The output not obtained, because of {exception} {str(self.message['warning'])}")



        def done(result):
            if not bool(self.n):
                self.n += 1
                QMessageBox.information(self, "Information", f"The output result obtained successfully")


        def started():
            self.progressbar.setRange(0, 0)


        def finished():
            self.progressbar.setRange(0, 1)
            pass


        self.progress.started.connect(started)
        self.progress.finished.connect(finished)
        self.progress.done.connect(done)
        self.progress.fail.connect(fail)
        self.progress.start()




    def onTextEdited(self):
        tfix = self.tfix.text()
        beta = self.beta.text()
        self.primary_data["tfix"] = tfix
        self.primary_data["beta"] = beta


    def mass_added(self, b):
        if b == "positive_pushbutton":
            value = self.positive_lineedit.text()
            value = "+" + value
            self.positive_listwidget.addItem(value)
            self.secondary_data["positive_formula"].append(value)

        if b == "negative_pushbutton":
            value = self.negative_lineedit.text()
            value = "-" + value
            self.negative_listwidget.addItem(value)
            self.secondary_data["negative_formula"].append(value)




    def primary_ion_text(self, s):

        self.primary_data["primary_ion"] = s


    def drift_gas_text(self, s):

        self.primary_data["drift_gas"] = s


    def mono_combobox_text(self, s):
        self.secondary_data["mono_combobox"] = s


    def c13_combobox_text(self, s):
        self.secondary_data["c13_combobox"] = s

    def abundance_combobox_text(self, s):
        self.secondary_data["abundance_combobox"] = s



    def file_dialog(self, msg, path):
        return QFileDialog.getOpenFileName(self, msg, path)[0]


    def remove1(self, item):
         value = self.positive_listwidget.currentRow()
         session_name = self.positive_listwidget.currentItem().text()
         self.positive_listwidget.takeItem(value)
         self.secondary_data["positive_formula"].remove(session_name)


    def remove2(self, item):
         value = self.negative_listwidget.currentRow()
         session_name = self.negative_listwidget.currentItem().text()
         self.negative_listwidget.takeItem(value)
         self.secondary_data["negative_formula"].remove(session_name)




    def open_file_formula(self):
        try:

            msg = '1) Select the formula tab (txt) or excel file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select formula file'
            formula = self.file_dialog(sel, ".")
            print(formula)
            obj1 = MF_Parser(formula)
            formula_list = obj1.run()
            print(formula_list)
            self.optional_data['formula'] = formula_list

        except:

            self.message['warning'] = "Please upload a correct file"
            pass


    def open_file_mzml(self):
        try:

            msg = '1) Select the mzml file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select mzml file'
            mzml = self.file_dialog(sel, ".")
            self.primary_data['mzml'] = mzml

        except:

            self.message['warning'] = "Please upload a correct file"
            pass




    def open_file_calibration(self):
        try:

            msg = '1) Select the calibration xml file \n'
            QMessageBox.information(self, 'Add input file', msg)
            sel = 'Select calibration xml file'
            calibration = self.file_dialog(sel, ".")
            obj1 = XML_parser(calibration)
            obj1.run()
            print(obj1.result)
            self.optional_data['calibration'] = obj1.result
            self.tfix.setText(obj1.result["TFix"])
            self.beta.setText(obj1.result["Beta"])
            self.primary_data["tfix"] = float(obj1.result["TFix"])
            self.primary_data["beta"] = float(obj1.result["Beta"])


        except:

            self.message['warning'] = "Please upload a correct file"
            pass




    def plaintextdata(self):
        text = self.formula.toPlainText()
        print(text)
        formulas = text.split('\n')
        regex = re.compile(r"[A-Z][a-z]?\d*|\((?:[^()]*(?:\(.*\))?[^()]*)+\)\d+")
        formulas = [x for x in formulas if bool(x)]
        formula_list = [x for x in formulas if re.match(regex, str(x))]
        self.primary_data['buffer_text'] = formula_list




    def find_ccs(self):
        self.n = 0
        self.result.clear()
        self._temp = tempfile.TemporaryDirectory(prefix = "drift_time_")
        self.drift = pathlib.Path(self._temp.name).as_posix()
        self.temp = tempfile.TemporaryDirectory(prefix = "rt_isotopic_")
        self.spectrum = pathlib.Path(self.temp.name).as_posix()
        self.launcher = src.Final(self.primary_data, self.secondary_data, self.optional_data, mass_values, function_values, self.drift, self.spectrum, self.message)
        self.ccs_table = self.launcher.run()
        self.optional_data = defaultdict(lambda: None)
        show_table = self.ccs_table[['molecular_formula',  'ion_type',  'mz_measured', "Error(PPM)", "#conformer", "drift_time", "ccs",  "resolving_power", "retention_time"]]
        self.write_df_to_qtable(show_table, self.result)




    def table_deleterow(self):
        row = self.result.currentRow()
        self.ccs_table = self.ccs_table.reset_index(drop = True)
        self.ccs_table = self.ccs_table.drop(row)
        self.ccs_table = self.ccs_table.reset_index(drop = True)
        self.write_df_to_qtable(self.ccs_table, self.result)

#
# os.path.join(self.temp_spectrum, f"{self.molecular_formula}_rt_overlay.html"))
# os.path.join(self.data.temp_spectrum, f"{self.data.molecular_formula}_{ion}_{rt}.png"
# os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{ion}_{rt}.html"))
# os.path.join(self.data.temp_drift, f"{self.data.molecular_formula}_{p}_overlay.html"))



    def table_im(self):
        row = self.result.currentRow()
        ff = self.ccs_table.loc[row, 'molecular_formula']
        ss = self.ccs_table.loc[row, 'rt']
        ss = str(ss)
        tt = self.ccs_table.loc[row, 'ion_type']
        tt = str(tt)
        import webbrowser
        new = 2
        xx= os.path.join(self.drift, f'{ff}_{tt}_{ss}_IM.html')
        url = f"file://{xx}"
        webbrowser.open(url,new=new)



    def table_xic(self):
        row = self.result.currentRow()
        ff = self.ccs_table.loc[row, 'molecular_formula']
        ss = self.ccs_table.loc[row, 'rt']
        ss = str(ss)
        import webbrowser
        new = 2
        xx= os.path.join(self.drift, f'{ff}_{ss}_IMoverlay.html')
        url = f"file://{xx}"
        webbrowser.open(url,new=new)



    def table_EIC(self):
        row = self.result.currentRow()
        ff = self.ccs_table.loc[row, 'molecular_formula']
        import webbrowser
        new = 2
        xx= os.path.join(self.spectrum, f'{ff}_rt_overlay.html')
        url = f"file://{xx}"
        webbrowser.open(url,new=new)



    def table_spectrum(self):
        row = self.result.currentRow()
        ff = self.ccs_table.loc[row, 'molecular_formula']
        ss = self.ccs_table.loc[row, 'rt']
        ss = str(ss)
        tt = self.ccs_table.loc[row, 'ion_type']
        tt = str(tt)
        import webbrowser
        new = 2
        xx= os.path.join(self.spectrum, f'{ff}_{tt}_{ss}.html')
        print(xx)
        url = f"file://{xx}"
        webbrowser.open(url,new=new)




    def download_result(self):

        msg = '1) Select the Output folder to save all files \n'
        QMessageBox.information(self, 'Add Output folder', msg)
        path_text = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        # dlg = QtWidgets.QFileDialog()
        # dlg.setDirectory(str(os.getcwd()))
        # dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        # dlg.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        #
        # proxy = ProxyModel(dlg)
        # dlg.setProxyModel(proxy)
        #
        # if dlg.exec_():
        #     filenames = dlg.selectedFiles()
        #     path_text= QDir.toNativeSeparators(str(filenames[0]))
        #     print(path_text)
        show_table = self.ccs_table[['molecular_formula',  'ion_type',  'mz_measured', "#conformer", "drift_time", "ccs",  "resolving_power", "retention_time"]]
        show_table.to_csv(os.path.join(path_text, 'Results.tab'), sep ='\t')
        import shutil
        shutil.move(self.drift, path_text)
        shutil.move(self.spectrum, path_text)






def main():
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
