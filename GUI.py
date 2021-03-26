# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:23:43 2019

@author: Thomas
"""
#from GUI_functions import *
from PyQt5.QtWidgets import (QLabel, QLineEdit, QWidget, QPushButton, QMessageBox,
                             QApplication, QFileDialog, QTableWidget, QTableWidgetItem,
                             QAbstractItemView)
from PyQt5 import QtCore

import experiment
import numpy as np
import sys

#from functools import partial

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.figure import Figure
#import matplotlib.pyplot as plt

#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# =============================================================================
# Graphical User Interface
# =============================================================================

class App(QWidget):
    def __init__(self):
        super().__init__()
        self.title = 'Calibration of Mössbauer Sources'
        self.left = 10
        self.top = 10
        self.width = 620
        self.height = 650
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.trig=False
        self.openFileNameDialog()
        
        
        self.name = QLabel("%s" %self.exp.title, self)
        self.name.move(5, 10)
        self.name.resize(600, 25)
        self.name.setStyleSheet("font: bold 12pt")
        
        
        self.dir_ = QLineEdit(self)
        self.dir_.setText(self.filename)
        self.dir_.move(5, 35)
        self.dir_.resize(600,20)
        self.dir_.returnPressed.connect(self.load_data)
        
        browse = QPushButton("Browse", self)
        browse.move(525 ,55)
        browse.clicked.connect(self.openFileNameDialog)
        
        save = QPushButton("Save", self)
        save.move(525 , 80)
        save.clicked.connect(self.save_cal)
        
        self.temp_txt = QLabel("Mean Temperature: %s C" %self.exp.temp, self)
        self.temp_txt.move(300, 70)
        self.temp_txt.resize(200,20)
        
        self.temp_txt = QLabel("Elapsed Time: %s h" %self.exp.time, self)
        self.temp_txt.move(300, 95)
        self.temp_txt.resize(200,20)

        
        self.cal_txt = QLabel("Calibration const.:", self)
        self.cal_txt.move(5, 70)
        self.cal_txt.resize(110, 20)
        
        self.cal = QLineEdit(self)
        self.cal.setText("%.6f" %self.exp.cal)
        self.cal.move(120, 70)
        self.cal.resize(90,20)
        
        
        self.zero_txt = QLabel("Zeroth Channel:", self)
        self.zero_txt.move(5, 95)
        self.zero_txt.resize(100, 20)
        
        self.zero = QLineEdit(self)
        self.zero.setText("%.4f" %self.exp.zero_ch)
        self.zero.move(120, 95)
        self.zero.resize(90,20)
        
        
        self.fold_txt = QLabel("Folding channel:", self)
        self.fold_txt.move(5, 120)
        self.fold_txt.resize(100, 20)
        
        self.fold = QLineEdit(self)
        self.fold.setText("%i" %self.exp.fold_ch)
        self.fold.move(120, 120)
        self.fold.resize(90,20)

        self.peak_txt = QLabel("No. peaks:", self)
        self.peak_txt.move(5, 170)
        self.peak_txt.resize(80, 20)
        
        self.peak = QLineEdit(self)
        self.peak.setText("6")
        self.peak.move(100, 170)
        self.peak.resize(15,20)        
        
        plot = QPushButton("Plot Spectrum", self)
        plot.move(120, 145)
        plot.resize(90, 20)
        plot.clicked.connect( self.plot_fold )
        
        
        cal = QPushButton("Calibrate", self)
        cal.move(120, 170)
        cal.resize(90, 20)
        cal.clicked.connect( self.calibrate )
        cal.setStyleSheet("font: bold 12pt")
        
        
        # Array of widths
        self.table = QTableWidget(self)
        self.table.move(5, 200)
        self.table.resize(600, 250)
        self.table.setColumnCount(3)
        self.table.setRowCount(6)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        
        self.table.setColumnWidth(0, 190)
        self.table.setColumnWidth(1, 190)
        self.table.setColumnWidth(2, 190)
        self.table.setHorizontalHeaderLabels(['Widths', 'Fitted Positions', 'Theoretical Positions [mm/s]'])
        self.update_table()
        
        
        # Weights in average input
        wTxt = QLabel("Set weights of 1st, 2nd, 3rd peak", self)
        wTxt.move(360, 125)
        wTxt.resize(200,20)
        wTxt.setStyleSheet("font: bold")

        
        self.wInput = QTableWidget(self)
        self.wInput.setGeometry(QtCore.QRect(230, 210, 311, 51))
#        self.wInput.setCornerButtonEnabled(True)
        self.wInput.setRowCount(1)
        self.wInput.setColumnCount(3)
#        self.wInput.setObjectName("wInput")
        self.w1 = QTableWidgetItem('50')
        self.w1.setTextAlignment(QtCore.Qt.AlignCenter)
        self.wInput.setItem(0, 0, self.w1)
        self.w2 = QTableWidgetItem('30')
        self.w2.setTextAlignment(QtCore.Qt.AlignCenter)
        self.wInput.setItem(0, 1, self.w2)
        self.w3 = QTableWidgetItem('20')
        self.w3.setTextAlignment(QtCore.Qt.AlignCenter)
        self.wInput.setItem(0, 2, self.w3)
        self.wInput.horizontalHeader().setVisible(True)
        self.wInput.horizontalHeader().setCascadingSectionResizes(False)
        self.wInput.horizontalHeader().setHighlightSections(True)
        self.wInput.horizontalHeader().setSortIndicatorShown(False)
        self.wInput.horizontalHeader().setStretchLastSection(True)
        self.wInput.verticalHeader().setVisible(False)
        self.wInput.verticalHeader().setSortIndicatorShown(False)
        self.wInput.verticalHeader().setStretchLastSection(False)
        
        self.wInput.move(295, 145)
        
        
        self.show()

        
    def openFileNameDialog(self):
        self.filename, _ = QFileDialog.getOpenFileName(self, "Select Experiment",
                                                  "C:\Md", "Experiment (*.exp)")
        if self.filename:
            if self.trig:
                self.exp = experiment.Experiment(self.filename)
                self.exp.fold()
                self.update_txt()
            else:
                self.exp = experiment.Experiment(self.filename)
                self.exp.fold()
                self.trig = True
            
            
    def load_data(self): 
        alert = QMessageBox()            
        try:
            self.exp = experiment.Experiment(self.dir_.text())
            self.exp.fold()
            self.update_txt()

            
        except FileNotFoundError:
            alert.setText("No such file exist!")
            alert.exec_()
        
        
    def update_txt(self):
        self.name.setText("%s" %self.exp.title)
        self.cal.setText("%.6f" %self.exp.cal)
        self.zero.setText("%.4f" %self.exp.zero_ch)
        self.fold.setText("%i" %self.exp.fold_ch)
        self.dir_.setText(self.filename)
        self.update_table()
        
        
        
    def update_values(self):
        self.exp.cal = float( self.cal.text() )
        self.exp.zero_ch = float( self.zero.text() )
        self.exp.fold_ch = int( self.fold.text() )
        
    def update_table(self):
        # Updating table
        N = int( self.peak.text() )
        popt = self.exp.fit_data(N, return_popt=True)
        
        self.table.setRowCount(N)
        self.table.setColumnCount(3)
        
        # Extracting widths and pos from popt
        widths = np.zeros(N)
        pos = np.zeros(N)
        for i in range(N):
            widths[i] = popt[4*i+1]*2
            pos[i]    = popt[4*i+2]
        
#        cal = np.array([-10.657/2, -6.167/2, -1.677/2,1.677/2, 6.167/2, 10.657/2 ])
        cal = np.array([1.677, 6.167, 10.657])
        P = int(N/2)
        cal = cal[:P]
        cal = np.append(cal, -cal)/2
        cal = np.sort(cal)
        
        for i in range(N):
            w = QTableWidgetItem("%.4f" %widths[i])
            w.setTextAlignment(QtCore.Qt.AlignCenter)
            pe = QTableWidgetItem("%.4f" %pos[i])
            pe.setTextAlignment(QtCore.Qt.AlignCenter)
            pt = QTableWidgetItem("%.4f" %cal[i])
            pt.setTextAlignment(QtCore.Qt.AlignCenter)
            
            self.table.setItem(i, 0,  w)
            self.table.setItem(i, 1,  pe)
            self.table.setItem(i, 2,  pt)
        
        
    def plot_fold(self):
        N = int( self.peak.text() )
        self.update_values()
        
        self.exp.fold()
        self.exp.plot(N)
#        self.canvas.ax.clear()
#        self.canvas1.plot(self.exp)
        
        
    def calibrate(self):
        N = int( self.peak.text() )
        self.update_values()
        
        weights = np.array([self.w1.text(), self.w2.text(), self.w3.text()]).astype(float)/100
        self.exp.calibration(N, weights)
        self.exp.optimal_fold(N)
        self.update_txt()
        
        self.update_table()
        
    def save_cal(self):
#        nr = self.exp.title.split()[0]
        filename=self.dir_.text()[:-4] + ".txt"
        with open(filename, 'w') as file:
            print("Calibration constant: %s \nZero Channel: %s \nFolding Channel: %s" %(self.cal.text(), self.zero.text(), self.fold.text()),
                  file=file)
        print("Calibration values have been saved to %s" %filename)
        

if __name__ == '__main__':
    app = QApplication([])
    app.setStyle('Fusion')
    
    ex = App()
    sys.exit(app.exec_())


