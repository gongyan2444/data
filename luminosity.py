'''
This small ui is designed to estimate the luminosity of any source. L= 4 pi d**2 Sv
This is written by Yan Gong. In case of problems, contact me via gongyan2444@gmail.com
'''
import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QPlainTextEdit, QGraphicsView)

from PyQt5.QtGui import QPalette, QFont

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class LEstimator(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        # Create widgets
        lbl_d = QLabel("Enter D (kpc):", self)
        self.txt_d = QLineEdit(self)
        lbl_d.setFont(QFont("Arial", 16))

        lbl_flux = QLabel("Enter ∫Sdv (Jy km s<sup>-1</sup>):", self)
        self.txt_flux = QLineEdit(self)
        #unit_label = QLabel(r"cm$^{-3}$")
        lbl_flux.setFont(QFont("Arial", 16))

        lbl_f = QLabel("Enter frequency (GHz):", self)
        self.txt_f = QLineEdit(self)
        lbl_f.setFont(QFont("Arial", 16))

        hbox_d = QHBoxLayout()
        hbox_d.addWidget(lbl_d)
        hbox_d.addWidget(self.txt_d)

        hbox_flux = QHBoxLayout()
        hbox_flux.addWidget(lbl_flux)
        hbox_flux.addWidget(self.txt_flux)
        #hbox_rho.addWidget(unit_label)

        hbox_f = QHBoxLayout()
        hbox_f.addWidget(lbl_f)
        hbox_f.addWidget(self.txt_f)

        
        self.btn_calculate = QPushButton("Calculate", self)
        self.btn_calculate.setFont(QFont("Arial", 16))

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_d)
        vbox.addLayout(hbox_flux)
        vbox.addLayout(hbox_f)
        vbox.addWidget(self.btn_calculate)

        self.result_label = QLabel(self)
        vbox.addWidget(self.result_label)
        self.result_label.setFont(QFont("Arial", 16))

        #self.result_box = QPlainTextEdit(self)
        #self.result_box.setReadOnly(True)
        #vbox.addWidget(self.result_box)

        
        # Connect signals and slots
        self.btn_calculate.clicked.connect(self.calLumi)

        # Set window properties
        self.setLayout(vbox)
        self.setWindowTitle('Luminosity Estimator')
        self.setGeometry(200, 200, 600, 200)

    def calLumi(self):
        try:
            d   = float(self.txt_d.text())
            flux = float(self.txt_flux.text())
            f    = float(self.txt_f.text())
            d = d*1000*3.086e18
            flux = flux*1e-23*1e5/2.99792458e10*f*1e9
            Ls = 4*np.pi*d**2*flux
            Ls_sun = Ls/3.9e33
            self.result_label.setText(f"The Luminosity is: {Ls:.2e} s<sup>-1</sup> or {Ls_sun:.2e} L<sub>☉</sub>")
            #self.result_box.setPlainText(r"The Luminosity is: {:.2e} s\u207B\u00B3".format(Ls))
            
        except ValueError:
            self.result_label.setText("Invalid input")
            #self.result_box.setPlainText("Invalid input")

if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = LEstimator()
    ex.show()
    sys.exit(app.exec_())
