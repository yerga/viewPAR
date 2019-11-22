# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import sys, os
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import transforms
from PyQt5.QtCore import Qt, QSettings
from PyQt5.QtWidgets import QFileDialog, QApplication, QMainWindow, QSizePolicy, QInputDialog, QLineEdit, QDialog, \
    QVBoxLayout, QComboBox, QLabel, QDialogButtonBox, QCheckBox, QPushButton
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import numpy as np
from scipy import stats
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import math
from numpy import diff, trapz
import peakutils
from itertools import zip_longest
import csv

progname = "ViewPAR"
progversion = "0.1"

from ui_viewwindow import Ui_ViewWindow
import extractdata
from tools import find_nearest


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=4, height=3, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        fig.tight_layout()

        self.x = 0
        self.y = 0
        self.xpoint = 0.0
        self.ypoint = 0.0

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""


    def update_figure(self, data, overlay):
        self.xdata, self.ydata, self.zdata, self.labels = data[0], data[1], data[2], data[3]
        if not overlay:
            self.axes.cla()


        self.theplot = self.axes.plot(np.array(self.xdata[0]).astype(float), np.array(self.ydata[0]).astype(float))

        # self.axes.set_title(self.title)
        self.axes.set_xlabel(self.labels[0], fontsize=14)
        self.axes.set_ylabel(self.labels[1], fontsize=14)
        self.figure.tight_layout()
        self.draw()

    def clear_figure(self):
        print ("clearing")
        self.axes.cla()
        self.draw()


class ViewWindow (QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("viewPAR")
        self.area = 1

        self.filename = ""

        self.settings = QSettings('org.dyerga', 'viewPAR')
        self.lastpath = str(self.settings.value("lastpath"))

        self.ui = Ui_ViewWindow()
        self.ui.setupUi(self)

        self.actionOpen = self.ui.actionOpen
        self.actionOpen.triggered.connect(lambda: self.openFile(False))

        self.actionOverlay = self.ui.actionOverlay
        self.actionOverlay.triggered.connect(lambda: self.openFile(True))

        self.actionClear = self.ui.actionClear
        self.actionClear.triggered.connect(self.clear_figure)

        self.actionDeriv = self.ui.actionDeriv
        self.actionDeriv.triggered.connect(self.derivativeCurve)

        self.actionSmooth = self.ui.actionSmooth
        self.actionSmooth.triggered.connect(self.showSmoothdialog)

        self.actioniR = self.ui.actioniR
        self.actioniR.triggered.connect(self.showIRdialog)

        self.actionREF = self.ui.actionREF
        self.actionREF.triggered.connect(self.showREFdialog)

        self.actionTafel = self.ui.actionTafel
        self.actionTafel.triggered.connect(self.showTafeldialog)

        self.actionPeaks = self.ui.actionPeaks
        self.actionPeaks.triggered.connect(self.showFitdialog)

        self.actionExport = self.ui.actionExport
        self.actionExport.triggered.connect(lambda: self.exportData(self.canvas, False))

        self.actionMulti = self.ui.actionMulti
        self.actionMulti.triggered.connect(self.showMultiCVdialog)

        self.actionArea = self.ui.actionArea
        self.actionArea.triggered.connect(self.showAreaDialog)

        self.actionRotate = self.ui.actionRotate
        self.actionRotate.triggered.connect(self.showRotateDialog)

        self.actionIntTime = self.ui.actionIntTime
        self.actionIntTime.triggered.connect(self.integrateOverTime)

        self.actionMultiPulse = self.ui.actionmultiPulse
        self.actionMultiPulse.triggered.connect(self.openMultiPulse)

        self.actionMultiPulseQ = self.ui.actionmultiPulseQ
        self.actionMultiPulseQ.triggered.connect(self.chargeMultiPulse)

        self.centralLayout = self.ui.verticalLayout

        self.statusbar = self.ui.statusbar

        self.canvas = MyStaticMplCanvas(self, width=4, height=3, dpi=100)

        #self.canvas.mpl_connect('button_press_event', self.onclickCanvas)

        toolbar = NavigationToolbar(self.canvas, self)
        self.centralLayout.addWidget(toolbar)
        self.centralLayout.addWidget(self.canvas)

        self.setFocus()

    def showAreaDialog(self):
        area, okPressed = QInputDialog.getText(self, "Electrode area", u"Area (cm^2): ", QLineEdit.Normal, str(self.area))
        if okPressed and area != '':
            self.area = float(area)


    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def clear_figure(self):
        self.canvas.clear_figure()

    def updateStatusBar(self):
        self.statusbar.showMessage('xdata=%f  ydata=%f' % (self.canvas.xdata, self.canvas.ydata))

    def onclickCanvas(self, event):
        self.canvas.x = event.x
        self.canvas.y = event.y
        self.canvas.xpoint = event.xdata
        self.canvas.ypoint = event.ydata
        self.updateStatusBar()


    def openMultiPulse(self):
        self.filenames, _ = QFileDialog.getOpenFileNames(self,"Open file", self.lastpath,"All Files (*);;PAR Files (*.cor)")

        print(self.filenames)

        if self.filenames:
            self.lastpath = os.path.dirname(self.filenames[0])

            xdataNall = []
            ydataNall = []
            zdata = []
            for i in range(len(self.filenames)):
                xdataN, ydataN, zdataN, labelsN = self.extractdata([self.filenames[i]])

                ydataNall.append(ydataN[0])
                xdataNall.append(xdataN[0])


            print("xdatanall: ", xdataNall)

            xdataNall.reverse()
            for i in range(1, len(xdataNall)):
                xdataNall[i] = np.asarray(xdataNall[i]).astype(float) + float(xdataNall[i-1][-1])

            xdata = np.asarray([item for sublist in xdataNall for item in sublist]).astype(float)

            ydataNall.reverse()
            ydata = np.asarray([item for sublist in ydataNall for item in sublist]).astype(float)

            labels = ["t / s", "i / mA cm$^{-2}$"]

            self.canvas.update_figure([[xdata], [ydata], zdata, labels], False)


    def chargeMultiPulse(self):
        self.filenames, _ = QFileDialog.getOpenFileNames(self,"Open file", self.lastpath,"All Files (*);;PAR Files (*.cor)")

        if self.filenames:
            self.lastpath = os.path.dirname(self.filenames[0])

            charges = []
            pulses = []
            zdata = []
            for i in range(len(self.filenames)):
                xdataN, ydataN, zdataN, labelsN = self.extractdata([self.filenames[i]])

                tsec = np.asarray(xdataN[0]).astype(float)
                ydata = np.asarray(ydataN[0]).astype(float)

                #area = math.pi * (self.diameter / 2) ** 2
                #print("Area (cm2): ", area)

                ycurrent = ydata

                #ycurrent = np.asarray(ydata).astype(float) * area  # Convert density to total current

                xnew = np.linspace(tsec[0], tsec[-1], num=int(len(tsec) - 1))
                interpdata = interp1d(tsec, ycurrent, 'cubic')
                ynew = interpdata(xnew)

                y_int = integrate.cumtrapz(ynew, xnew, initial=0)

                chargebypulse = y_int[-1]
                charges.append(chargebypulse)
                pulses.append(i)


            charges.reverse()

            labels = ["pulse number", "Q / mC"]

            self.canvas.update_figure([[pulses], [charges], zdata, labels], False)





    def openFile(self, overlay):
        self.filename, _ = QFileDialog.getOpenFileName(self,"Open file", self.lastpath,"All Files (*);;PAR Files (*.cor)")

        if self.filename:
            self.lastpath = os.path.dirname(self.filename)
            xdata, ydata, zdata, labels = self.extractdata([self.filename])
            #TODO: detect labels
            labels = ["E / V", "i / mA cm$^{-2}$"]

            print ("overlay: ", overlay)
            #ydata = savgol_filter(ydata, 21, 3)

            self.canvas.update_figure([xdata, ydata, zdata, labels], overlay)


#            fig, ax = plt.subplots()
#            line1, = ax.plot(np.array(xdata[0]).astype(float), np.array(ydata[0]).astype(float))
#            plt.show()

            self.settings.setValue("lastpath", self.lastpath)

    def extractdata(self, filename):
        milliamps = True
        area = self.area
        exdata = extractdata.ExtractData(filename, milliamps, area)
        xdata, ydata, zdata, labels = exdata.get_data()

        return xdata, ydata, zdata, labels

    def showIRdialog(self):
        resistance, okPressed = QInputDialog.getText(self, "Solution resistance", u"R (\u03A9): ", QLineEdit.Normal, "")
        if okPressed and resistance != '':
            self.correctIR(float(resistance))

    def correctIR(self, resistance):
        numlines = len(self.canvas.axes.lines)
        newxs = []
        oldys = []
        for i in range(numlines):
            print ("line: ", i)
            xdata, ydata = self.canvas.axes.lines[i].get_data()
            #FIXME: detect current range
            ampohms = resistance/1000
            ytocur = np.array(ydata).astype(float) * self.area
            newxdata = np.array(xdata).astype(float) - ampohms*np.array(ytocur).astype(float)
            newxdata = newxdata.tolist()
            newxs.append(newxdata)
            oldys.append(ydata)

        self.canvas.axes.cla()
        for i in range(len(newxs)):
            self.canvas.update_figure([[newxs[i]], [oldys[i]], self.canvas.zdata, self.canvas.labels], True)


    def showRotateDialog(self):
        degrotate, okPressed = QInputDialog.getText(self, "Rotate", u"deg (\u03A9): ", QLineEdit.Normal, "")
        if okPressed and degrotate != '':
            self.Rotate(float(degrotate))

    def Rotate(self):

        data = np.random.randn(100)

        # first of all, the base transformation of the data points is needed
        base = plt.gca().transData
        rot = transforms.Affine2D().rotate_deg(45)

        numlines = len(self.canvas.axes.lines)
        for i in range(numlines):
            line = self.canvas.axes.lines[i]
            line.set_transform(rot + base)

        self.canvas.draw()


        # define transformed line
        #line = pyplot.plot(data, 'r--', transform=rot + base)
        # or alternatively, use:
        # line.set_transform(rot + base)

        #pyplot.show()



        # self.canvas.axes.cla()
        # for i in range(len(newxs)):
        #     self.canvas.update_figure([[newxs[i]], [oldys[i]], self.canvas.zdata, self.canvas.labels], True)

        #
        # numlines = len(self.canvas.axes.lines)
        # newxs = []
        # oldys = []
        # for i in range(numlines):
        #     print ("line: ", i)
        #     xdata, ydata = self.canvas.axes.lines[i].get_data()
        #
        #
        #
        #
        #     if yaxiscbox:
        #         datachanged = ydata
        #     else:
        #         datachanged = xdata
        #
        #     newdata = np.array(datachanged).astype(float) + addpH[i] + addref - overpotential
        #     newdata = newdata.tolist()
        #
        #     if yaxiscbox:
        #         newxs.append(xdata)
        #         oldys.append(newdata)
        #     else:
        #         newxs.append(newdata)
        #         oldys.append(ydata)





    def showREFdialog(self):

        refdialog = QDialog(self)
        refdialog.setWindowTitle("Convert Reference Potential")
        layout = self.setupUI_REFdialog(refdialog)

        refdialog.setLayout(layout)
        refdialog.setGeometry(300, 300, 300, 250)
        refdialog.show()

    def setupUI_REFdialog(self, refdialog):
        electrodes = ["Hg/HgO", "RHE", "SHE", "Ag/AgCl", "SCE"]
        electrodes2 = ["RHE", "Hg/HgO", "SHE", "Ag/AgCl", "SCE"]
        layout = QVBoxLayout()
        lbl = QLabel("Original Reference Electrode")
        combo = QComboBox()

        lbl2 = QLabel("New Reference Electrode")
        combo2 = QComboBox()

        for electrode in electrodes:
            combo.addItem(electrode)

        for electrode in electrodes2:
            combo2.addItem(electrode)

        combo.setCurrentIndex(2)

        layout.addWidget(lbl)
        layout.addWidget(combo)
        layout.addStretch()

        layout.addWidget(lbl2)
        layout.addWidget(combo2)
        layout.addStretch()

        lbl3 = QLabel("pH")
        pHedit = QLineEdit("13")

        lbl4 = QLabel("Overpotential (V)")
        overpotedit = QLineEdit("1.22")

        layout.addWidget(lbl3)
        layout.addWidget(pHedit)

        layout.addWidget(lbl4)
        layout.addWidget(overpotedit)

        yaxiscbox = QCheckBox("Y axis Potential")
        yaxiscbox.setChecked(False)
        layout.addWidget(yaxiscbox)

        layout.addStretch()

        buttonBox = QDialogButtonBox(refdialog)
        buttonBox.setGeometry(50, 240, 341, 32)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        buttonBox.accepted.connect(lambda: self.convertREFpotential(combo.currentText(), combo2.currentText(), pHedit.text(), overpotedit.text(), yaxiscbox.isChecked()))
        buttonBox.accepted.connect(refdialog.accept)
        buttonBox.rejected.connect(refdialog.reject)
        layout.addWidget(buttonBox)

        return layout


    def showSmoothdialog(self):

        dialog = QDialog(self)
        dialog.setWindowTitle("Smooth plot")
        layout = self.setupUI_Smoothdialog(dialog)

        dialog.setLayout(layout)
        dialog.setGeometry(300, 300, 300, 250)
        dialog.show()

    def setupUI_Smoothdialog(self, dialog):
        layout = QVBoxLayout()
        lbl = QLabel("Number points")
        polyorder = QLineEdit()

        lbl2 = QLabel("Polyorder")
        numpoints = QLineEdit()

        layout.addWidget(lbl)
        layout.addWidget(numpoints)
        layout.addStretch()

        layout.addWidget(lbl2)
        layout.addWidget(polyorder)
        layout.addStretch()

        buttonBox = QDialogButtonBox(dialog)
        buttonBox.setGeometry(50, 240, 341, 32)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        buttonBox.accepted.connect(lambda: self.smoothCurves(int(numpoints.text()), int(polyorder.text())))
        buttonBox.accepted.connect(dialog.accept)
        buttonBox.rejected.connect(dialog.reject)
        layout.addWidget(buttonBox)

        return layout

    def smoothCurves(self, numpoints, polyorder):
        numlines = len(self.canvas.axes.lines)
        newxs = []
        oldys = []
        for i in range(numlines):
            print ("line: ", i)
            xdata, ydata = self.canvas.axes.lines[i].get_data()

            newydata = savgol_filter(np.array(ydata).astype(float), numpoints, polyorder)

            newdata = newydata.tolist()

            newxs.append(xdata)
            oldys.append(newdata)


        self.canvas.axes.cla()
        for i in range(len(newxs)):
            self.canvas.update_figure([[newxs[i]], [oldys[i]], self.canvas.zdata, self.canvas.labels], True)



    def derivativeCurve(self):
        numlines = len(self.canvas.axes.lines)
        newxs = []
        oldys = []
        for i in range(numlines):
            print ("line: ", i)
            xdata, ydata = self.canvas.axes.lines[i].get_data()

            dydx = diff(np.array(ydata).astype(float)) / diff(np.array(xdata).astype(float))

            newdata = dydx.tolist()

            newxs.append(xdata[:-1])
            oldys.append(newdata)


        self.canvas.axes.cla()
        for i in range(len(newxs)):
            self.canvas.update_figure([[newxs[i]], [oldys[i]], self.canvas.zdata, self.canvas.labels], True)



    def convertREFpotential(self, ref1, ref2, pH, overpotential, yaxiscbox):
        print("convert")
        print("Ref1: ", ref1)
        print("Ref2: ", ref2)
        print("pH: ", pH)

        pH = pH.split(";")

        if ref1 == ref2:
            print ("same electrode, same potential")
        elif ref1 == "Ag/AgCl" and ref2 == "RHE":
            addpH = 0.0591*np.array(pH).astype(float)
            addref = 0.197
        elif ref1 == "Hg/HgO" and ref2 == "RHE":
            addpH = 0.0591*np.array(pH).astype(float)
            addref = 0.140
        elif ref1 == "Hg/HgO" and ref2 == "Ag/AgCl":
            addpH = 0*np.array(pH).astype(float)
            addref = 0.057

        if overpotential:
            overpotential = float(overpotential)
        else:
            overpotential = 0


        numlines = len(self.canvas.axes.lines)
        newxs = []
        oldys = []
        for i in range(numlines):
            print ("line: ", i)
            xdata, ydata = self.canvas.axes.lines[i].get_data()
            if yaxiscbox:
                datachanged = ydata
            else:
                datachanged = xdata

            newdata = np.array(datachanged).astype(float) + addpH[i] + addref - overpotential
            newdata = newdata.tolist()

            if yaxiscbox:
                newxs.append(xdata)
                oldys.append(newdata)
            else:
                newxs.append(newdata)
                oldys.append(ydata)


        self.canvas.axes.cla()
        for i in range(len(newxs)):
            self.canvas.update_figure([[newxs[i]], [oldys[i]], self.canvas.zdata, self.canvas.labels], True)



    def showTafeldialog(self):

        tafeldialog = QDialog(self)
        tafeldialog.setWindowTitle("Tafel plot")
        layout = self.setupUI_Tafeldialog(tafeldialog)

        tafeldialog.setLayout(layout)
        tafeldialog.setGeometry(300, 300, 300, 250)
        tafeldialog.show()

    def setupUI_Tafeldialog(self, tafeldialog):
        layout = QVBoxLayout()
        lbl = QLabel("Initial X value")
        startpot = QLineEdit()

        lbl2 = QLabel("Final X value")
        endpot = QLineEdit()

        layout.addWidget(lbl)
        layout.addWidget(startpot)
        layout.addStretch()

        layout.addWidget(lbl2)
        layout.addWidget(endpot)
        layout.addStretch()

        tafelcurrcb = QCheckBox("Current is X data")
        tafelcurrcb.setChecked(True)
        layout.addWidget(tafelcurrcb)
        layout.addStretch()

        buttonBox = QDialogButtonBox(tafeldialog)
        buttonBox.setGeometry(50, 240, 341, 32)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        buttonBox.accepted.connect(lambda: self.calculateTafel(startpot.text(), endpot.text(), tafelcurrcb.isChecked()))
        buttonBox.accepted.connect(tafeldialog.accept)
        buttonBox.rejected.connect(tafeldialog.reject)
        layout.addWidget(buttonBox)

        return layout


    def calculateTafel(self, startvalue, endvalue, tafelcurr):
        #FIXME: get all lines
        #line = self.canvas.theplot[0]
        xdata, ydata = self.canvas.theplot[0].get_data()

        if tafelcurr:
            dataregion = xdata
        else:
            dataregion = xdata[0:int(len(xdata) / 2)]

        ind1, indvalue1 = find_nearest(np.array(dataregion).astype(float), float(startvalue))
        ind2, indvalue2 = find_nearest(np.array(dataregion).astype(float), float(endvalue))

        print(ind1, indvalue1)
        print(ind2, indvalue2)

        xdata = xdata[ind1:ind2]
        ydata = ydata[ind1:ind2]

        print("xdata: ", xdata)
        print("ydata: ", ydata)

        if tafelcurr:
            lognumber = np.log10(np.array(xdata).astype(float))
            yplot = ydata
        else:
            lognumber = np.log10(np.array(ydata).astype(float))
            yplot = xdata

        fig, axis = plt.subplots(1,1)
        line2, = axis.plot(lognumber, yplot, ls="-")
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(lognumber, np.array(yplot).astype(float))
        linefit = slope1 * np.array(lognumber) + intercept1
        axis.plot(lognumber, linefit, ':', label=str(int(round(slope1, 3) * 1000)) + " mV/dec")
        # print ('slope: ', slope1, '; intercept: ', intercept1)
        axis.set_xlabel("Log i")
        axis.set_ylabel("E")
        leg = axis.legend(loc='best', shadow=False)
        fig.tight_layout()
        #print("lines: ", axis.lines)
        savebutton = QPushButton("CSV")
        savebutton.clicked.connect(lambda: self.exportData(axis, True))
        fig.canvas.manager.toolbar.addWidget(savebutton)
        fig.show()


    def showFitdialog(self):

        fitdialog = QDialog(self)
        fitdialog.setWindowTitle("Peak search")
        layout = self.setupUI_Fitdialog(fitdialog)

        fitdialog.setLayout(layout)
        fitdialog.setGeometry(300, 300, 300, 250)
        fitdialog.show()

    def setupUI_Fitdialog(self, fitdialog):
        layout = QVBoxLayout()
        lbl = QLabel("Start potential (V)")
        startpot = QLineEdit()

        lbl2 = QLabel("End potential (V)")
        endpot = QLineEdit()

        forwcbox = QCheckBox("Forward curve")
        forwcbox.setChecked(True)

        layout.addWidget(lbl)
        layout.addWidget(startpot)
        layout.addStretch()

        layout.addWidget(lbl2)
        layout.addWidget(endpot)
        layout.addStretch()

        layout.addWidget(forwcbox)
        layout.addStretch()


        lbl3 = QLabel("Scan rate (V)")
        scanrate = QLineEdit()
        layout.addWidget(lbl3)
        layout.addWidget(scanrate)
        layout.addStretch()

        buttonBox = QDialogButtonBox(fitdialog)
        buttonBox.setGeometry(50, 240, 341, 32)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        buttonBox.accepted.connect(lambda: self.fitCurve(startpot.text(), endpot.text(), forwcbox.isChecked(), scanrate.text()))
        buttonBox.accepted.connect(fitdialog.accept)
        buttonBox.rejected.connect(fitdialog.reject)
        layout.addWidget(buttonBox)

        return layout


    def fitCurve(self, startpot, endpot, forward, scanrate):
        #FIXME: get all lines
        #line = self.canvas.theplot[0]
        xdata, ydata = self.canvas.theplot[0].get_data()

        ydata = savgol_filter(np.array(ydata).astype(np.float), 11, 1)


        #FIXME: when not cyclic voltammetry (same number of forward/backward data points)
        if forward:
            limit0 = 0
            limit1 = int(len(xdata)/2)
            offset = len(xdata) - limit1
        else:
            limit0 = int(len(xdata)/2)
            limit1 = len(xdata)
            offset = limit1-limit0


        #limit0 = 0
        #limit1 = len(xdata)
        #offset = limit1-limit0

        print("limits: ", limit0, "  ", limit1)
        print("potentials: ", startpot, "  ", endpot)

        ind1, indvalue1 = find_nearest(np.array(xdata[limit0:limit1]).astype(float), float(startpot))
        ind2, indvalue2 = find_nearest(np.array(xdata[limit0:limit1]).astype(float), float(endpot))

        print (ind1, indvalue1)
        print (ind2, indvalue2)

        if not forward:
            ind1 += offset
            ind2 += offset

        print("ind1, ind2: ", ind1, ind2)

        xdata = xdata[ind1:ind2]
        ydata = ydata[ind1:ind2]

        xdata = np.array(xdata).astype(float)
        ydata = np.array(ydata).astype(float)

        if not forward:
            ydata = -ydata

        baseline = peakutils.baseline(ydata, 1)
        diff = np.array(ydata) - np.array(baseline)

        from scipy.integrate import simps


        # Compute the area using the composite trapezoidal rule.
        xdatatime = np.array(xdata).astype(float) / float(scanrate)  # convert to s
        print ("potentia range: ", xdata)
        print("scanrate: ", scanrate)
        print("time: ", xdatatime)
        ydataamps = np.array(diff).astype(float) / 1000  # convert to A
        print ("currents: ", ydataamps)
        area = trapz(y=ydataamps, x=xdatatime, dx=100)
        print("area = %.6f C" % area)


        indexes = peakutils.indexes(ydata, thres=0.01, min_dist=0.01)
        print(indexes)

        peaks_x = peakutils.interpolate(xdata, ydata, ind=indexes)
        print("peaks: ", peaks_x)

        peakheights = []

        for i in range(len(indexes)):
            height = ydata[indexes[i]] - baseline[indexes[i]]
            peakheights.append(height)

        peakpotentials = xdata[indexes]

        print("peak potentials: ", peakpotentials)
        print("peak heights: ", peakheights)

        #
        # from scipy import optimize
        #
        # def gaussian(x, height, center, width, offset):
        #     #return height * np.exp(-(x - center) ** 2 / (2 * width ** 2)) + offset
        #
        #     return height * width ** 2 / ((x - center) ** 2 + width ** 2)
        #
        # def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3, offset):
        #     return (gaussian(x, h1, c1, w1, offset=0) +
        #             gaussian(x, h2, c2, w2, offset=0) +
        #             gaussian(x, h3, c3, w3, offset=0) + offset)
        #
        # def two_gaussians(x, h1, c1, w1, h2, c2, w2, offset):
        #     return three_gaussians(x, h1, c1, w1, h2, c2, w2, 0, 0, 1, offset)
        #
        # #errfunc3 = lambda p, x, y: (three_gaussians(x, *p) - y) ** 2
        # errfunc2 = lambda p, x, y: (two_gaussians(x, *p) - y) ** 2
        #
        # #guess3 = [0.49, 0.55, 0.01, 0.6, 0.61, 0.01, 1, 0.64, 0.01, 0]
        # # I guess there are 3 peaks, 2 are clear, but between them there seems to be another one, based on the change in slope smoothness there
        # guess2 = [8, 0.07, 0.01, 2, 0.27, 0.01, 0]  # I removed the peak I'm not too sure about
        # #optim3, success = optimize.leastsq(errfunc3, guess3[:], args=(xdata, diff))
        # optim2, success = optimize.leastsq(errfunc2, guess2[:], args=(xdata, diff))
        #
        # plt.plot(xdata, diff, lw=5, c='g', label='measurement')
        # #plt.plot(xdata, three_gaussians(xdata, *optim3), lw=3, c='b', label='fit of 3 Gaussians')
        # plt.plot(xdata, two_gaussians(xdata, *optim2), lw=1, c='r', ls='--', label='fit of 2 Gaussians')
        # plt.legend(loc='best')
        # #plt.savefig('result.png')
        # plt.show()

        # from scipy.optimize import curve_fit
        #
        # from scipy.special import erf
        #
        # def asym_peak(t, pars):
        #     'from Anal. Chem. 1994, 66, 1294-1301'
        #     a0 = pars[0]  # peak area
        #     a1 = pars[1]  # elution time
        #     a2 = pars[2]  # width of gaussian
        #     a3 = pars[3]  # exponential damping term
        #     f = (a0 / 2 / a3 * np.exp(a2 ** 2 / 2.0 / a3 ** 2 + (a1 - t) / a3)
        #          * (erf((t - a1) / (np.sqrt(2.0) * a2) - a2 / np.sqrt(2.0) / a3) + 1.0))
        #     return f
        #
        # def two_peaks(t, *pars):
        #     'function of two overlapping peaks'
        #     a10 = pars[0]  # peak area
        #     a11 = pars[1]  # elution time
        #     a12 = pars[2]  # width of gaussian
        #     a13 = pars[3]  # exponential damping term
        #     a20 = pars[4]  # peak area
        #     a21 = pars[5]  # elution time
        #     a22 = pars[6]  # width of gaussian
        #     a23 = pars[7]  # exponential damping term
        #     p1 = asym_peak(t, [a10, a11, a12, a13])
        #     p2 = asym_peak(t, [a20, a21, a22, a23])
        #     return p1 + p2
        #
        # parguess = (50, 0.07, 0.05, 0.1,    50, 0.27, 0.05, 0.1)
        # popt, pcov = curve_fit(two_peaks, xdata, diff, parguess)
        #
        # pars1 = popt[0:4]
        # pars2 = popt[4:8]
        #
        # peak1 = asym_peak(xdata, pars1)
        # peak2 = asym_peak(xdata, pars2)
        #
        # plt.figure()
        # plt.plot(xdata, diff)
        # plt.plot(xdata, peak1, 'r-')
        # plt.plot(xdata, peak2, 'g-')
        # plt.show()

        #a,b,c = peakutils.gaussian_fit(xdata, ydata, center_only=False)
        #print("gauss: ", a, " ", b, " ", c)
        #ygauss = peakutils.gaussian(xdata, a,b,c)


        # from scipy.optimize import curve_fit
        # from scipy import asarray as ar, exp
        #
        # n = len(xdata)  # the number of data
        # mean = sum(xdata * ydata) / n  # note this correction
        # sigma = sum(ydata * (xdata - mean) ** 2) / n  # note this correction
        #
        # def gaus(x, a, x0, sigma):
        #     #return a * sigma ** 2 / ((x - x0) ** 2 + sigma ** 2)
        #     return a * exp(-(x - x0) ** 2 / (2 * sigma ** 2))
        #
        # popt, pcov = curve_fit(gaus, xdata, ydata, p0=[1, mean, sigma])
        #
        # plt.plot(xdata, ydata, 'b+:', label='data')
        # plt.plot(xdata, gaus(xdata, *popt), 'ro:', label='fit')
        # plt.legend()
        # plt.title('Fig. 3 - Fit for Time Constant')
        # plt.xlabel('Time (s)')
        # plt.ylabel('Voltage (V)')
        # plt.show()



        # import fitraman
        # # TODO: adjust number of peaks to find, initialwidth, curve type, sigmavalue
        # peaks_to_find = 2
        # initialwidth = 0.05
        # fitraman.CURVE = "Gaussian"
        # # fitraman.SIGMAVALUE = np.full(len(subtract), 5)
        # params, fit, ys, n_peaks = fitraman.predict_and_plot_lorentzians(xdata, diff, peaks_to_find, initialwidth)
        # print ('params: ', params)
        #
        # peakdata = []
        # for j in range(0, len(params), 3):
        #     ctr = params[j]
        #     amp = params[j + 1]
        #     width = params[j + 2]
        #     peakdata.append(["%.2f" % ctr, "%.2f" % amp, "%.2f" % width])
        #     ysignal = fitraman.lorentzian(xdata, amp, ctr, width)
        #     ymax = np.max(ysignal)
        #     idxmax = np.argmax(ysignal)
        #     # plot max points in fitted curves
        #     # plt.plot(xdata[idxmax], ymax, ls='', marker='x')
        #     #Plot max points in experimental curve
        #     plt.plot(xdata[idxmax], diff[idxmax], ls='', marker='x')
        #     plt.plot(xdata, ysignal, ls='-', label="ysignal")
        #
        # #self.printpeakdata(peakdata)
        # #plt.plot(xdata, fit, 'r-', label='fit', c='red', lw=1.2, ls='--')
        # plt.plot(xdata, diff, label="data")
        #
        # plt.legend(loc='upper right', fontsize=10, shadow=True)
        # plt.show()



        if not forward:
            ydata = -ydata
            diff = -diff
            baseline = -baseline




        fig, axis = plt.subplots(1, 1)
        #lognumber = np.log10(np.array(ydata).astype(float))
        line2, = axis.plot(xdata, ydata, ls="-", label="data")
        line3, = axis.plot(xdata, baseline, ls=":", label="baseline")
        line3, = axis.plot(xdata, diff, ls=":", label="corrected")
        #line4, = axis.plot(xdata[indexes], ydata[indexes], "r+")
        #line5, = axis.plot(xdata, ygauss, ls=":", label="gauss")
        #slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(lognumber, np.array(xdata).astype(float))
        #linefit = slope1 * np.array(lognumber) + intercept1
        #axis.plot(lognumber, linefit, ':', label=str(int(round(slope1, 3) * 1000)) + " mV/dec")
        # print ('slope: ', slope1, '; intercept: ', intercept1)
        axis.set_xlabel("E")
        axis.set_ylabel("i")
        plt.title("Charge = %.6f C" % area, fontsize=12)
        leg = axis.legend(loc='best', shadow=False)
        fig.tight_layout()
        fig.show()


    def exportData(self, canvas, tafel):
        if tafel:
            axes = canvas
        else:
            axes = canvas.axes

        numlines = len(axes.lines)
        datas = []
        for i in range(numlines):
            x, y = axes.lines[i].get_data()
            datas.append(x)
            datas.append(y)

        export_data = zip_longest(*datas, fillvalue='')

        self.saveFile(export_data)

    def saveFile(self, export_data):
        savefilename, _ = QFileDialog.getSaveFileName(self,"Save file", self.lastpath,"All Files (*);;CSV Files (*.csv)")

        if savefilename:
            with open(savefilename, 'w', encoding="ISO-8859-1", newline='') as f:
                wr = csv.writer(f, delimiter=";")
                #FIXME: add two initial rows? -> check plotscifigs
                #wr.writerow(("List1", "List2", "List3", "List4"))
                wr.writerows(export_data)

            f.close()


    def showMultiCVdialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Show multiCV")
        layout = self.setupUI_Multidialog(dialog)

        dialog.setLayout(layout)
        dialog.setGeometry(300, 300, 300, 250)
        dialog.show()

    def setupUI_Multidialog(self, dialog):
        layout = QVBoxLayout()

        lbl3 = QLabel("Number of CVs")
        numcvs = QLineEdit("")

        lbl4 = QLabel("Show CVs")
        showcvs = QLineEdit("")

        layout.addWidget(lbl3)
        layout.addWidget(numcvs)

        layout.addWidget(lbl4)
        layout.addWidget(showcvs)

        layout.addStretch()

        buttonBox = QDialogButtonBox(dialog)
        buttonBox.setGeometry(50, 240, 341, 32)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        buttonBox.accepted.connect(
            lambda: self.showMultiCVs(numcvs.text(), showcvs.text()))
        buttonBox.accepted.connect(dialog.accept)
        buttonBox.rejected.connect(dialog.reject)
        layout.addWidget(buttonBox)

        return layout


    def showMultiCVs(self, numcvs, showcvs):
        x, y = self.canvas.axes.lines[0].get_data()

        lenx = len(x)
        print("lenx: ", lenx)

        leneachcv = int(lenx/int(numcvs))

        xdata = []
        ydata = []

        for i in range(0, lenx, leneachcv):
            print("indexcv: ", i, "leneachcv: ", leneachcv+i-1)
            onecvx = x[i:leneachcv+i-1]
            onecvy = y[i:leneachcv+i-1]
            xdata.append(onecvx)
            ydata.append(onecvy)


        showcvs = showcvs.split(";")

        newcvxdata = []
        newcvydata = []

        for i in range(len(showcvs)):
            newcvx = xdata[int(showcvs[i])]
            newcvy = ydata[int(showcvs[i])]
            newcvxdata.append(newcvx)
            newcvydata.append(newcvy)

        self.canvas.axes.cla()
        for i in range(len(newcvxdata)):
            self.canvas.update_figure([[newcvxdata[i]], [newcvydata[i]], self.canvas.zdata, self.canvas.labels], True)


    def integrateOverTime(self):
        tsec, ydata = self.canvas.theplot[0].get_data()

        print("Area (cm2): ", self.area)

        ycurrent = np.asarray(ydata).astype(float)*self.area # Convert density to total current

        xnew = np.linspace(tsec[0], tsec[-1], num=int(len(tsec)-1))
        interpdata = interp1d(tsec, ycurrent, 'cubic')
        ynew = interpdata(xnew)

        y_int = integrate.cumtrapz(ynew, xnew, initial=0)
        #plt.plot(xnew, y_int, 'ro')
        #plt.show()

        self.canvas.update_figure([[xnew], [y_int], self.canvas.zdata, self.canvas.labels], False)


if __name__ == "__main__":
    qApp = QApplication(sys.argv)

    aw = ViewWindow()
    aw.setWindowTitle("%s" % progname)
    aw.show()
    sys.exit(qApp.exec_())
