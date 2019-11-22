# -*- coding: utf-8 -*-
import numpy as np
import os
import csv
import pandas as pd
from scipy.signal import savgol_filter
import math

class ExtractData:
    def __init__(self, filenames, milliamps, area):

        detection = DetectType(filenames)
        filetypes = detection.get_filetypes()
        print(filetypes)

        if filetypes[0] == ".cor":
            format = "PAR"
            parfiles = PARfiles(filenames)
            self.xdata, self.ydata, self.zdata, self.labels = parfiles.get_data(milliamps, area)
        elif filetypes[0] == ".csv" or ".CSV":
            format = "CSV"
            csvfiles = CSVfiles(filenames)
            #TODO: labels
            self.labels = []
            print("filenames: ", filenames[0])
            if "HPLC" in filenames[0]:
                self.labels = ["t / min", "Counts"]
                self.xdata, self.ydata, self.zdata = csvfiles.get_hplc_data()
            else:
                self.xdata, self.ydata, self.zdata = csvfiles.get_data(milliamps, area)
        elif filetypes[0] == ".txt":
            txtfiles = TXTfiles(filenames)
            self.xdata, self.ydata, self.zdata, self.labels = txtfiles.get_data()

    def get_data(self):
        return self.xdata, self.ydata, self.zdata, self.labels


class TXTfiles:
    def __init__(self, filenames):
        self.filenames = filenames

    def get_data(self):
        potentials = []
        currents = []
        exptimes = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            datafile = textfile.readlines()
            datafilelen = len(datafile)
            textfile.close()

            xdata = []
            ydata = []
            zdata = []
            labels = []



            # SIGNALTYPE  : EDS

            # SPECTRUM    : Spectral Data Starts Here

            import re


            for i in range(0, datafilelen):
                datafile[i] = re.sub(r"\s+", "", datafile[i], flags=re.UNICODE)


            for i in range(0, 100):
                print (datafile[i])
                data = datafile[i].split(":")
                if data[1] == "SpectralDataStartsHere":
                    print("index: ", i)
                    indexvalues = i
                    break


            xdata1 = []
            ydata1 = []
            for i in range(indexvalues+1, datafilelen-1):
                #print(datafile[i])
                #print (datafile[i].split(","))

                splitteddata = datafile[i].split(",")
                x1, y1 = splitteddata[0], splitteddata[1]
                xdata1.append(float(x1))
                ydata1.append(float(y1))


        xdata.append(xdata1)
        ydata.append(ydata1)

        #labels = [xlabel, ylabel, zlabel]



        return xdata, ydata, zdata, labels


class ExtractSpectraADData:
    def __init__(self, filenames, milliamps):
            format = "CSV"
            csvfiles = CSVfiles(filenames)
            # TODO: labels
            self.labels = []
            self.xdata, self.ydata, self.xdata2, self.ydata2 = csvfiles.get_data_spectraAD(milliamps)

    def get_data(self):
        return self.xdata, self.ydata, self.xdata2, self.ydata2, self.labels



class ExtractMultiSpectra:
    def __init__(self, filenames):
            format = "CSV"
            csvfiles = CSVfiles(filenames)
            # TODO: labels
            self.labels = []
            self.xdata, self.ydata = csvfiles.get_data_multiSpectra()

    def get_data(self):
        return self.xdata, self.ydata, self.labels

class ExtractSpectraEvo:
    def __init__(self, filenames):
            format = "CSV"
            csvfiles = CSVfiles(filenames)
            # TODO: labels
            self.labels = []
            self.xdata, self.ydata = csvfiles.get_data_SpectraEvo()

    def get_data(self):
        return self.xdata, self.ydata, self.labels


class ExtractDotplot:
    def __init__(self, filenames):
            format = "Excel"
            csvfiles = CSVfiles(filenames)
            # TODO: labels
            self.labels = []
            self.xdata, self.ydata = csvfiles.get_data_dotplot()

    def get_data(self):
        return self.xdata, self.ydata, self.labels

class ExtractCalplot:
    def __init__(self, filenames):
            format = "Excel"
            excelfiles = EXCELfiles(filenames)
            # TODO: labels
            self.labels = []
            self.xdata, self.ydata, self.desvstd = excelfiles.get_data_calplot()

    def get_data(self):
        return self.xdata, self.ydata, self.desvstd, self.labels

class DetectType:
    def __init__(self, filenames):
        self.filenames = filenames

    def get_filetypes(self):
        extensions = []
        for afile in self.filenames:
            filename, file_extension = os.path.splitext(afile)
            extensions.append(file_extension)
        return extensions

class EXCELfiles:
    def __init__(self, filenames):

        self.filenames = filenames

    def get_data_calplot(self):
        excelfile = pd.ExcelFile(self.filenames[0])
        csvrows1 = excelfile.parse()
        csvrows = csvrows1.get_values()
        numrows = len(csvrows)
        numcolumns = len(csvrows[0])

        xdata = [item[0] for item in csvrows]

        medias = []
        desvstd = []

        for i in range(numrows):
            items = csvrows[i][1:numcolumns]
            media1 = np.nanmean(items)
            desvstd1 = np.nanstd(items)

            medias.append(media1)
            desvstd.append(desvstd1)

        return xdata, medias, desvstd

class PARfiles:
    def __init__(self, filenames):

        self.filenames = filenames


    def get_data(self, milliamps, area):
        potentials = []
        currents = []
        exptimes = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            datafile = textfile.readlines()
            datafilelen = len(datafile)
            textfile.close()

            potential = []
            current = []
            exptime = []

            for i in range(0, 110):
                print(i)
                print (datafile[i].replace(" ", ""))
                if datafile[i] == "End Comments\n":
                    indexvalues = i
                elif datafile[i].replace(" ", "") == "EndComments:1\n":
                    indexvalues = i-1

            technique = datafile[indexvalues-3].split()
            print("technique: ", technique[2:4])

            if "Voltammogram" in technique[2:4]:
                print ("Voltammetry")
                technique = "Voltammetry"
            elif "Galvanostatic" in technique[2:4]:
                print ("Galvanostatic")
                technique = "Galvanostatic"
            elif "Galvanic" and "Square-Wave" in technique[2:4]:
                print("Galvanic Square-Wave")
                technique = "GalvanicSW"
            elif "Potentiostatic" in technique[2:4]:
                print("Potentiostatic")
                technique = "Potentiostatic"
            elif ("Potential" and "Scan/Hold") in technique[2:4]:
                print("Potential Scan/Hold")
                technique = "Potentiostatic"

            print ("Area (cm2): ", area)
            normalizearea = True

            for i in range(indexvalues+2, datafilelen):
                pot, cur, tim = datafile[i].split()[0:3]
                potential.append(pot)
                if milliamps:
                    cur = np.array(cur).astype(float) * 1000

                if normalizearea:
                    cur = cur/area
                    cur = cur

                current.append(cur)
                exptime.append(tim)

            xlabel, ylabel, zlabel = datafile[indexvalues-1].split()

            potentials.append(potential)
            currents.append(current)
            exptimes.append(exptime)

            if technique == "Galvanostatic":
                xdata = exptimes
                ydata = potentials
                zdata = currents
                labels = [zlabel, xlabel, ylabel]
            elif technique == "GalvanicSW":
                xdata = exptimes
                ydata = potentials
                zdata = currents
                labels = [zlabel, xlabel, ylabel]
            elif technique == "Potentiostatic":
                xdata = exptimes
                ydata = currents
                zdata = potentials
                labels = [zlabel, ylabel, xlabel]
            else:
                xdata = potentials
                ydata = currents
                zdata = exptimes
                labels = [xlabel, ylabel, zlabel]

        # indexes = []
        # for i in range(1, len(xdata[0])-2):
        #     print ("i: ", i)
        #     if xdata[0][i] == xdata[0][i-1]:
        #         indexes.append(i)
        #
        # def remove_by_indices(iter, idxs):
        #     return [e for i, e in enumerate(iter) if i not in idxs]
        #
        # xdata[0] = remove_by_indices(xdata[0], indexes)
        # ydata[0] = remove_by_indices(ydata[0], indexes)
        #
        # print ("xdata: ", xdata)
        # print ("ydata: ", ydata)

        print (xdata)

        return xdata, ydata, zdata, labels

class CSVfiles:
    def __init__(self, filenames):

        self.filenames = filenames
        self.datacsv = []
        #FIXME: detct dots or commas
        dots = True
        if dots:
            self.dotsreplace1 = ''
            self.dotsreplace2 = ''
            self.dotsreplace3 = ','
        else:
            self.dotsreplace1 = ','
            self.dotsreplace2 = '.'
            self.dotsreplace3 = '.'



    def get_hplc_data(self):
        print("Getting HPLC data")
        textfile = open(self.filenames[0], newline='', encoding='utf-16')
        reader = csv.reader(textfile, delimiter='\t')
        for row in reader:
            self.datacsv.append(row)
        textfile.close()
        self.numrows = len(self.datacsv)

        exptime = [np.asarray(self.get_column_data(0)).astype(float)]
        print("exptime: ", exptime)
        counts = [np.asarray(self.get_column_data(1)).astype(float)]
        print("counts: ", counts)
        zdata = []

        return exptime, counts, zdata


    def get_data(self, milliamps, area):
        potentials = []
        currents = []
        exptimes = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            reader = csv.reader(textfile, delimiter=';')
            for row in reader:
                self.datacsv.append(row)
            textfile.close()
            self.numrows = len(self.datacsv)
            numcolumns = len(self.datacsv[1])
            print("columns: ", numcolumns)

            potential = []
            current = []
            exptime = []

            print ("Area (cm2): ", area)
            normalizearea = True

            for i in range(0, numcolumns, 2):
                pot = self.get_column_data(i)
                potential.append(pot)
                cur = self.get_column_data(i+1)
                # FIXME: current range
                milliamps = False
                if milliamps:
                    cur = np.array(cur).astype(float) / 1000

                if normalizearea:
                    cur = np.asarray(cur).astype(float)/area
                    cur = cur

                current.append(cur)

            potentials.append(potential)
            currents.append(current)
            exptimes.append(exptime)

        return potentials[0], currents[0], exptimes


    def get_data_spectraAD(self, milliamps):
        times = []
        currents = []
        wavelenghts = []
        intcounts = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            reader = csv.reader(textfile, delimiter=';')
            for row in reader:
                self.datacsv.append(row)
            textfile.close()
            self.numrows = len(self.datacsv)
            numcolumns = len(self.datacsv[1])

            atime = []
            current = []
            alenght = []
            anint = []

            for i in range(0, numcolumns-1, 4):
                print (i)
                # if i == 4:
                #     continue
                btime = self.get_column_data(i)
                atime.append(btime)
                cur = self.get_column_data(i+1)
                # FIXME: current range
                milliamps = False
                if milliamps:
                    cur = np.array(cur).astype(float) / 1000
                current.append(cur)
                blenght = self.get_column_data(i+2)
                alenght.append(blenght)
                bint = self.get_column_data(i+3)
                anint.append(bint)

            times.append(atime)
            currents.append(current)
            wavelenghts.append(alenght)
            intcounts.append(anint)

        return times[0], currents[0], wavelenghts[0], intcounts[0]


    def get_data_multiSpectra(self):
        wavelenghts = []
        intcounts = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            reader = csv.reader(textfile, delimiter=';')
            for row in reader:
                self.datacsv.append(row)
            textfile.close()
            self.numrows = len(self.datacsv)
            numcolumns = len(self.datacsv[1])

            alenght = []
            anint = []

            for i in range(0, numcolumns, 2):
                print(i)
                blenght = self.get_column_data(i)
                alenght.append(blenght)
                bint = self.get_column_data(i+1)
                anint.append(bint)

            wavelenghts.append(alenght)
            intcounts.append(anint)

        return wavelenghts[0], intcounts[0]


    def get_data_dotplot(self):
        xdatas = []
        ydatas = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            reader = csv.reader(textfile, delimiter=';')
            for row in reader:
                self.datacsv.append(row)
            textfile.close()
            self.numrows = len(self.datacsv)
            numcolumns = len(self.datacsv[1])


            xdatafiles = self.get_column_data(0)

            ydatasfile = []

            for i in range(1, numcolumns):
                ydata = self.get_column_data(i)
                ydatasfile.append(ydata)

            ydatas.append(ydatasfile)
            xdatas.append(xdatafiles)

        return xdatas[0], ydatas[0]


    def get_data_SpectraEvo(self):
        wavelenghts = []
        intcounts = []
        for i in range(len(self.filenames)):
            textfile = open(self.filenames[i])
            reader = csv.reader(textfile, delimiter=';')
            for row in reader:
                self.datacsv.append(row)
            textfile.close()
            self.numrows = len(self.datacsv)
            numcolumns = len(self.datacsv[1])

            xdata = self.datacsv[1][10:numcolumns - 10]
            for i in range(len(xdata)):
                if type(xdata[i]) is not float:
                    xdata[i] = float(xdata[i].replace(self.dotsreplace3, '').replace(self.dotsreplace1, self.dotsreplace2))


            ydatatotal = []
            for k in range(self.numrows - 2):
                # print k
                if k < self.numrows - 3:
                    ydata = self.datacsv[k + 2][10:numcolumns - 10]
                    # ydata = ydata[20:numcolumns-80]
                    ydatatotal.append(ydata)

            for ydataA in ydatatotal:
                for j in range(len(ydataA)):
                    ydata = ydataA[j].replace(self.dotsreplace3, '').replace(self.dotsreplace1, self.dotsreplace2)
                    ydata = float(ydata)
                    ydataA[j] = ydata



            def find_nearest(array, value):
                idx = (np.abs(array - value)).argmin()
                return idx, array[idx]

            indstart, startvalue = find_nearest(np.array(xdata), 460)
            indend, endvalue = find_nearest(np.array(xdata), 700)


            for i in range(len(ydatatotal)):
                ydatanew = ydatatotal[i][indstart:indend]
                ydatax2 = savgol_filter(np.array(ydatanew).astype(np.float), 35, 1)
                ydatatotal[i] = ydatax2

            xdata = xdata[indstart:indend]

            wavelenghts.append(xdata)
            intcounts.append(ydatatotal)

            return wavelenghts[0], intcounts[0]





    def get_column_data(self, numcol):
        datatotal = []
        for i in range(0, self.numrows):
            data = self.datacsv[i][numcol].replace(self.dotsreplace3, '').replace(self.dotsreplace1, self.dotsreplace2)
            if data != '':
                data = float(data)
                datatotal.append(data)
        return datatotal

