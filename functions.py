from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import time

from setuptools.command.rotate import rotate

# ===================================== Set Options ======================================
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 300)
pd.set_option('display.float_format', '{:,.5f}'.format)

# Colors: Console
greyDark = '\033[38;2;144;144;144m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
blue = '\033[38;5;51m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenLightB = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'



# =================================== Define Functions ===================================
def pressKey(event):
    if event.key == 'escape':
        plt.close()
    elif event.key == 'backspace':
        sys.exit()



class BrainData:
    def __init__(self, pathFolder, perAT8Cutoff, plotAT8, plotAT8Cutoff,
                 printData=False, NBars=None):
        # Parameters: Files
        self.pathFolder = pathFolder

        # Parameters: Data
        self.patientsOfInterest = [] # Selected donors with localized AT8
        self.dataPOI = pd.DataFrame()
        self.perAT8 = None
        self.perAT8High = None

        # Parameters: Figures
        self.plotAT8 = plotAT8
        self.plotAT8Cutoff = plotAT8Cutoff
        self.figSize = (9.5, 8)
        self.figSizeTall = (12, 9.5)
        self.figSizeW = (16, 8)
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.NBars = NBars

        # Parameters: Variables
        self.perAT8Cutoff = perAT8Cutoff
        self.printData = printData
        


    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'bars':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF','#ABFF9B','#39FF14','#2E9418','#2E9418',
                          '#005000']
            else:
                # Orange
                colors = ['white','white','#FF76FA','#FF50F9','#FF00F2',
                          '#CA00DF','#BD16FF']
        elif colorType == 'stdev':
            colors = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
        elif colorType == 'word cloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#CC5500','#F79620','#FAA338','#00C01E','#003000','black']
            # colors = ['#008631','#39E75F','#CC5500','#F79620','black']
        elif colorType == 'em':
            colors = ['navy','royalblue','dodgerblue','lightskyblue','white','white',
                      'lightcoral','red','firebrick','darkred']
        else:
            print(f'{orange}ERROR: Cannot create colormap. '
                  f'Unrecognized colorType parameter: {colorType}{resetColor}\n')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)



    @staticmethod
    def compairDF(data1, name1, data2, name2):
        print('=============================== Compare Datasets '
              '================================')
        print(f'Datasets:\n'
              f'     {greenLight}{name1}\n'
              f'     {name2}{resetColor}\n')
        matchID = list(data2.index)
        missingID, missingIDCount = [], 0
        for index, donorID in enumerate(data1.index):
            if donorID not in matchID:
                missingIDCount += 1
        if missingIDCount == 0:
            print(f'All Donor IDs match\n\n')
        else:
            print(f'Missing Donor IDs: {missingIDCount}\n{missingID}\n\n')



    def loadData(self, fileName):
        data = None
        fileLocation = os.path.join(self.pathFolder, fileName)
        if not os.path.exists(fileLocation):
            print(f'\n{orange}ERROR: File not found\n'
                  f'     {cyan}{fileLocation}\n')
            sys.exit(1)

        # Evaluate: File extension
        fileExt = fileName[-7:]
        if '.csv' in fileExt:
            print('================================= Load CSV File '
                  '=================================')
            print(f'Loading File: {greenLight}{fileName}\n'
                  f'    {greenDark}{fileLocation}{resetColor}')
            data = pd.read_csv(fileLocation, index_col=1)
            if 'Unnamed: 0' in data.columns:
                data = data.drop('Unnamed: 0', axis=1)
        elif '.fastq' in fileExt:
            print('================================ Load Fastq File '
                  '================================')
            print(f'Loading File: {greenLight}{fileName}\n'
                  f'    {greenDark}{fileLocation}{resetColor}')
            sys.exit(1)
        elif '.xlsx' in fileExt:
            print('================================ Load Execl File '
                  '================================')
            print(f'Loading File: {greenLight}{fileName}\n'
                  f'    {greenDark}{fileLocation}{resetColor}')
            if os.path.exists(fileLocation):
                data = pd.read_excel(fileLocation, index_col=0)
        else:
            print(f'\n{orange}ERROR: Unknown File Extension: {cyan}{fileName}\n')
            sys.exit(1)

        print(f'\nTotal Patients: {red}{len(data.iloc[:, 0])}{resetColor}\n')
        if self.printData:
            print(f'Loaded Data: {greenLight}{fileName}{resetColor}\n'
                  f'{data}\n\n')
        else:
            print()

        return data



    def plotBarGraph(self, data, dataType, barColors='#CC5500',
                     barWidth=0.75, cutoff=False):
        if self.NBars is not None:
            data = len(data.iloc[:, 0:self.NBars])
        title = ''
        figSize = self.figSize
        yValues = []
        if '% at8' in dataType.lower():
            if cutoff:
                title = f'Localized AT8 Distribution'
            else:
                title = 'AT8 Distribution'
            figSize = self.figSizeW
            edgeWidth = 0

            # Evaluate: Y axis
            maxVal = data.max().max()
            maxValue = math.ceil(maxVal)
            magnitude = math.floor(math.log10(maxValue))
            unit = 10**(magnitude-1)
            yMax = math.ceil(maxValue / unit) * unit
            if yMax < maxVal:
                increaseValue = unit / 2
                while yMax < max(yValues):
                    yMax += increaseValue
            yMin = 0
        else:
            print(f'Unknown data type for bar graph: {dataType}')
            sys.exit(1)


        # Params: x ticks
        xLabels = data.index
        spacingFactor = 1.5  # or try 2.0, 2.5, etc.
        xTicks = np.arange(len(xLabels)) * spacingFactor
        totalGroupWidth = (len(data.columns.tolist()) - 1) * barWidth
        xMin = xTicks[0] - totalGroupWidth / 2
        xMax = xTicks[-1] + totalGroupWidth / 2

        # Plot the figure
        columns = data.columns.tolist()
        nColumns = len(columns)  # Number of layers = 5
        fig, ax = plt.subplots(figsize=figSize)
        for index, column in enumerate(columns):
            # Plot each layer with an offset
            offsets = xTicks + (index - nColumns / 2) * barWidth + barWidth / 2
            ax.bar(offsets, data[column], width=barWidth,
                   label=column.replace('% AT8 positive ', ''),
                   color=barColors[index % len(barColors)],
                   edgecolor='black', linewidth=edgeWidth)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        ax.set_ylabel(dataType, fontsize=self.labelSizeAxis)
        ax.axhline(y=0, color='black', linewidth=self.lineThickness)
        ax.legend(title='Layer', fontsize=8, title_fontsize=9)
        ax.set_xlim(xMin, xMax)
        ax.set_ylim(yMin, yMax)

        # Set: xticks
        padding = barWidth * 2
        ax.set_xlim(xMin - padding, xMax + padding)
        ax.set_xticks(xTicks)
        ax.set_xticklabels([str(label) for label in xLabels], rotation=90, fontsize=8)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        plt.xticks(rotation=90, ha='center')

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        fig.canvas.mpl_connect('key_press_event', pressKey)
        fig.tight_layout()
        plt.show()



    def processAT8(self, data, name, header, divisorHeader):
        print('============================== Evaluating Datasets '
              '==============================')
        print(f'Datasets:\n'
              f'     {greenLight}{name}{resetColor}\n\n'
              f'Evaluating Headers: {purple}{header}{resetColor}')

        # Determent Relevant Columns
        columns = []
        columnsNew = []
        indexColumns = []
        for index, column in enumerate(data.columns):
            if header in column:
                print(index, column)
                if divisorHeader in column:
                    divisorHeader = column
                    print(f'Divisor: {pink}{divisorHeader}{resetColor}')
                else:
                    columns.append(column)
                    columnsNew.append(f'% AT8 positive {column.split("_")[1]}')
                    indexColumns.append(index)
                    print(f'     Column: {pink}{column}{resetColor}')
        print('\n')

        # Initialize DFs
        colStart, colEnd = indexColumns[0], indexColumns[-1] + 1
        self.perAT8 = pd.DataFrame(data.loc[:, data.columns[colStart]:
                                               data.columns[colEnd]],
                                   columns=columnsNew,
                                   index=list(data.index))
        self.perAT8High = pd.DataFrame(data.loc[:, data.columns[colStart]:
                                                   data.columns[colEnd]],
                                       columns=columnsNew, index=[])

        # Evaluate Data
        for index, donorID in enumerate(data.index):
            totalSignal = data.iloc[index, colStart:colEnd]
            totalSignal.astype(float) # Convert to floats
            totalSignal = totalSignal.sum()
            for indexCol, column in enumerate(self.perAT8.columns):
                self.perAT8.loc[donorID, column] = (
                        (data.loc[donorID, columns[indexCol]] / totalSignal) * 100)
            
            # Collect samples with localized AT8 distributions
            maxVal = self.perAT8.iloc[index, :].max()
            if maxVal > self.perAT8Cutoff:
                self.patientsOfInterest.append(donorID)
                self.perAT8High.loc[donorID] = self.perAT8.iloc[index, :]
        print(f'Percent AT8 positive in each layer:\n{self.perAT8}\n\n')
        if self.patientsOfInterest:
            print('Donors localized AT8:')
            for donorID in self.patientsOfInterest:
                print(f'{purple}{donorID}{resetColor}:\n'
                      f'{self.perAT8High.loc[donorID, :]}\n')
            self.dataPOI.index = self.patientsOfInterest
        print(f'{resetColor}')

        # Plot the data
        barColors = ['#FF8800', '#2E9418', '#56F1FF', '#FF0080', '#A800FF']
        if self.plotAT8:
            self.plotBarGraph(data=self.perAT8, dataType='% AT8',
                              barColors=barColors, barWidth=0.2)
        if self.plotAT8Cutoff:
            self.plotBarGraph(data=self.perAT8High, dataType='% AT8',
                              barColors=barColors, barWidth=0.2, cutoff=True)

        return self.patientsOfInterest



    # def DOI(self, data):
    #     print('============================== Donors Of Interest '
    #           '===============================')
    #     print(f'Selected donors: {pink}AT8{resetColor} > '
    #           f'{red}{self.perAT8Cutoff} %{resetColor}')
    #     numColumns = len(data.columns) - 1
    #     print(f'Num Columns: {red}{numColumns}{resetColor}')
    #     selectColumns = 5
    #     print(f'Find step size with even division for consistent '
    #           f'numbers of columns: {red}{numColumns}{resetColor}')
    #
    #     for donorID in self.patientsOfInterest:
    #         print(f'     Donor ID: {purple}{donorID}{resetColor}\n')
    #         for index in range(0, numColumns, selectColumns + 1):
    #             if index >= numColumns:
    #                 break
    #             indexLastCol = index + selectColumns
    #             print(f'Index: {red}{index}{resetColor}\n'
    #                   f'     End: {numColumns}\n'
    #                   f'     Last: {indexLastCol}\n')
    #             colStart = data.columns[index]
    #             if indexLastCol > numColumns:
    #                 colEnd = data.columns[numColumns]
    #             else:
    #                 colEnd = data.columns[index + selectColumns]
    #
    #             print(f'Select: {purple}{colStart}{resetColor} - '
    #                   f'{purple}{colEnd}{resetColor}')
    #             print(f'{data.loc[donorID, colStart:colEnd]}\n')
    #             print(f'{red}{index} - {index + selectColumns}{resetColor}\n\n')
    #             time.sleep(0)
    #         sys.exit()
    #     print()


    def DOI(self, data):
        print('============================== Donors Of Interest '
              '===============================')
        print(f'Selected donors: {pink}AT8{resetColor} > '
              f'{red}{self.perAT8Cutoff} %{resetColor}')
        # dataPOI = pd.DataFrame(None, index=self.patientsOfInterest)
        for donorID in self.patientsOfInterest:
            print(f'     Donor ID: {purple}{donorID}{resetColor}')
            self.dataPOI.loc[donorID, :] += data.loc[donorID, :]
        print('\n')

        dataTag = f'Localized AT8 > {self.perAT8Cutoff} %'
        self.getMetadata(data, dataTag=dataTag)

    def getMetadata(self, data, dataTag):
        print('=================================== Metadata '
              '====================================')
        # Select data for POI
        dataFields = ['Highest level of education', 'APOE Genotype', 'Cognitive Status',
                      'Age of onset cognitive symptoms', 'Age of Dementia diagnosis',
                      'Last CASI Score', 'Interval from last CASI in months',
                      'Last MMSE Score', 'Interval from last MMSE in months',
                      'Last MOCA Score', 'Interval from last MOCA in months',
                      'PMI', 'Brain pH', 'Overall AD neuropathological Change',
                      'Thal', 'Braak', 'CERAD score', 'Overall CAA Score',
                      'Highest Lewy Body Disease', 'Atherosclerosis',
                      'Arteriolosclerosis', 'LATE', 'RIN']
        for donorID in self.dataPOI.index:
            for field in dataFields:
                self.dataPOI.loc[donorID, field] = data.loc[donorID, field]


        # Evaluate metadata
        metaData = {}
        for field in self.dataPOI.columns:
            datapoints = {}
            metaData[field] = {}
            for donorID in self.dataPOI.index:
                datapoint = self.dataPOI.loc[donorID, field]
                if isinstance(datapoint, (int, float)) and np.isnan(datapoint):
                    datapoint = 'NaN'
                if datapoint in datapoints.keys():
                    datapoints[datapoint] += 1
                else:
                    datapoints[datapoint] = 1

            # Sort dictionaries
            try:
                # Numerical sort
                datapoints = dict(sorted(datapoints.items(), key=lambda x: float(
                    x) if x == x and x is not None else float('inf')))
            except:
                # Alphabetical sort
                datapoints = dict(sorted(datapoints.items(), key=lambda x: str(x)))
            metaData[field] = datapoints


        # Print data
        for field in list(metaData.keys()):
            print(f'Category: {greenLight}{field}{resetColor}')
            values = metaData[field]
            for key, value in values.items():
                print(f'     {pink}{key}{resetColor}, Count: {red}{value}{resetColor}')
            print()
        print()

        # Plot the data
        self.plotMetadata(data=metaData, dataTag=dataTag)



    def plotMetadata(self, data, dataTag):
        barColors = plt.cm.Accent.colors
        title = dataTag

        # Get: Dataset fields
        fields = list(data.keys())

        # Get: Max value
        xMax = 0
        labelsBar = []
        for field in fields:
            values = data[field]
            for key, value in values.items():
                labelsBar.append(key)
                if value > xMax:
                    xMax = value
        xMax += 1


        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSizeTall)

        yTicks = []
        yLabels = []
        maxCategories = max(len(data[field]) for field in fields)
        barHeight = 0.9 / maxCategories  # Shrink bar height to fit all bars
        for indexYCluster, field in enumerate(fields):
            values = data[field]
            for j, (label, count) in enumerate(values.items()):
                yPos = indexYCluster + j * barHeight - (barHeight * len(values)) / 2
                color = barColors[j % len(barColors)]
                ax.barh(yPos, count, height=barHeight, color=color)

                # Add count label to the end of each bar
                ax.text(count + 0.1, yPos, f'{label}', va='center', color='black',
                        fontsize=8, fontweight='bold')
            yTicks.append(indexYCluster)
            yLabels.append(field)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        ax.set_xlabel('Counts', fontsize=self.labelSizeAxis)

        # Set: xAxis
        ax.set_xlim(0, xMax)

        # Set: yAxis
        ax.set_ylim(-1, len(yTicks))
        ax.set_yticks(yTicks)
        ax.set_yticklabels([str(label) for label in yLabels])

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        ax.grid(axis='x', linestyle='--', color='black')
        fig.canvas.mpl_connect('key_press_event', pressKey)
        plt.tight_layout()
        plt.show()
