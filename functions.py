import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import time



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
    elif event.key == 'e':
        sys.exit()
    elif event.key == 'r':
        python = sys.executable
        os.execl(python, python, *sys.argv)



class Brains:
    def __init__(self, pathFolder, files, perAT8Cutoff, plotAT8, plotAT8Cutoff,
                 printData=False, NBars=None):
        # Parameters: Files
        self.pathFolder = pathFolder
        self.files = files

        # Parameters: Datasets
        self.neuropathy = None
        self.metadata = None
        self.biomarkers = None

        # Parameters: Data
        self.biomarkerExtractions = None
        self.patientsOfInterest = [] # Selected donors with localized AT8
        self.numPatients = 0
        self.dataPOI = pd.DataFrame()
        self.perAT8 = None
        self.perAT8Select = None
        self.metadataPOI = {}
        self.biomarkersPOI = {}

        # Parameters: Figures
        self.plotAT8 = plotAT8
        self.plotAT8Cutoff = plotAT8Cutoff
        self.figSize = (9.5, 8)
        self.figSizeLarge = (16, 9.5)
        self.figSizeW = (16, 8)
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.NBars = NBars

        # Parameters: Variables
        self.perAT8Cutoff = perAT8Cutoff
        self.selectionType = ''
        self.printData = printData

        # Parameters: Miscellaneous
        self.roundDeci = 3



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



    def compairDF(self):
        print('=============================== Compare Datasets '
              '================================')
        fileBase = self.files[0]
        print(f'Comparing Datasets:{greenLight}')
        for file in self.files:
            print(f'     {file}')
        print(f'{resetColor}\n')


        matchID = list(self.metadata.index)
        missingID, missingIDCount = [], 0
        for index, donorID in enumerate(self.neuropathy.index):
            if donorID not in matchID:
                missingIDCount += 1
                missingID.append(donorID)
        if missingIDCount == 0:
            print(f'All Donor IDs match\n\n')
        else:
            print(f'Missing Donor IDs: {missingIDCount}\n{missingID}\n\n')

        sys.exit()



    def loadData(self):
        for fileName in self.files:
            numRows = 0

            # Inspect file
            fileLocation = os.path.join(self.pathFolder, fileName)
            if not os.path.exists(fileLocation):
                print(f'\n{orange}ERROR: File not found\n'
                      f'     {cyan}{fileLocation}\n')
                sys.exit(1)

            # Evaluate: File extension
            fileExt = fileName[-7:]
            if '.csv' in fileExt:
                # Load: CSV
                print('================================= Load CSV File '
                      '=================================')
                print(f'Loading File: {greenLight}{fileName}\n'
                      f'    {greenDark}{fileLocation}{resetColor}')
                if fileName == 'sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv':
                    self.neuropathy = pd.read_csv(fileLocation, index_col=1)
                    if 'Unnamed: 0' in self.neuropathy.columns:
                        self.neuropathy = self.neuropathy.drop('Unnamed: 0', axis=1)
                    numRows = len(self.neuropathy.index)
                else:
                    print(f'\n{orange}ERROR: File not found: {cyan}{fileName}\n')
                    sys.exit(1)

            # Load: Excel
            elif '.xlsx' in fileExt:
                print('================================ Load Execl File '
                      '================================')
                print(f'Loading File: {greenLight}{fileName}\n'
                      f'    {greenDark}{fileLocation}{resetColor}')
                if os.path.exists(fileLocation):
                    if fileName == 'sea-ad_cohort_donor_metadata_072524.xlsx':
                        self.metadata = pd.read_excel(fileLocation, index_col=0)

                        # Drop rows where the index is NaN
                        self.metadata = self.metadata[self.metadata.index.notna()]
                        numRows = len(self.metadata.index)
                    elif fileName == 'sea-ad_cohort_mtg-tissue_extractions-luminex_data.xlsx':
                        self.biomarkers = pd.read_excel(fileLocation, header=[0, 1],
                                                        index_col=0)
                        self.biomarkerExtractions = (
                            self.biomarkers.columns.get_level_values(0))
                        self.biomarkers.columns = (
                            self.biomarkers.columns.get_level_values(1))
                        numRows = len(self.biomarkers.index)
                    else:
                        print(f'\n{orange}ERROR: File not found: {cyan}{fileName}\n')
                        sys.exit(1)
            else:
                print(f'\n{orange}ERROR: Unknown File Extension: {cyan}{fileName}\n')
                sys.exit(1)
            print(f'\nTotal Rows: {red}{numRows}{resetColor}\n\n')



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



    def processAT8(self, name, header, divisorHeader):
        print('============================== Evaluating Datasets '
              '==============================')
        print(f'Datasets:\n'
              f'     {greenLight}{name}{resetColor}\n\n'
              f'Evaluating Headers: {purple}{header}{resetColor}')

        # Determent Relevant Columns
        columns = []
        columnsNew = []
        indexColumns = []
        for index, column in enumerate(self.neuropathy.columns):
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
        colStart = self.neuropathy.columns[indexColumns[0]]
        colEnd = self.neuropathy.columns[indexColumns[-1] + 1]
        self.perAT8 = self.neuropathy.loc[:, colStart:colEnd]
        print(f'AT8 Levels: {pink}Neuropathy{resetColor}\n{self.perAT8}\n\n')

        sys.exit()

        self.perAT8Select = pd.DataFrame(self.neuropathy.loc[:, colStart:colEnd],
                                         columns=columnsNew, index=[])
        # print(f'Neurpoathy:\n{self.neuropathy}\n\nAT8:\n{self.perAT8}\n\n')

        def evaluateAT8():
            # Evaluate AT8 levels
            if self.selectionType == '>':
                if maxVal > self.perAT8Cutoff:
                    self.patientsOfInterest.append(donorID)
                    self.perAT8Select.loc[donorID] = self.perAT8.iloc[index, :]
            elif self.selectionType == '<':
                if maxVal < self.perAT8Cutoff:
                    self.patientsOfInterest.append(donorID)
                    self.perAT8Select.loc[donorID] = self.perAT8.iloc[index, :]
            elif self.selectionType == 'Range':
                if self.perAT8Cutoff[0] > maxVal > self.perAT8Cutoff[1]:
                    self.patientsOfInterest.append(donorID)
                    self.perAT8Select.loc[donorID] = self.perAT8.iloc[index, :]
            else:
                print(f'{orange}ERROR: I have no use for this selectionType {cyan}'
                      f'{self.selectionType}\n')
                sys.exit()


        # Define: Selection type
        if (len(self.perAT8Cutoff)) == 1:
            self.perAT8Cutoff = self.perAT8Cutoff[0]
            if self.perAT8Cutoff > 0:
                self.selectionType = '>'
            else:
                self.perAT8Cutoff = abs(self.perAT8Cutoff)
                self.selectionType = '<'
        else:
            self.perAT8Cutoff = sorted(self.perAT8Cutoff, reverse=True)
            self.selectionType = 'Range'


        # Evaluate Data
        for index, donorID in enumerate(self.neuropathy.index):
            print(index, donorID)
            totalSignal = self.neuropathy.iloc[index, colStart:colEnd]
            totalSignal.astype(float) # Convert to floats
            totalSignal = totalSignal.sum()
            for indexCol, column in enumerate(self.perAT8.columns):
                self.perAT8.loc[donorID, column] = (
                        self.neuropathy.loc[donorID, columns[indexCol]] / totalSignal)
                self.perAT8.loc[donorID, column] *= 100
            # Collect samples with localized AT8 distributions
            maxVal = self.perAT8.iloc[index, :].max()
            evaluateAT8()
        print(f'Percent AT8 positive in each layer:\n{self.perAT8}\n\n')
        self.dataPOI.index = self.patientsOfInterest
        self.numPatients = len(self.patientsOfInterest)

        print(f'AT8: Select\n{self.perAT8Select}\n\n')


        # Plot the data
        barColors = ['#FF8800', '#2E9418', '#56F1FF', '#FF0080', '#A800FF']
        if self.plotAT8:
            self.plotBarGraph(data=self.perAT8, dataType='% AT8',
                              barColors=barColors, barWidth=0.2)
        if self.plotAT8Cutoff:
            self.plotBarGraph(data=self.perAT8Select, dataType='% AT8',
                              barColors=barColors, barWidth=0.2, cutoff=True)



    def DOI(self):
        print('============================== Donors Of Interest '
              '===============================')
        print(f'Selected donors: {pink}AT8{resetColor} > '
              f'{red}{self.perAT8Cutoff} %{resetColor}\n')
        print(f'Donors Of Interest: {pink}Localized AT8{resetColor}')
        for donorID in self.patientsOfInterest:
            print(f'{purple}{donorID}{resetColor}:\n'
                  f'{self.perAT8Select.loc[donorID, :]}\n')
        print(f'Number of selected donors: {red}{self.numPatients}{resetColor}\n\n')
        if self.numPatients == 0:
            print(f'No donor were selected')
            sys.exit()
        sys.exit()


        # Define: Selection type
        if self.selectionType == 'Range':
            dataTag = (f'{self.perAT8Cutoff[0]} %  > Maximum AT8 Signal > '
                       f'{self.perAT8Cutoff[1]} %')
        else:
            dataTag = f'Maximum AT8 Signal {self.selectionType} {self.perAT8Cutoff} %'

        # Process data
        self.getMetadata(dataTag=dataTag)
        self.biomarkers()



    def getMetadata(self, dataTag):
        print('=================================== Metadata '
              '====================================')
        # Select data for POI
        dataFields = ['Highest level of education', 'APOE Genotype', 'Cognitive Status',
                      'Age of onset cognitive symptoms', 'Age of Dementia diagnosis',
                      'Last CASI Score', 'Last MMSE Score', 'Last MOCA Score',
                      'PMI', 'Brain pH', 'Overall AD neuropathological Change',
                      'Thal', 'Braak', 'CERAD score', 'Overall CAA Score',
                      'Highest Lewy Body Disease', 'Atherosclerosis',
                      'Arteriolosclerosis', 'LATE', 'RIN']
        # 'Interval from last CASI in months', 'Interval from last MMSE in months',
        # 'Interval from last MOCA in months',

        # Get: Metadata
        for donorID in self.patientsOfInterest:
            for field in dataFields:
                print(f'Field: {field}\n\n'
                      f'{self.metadata}')
                self.dataPOI.loc[donorID, field] = self.metadata.loc[donorID, field]

        # Evaluate metadata
        for field in self.dataPOI.columns:
            print(field)
            datapoints = {}
            self.metadataPOI[field] = {}
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
            self.metadataPOI[field] = datapoints


        # Print data
        for field in list(self.metadataPOI.keys()):
            print(f'Category: {greenLight}{field}{resetColor}')
            values = self.metadataPOI[field]
            for key, value in values.items():
                print(f'     {pink}{key}{resetColor}, Count: {red}{value}{resetColor}')
            print()
        print()

        # Evaluate: Dementia stats
        dementiaCounts = 0
        for key, count in self.metadataPOI['Cognitive Status'].items():
            if 'dementia' in key.lower():
                dementiaCounts += count
                break
        dementiaPercent = (dementiaCounts / self.numPatients) * 100
        print(f'Dementia cases in subset: {purple}{dataTag}{resetColor}\n'
              f'     Total Patients: {red}{self.numPatients}{resetColor}\n'
              f'     Dementia Cases: {red}{dementiaCounts}{resetColor}\n'
              f'     Prevalence: {red}{round(dementiaPercent, self.roundDeci)} %'
              f'{resetColor}\n')

        # Plot the data
        self.plotMetadata(dataTag=dataTag, N=self.numPatients)



    def plotMetadata(self, dataTag, N):
        barColors = plt.cm.Accent.colors
        title = f'{dataTag}\nN = {N} Patients'

        # Get: Dataset fields
        fields = list(self.metadataPOI.keys())

        # Get: Max value
        xMax = 0
        labelsBar = []
        for field in fields:
            values = self.metadataPOI[field]
            for key, value in values.items():
                labelsBar.append(key)
                if value > xMax:
                    xMax = value
        xMax += 1


        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSizeLarge)

        yTicks = []
        yLabels = []
        maxCategories = max(len(self.metadataPOI[field]) for field in fields)
        barHeight = 0.9 / maxCategories  # Shrink bar height to fit all bars
        for indexYCluster, field in enumerate(fields):
            values = self.metadataPOI[field]
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
        ax.invert_yaxis()

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



    def processBioMarkers(self):
        print('============================== Evaluating Datasets '
              '==============================')

        # Get DOI biomarkers
        for donorID in self.patientsOfInterest:
            self.biomarkersDOI = 0

        print(f'Biomarkers:')


        sys.exit()
