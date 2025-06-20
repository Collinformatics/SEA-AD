import math
from unicodedata import category

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
import pandas as pd
import sys



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
                 plotMetadata, plotBiomarkers, printData=False, NBars=None):
        # Parameters: Files
        self.pathFolder = pathFolder
        self.files = files

        # Parameters: Datasets
        self.biomarkers = None
        self.biomarkersPOI = {}
        self.metadata = None
        self.metadataPOI = {}
        self.neuropathy = None
        self.volume = None

        # Parameters: Data
        self.biomarkerExtractions = None
        self.dataPOI = pd.DataFrame()
        self.dataTag = None
        self.dementiaPrevalencePOI = None
        self.numPatients = 0
        self.patientsOfInterest = [] # Selected donors with localized AT8
        self.perAT8 = None
        self.perAT8Select = None

        # Parameters: Figures
        self.plotAT8 = plotAT8
        self.plotAT8Cutoff = plotAT8Cutoff
        self.plotMetadata = plotMetadata
        self.plotBiomarkers = plotBiomarkers
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
        self.roundDecimal = 3



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
            print(f'{orange}ERROR: Cannot Create Colormap. '
                  f'Unrecognized colorType parameter: {colorType}{resetColor}\n')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)



    def compairDatasetDonors(self):
        print('=============================== Compare Datasets '
              '================================')
        for index in range(1, len(self.files)):
            print(f'Comparing Datasets:{greenLight}\n'
                  f'     {self.files[0]}\n'
                  f'     {self.files[index]}{resetColor}')
            if index == 1:
                matchData = self.metadata
            elif index == 2:
                matchData = self.biomarkers
            else:
                print(f'\n{orange}ERROR: Write A Script To Check The {cyan}Donor IDs'
                      f'{orange} For\n'
                      f'    {cyan}{self.files[index]}\n')
                sys.exit(1)

            # Compair: Donor IDs
            matchID = list(matchData.index)
            missingID, missingIDCount = [], 0
            for index, donorID in enumerate(self.neuropathy.index):
                if donorID not in matchID:
                    missingIDCount += 1
                    missingID.append(donorID)
            if missingIDCount == 0:
                print(f'All Donor IDs match\n')
            else:
                print(f'Missing Donor IDs: {missingIDCount}\n{missingID}\n')
                sys.exit(1)
        print()



    def loadData(self):
        for fileName in self.files:
            numRows = 0
            # Inspect file
            fileLocation = os.path.join(self.pathFolder, fileName)
            if not os.path.exists(fileLocation):
                print(f'\n{orange}ERROR: File Not Found\n'
                      f'     {cyan}{fileLocation}\n')
                sys.exit(1)

            # Evaluate: File extension
            fileExt = fileName[-7:]
            if '.csv' in fileExt:
                # Load: CSV
                print('================================= Load CSV File '
                      '=================================')
                print(f'Loading File: {greenLight}{fileName}\n'
                      f'    {greenDark}{fileLocation}{resetColor}\n')
                if fileName == 'sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv':
                    self.neuropathy = pd.read_csv(fileLocation, index_col=1)
                    if 'Unnamed: 0' in self.neuropathy.columns:
                        self.neuropathy = self.neuropathy.drop('Unnamed: 0', axis=1)
                    numRows = len(self.neuropathy.index)
                    print(f'Dataset: {pink}Neuropathy{resetColor}')
                else:
                    print(f'\n{orange}ERROR: Add Code To Load File\n'
                          f'     {cyan}{fileName}\n')
                    sys.exit(1)

            # Load: Excel
            elif '.xlsx' in fileExt:
                print('================================ Load Execl File '
                      '================================')
                print(f'Loading File: {greenLight}{fileName}\n'
                      f'    {greenDark}{fileLocation}{resetColor}\n')
                if os.path.exists(fileLocation):
                    if fileName == 'sea-ad_cohort_donor_metadata_072524.xlsx':
                        self.metadata = pd.read_excel(fileLocation, index_col=0)

                        # Drop rows where the index is NaN
                        self.metadata = self.metadata[self.metadata.index.notna()]
                        numRows = len(self.metadata.index)
                        print(f'Dataset: {pink}Metadata{resetColor}')
                    elif (fileName ==
                          'sea-ad_cohort_mtg-tissue_extractions-luminex_data.xlsx'):
                        self.biomarkers = pd.read_excel(fileLocation, header=[0, 1],
                                                        index_col=0)
                        self.biomarkerExtractions = (
                            self.biomarkers.columns.get_level_values(0)).tolist()
                        self.biomarkers.columns = (
                            self.biomarkers.columns.get_level_values(1))
                        numRows = len(self.biomarkers.index)
                        print(f'Dataset: {pink}Biomarkers{resetColor}')
                    elif fileName == 'sea-ad_cohort_mri_volumetrics.xlsx':
                        self.volume = pd.read_excel(fileLocation, index_col=0)
                    else:
                        print(f'\n{orange}ERROR: Add Code To Load File\n'
                              f'     {cyan}{fileName}\n')
                        sys.exit(1)
            else:
                print(f'{orange}ERROR: Unknown File Extension {cyan}{fileName}\n')
                sys.exit(1)
            print(f'Total Rows: {red}{numRows}{resetColor}\n\n')

        self.compairDatasetDonors()



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



    def processNeuropathy(self, header):
        print('========================== Evaluating Neuropathy Data '
              '===========================')
        print(f'Dataset: {pink}Neuropathy{resetColor}\n'
              f'Extracting Data: {purple}{header}{resetColor}')

        # Determent Relevant Columns
        columns = []
        columnsNew = []
        indexColumns = []
        for index, column in enumerate(self.neuropathy.columns):
            if 'grey matter' in column.lower():
                continue
            if header in column:
                columns.append(column)
                columnsNew.append(f'% AT8 positive {column.split("_")[1]}')
                indexColumns.append(index)
        print('\n')


        # Initialize DFs
        colStart = self.neuropathy.columns[indexColumns[0]]
        colEnd = self.neuropathy.columns[indexColumns[-1]]
        self.perAT8 = self.neuropathy.loc[:, colStart:colEnd]
        self.perAT8Select = pd.DataFrame(self.neuropathy.loc[:, colStart:colEnd],
                                         columns=self.perAT8.columns, index=[])
        print(f'AT8 Levels: {pink}Neuropathy{resetColor}\n{self.perAT8}\n\n')
        print(f'AT8 Distribution: {pink}Percent AT8 signal in each layer{resetColor}\n'
              f'{self.perAT8}\n\n')


        def evaluateSignal():
            # Evaluate AT8 levels
            if self.selectionType == '>':
                if maxVal > self.perAT8Cutoff:
                    self.perAT8Select.loc[donorID, :] = self.perAT8.loc[donorID, :]
            elif self.selectionType == '<':
                if maxVal < self.perAT8Cutoff:
                    self.perAT8Select.loc[donorID, :] = self.perAT8.iloc[index, :]
            elif self.selectionType == 'Range':
                if self.perAT8Cutoff[0] > maxVal > self.perAT8Cutoff[1]:
                    self.perAT8Select.loc[donorID, :] = self.perAT8.iloc[index, :]
            else:
                print(f'{orange}ERROR: I have No Use For This selectionType {cyan}'
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
            signal = self.neuropathy.loc[donorID, colStart:colEnd]
            signal.astype(float) # Convert to floats
            totalSignal = signal.sum()
            for column in columns:
                self.perAT8.loc[donorID, column] = (
                        self.neuropathy.loc[donorID, column] / totalSignal)
                self.perAT8.loc[donorID, column] *= 100
            maxVal = self.perAT8.iloc[index, :].max()
            evaluateSignal() # Collect samples based on AT8 distribution
        self.perAT8Select = self.perAT8Select.sort_index()
        self.patientsOfInterest = self.perAT8Select.index.tolist()
        self.dataPOI.index = self.patientsOfInterest
        self.numPatients = len(self.patientsOfInterest)


        # Plot the data
        barColors = ['#FF8800', '#2E9418', '#56F1FF', '#FF0080', '#A800FF']
        if self.plotAT8:
            self.plotBarGraph(data=self.perAT8, dataType='% AT8',
                              barColors=barColors, barWidth=0.2)
        if self.plotAT8Cutoff:
            self.plotBarGraph(data=self.perAT8Select, dataType='% AT8',
                              barColors=barColors, barWidth=0.2, cutoff=True)

        # Scan Donors Of Interest
        self.DOI()



    def DOI(self):
        print('============================== Donors Of Interest '
              '===============================')
        print(f'Selected donors: {pink}AT8{resetColor} > '
              f'{red}{self.perAT8Cutoff} %{resetColor}\n')
        print(f'Patients Of Interest: {pink}AT8 Distribution{resetColor}\n'
              f'{self.perAT8Select}\n\n'
              f'Total Patients Of Interest: {red}{self.numPatients}{resetColor}\n\n')
        if self.numPatients == 0:
            print(f'No donor were selected')
            sys.exit()

        # Define: Selection type
        if self.selectionType == 'Range':
            self.dataTag = (f'{self.perAT8Cutoff[0]} %  > Maximum AT8 Signal > '
                       f'{self.perAT8Cutoff[1]} %')
        else:
            self.dataTag = (f'Maximum AT8 Signal {self.selectionType} '
                            f'{self.perAT8Cutoff} %')

        # Process data
        self.getMetadata()
        self.processBiomarkers()



    def inspectMetadataPOI(self, category):
        print('=============================== Inspect Metadata '
              '================================')
        print(f'Category: {pink}{category}{resetColor}')
        if not self.patientsOfInterest:
            return

        for donorID in self.patientsOfInterest:
            if self.metadata.loc[donorID, category] == 'Dementia':
                print(f'Donor: {purple}{donorID} '
                      f'{red}{self.metadata.loc[donorID, category]}{resetColor} ')
            else:
                print(f'Donor: {donorID} {self.metadata.loc[donorID, category]}')
        print('\n')



    def plotDictionary(self, data, datasetType, inspectMetadata=True):
        barColors = plt.cm.Accent.colors
        title = (f'{datasetType}\n'
                 f'N = {self.numPatients} Patients\n'
                 f'{self.dataTag}\n'
                 f'Dementia Prevalence: {self.dementiaPrevalencePOI} %')
        datasetType = datasetType.lower()

        if inspectMetadata:
            # Evaluate: Cognitive status of POI
            self.inspectMetadataPOI(category='Cognitive Status')

        # Get: Dataset fields
        fields = list(data.keys())

        # Get: Max value
        if 'biomarkers' in datasetType or 'biomarker' in datasetType:
            xMax = 0
            labelsBar = []
            for field in fields:
                values = data[field]
                for key, value in values.items():
                    labelsBar.append(key)
                    if value > xMax:
                        xMax = value
            magnitude = np.floor((np.log10(xMax)) - 1)
            xMax += 10**magnitude
        elif datasetType == 'metadata':
            xMax = 0
            labelsBar = []
            for field in fields:
                values = data[field]
                for key, value in values.items():
                    labelsBar.append(key)
                    if value > xMax:
                        xMax = value
            xMax += 1
        else:
            print(f'{orange}ERROR: Add Code To Plot Dataset Type\n'
                  f'     {cyan}{datasetType}\n')
            sys.exit(1)


        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSizeLarge)

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



    def getMetadata(self):
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
                self.dataPOI.loc[donorID, field] = self.metadata.loc[donorID, field]

        # Evaluate metadata
        for field in self.dataPOI.columns:
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
        self.dementiaPrevalencePOI = (dementiaCounts / self.numPatients) * 100
        self.dementiaPrevalencePOI = round(self.dementiaPrevalencePOI, self.roundDecimal)
        print(f'Dementia cases in subset: {purple}{self.dataTag}{resetColor}\n'
              f'     Total Patients: {red}{self.numPatients}{resetColor}\n'
              f'     Dementia Cases: {red}{dementiaCounts}{resetColor}\n'
              f'     Prevalence: {red}{self.dementiaPrevalencePOI} %'
              f'{resetColor}\n\n')

        # Plot the data
        if self.plotMetadata:
            self.plotDictionary(data=self.metadataPOI, datasetType='Metadata')



    def processBiomarkers(self):
        print('=========================== Evaluating Biomarker Data '
              '===========================')
        methods = []
        indices = []
        for index, extractionMethod in enumerate(self.biomarkerExtractions):
            if extractionMethod not in methods:
                methods.append(extractionMethod)
                indices.append(index)
        print(f'Biomarker Extractions: {pink}{methods[0]}{resetColor}, {pink}{methods[1]}'
              f'{resetColor}\n{self.biomarkers}\n\n')

        # Get POI biomarkers
        biomarkerLevels = pd.DataFrame(index=self.patientsOfInterest,
                                       columns=self.biomarkers.columns)
        extractions, extractionA, extractionB = [], {}, {}
        for donorID in self.patientsOfInterest:
            # Split the biomarkers by the extraction type
            index = self.biomarkers.index.get_loc(donorID)
            extractionA[donorID] = self.biomarkers.iloc[index, indices[0]:indices[1]]
            extractionB[donorID] = self.biomarkers.iloc[index, indices[1]:]

            biomarkerLevels.loc[donorID, :] = self.biomarkers.loc[donorID, :]
            # self.biomarkersPOI[donorID] = self.biomarkers.loc[donorID, :]

        print(f'Biomarker Levels: {pink}Patients Of Interest{resetColor}\n'
              f'{biomarkerLevels}\n\n')

        for extractionMethod in methods:
            extractionMethod = extractionMethod.replace('extractions', 'Extractions')
            if extractionMethod == 'RIPA Buffer Tissue Extractions':
                self.biomarkersPOI[extractionMethod] = extractionA
            elif (extractionMethod ==
                  'GuHCl (Guanidine Hydrochloride) Buffer Tissue Extractions'):
                self.biomarkersPOI[extractionMethod] = extractionB
            else:
                print(f'{orange}ERROR: Add Code To Plot Extraction Type\n'
                      f'     {cyan}{extractionMethod}\n')


        if self.plotBiomarkers:
            firstFig = True
            for extractionMethod, values in self.biomarkersPOI.items():
                self.plotDictionary(data=values,
                                    datasetType=f'Biomarkers\n{extractionMethod}',
                                    inspectMetadata=firstFig)
                firstFig = False
