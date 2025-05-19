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

from matplotlib.pyplot import xlabel

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
    def __init__(self, pathFolder, printData, NBars=None):
        # Parameters: Files
        self.pathFolder = pathFolder

        # Parameters: Figures
        self.figSize = (9.5, 8)
        self.figSizeW = (16, 8)
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.NBars = NBars

        # Parameters: Variables
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
        matchID = list(data2.loc[:, 'Donor ID'])
        missingID, missingIDCount = [], 0
        for index, donorID in enumerate(data1.loc[:, 'Donor ID']):
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
            data = pd.read_csv(fileLocation, index_col=0)
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
                data = pd.read_excel(fileLocation)

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

        # Initialize DF
        colStart, colEnd = indexColumns[0], indexColumns[-1] + 1
        df = pd.DataFrame(data.loc[:, data.columns[colStart]:data.columns[colEnd]],
                          columns=columnsNew, index=list(data.loc[:, 'Donor ID']))

        # Evaluate Data
        sumAT8 = 0
        lowDiv = False
        for index, donorID in enumerate(data.loc[:, 'Donor ID']):
            totalSignal = data.iloc[index, colStart:colEnd]
            totalSignal.astype(float) # Convert to floats
            totalSignal = totalSignal.sum()
            for indexCol, column in enumerate(df.columns):
                sumAT8 += data.loc[index, columns[indexCol]]
                df.loc[donorID, column] = (data.loc[index, columns[indexCol]] /
                                           totalSignal) * 100
            if not lowDiv and sumAT8 > totalSignal:
                lowDiv = True
                print(f'{yellow}Warning: In at least 1 sample the {cyan}{header}{yellow} '
                      f'total signal from the layers is > the total signal used for '
                      f'the {cyan}Divisor{yellow}\n'
                      f'Sample: {pink}{donorID}\n'
                      f'     {cyan}{header}: {red}{sumAT8:,}\n'
                      f'     {cyan}Divisor: {red}{totalSignal:,}{resetColor}\n')
            sumAT8 = 0
        print(f'Percent AT8 positive in each layer:\n{df}\n\n')

        barColors = ['#FF8800', '#2E9418', '#56F1FF', '#FF0080', '#A800FF']
        self.plotBarGraph(data=df, dataType='% AT8', barColors=barColors, barWidth=0.2)

        return df



    def plotBarGraph(self, data, dataType, barColors='#CC5500', barWidth=0.75):
        print('================================ Plot: Bar Graph '
              '================================')
        if self.NBars is None:
            self.NBars = len(data.iloc[:, 0])
        else:
            data = len(data.iloc[:, 0:self.NBars])
        title = ''
        figSize = self.figSize
        yValues = []
        if '% at8' in dataType.lower():

            title = 'AT8 Distribution'
            figSize = self.figSizeW
            edgeWidth = 0
            print(f'Plotting: {purple}{title}{resetColor}\n'
                  f'{data}\n\n')

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
            print(f'Unknown Data Type: {dataType}')
            sys.exit(1)


        # Params: xticks
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
