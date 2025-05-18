from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import math
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
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
    elif event.key == 'backspace':
        sys.exit()



class BrainData:
    def __init__(self, pathFolder, printData, NBars=None):
        # Parameters: Files
        self.pathFolder = pathFolder

        # Parameters: Figures
        self.figSize = (9.5, 8)
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.NBars = NBars

        # Parameters: Variables
        self.printData = printData



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

        print(f'\nDatapoints: {red}{len(data.iloc[:, 0])}{resetColor}\n')
        if self.printData:
            print(f'Loaded Data: {greenLight}{fileName}{resetColor}\n'
                  f'{data}\n\n')
        else:
            print()

        return data



    def plotBarGraph(self, data, dataType, barColor='#CC5500', barWidth=0.75):
        print('================================ Plot: Bar Graph '
              '================================')
        if self.NBars is None:
            self.NBars = len(data.iloc[:, 0])
        xValues = []
        yValues = []

        # Collect substrates
        iteration = 0
        countsTotal = 0
        for count in substrates.values():
            countsTotal += count
        print(f'Total Substrates: {red}{countsTotal:,}{resetColor}')

        if 'counts' in dataType.lower():
            # Evaluate: Substrates
            for substrate, count in substrates.items():
                xValues.append(str(substrate))
                yValues.append(count)
                iteration += 1
                if iteration == self.NSubBars:
                    break

            # Evaluate: Y axis
            maxValue = math.ceil(max(yValues))
            magnitude = math.floor(math.log10(maxValue))
            unit = 10**(magnitude-1)
            yMax = math.ceil(maxValue / unit) * unit
            if yMax < max(yValues):
                increaseValue = unit / 2
                while yMax < max(yValues):
                    print(f'Increase yMax by:{yellow} {increaseValue}{resetColor}')
                    yMax += increaseValue
                print('\n')
            yMin = 0 # math.floor(min(yValues) / unit) * unit - spacer
        elif 'probability' in dataType.lower():
            # Evaluate: Substrates
            for substrate, count in substrates.items():
                xValues.append(str(substrate))
                yValues.append(count / countsTotal)
                iteration += 1
                if iteration == self.NSubBars:
                    break

            # Evaluate: Y axis
            maxValue = max(yValues)
            magnitude = math.floor(math.log10(maxValue))
            adjustedMax = maxValue * 10**abs(magnitude)
            yMax = math.ceil(adjustedMax) * 10**magnitude
            adjVal = 5 * 10**(magnitude-1)
            yMaxAdjusted = yMax - adjVal
            if yMaxAdjusted > maxValue:
                yMax = yMaxAdjusted
            yMin = 0
        else:
            # Evaluate: Substrates
            for substrate, count in substrates.items():
                xValues.append(str(substrate))
                yValues.append(count)
                iteration += 1
                if iteration == self.NSubBars:
                    break

            # Evaluate: Y axis
            spacer = 0.2
            yMax = math.ceil(max(yValues)) + spacer
            yMin = math.floor(min(yValues))
        NSubs = len(xValues)
        print(f'Number of plotted sequences: {red}{NSubs}{resetColor}\n\n')

        # Plot the data
        fig, ax = plt.subplots(figsize=self.figSize)
        bars = plt.bar(xValues, yValues, color=barColor, width=barWidth)
        plt.ylabel(dataType, fontsize=self.labelSizeAxis)
        plt.title(f'\n{self.enzymeName}\n{self.datasetTag}\n'
                 f'Top {NSubs} Substrates',
                  fontsize=self.labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=self.lineThickness)
        plt.ylim(yMin, yMax)
        # plt.subplots_adjust(top=0.873, bottom=0.12, left=0.101, right=0.979)


        # Set the edge color
        for bar in bars:
            bar.set_edgecolor('black')

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
