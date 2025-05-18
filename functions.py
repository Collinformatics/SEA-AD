from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
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
    def __init__(self, pathFolder, fileName):
        # Define: Variables
        self.counts = pd.DataFrame()

        # Define: Files
        self.pathFolder = pathFolder
        self.fileName = fileName
        self.filePath = os.path.join(self.pathFolder, self.fileName)
        self.loadData()



    def loadData(self):
        # Evaluate: File extension
        fileExt = self.fileName[-5:]
        if '.csv' in fileExt:
            self.loadCSV()
        elif '.xlsx' in fileExt:
            self.loadExcel()



    def loadCSV(self):
        print('================================= Load CSV File '
              '=================================')
        # Load data
        print(f'Loading File:\n'
              f'    {greenDark}{self.filePath}{resetColor}\n\n')
        if os.path.exists(self.filePath):
            data = pd.read_csv(self.filePath, index_col=0)
        else:
            print(f'\n{orange}ERROR: File not found\n'
                  f'     {cyan}{self.filePath}\n')
            sys.exit(1)
        print(f'Loaded Data: {greenLight}{self.fileName}{resetColor}\n'
              f'{data}\n\n')


    def loadExcel(self):
        print('================================ Load Execl File '
              '================================')
        # Load data
        print(f'Loading File:\n'
              f'    {greenDark}{self.filePath}{resetColor}\n\n')
        if os.path.exists(self.filePath):
            data = pd.read_excel(self.filePath)
        else:
            print(f'\n{orange}ERROR: File not found\n'
                  f'     {cyan}{self.filePath}\n')
            sys.exit(1)
        print(f'Loaded Data: {greenLight}{self.fileName}{resetColor}\n'
              f'{data}\n\n')

