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
    def __init__(self, pathFolder):
        # Define: Variables
        self.counts = pd.DataFrame()

        # Define: Files
        self.pathFolder = pathFolder



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
            print(f'Loading File:\n'
                  f'    {greenDark}{fileLocation}{resetColor}\n\n')
            data = pd.read_csv(fileLocation, index_col=0)
        elif '.fastq' in fileExt:
            print('================================ Load Fastq File '
                  '================================')
            print(f'Loading File:\n'
                  f'    {greenDark}{fileLocation}{resetColor}\n\n')
            sys.exit(1)
        elif '.xlsx' in fileExt:
            print('================================ Load Execl File '
                  '================================')
            print(f'Loading File:\n'
                  f'    {greenDark}{fileLocation}{resetColor}\n\n')
            if os.path.exists(fileLocation):
                data = pd.read_excel(fileLocation)

        else:
            print(f'\n{orange}ERROR: Unknown File Extension: {cyan}{fileName}\n')
            sys.exit(1)

        print(f'Loaded Data: {greenLight}{fileName}{resetColor}\n'
              f'{data}\n\n')

        return data
