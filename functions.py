from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import gzip
import math
from itertools import combinations, product
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pandas as pd
import pickle as pk
import seaborn as sns
import sys
import threading



# ===================================== Set Options ======================================
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
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
class BrainData:
    def __init__(self, pathFolder, fileName):
        # Define: Variables
        self.counts = pd.DataFrame()

        # Define: Files
        self.pathFolder = pathFolder
        self.fileName = fileName
        self.filePath = []
        self.filePath.append(os.path.join(self.pathFolder, self.fileName))

        # Load Data
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
            print('')
            # data = pd.read_csv(filePath, index_col=0)
        else:
            print(f'\n{orange}ERROR: File not found\n'
                  f'     {cyan}{self.filePath}\n')
            sys.exit(1)
        print(f'{resetColor}')


    def loadExcel(self):
        print('================================ Load Execl File '
              '================================')
        # Load data
        print(f'Loading File:\n'
              f'    {greenDark}{self.filePath}{resetColor}\n\n')
        if os.path.exists(self.filePath):
            print('')
            # data = pd.read_csv(filePath, index_col=0)
        else:
            print(f'\n{orange}ERROR: File not found\n'
                  f'     {cyan}{self.filePath}\n')
            sys.exit(1)
        print(f'{resetColor}')

