# MTG: Middle temporal gyrus
    # Higher level cognitive function: Memory, verbal memory, visual processing
    # Buildup of Tau in preclinical AD, increases as the disease progresses


from functions import BrainData
import pandas as pd
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inPathFolder = '/Users/ca34522/Documents/Bioinformatics/SEA-AD/'
inFiles = ['sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv',
           'sea-ad_cohort_donor_metadata_072524.xlsx']
inLoadFile = inFiles[0]

# Input 2: Data Inspection
inPrintLoadedData = False


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



# ===================================== Import Data ======================================
# Initialize class
brains = BrainData(pathFolder=inPathFolder, printData=inPrintLoadedData)

# Load Data
quantNeuropathy = brains.loadData(fileName=inFiles[0])
metaData = brains.loadData(fileName=inFiles[1])



# ==================================== Evaluate Data =====================================
print('=============================== Compare Datasets ================================')
matchID = list(metaData.loc[:, 'Donor ID'])
missingID, missingIDCount  = [], 0
for index, donorID in enumerate(quantNeuropathy.loc[:, 'Donor ID']):
    if donorID not in matchID:
        missingIDCount += 1
if missingIDCount == 0:
    print(f'All Donor IDs match\n')
else:
    print(f'Missing Donor IDs: {missingIDCount}\n{missingID}\n')



