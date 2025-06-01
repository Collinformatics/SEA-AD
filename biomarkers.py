# MTG: Middle temporal gyrus
    # Higher level cognitive function: Memory, verbal memory, visual processing
    # Buildup of Tau in preclinical AD, increases as the disease progresses

# AT8: Anti-Phospho-Tau (Ser202, Thr205) Monoclonal Antibody


from functions import Brains
import pandas as pd
import sys



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inPathFolder = '/path/SEA-AD/'
inLoadFiles = ['sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv',
               'sea-ad_cohort_donor_metadata_072524.xlsx',
               'sea-ad_cohort_mtg-tissue_extractions-luminex_data.xlsx']

# Input 2: Data Inspection
inSelectDataColumns = 'total AT8 positive'
inAT8Cutoff = [52]


# Input 3: Plotting Data
inPlotAT8 = False
inPlotAT8Cutoff = False
inPlotMetadata = True
inPlotBiomarkers = True



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



# =================================== Initialize Class ===================================
brains = Brains(pathFolder=inPathFolder, files=inLoadFiles, perAT8Cutoff=inAT8Cutoff,
                plotAT8=inPlotAT8, plotAT8Cutoff=inPlotAT8Cutoff,
                plotMetadata=inPlotMetadata, plotBiomarkers=inPlotBiomarkers)



# ==================================== Evaluate Data =====================================
# Load Data
brains.loadData()

# Tau distributions
brains.processNeuropathy(header=inSelectDataColumns)
