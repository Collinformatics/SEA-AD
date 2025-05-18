from functions import BrainData


# User Inputs
inPathFolder = '/Users/ca34522/Documents/Bioinformatics/SEA-AD/'
inFiles = ['sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv',
           'sea-ad_cohort_donor_metadata_072524.xlsx']
inLoadFile = inFiles[0]

# Initialize class
brains = BrainData(pathFolder=inPathFolder, fileName=inLoadFile)
