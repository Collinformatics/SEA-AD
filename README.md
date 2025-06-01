# Install Modules:

    pip install numpy
    pip install pandas
    pip install openpyxl
    pip install seaborn


# Keyboard Shortcuts

- Esc: Close the current figure

- E: End the script that is currently running

- R: Rerun the script

# Loading Datasets

- User Inputs:

  - Define: File Path
 
    Define the path to the folder with the SEA-AD datasets

        # Input 1: Select Dataset
        inPathFolder = '/path/SEA-AD/'

    - Define: Files
 
      Define the names of the SEA-AD datasets

          # Input 1: Select Dataset
          inLoadFiles = ['sea-ad_all_mtg_quant_neuropath_bydonorid_081122.csv',
                         'sea-ad_cohort_donor_metadata_072524.xlsx',
                         'sea-ad_cohort_mtg-tissue_extractions-luminex_data.xlsx']

# Selecting Patients of Interest (POI)

Distribution of AT8 Signal:

- User Inputs:

  - Define: Columns With AT8 Data

        # Input 2: Data Inspection
        inSelectDataColumns = 'total AT8 positive'

  - Define: Threshold For AT8 Signal For POI

        # Input 2: Data Inspection
        inAT8Cutoff = []
  
    This input will allow you to select patients of interest based upon the distributions of AT8 signal in layers: 1, 2, 3, 4, 5-6

      - To select donors with: Maximum % of AT8 Signal > 45%
    
              inAT8Cutoff = [45]
    
      - To select donors with: Maximum % of AT8 Signal < 45%
    
              inAT8Cutoff = [-45]
    
      - To select donors with: 45 %  > Maximum % of AT8 Signal > 40 %
    
              inAT8Cutoff = [45, 40]
