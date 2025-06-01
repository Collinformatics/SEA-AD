# Install Modules:

    pip install numpy
    pip install pandas
    pip install openpyxl
    pip install seaborn


# Keyboard Shortcuts

- Esc: Close the current figure

- E: End the script that is currently running

- R: Rerun the script


# Selecting Patients of Interest (POI)

Distribution of AT8 Signal:

- Input:

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
