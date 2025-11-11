# Detection-Limits
Code to calculate detection limits for instruments in Bergstrom Watershed Biogeochemistry Lab at Boise State University

The IC detection limit has 2 files.
1. read_drive_csv.R - a function script that allows IC data to be pulled directly from google drive, temporarily storing in the memory to create variable
2. IC_Error_Final - the main detection limit and error script for the IC 
    - IC sample data needs the Determination.Start, Info.1 , and all ion concentrations (also exporting Ident column is preferred)
    - IC standards run data needs determination, sample type, ion concentrations, and ion areas exported. 
    ** There are export templates already created: Anions, Cations, Anion Calibration, and Cation Calibration, using these templates will ensure that you have the necessary columns exported. 

The AA500 detection limit code has 2 files 
1. AA500_DL_STD.R - Regular code to run detection limit for the AA500. This will only run of sample names are fully numeric
2. AA500_DL_STD_StringID.R - the same code but works if the sample ID is a string and non-numeric

The archive folder contains code that was previously used to calculate detection limits but is no longer used.
1. Detection Limit.qmd - a function for calculating detection limits based on the export tempelate labeled HR_IC_etc. Run this function to add it to your environment.
2. IC_Detectionlimit_Functionintegrated.qmd - prompts the user to load their own calibration data and run the function.
3. IC_Detectionlimit.qmd - breaks the code in the function into individual chunks and explains the code and the decisions I(Hannah Richardson) made. Should be the more human-friendly version of function itself.

4. IC_Error_DL.R - the main script to run to get detection limits for IC data. This was turned into a quarto document (IC_Error_Final.qmd) for easier use. The quarto document was then modified further. 
