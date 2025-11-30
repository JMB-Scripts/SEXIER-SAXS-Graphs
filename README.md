# Update 
Major update 7.x

3-Panel Layout**: Central control panel, 4-panels: Guinier analysis plot, Kratkt, Volume of correlation and a P(r) distribution plot.
-   **Internal BIFT for P(r) Analysis**: Pure Python implementation of the BIFT algorithm (Thx a lot BioXTAS RAW).
-   **Data Loading & Unit Conversion**: Automatically loads a .dat file and converts q-units from nm⁻¹ to Å⁻¹.
-   **Automated & Manual Guinier Analysis**: For determining Rg and I(0).
-   **Molecular Weight (MW) & Oligomeric State Estimation**: Based on the Volume of Correlation.
-   **Comprehensive Plotting & Automatic Saving**: Automatic saving of plots (PNG/SVG) and data files.


# SEXIER-SAXS-Graphs

Regain control of your SAXS data.

This Python script is designed to process SAXS (Small Angle X-ray Scattering) data and to generate results such as the Guinier approximation, the Kratky graph and the volume of correlation in a text file format.

Then you can make YOUR OWN representations for these graphs, with YOUR preferred software (Excel, Prism, Origin etc.). 

If you found this project useful, used it, or needed to customize it, please let me know! 

Your feedback is essential to help me improve and continue this project. You can reach out to me directly at [reach out to me via email](jean-marie.bourhis@univ-grenoble-alpes.fr).


## User manual

# Prerequisites
- Python 3.x
- Necessary packages: numpy, matplotlib, scipy, Qt5
- SAXS Data need to be in Å-1

## Command syntax

v4 introduce a gui where you can :

1. Browse to your `filename.dat` : the name of the .dat file containing the SAXS experimental data.
    
2. You can estimate qmin and qmax automatically (not as good as Raw or Primus), or enter the known values.
    
    - `qmin_offset` : the offset (in number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS or RAW.
    
    - `qmax_offset`: the offset (in number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS or RAW.
    
    3. Process the data 

<img width="1440" height="962" alt="image" src="https://github.com/user-attachments/assets/99940482-d1ac-489d-bc0b-d9d17d1b089a" />


## Features

 1. Guinier approximation :
    
 - Read .dat file and determine a first usable line.
   
 - Data extraction for q and I(q) in the selected range.
   
 - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
   
 - Display graph with experimental points and theoretical curve.

 2. Kratky 2:
    
 - Extract data in the selected range for Kratky.
   
 - Calculation and normalisation of values for Kratky (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
   
 - Display Kratky graph.

 3. Correlation volume (CV):
    
 - Extract data up to q=0.3 or up to 8/Rg.
   
 - Calculate the integral of the product I(q)*q:
   
	 - Calculation of VC, QR (quality factor) and MW (molecular weight).
    
	 - Calculation of the Mw from Porod using SAXSMOW approach with a cut at 8/Rg
    
 - If you enter the Mw estimated from the sequence an estimation (+/- 10%) of the oligomeric state will be perform 

<img width="283" height="154" alt="image" src="https://github.com/user-attachments/assets/bd6df925-fbd0-4c04-947d-c2cba1a42d44" />
    

 4. BIFT to estimate Dmax

   <img width="1440" height="962" alt="image" src="https://github.com/user-attachments/assets/ded85805-a03a-46fd-8d5b-5c93e43ec4e6" />



## Output

 The script will generate the following files:
 
 1- filename_01_Rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalized residuals)
 
 2- file_name_02_Norm-Kratky.txt: data for normalized Kratky (x, y)
 
 3- file_name_03_VC.txt: data for VC (q, I(q)*q)
 
 4- file_name_04_VC_integral.txt: data for VC
 
 5- file_name_05_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW, Dmax)
 
 6- A summary Graphs as png or svg ready to put in your labbook

<img width="288" height="225" alt="image" src="https://github.com/user-attachments/assets/bdc92d55-ca3c-48d6-8845-af79ddf52a46" />


## Plots:

  1- Form factor (can be display as log-log by checking the box)
  
  2- Guinier Fit (Rg, I(0), qmin*Rg qmax*Rg,nbeg, nend) & Residuals of the fit (check aggregation and or repulsion, here looks nice)
  
  3- Normalised Kratky plot (here presence of disordered regions)
  
  4- Plot of the Volume of correlation and cumulative integral for MW determination (severalat 0.3, 8/Rg and from Porod cut at 8/Rg)

  5- P(r) and the fit 


## Stand alone :


Macapp stand alone version can be downloaded here :



If the app doesn't start go to 
	1.	Open System Settings (or System Preferences on older macOS versions).
 
	2.	Go to Privacy & Security.
 
	3.	Scroll down to the Security section.
 
	4.	If macOS has blocked the app, you’ll see a message saying it was prevented from opening.
 
	5.	Click “Open Anyway” to allow the app to run.
 
 Windows stand alone can be found:
 


 
 Linux  stand alone can be found here (coming soon):






 
