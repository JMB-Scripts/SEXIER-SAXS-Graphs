# Update 

Qt.4.0.2
- Add a Qt interface to see directly your changes 
- Add a check box to plot for the form factor log(I(q) vs log(q)
  
4.0.1
Minor:

- Correction on the label of the axis (bad copy/paste) for the Correlation volume the right axis is I(q)*q).


# SEXIER-SAXS-Graphs
Regain control of your SAXS data.
This Python script is designed to process SAXS (Small Angle X-ray Scattering) data and to generate results such as the Guinier approximation, the Kratky graph and the volume of correlation in a text file format.
Then you can make YOUR OWN representations for these graphs, with YOUR preferred software (Excel, Prism, Origin etc.). 

If you found this project useful, used it, or needed to customize it, please let me know! 
Your feedback is essential to help me improve and continue this project. You can reach out to me directly at [reach out to me via email](jean-marie.bourhis@univ-grenoble-alpes.fr).


## User manual

# Prerequisites
- Python 3.x
- Necessary packages: numpy, matplotlib, scipy
- SAXS Data need to be in Å-1

## Command syntax

v4 introduce a gui where you can :
    1. Browse to your `filename.dat` : the name of the .dat file containing the SAXS experimental data.
    2. You can estimate qmin and qmax automatically (not as good as Raw or Primus), or enter the known values.
    - `qmin_offset` : the offset (in number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS or RAW.
    - `qmax_offset`: the offset (in number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS or RAW.
    3. Process the data 

![image](https://github.com/user-attachments/assets/b04afb01-bf6f-4e54-a428-6684e08b3acd)



## Features

 1. Guinier approximation :
 - Read .dat file and determine a first usable line.
 - Data extraction for q and I(q) in the selected range.
 - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
 - Write data to text file.
 - Display graph with experimental points and theoretical curve.

 2. Kratky 2:
 - Extract data in the selected range for Kratky.
 - Calculation and normalisation of values for Kratky (𝑞𝑅𝑔)^2.𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
 - Write data to a text file.
 - Display Kratky graph.

 3. Correlation volume (CV):
 - Extract data up to q=0.3 or up to 8/Rg.
 Calculate the integral of the product I(q)*q.
 - Calculation of VC, QR (quality factor) and MW (molecular weight).
 - Calculation of the Mw from Porod using SAXSMOW approach with a cut at 8/Rg
 - Write data to text file.

 4. Summary file:
 - Write values for Rg, I0, qmin_Rg, qmax_Rg, MW to a text file.

## Output
 The script will generate the following files:
 
 1- filename_01_Rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalized residuals)
 
 2- file_name_02_Norm-Kratky.txt: data for normalized Kratky (x, y)
 
 3- file_name_03_VC.txt: data for VC (q, I(q)*q)
 
 4- file_name_04_VC_integral.txt: data for VC
 
 5- file_name_05_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW)
 
 6- A summary Graphs as png or svg ready to put in your labbook
 

![image](https://github.com/user-attachments/assets/834edc29-8e5a-4ac1-9952-096e2127d903)

   
## Plots:
  1- Form factor
  2- Guinier Fit (Rg, I(0), qmin*Rg qmax*Rg,nbeg, nend) & Residuals of the fit (check aggregation and or repulsion, here looks nice)
  3- Normalised Kratky plot (here presence of disordered regions)
  4- Plot of the Volume of correlation and cumulative integral (MW determination at 0.3, 8/Rg and from Porod cut at 8/Rg)

![Capture d’écran 2024-10-25 à 13 39 00](https://github.com/user-attachments/assets/3e07d5a6-09f5-4b3b-9b10-cd18248121e5)

## Notes:

/!\ The first start is always very long, it seems that nothing append but be patient. This is linked to matplot that build its font /!\

Macapp stand alone version can be downloaded here :

[Here](https://cloud.univ-grenoble-alpes.fr/s/DowBSzr9TzC7Qkg)

If the app doesn't start go to 
	1.	Open System Settings (or System Preferences on older macOS versions).
	2.	Go to Privacy & Security.
	3.	Scroll down to the Security section.
	4.	If macOS has blocked the app, you’ll see a message saying it was prevented from opening.
	5.	Click “Open Anyway” to allow the app to run.
 
 Windows stand alone can be found:
 
 [Here](https://cloud.univ-grenoble-alpes.fr/s/TC2WecCtjyBz3ii)

 
 Linux  stand alone can be found here:

 [Here](https://cloud.univ-grenoble-alpes.fr/s/HNrTmHtDrESr7Cg)

 
## References :

Putnam, C. D., Hammel, M., Hura, G. L. & Tainer, J. a (2007). Q. Rev. Biophys. 40, 191–285. DOI: 10.1017/S0033583507004635

Rambo, R. P. & Tainer, J. A. (2013). Nature. 496, 477–481. DOI: 10.1038/nature12070

J. Perez, P. Vachette, D. Russo, M. Desmadril, D. Durand (2001) J. Mol. Biol., 308 , pp. 721-743. DOI: 10.1006/jmbi.2001.4611
