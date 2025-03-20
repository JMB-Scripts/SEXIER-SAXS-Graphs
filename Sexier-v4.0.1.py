"""
User manual
Jean-Marie Bourhis (jean-marie.bourhis@ibs.fr) and GPT via mac-gpt (because I'm a "tanche" (kind of fish) with python)
2025
 This Python script is designed to asses the quality of your SAXS data and keep a trace of it :
 1- From Guinier approximation, I(0), Rg
 2- the Kratky graph (presence or not of disordered region)
 3- the correlation volume to get the MW.

 It will allow you to keep all these information in several simple text files, allowing you to use your own graphical software (Excel, Graphpad, Origin, ...).
 
 Prerequisites
- Python 3.x
- Necessary packages: numpy, matplotlib, scipy
- SAXS Data need to be in Å-1


 v4 introduce a gui where you can 
    1.browse to you `filename.dat` : the name of the .dat file containing the SAXS experimental data.
    2. you can estimate qmin and qmax automatically (not as good as Raw or Primus), or enter the values
    - `qmin_offset` : the offset (in number of lines) to be added to the first usable line to determine qmin, use the value from PRIMUS or RAW.
    - `qmax_offset`: the offset (in number of lines) to be added to the first usable line to determine qmax use the value from PRIMUS or RAW.
 
 Features:

 1. Guinier approximation :
 - Read .dat file and determine first usable line.
 - Extraction of data for q and I(q) in the selected range.
 - Linear regression to calculate Rg (radius of gyration) and I0 (intensity at q=0).
 - Write data to text file.
 - Display graph with experimental points and theoretical curve.
 2. Kratky 2:
 - Extract data in selected range for Kratky.
 - Calculation and normalization of values for Kratky ((𝑞𝑅𝑔)^2).𝐼(𝑞)/𝐼(0) vs 𝑞𝑅𝑔.
 - Write data to text file.
 - Display Kratky graph.
 3. Correlation volume (CV):
 - Extract data up to q=0.3.
 Calculate the integral of the product I(q)*q.
 - Calculation of VC, QR (quality factor) and MW (molecular weight).
 - Write data to text file.
 4. Summary file:
 - Write values for Rg, I0, qmin_Rg, qmax_Rg, MW to a text file.

 5. Summary Graphs as png or svg ready to put in your labbook

 Output
 The script will generate the following files:
   - filename_rg.txt: Guinier approximation data (q^2, ln(I_exp), ln(I_theo), normalised residuals)
   - file_name_Normalized-Kratky.txt: data for normalised Kratky (x, y)
   - nom_du_fichier_VC.txt: data for VC (q, I(q)*q)
   - file_name_Summary.txt: summary file (Rg, I0, qmin_Rg, qmax_Rg, MW)
4.1 
- Change the Y axis label of the VC plot 

"""
import tkinter as tk
from tkinter import filedialog, messagebox
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.backends.backend_svg import FigureCanvasSVG
from scipy import integrate
from scipy.integrate import simpson 
from scipy.stats import linregress


# def load data 
def load_data(data_file):  # function to load and pre-process data

    first_usable_line = None
    last_usable_line = None

    try:
        with open(data_file, 'r') as file:
            lines = file.readlines()
            for line_number, line in enumerate(lines):
                line = line.strip()
                if line and len(line.split()) >= 3:
                    try:
                        values = [float(x) for x in line.split()[:3]]
                        if not np.any(np.isnan(values)):
                            if first_usable_line is None:
                                first_usable_line = line_number + 1
                            last_usable_line = line_number + 1
                    except ValueError:
                        continue

        if first_usable_line is None:
            raise ValueError("No usable lines found in the file.")  # Raise an exception

        nbr_lines = last_usable_line - first_usable_line
        data = np.loadtxt(data_file, skiprows=first_usable_line - 1, max_rows=nbr_lines)
        return data, first_usable_line, nbr_lines # Return all values

    except FileNotFoundError:
        raise  # Re-raise the FileNotFoundError
    except ValueError as e:
        raise  # Re-raise the ValueError
    except Exception as e:
        raise  # Re-raise other exceptions
#
# AUTO GUINIER so difficult to do 
#
def estimate_guinier_region():
    data_file = data_file_entry.get()

    # Create a new top-level window to show calculation progress
    calculation_window = tk.Toplevel(root)  
    calculation_window.title("Calculating...")
    calculation_label = tk.Label(calculation_window, text="Finding Guinier region (please wait...)It's not so trivial)")
    calculation_label.pack(padx=10, pady=10)
    calculation_window.update()
    
    try:
        # Load data
        data, _, _ = load_data(data_file)
        q = data[:, 0]
        I = data[:, 1]
        q_squared = q**2
        ln_I = np.log(I)
        q_squared_valid = q_squared
        ln_I_valid = ln_I

        # Set constraints for the search
        min_data_points = 5  # Minimum points for a valid fit
        max_data_points = 160  # Maximum points to consider
        max_qmin_index = 40  # Restrict qmin to the first 30 points
        
        best_score = -np.inf
        best_qmin_idx = 0
        best_qmax_idx = min_data_points - 1
        best_rg = 0
        best_i0 = 0
        best_r2 = 0
        
        # Try different window sizes to find the best linear region
        for window_size in range(min_data_points, min(max_data_points, len(q_squared_valid))):
            for start_idx in range(0, min(max_qmin_index, len(q_squared_valid) - window_size)):
                end_idx = start_idx + window_size
                
                # Get window data
                q_window = q_squared_valid[start_idx:end_idx]
                ln_I_window = ln_I_valid[start_idx:end_idx]
                
                # Calculate linear regression
                slope, intercept, r_value, p_value, std_err = linregress(q_window, ln_I_window)
                r2 = r_value**2
                
                # Only consider negative slopes (Guinier should have negative slope)
                if slope >= 0:
                    continue
                    
                # Calculate Rg from slope and I0 from intercept
                rg = np.sqrt(-3 * slope)
                i0 = np.exp(intercept)
                
                # Calculate qRg values
                qmin_rg = q[start_idx] * rg
                qmax_rg = q[end_idx - 1] * rg
                
                # Skip if outside desired qRg range
                if qmax_rg < 1.1 or qmax_rg > 1.5:
                    continue
                
                # Skip if qmin_rg is too large (prioritize smaller qmin_rg)
                if qmin_rg > 0.8:
                    continue
                    
                # Calculate comprehensive score based on all criteria
                score = 0
                
                # Criterion 1: qmin_rg as small as possible
                # Smaller values get higher scores, with diminishing returns below 0.3
                if qmin_rg < 0.3:
                    score += 50  # Bonus for very small qmin_rg
                else:
                    score += 40 * (1 - qmin_rg/0.8)  # Linear decrease as qmin_rg approaches 0.8
                
                # Criterion 2: qmax_rg as close to 1.3 as possible
                # Perfect score at 1.3, decreasing as we move away
                score += 40 * (1 - abs(qmax_rg - 1.3)/0.4)
                
                # Criterion 3: Smallest fit error in Rg
                # Normalize by the average value of Rg to make this criterion scale-independent
                rg_error = std_err * rg / (-3 * slope)
                score += 30 * max(0, 1 - 20 * rg_error)  # Higher penalty for large errors
                
                # Criterion 4: Smallest fit error in I(0)
                # Normalize by the value of I(0)
                i0_error = std_err / abs(intercept)
                score += 30 * max(0, 1 - 20 * i0_error)
                
                # Criterion 5: R² as close to 1 as possible
                # Scale exponentially - small differences near 1.0 matter more
                if r2 > 0.98:
                    score += 50  # Bonus for excellent fits
                else:
                    score += 40 * (r2**4)  # Strongly favor high R²
                
                # Criterion 6: Largest range of q values used as possible
                # Normalize by the total range to make this scale-independent
                q_range = q[end_idx - 1] - q[start_idx]
                q_total_range = q[-1] - q[0]
                score += 20 * (q_range / q_total_range) * window_size / max_data_points
                
                # Update best values if this score is better
                if score > best_score:
                    best_score = score
                    best_qmin_idx = start_idx
                    best_qmax_idx = end_idx - 1
                    best_rg = rg
                    best_i0 = i0
                    best_r2 = r2
        
        if best_score > -np.inf:
            # Update the entry fields with the estimated values
            qmin_offset_entry.delete(0, tk.END)
            qmin_offset_entry.insert(0, str(best_qmin_idx))
            
            qmax_offset_entry.delete(0, tk.END)
            qmax_offset_entry.insert(0, str(best_qmax_idx))
            
            qmin_rg = q[best_qmin_idx] * best_rg
            qmax_rg = q[best_qmax_idx] * best_rg
            q_range = q[best_qmax_idx] - q[best_qmin_idx]
            
            # Calculate the Rg error properly
            slope, intercept, _, _, std_err = linregress(
                q_squared_valid[best_qmin_idx:best_qmax_idx+1], 
                ln_I_valid[best_qmin_idx:best_qmax_idx+1]
            )
            rg_error = np.sqrt((-3 * std_err / (2 * slope))**2)
            
            messagebox.showinfo("Estimation Complete",
                                f"Estimated range: indices {best_qmin_idx}-{best_qmax_idx}\n"
                                f"Corresponding to q range: {q[best_qmin_idx]:.4f}-{q[best_qmax_idx]:.4f} ({q_range:.4f} Å⁻¹)\n"
                                f"Estimated Rg: {best_rg:.2f} ± {rg_error:.2f} Å\n"
                                f"Estimated I(0): {best_i0:.2f}\n"
                                f"qmin·Rg: {qmin_rg:.2f}, qmax·Rg: {qmax_rg:.2f}\n"
                                f"R² value: {best_r2:.4f}\n"
                                f"Score: {best_score:.1f}")
        else:
            messagebox.showwarning("Estimation Result", 
                                  "Could not find a satisfactory Guinier region.\n"
                                  "Please set qmin and qmax manually.")
    
    except FileNotFoundError as e:
        messagebox.showerror("Error", str(e))
    except ValueError as e:
        messagebox.showerror("Error", str(e))
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {e}")
    finally:
        calculation_window.destroy()

#######
#
# SEXIER
#
########

def process_data():

    try:
        data_file = data_file_entry.get()
        data, first_usable_line, nbr_lines = load_data(data_file)  # Call the new function
        
        qmin_offset = int(qmin_offset_entry.get())
        qmax_offset = int(qmax_offset_entry.get())

        if not os.path.exists(data_file):
            raise FileNotFoundError(f"File not found: {data_file}")
        source_file_prefix = os.path.splitext(os.path.basename(data_file))[0]

# define here all the output file for processing
        # Get the directory of the input data file
        data_dir = os.path.dirname(data_file)  # to put files in the same directory than the dat file
        output_file = os.path.join(data_dir,f'{source_file_prefix}_01_Rg.txt')
        output_file2 = os.path.join(data_dir,f'{source_file_prefix}_02_Norm-Kratky.txt')
        output_file3 = os.path.join(data_dir,f'{source_file_prefix}_03_VC.txt')
        output_file4 = os.path.join(data_dir,f'{source_file_prefix}_04_VC_integral.txt')
        output_file5 = os.path.join(data_dir,f'{source_file_prefix}_05_Summary.txt')
        output_file6 = os.path.join(data_dir,f'{source_file_prefix}_Graphs.png')
        output_file7 = os.path.join(data_dir,f'{source_file_prefix}_Graphs.svg')

# Define what is on each colomun
        q = data[:, 0]
        I = data[:, 1]
        E = data[:, 2]

# Define q and I based on the line 
        qmin = q[qmin_offset]
        qmax = q[qmax_offset]
        I_qmin = I[qmin_offset]
        I_qmax = I[qmax_offset]

# Define qrange
        q_range = q[(q >= qmin) & (q <= qmax)]
        I_range = I[(q >= qmin) & (q <= qmax)]
        E_range = E[(q >= qmin) & (q <= qmax)]

# make the linear fit 
        m, c = np.polyfit(q_range**2, np.log(I_range), deg=1)

#Extract Rg and I0
        Rg = np.sqrt(-3 * m)
        I0 = np.exp(c)

#Extract qminRg and qmaxRg
        qmin_Rg = qmin * Rg
        qmax_Rg = qmax * Rg

#Prepare residual data 
        q_carre = q_range * q_range #not nice but works
        ln_intensity_theoretical = m * q_range**2 + c
        ln_intensity_exp = np.log(I_range)
        residuals = (np.log(I_range) - ln_intensity_theoretical) / (E_range * 10)

        with open(output_file, 'w') as file:
            file.write("q**2\tln(I_exp)\tln(I_theo)\tnormalized-residuals\n")
            for i in range(len(q_range)):
                file.write(f"{q_carre[i]}\t{ln_intensity_exp[i]}\t{ln_intensity_theoretical[i]}\t{residuals[i]}\n")

        firstqmin = q[1]
        firstqmax = q[nbr_lines - 1]

        q_full = q[(q >= firstqmin) & (q <= firstqmax)]
        I_full = I[(q >= firstqmin) & (q <= firstqmax)]

# Prepare kratky data 
        xkrat = q_full * Rg
        zkrat = xkrat * xkrat * I_full
        ykrat = np.divide(zkrat, I0)

        with open(output_file2, 'w') as file:
            file.write("q.Rg\t(q.Rg)^2.(I(q)/I(0)\n")
            for j in range(len(q_full)):
                file.write(f"{xkrat[j]}\t{ykrat[j]}\n")

        q_filtered = q[q <= 0.3]
        I_filtered = I[q <= 0.3]

# Prepare VC plot
        yvc = I_filtered * q_filtered
        Intgr = integrate.simpson(yvc,x=q_filtered)

        VC = I0 / Intgr
        QR = VC**2 / Rg
        MW1 = QR / 0.1231
        
# Prepare VC plot if the limit of integration is 8//rg
        q_vc = 8 / Rg
        q_alt = q_full[q_full <= q_vc]
        I_alt = I_full[q_full <= q_vc]

        yvc_alt = I_alt * q_alt
        Intgr_alt = integrate.simpson(yvc_alt,x=q_alt)

        VC_alt = I0 / Intgr_alt
        QR_alt = VC_alt**2 / Rg
        MW_alt = QR_alt / 0.1231

        q_integral_cumulative = [integrate.simpson(I_full[:i] * q_full[:i], x= q_full[:i]) for i in range(1, len(q_full) + 1)]

        with open(output_file4, 'w') as file:
            file.write("q\tCumulative Integral\n")
            for i in range(len(q_full)):
                file.write(f"{q_full[i]}\t{q_integral_cumulative[i]:.6f}\n")

        yvc_plot = I_full * q_full
        Intgr_plot = integrate.simpson(yvc_plot,x=q_full)

        with open(output_file3, 'w') as file:
            file.write("q\tI(q)*q\n")
            for i in range(len(q_full)):
                file.write(f"{q_full[i]}\t{yvc_plot[i]}\n")
# Prepare porod plot
        q_porod = q_full[q_full <= q_vc]
        I_porod = I_full[q_full <= q_vc]

        q_squared_I = q_porod**2 * I_porod
        Q_p = integrate.simpson(q_squared_I, x=q_porod)
        V_p = (2 * np.pi**2 * I0) / Q_p

# Calculate MW from VC plot
        A = (-2.114e6 * (q_vc**4) + 2.920e6 * (q_vc**3) - 1.472e6 * (q_vc**2) + 3.349e5 * (q_vc) - 3.577e4)
        B = (12.09 * (q_vc**3) - 9.39 * (q_vc**2) + 3.03 * (q_vc) + 0.29)
        V_prime = A + V_p * B
        MV_P = (V_prime / 1.2)

#####################################
#
# Summary
#
#####################################
# Output summary 
        with open(output_file5, 'w') as file:
            file.write("Rg\n")
            file.write(f"{Rg:.3f}\n")
            file.write("I0\n")
            file.write(f"{I0:.3f}\n")
            file.write("qmin_Rg\n")
            file.write(f"{qmin_Rg:.3f}\n")
            file.write("qmax_Rg\n")
            file.write(f"{qmax_Rg:.3f}\n")
            file.write("nbeg\n")
            file.write(f"{qmin_offset}\n")
            file.write("nend\n")
            file.write(f"{qmax_offset}\n")  
            file.write("MW from Vc cut q=0.3\n")
            file.write(f"{MW1:3f}\n")
            file.write("MW from Vc cut q=8/Rg\n")
            file.write(f"{MW_alt:3f}\n")
            file.write("Vp Porod\n")
            file.write(f"{V_p:2f}\n")
            file.write("MW Porod\n")
            file.write(f"{MV_P:2f}\n")
            
#####################################
#
# PLOT
#
#####################################

        col1='#0D92F4'#plot guinier and residual
        col2='#77CDFF'
        col3='#F95454'
        col4='#C62E2E'
        col5='#F3F3E0'#boxbackground
        col6='#BC7C7C' #cumulative 
        # Create a figure with 2x2 subplots, four panels A,B,C,D

        fig, axs = plt.subplots(2, 2, figsize=(13, 11))

        # A- Upper-left: Form Factor plot (log(I) vs q)
        # Filter out negative I_full values
        positive_indices = I_full > 0
        q_full_positive = q_full[positive_indices]
        I_full_positive = I_full[positive_indices]
        # Plot only positive values
        axs[0, 0].plot(q_full_positive, np.log(I_full_positive), label='Log(I) vs q', color=col1)
        axs[0, 0].set_xlabel('q')
        axs[0, 0].set_ylabel('log(I(q))')
        axs[0, 0].set_title(str(source_file_prefix) +' Form Factor' ) #add name of the file
        axs[0, 0].legend(loc='upper center')
        axs[0, 0].grid(True)

        # B-1-Lower-left: Guinier plot
        axs[1, 0].errorbar(q_range**2, np.log(I_range), yerr=E_range, fmt='o', markersize=4, label='Experimental', color=col1)
        axs[1, 0].plot(q_range**2, ln_intensity_theoretical, label='Fit', color=col4)
        axs[1, 0].set_xlabel('q**2')
        axs[1, 0].set_ylabel('ln(I(q))')
        axs[1, 0].set_title('Guinier Plot')
        axs[1, 0].ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True) #sci test
        axs[1, 0].legend(loc='upper center')
        axs[1, 0].grid(True)

        # Add Guinier details as text
        guinier_text = f'Rg = {Rg:.3f}\nI0 = {I0:.3f}\nqmin * Rg = {qmin_Rg:.3f}\nqmax * Rg = {qmax_Rg:.3f}\nnbeg = {qmin_offset}\nnend = {qmax_offset}'
        axs[1, 0].text(0.72, 0.6, guinier_text, transform=axs[1, 0].transAxes, fontsize=10,
                    verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor=col5, alpha=0.5))# Adjusted position and size: left, bottom,

        # B-2-Create inset for residuals (smaller and repositioned)
        inset_residuals = fig.add_axes([0.09, 0.1, 0.21, 0.09])  # Adjusted position and size: left, bottom, width, height  
        inset_residuals.plot(q_range**2, residuals, 'o', markersize=3, color=col1, label='Residuals')
        inset_residuals.axhline(y=0, color=col4, linewidth=1)
        #inset_residuals.set_title('Residuals', fontsize=9)
        inset_residuals.set_xlabel('q**2', fontsize=7)
        inset_residuals.set_ylabel('Residuals/sigma', fontsize=7)
        inset_residuals.ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True) #sci test
        inset_residuals.xaxis.get_offset_text().set_fontsize(7)
        inset_residuals.tick_params(axis='both', which='major', labelsize=7)
        #inset_residuals.grid(True)


        # C- Plot Kratky curve (top-right)
        axs[0, 1].plot(xkrat, ykrat, color=col3)
        axs[0, 1].set_xlabel('(qRg)^2')
        axs[0, 1].set_ylabel('(qRg)^2 * I(q) / I(0)')
        axs[0, 1].set_title('Normalized Kratky Plot')
        axs[0, 1].axvline(x=1.73, color='grey', linestyle='--')
        axs[0, 1].axhline(y=1.1, color='grey', linestyle='--')
        axs[0, 1].grid(True)

        # D- Plot the Volume of Correlation (VC) and Cumulative Integral
        ax_vc = axs[1, 1]  # Left y-axis for Volume of Correlation
        ax_integral = ax_vc.twinx()  # Create a twin axes for the right y-axis

        # Plot the Volume of Correlation
        ax_vc.plot(q_full, yvc_plot, label='Volume of Correlation', color=col1)
        ax_vc.set_xlabel('q')
        ax_vc.set_ylabel('(q*I(q) / I(0)', color=col1)
        ax_vc.tick_params(axis='y', labelcolor=col1)
        ax_vc.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)#sci test
        ax_vc.set_title('Volume of Correlation and Cumulative Integral')
        ax_vc.axvline(x=0.3,color=col1, linestyle='--',alpha=0.8)
        ax_vc.axvline(x=8/Rg,color=col2, linestyle='--',alpha=0.8)
        # Add annotation for molecular weight (MW) in the VC plot

        # Define text parts with different colors and formatted with a space for thousands
        vc_text1 = f'MW Vc(q<0.3) = {MW1:,.0f}'.replace(',', ' ')  # Blue, formatted with space
        vc_text2 = f'MW Vc(q<8/Rg) = {MW_alt:,.0f}'.replace(',', ' ')  # Red, formatted with space
        vc_text3 = f'MW Porod (q<8/Rg) = {MV_P:,.0f}'.replace(',', ' ')  # Red, formatted with space
        # Create a background box for the text
        bbox_props = dict(boxstyle='round,pad=0.3', facecolor=col5, alpha=0.5)

        # Add text with colored box for Volume of Correlation
        ax_vc.text(0.17, 0.29, vc_text1,
                transform=ax_vc.transAxes, fontsize=10,
                verticalalignment='top', color=col1,
                bbox=bbox_props)  # Color for first mass

        ax_vc.text(0.17, 0.22, vc_text2,
                transform=ax_vc.transAxes, fontsize=10,
                verticalalignment='top', color=col2,
                bbox=bbox_props)  # Color for second mass
        ax_vc.text(0.17, 0.15, vc_text3,
                transform=ax_vc.transAxes, fontsize=10,
                verticalalignment='top', color=col4,
                bbox=bbox_props)  # Color for Porod

        # Plot the cumulative integral on the right y-axis
        ax_integral.plot(q_full, q_integral_cumulative, color=col3, label='Cumulative Integral')
        ax_integral.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)#sci test
        ax_integral.set_ylabel('Integral of qI(q) dq', color=col3)
        ax_integral.tick_params(axis='y', labelcolor=col3)

        # Adding legends for both axes
        ax_vc.legend(loc='lower left')
        ax_integral.legend(loc='lower right')
        ax_vc.grid(True)

        ###############################################
        # Adjust layout to avoid overlap
        #plt.tight_layout()
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05) #marge on the figure  
        plt.savefig(output_file6)  # Save the graph as a PNG image
        plt.savefig(output_file7,format="svg")  # Save the graph as a svg image
        # Show the figure
        plt.show()
        plt.close()  # VERY IMPORTANT: Close the plot after saving!


# If try fails
    except FileNotFoundError as e:
        messagebox.showerror("Error", str(e))
    except ValueError as e:
        messagebox.showerror("Error", "Invalid input. qmin and qmax offsets must be integers.")
    except OSError as e:
        messagebox.showerror("Error", f"Error writing file: {e}")  # Handle file writing errors
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

######
# Action under browse button
######

def browse_file():
    filename = filedialog.askopenfilename(title="Select Data File", filetypes=[("Data Files", "*.dat"), ("All Files", "*.*")])
    if filename:
        data_file_entry.delete(0, tk.END)
        data_file_entry.insert(0, filename)

######
# Action under Reset button
#####        

def reset_entries():
    """Clears the file entry and offset entries."""
    data_file_entry.delete(0, tk.END)
    qmin_offset_entry.delete(0, tk.END)
    qmax_offset_entry.delete(0, tk.END)

######
# GUI setup
######

root = tk.Tk()
root.title("SEXIER Data Processing")

# File Selection
data_file_label = tk.Label(root, text="Data File:")
data_file_label.grid(row=0, column=0, padx=5, pady=5)

data_file_entry = tk.Entry(root, width=50)
data_file_entry.grid(row=0, column=1, padx=5, pady=5)

browse_button = tk.Button(root, text="Browse", command=browse_file)
browse_button.grid(row=0, column=2, padx=5, pady=5)

# Offset Inputs
qmin_offset_label = tk.Label(root, text="qmin Offset:")
qmin_offset_label.grid(row=1, column=0, padx=5, pady=5)

qmin_offset_entry = tk.Entry(root)
qmin_offset_entry.grid(row=1, column=1, padx=5, pady=5)

qmax_offset_label = tk.Label(root, text="qmax Offset:")
qmax_offset_label.grid(row=2, column=0, padx=5, pady=5)

qmax_offset_entry = tk.Entry(root)
qmax_offset_entry.grid(row=2, column=1, padx=5, pady=5)

# Add Auto-Estimate button
estimate_button = tk.Button(root, text="Auto (roughly) Estimate Guinier Region", command=estimate_guinier_region)
estimate_button.grid(row=3, column=0, columnspan=3, pady=5)

# Process button
process_button = tk.Button(root, text="Process Data", command=process_data)
process_button.grid(row=4, column=0, columnspan=3, pady=5)

# Reset Button
reset_button = tk.Button(root, text="Reset", command=reset_entries)
reset_button.grid(row=5, column=0, columnspan=3, pady=5)  # Place it below the other buttons

# Quit button
quit_button = tk.Button(root, text="Quit", command=root.quit)
quit_button.grid(row=6, column=0, columnspan=3, pady=5)

# Version
made_by_label = tk.Label(root, text="JMB-Scripts - 4.0.1 -")
made_by_label.grid(row=7, column=0, columnspan=3, pady=5)  # Spanning 2 columns for centering

root.mainloop()