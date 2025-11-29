"""
Sexier-PyS v7.2 : A tool for Small-Angle X-ray Scattering (SAXS) data analysis.

This script provides a graphical user interface (GUI) built with PySide6 and Matplotlib
to perform a comprehensive analysis of a single SAXS profile, including Guinier
analysis and Pair-Distance Distribution Function (P(r)) analysis.

---
SUMMARY OF FUNCTIONALITY (v7.2):
---
1.  **P(r) Analysis Enhancements**:
    -   Added an input field for `q_min` for more control over the BIFT calculation range (thx to BioXTAS RAW).
    -   Saving of P(r) plots and data is now fully automatic after a successful run.
    -   The BIFT fit to the experimental data is now saved as a separate text file.
    -   The calculated Dmax and Rg from P(r) are appended to the main summary file.
    -   Removed manual save buttons for P(r) analysis for a cleaner interface.

---
PREVIOUS FUNCTIONALITY (v7.0):
---
-   **New 3-Panel Layout**: Central control panel, 4-panel Guinier analysis plot, and a P(r) distribution plot.
-   **Internal BIFT for P(r) Analysis**: Pure Python implementation of the BIFT algorithm.
-   **Data Loading & Unit Conversion**: Automatically loads a .dat file and converts q-units from nm⁻¹ to Å⁻¹.
-   **Automated & Manual Guinier Analysis**: For determining Rg and I(0).
-   **Molecular Weight (MW) & Oligomeric State Estimation**: Based on the Volume of Correlation.
-   **Comprehensive Plotting & Automatic Saving**: Automatic saving of plots (PNG/SVG) and data files.

"""
import os
import numpy as np
import sys
import matplotlib
matplotlib.use('qtagg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy import integrate
from scipy.integrate import simpson
from scipy.stats import linregress
import shutil
import multiprocessing
import functools
import threading
import scipy.optimize
import scipy.interpolate
from numba import jit

from PySide6 import QtWidgets
from PySide6.QtWidgets import (QMainWindow, QFileDialog, QMessageBox, QVBoxLayout, QWidget,
                             QHBoxLayout, QGridLayout, QCheckBox, QFrame, QPushButton,
                             QLineEdit, QLabel, QProgressDialog)
from PySide6.QtCore import Qt, QThread, Signal, Slot
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.ticker import FormatStrFormatter

#
# BIFT PYTHON IMPLEMENTATION & HELPER FUNCTIONS
#

def load_data(data_file):
    """Loads SAXS data from a .dat file, automatically detecting the data block."""
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
                            if first_usable_line is None: first_usable_line = line_number + 1
                            last_usable_line = line_number + 1
                    except ValueError:
                        continue
        if first_usable_line is None:
            raise ValueError("No usable data lines found in the file.")
        nbr_lines = last_usable_line - first_usable_line
        data = np.loadtxt(data_file, skiprows=first_usable_line - 1, max_rows=nbr_lines)
        return data
    except Exception as e:
        raise

@jit(nopython=True, cache=True)
def createTransMatrix(q, r):
    """
    Matrix such that when you take T dot P(r) you get I(q),
    The A_ij matrix in equation (2) of Hansen 2000.
    """
    T = np.outer(q, r)
    T = 4*np.pi*(r[1]-r[0])*np.where(T==0, 1, np.sin(T)/T)
    return T

@jit(nopython=True, cache=True)
def distDistribution_Sphere(i0_meas, N, dmax):
    """Creates the initial P(r) function for the prior as a sphere."""
    r = np.linspace(0, dmax, N+1)
    p = r**2 * (1 - 1.5*(r/dmax) + 0.5*(r/dmax)**3)
    p = p * i0_meas/(4*np.pi*dmax**3/24.)
    return p, r

@jit(nopython=True, cache=True)
def makePriorDistribution(i0, N, dmax, dist_type='sphere'):
    if dist_type == 'sphere':
        p, r = distDistribution_Sphere(i0, N, dmax)
    else:
        r = np.linspace(0, dmax, N+1)
        p = np.zeros_like(r)
    return p, r

@jit(nopython=True, cache=True)
def bift_inner_loop(f, p, B, alpha, N, sum_dia):
    ite, maxit, minit, xprec, dotsp, omega = 0, 2000, 100, 0.999, 0, 0.5
    sigma = np.zeros_like(p)
    while ite < maxit and not (ite > minit and dotsp > xprec):
        ite += 1
        sigma[1:-1] = np.abs(p[1:-1]+1e-10)
        p_neg_idx = p[1:-1]<=0
        f_neg_idx = f[1:-1]<=0
        p[1:-1][p_neg_idx] = p[1:-1][p_neg_idx]*-1+1e-10
        f[1:-1][f_neg_idx] = f[1:-1][f_neg_idx]*-1+1e-10
        for k in range(2, N-1): p[k] = (f[k-1] + f[k+1])/2.
        p[1], p[-2] = f[2]/2., p[-3]/2.
        p[0], p[-1] = f[0], f[-1]
        sigma[0] = 10
        for k in range(1, N):
            fsumi = 0
            for j in range(1, N): fsumi += B[k, j]*f[j]
            fsumi -= B[k, k]*f[k]
            fx = (2*alpha*p[k]/sigma[k]+sum_dia[k]-fsumi)/(2*alpha/sigma[k]+B[k,k])
            f[k] = (1-omega)*f[k]+omega*fx
        gradsi = -2*(f[1:-1]-p[1:-1])/sigma[1:-1]
        gradci = 2*(np.sum(B[1:-1,1:-1]*f[1:-1], axis=1)-sum_dia[1:-1])
        wgrads = np.sqrt(np.abs(np.sum(gradsi**2)))
        wgradc = np.sqrt(np.abs(np.sum(gradci**2)))
        if wgrads*wgradc == 0: dotsp = 1
        else: dotsp = np.sum(gradsi*gradci)/(wgrads*wgradc)
    return f, p, sigma, dotsp, xprec

@jit(nopython=True, cache=True)
def getEvidence(params, q, i, err, N):
    alpha, dmax = params
    alpha = np.exp(alpha)
    err = err**2
    p, r = makePriorDistribution(i[0], N, dmax, 'sphere')
    T = createTransMatrix(q, r)
    p[0] = 0
    f = np.zeros_like(p)
    norm_T = T/err.reshape((err.size, 1))
    sum_dia = np.sum(norm_T*i.reshape((i.size, 1)), axis=0)
    sum_dia[0] = 0
    B = np.dot(T.T, norm_T)
    B[0,:], B[:,0] = 0, 0
    c1 = np.sum(np.sum(T[1:4,1:-1]*p[1:-1], axis=1)/err[1:4])
    c2 = np.sum(i[1:4]/err[1:4])
    p[1:-1] = p[1:-1]*(c2/c1)
    f[1:-1] = p[1:-1]*1.001
    f, p, sigma, dotsp, xprec = bift_inner_loop(f, p, B, alpha, N, sum_dia)
    s = np.sum(-(f[1:-1]-p[1:-1])**2/sigma[1:-1])
    c = np.sum((i[1:-1]-np.sum(T[1:-1,1:-1]*f[1:-1], axis=1))**2/err[1:-1])/(i.size - p.size)
    u = np.sqrt(np.abs(np.outer(f[1:-1], f[1:-1])))*B[1:-1, 1:-1]/alpha
    for j in range(0, u.shape[0]): u[j, j] += 1
    rlogdet = np.log(np.abs(np.linalg.det(u)))
    evidence = -np.log(abs(dmax))+(alpha*s-0.5*c*i.size)-0.5*rlogdet-np.log(abs(alpha))
    if evidence <= 0 and dotsp < xprec: evidence *= 30
    elif dotsp < xprec: evidence /= 30.
    return evidence, c, f, r

def getEvidenceOptimize(params, q, i, err, N):
    evidence, _, _, _ = getEvidence(params, q, i, err, N)
    return -evidence

def calc_bift_errors(opt_params, q, i, err, N, mc_runs=300, abort_check=False, single_proc=False, nprocs=0):
    alpha_opt, dmax_opt = opt_params
    mult = 3.0
    ev_array, c_array = np.zeros(mc_runs), np.zeros(mc_runs)
    f_array = np.zeros((mc_runs, N+1)); r_array = np.zeros((mc_runs, N+1))
    max_dmax = dmax_opt+0.1*dmax_opt*0.5*mult
    _, ref_r = makePriorDistribution(i[0], N, max_dmax, 'sphere')
    run_mc = True
    while run_mc:
        alpha_array = alpha_opt+0.1*alpha_opt*(np.random.random(mc_runs)-0.5)*mult
        dmax_array = dmax_opt+0.1*dmax_opt*(np.random.random(mc_runs)-0.5)*mult
        alpha_array[0], dmax_array[0] = alpha_opt, dmax_opt
        pts = list(zip(alpha_array, dmax_array))
        results = [getEvidence(params, q, i, err, N) for params in pts]
        for res_idx, res in enumerate(results):
            evidence, c, f, r = res
            dmax = dmax_array[res_idx]
            interp = scipy.interpolate.interp1d(r, f, copy=False)
            f_interp = np.zeros_like(ref_r)
            f_interp[ref_r<dmax] = interp(ref_r[ref_r<dmax])
            ev_array[res_idx], c_array[res_idx] = evidence, c
            f_array[res_idx,:], r_array[res_idx,:] = f_interp, ref_r
        if np.abs(ev_array).max() >= 9e8:
            mult /= 2.
            if mult < 0.001 or (hasattr(abort_check, 'is_set') and abort_check.is_set()): run_mc = False
        else: run_mc = False
    ev_max = ev_array.max()
    prob = np.exp(ev_array - ev_max)**(1./c_array.min()); prob /= prob.sum()
    p_avg = np.sum(f_array*prob[:,None], axis=0)
    err = np.sqrt(np.abs(np.sum((f_array-p_avg)**2*prob[:,None], axis=0)))
    alpha, dmax, c, evidence = np.sum(alpha_array*prob), np.sum(dmax_array*prob), np.sum(c_array*prob), np.sum(ev_array*prob)
    area = np.trapz(f_array, r_array, axis=1)
    area2 = np.trapz(f_array*r_array**2, r_array, axis=1)
    rg_array = np.sqrt(abs(area2/(2.*area))); i0_array = area*4*np.pi
    rg, i0 = np.sum(rg_array*prob), np.sum(i0_array*prob)
    sd_alpha, sd_dmax = np.sqrt(np.sum((alpha_array - alpha)**2*prob)), np.sqrt(np.sum((dmax_array - dmax)**2*prob))
    sd_c, sd_ev = np.sqrt(np.sum((c_array - c)**2*prob)), np.sqrt(np.sum((ev_array - evidence)**2*prob))
    sd_rg, sd_i0 = np.sqrt(np.sum((rg_array - rg)**2*prob)), np.sqrt(np.sum((i0_array - i0)**2*prob))
    return ref_r, p_avg, err, (alpha, sd_alpha), (dmax, sd_dmax), (c, sd_c), (evidence, sd_ev), (rg, sd_rg), (i0, sd_i0)

def make_fit(q, r, pr):
    qr = np.outer(q, r)
    sinc_qr = np.where(qr==0, 1, np.sin(qr)/qr)
    i = 4*np.pi*np.trapz(pr*sinc_qr, r, axis=1)
    return i

def doBift(q, i, err, filename, npts, alpha_min, alpha_max, alpha_n, dmax_min, dmax_max, dmax_n, mc_runs, queue=None, abort_check=threading.Event(), single_proc=True, nprocs=0):
    start_idx = 0
    for j in range(i.size):
        if i[j] == 0 or err[j] == 0: start_idx = j+1
        else: break
    end_idx = i.size
    for j in range(i.size, 0, -1):
        if i[j-1] == 0 or err[j-1] == 0: end_idx = j-1
        else: break
    q, i, err = q[start_idx:end_idx], i[start_idx:end_idx], err[start_idx:end_idx]
    i_zeros, err_zeros = np.argwhere(i==0), np.argwhere(err==0)
    both_zeros = np.intersect1d(i_zeros, err_zeros)
    q, i, err = np.delete(q, both_zeros), np.delete(i, both_zeros), np.delete(err, both_zeros)
    if npts > len(q)//2: npts = len(q)//2
    alpha_min, alpha_max = np.log(alpha_min), np.log(alpha_max)
    alpha_points = np.linspace(alpha_min, alpha_max, alpha_n)
    dmax_points = np.linspace(dmax_min, dmax_max, dmax_n)
    all_posteriors = np.zeros((dmax_points.size, alpha_points.size))
    N = npts - 1
    if abort_check.is_set(): return None
    for d_idx, dmax in enumerate(dmax_points):
        pts = [(alpha, dmax) for alpha in alpha_points]
        results = [getEvidence(params, q, i, err, N) for params in pts]
        for res_idx, res in enumerate(results): all_posteriors[d_idx, res_idx] = res[0]
        if queue is not None: queue.put({'update': {'spoint': (d_idx+1)*len(alpha_points), 'tpoint': len(alpha_points)*len(dmax_points)}})
        if abort_check.is_set(): return None
    if queue is not None: queue.put({'update': {'status': 'Running minimization'}})
    min_idx = np.unravel_index(np.argmax(all_posteriors, axis=None), all_posteriors.shape)
    min_dmax, min_alpha = dmax_points[min_idx[0]], alpha_points[min_idx[1]]
    opt_res = scipy.optimize.minimize(getEvidenceOptimize, (min_alpha, min_dmax), (q, i, err, N), method='Powell')
    if abort_check.is_set(): return None
    if opt_res.get('success'):
        alpha, dmax = opt_res.get('x')
        if dmax > 0:
            evidence, c, f, r = getEvidence((alpha, dmax), q, i, err, N)
            if queue is not None: queue.put({'update': {'status': 'Calculating Monte Carlo errors'}})
            pr = f
            area, area2 = np.trapz(pr, r), np.trapz(np.array(pr)*np.array(r)**2, r)
            rg, i0 = np.sqrt(abs(area2/(2.*area))), area*4*np.pi
            fit, q_extrap = make_fit(q, r, pr), np.concatenate((np.arange(0, q[1]-q[0], q[1]), q))
            fit_extrap = make_fit(q_extrap, r, pr)
            err_calc = calc_bift_errors((alpha, dmax), q, i, err, N, mc_runs, abort_check, single_proc, nprocs)
            if abort_check.is_set(): return None
            r_err, _, pr_err, a_res, d_res, c_res, ev_res, rg_res, i0_res = err_calc
            interp = scipy.interpolate.interp1d(r_err, pr_err, copy=False)
            err_interp = interp(r)
            results = {'dmax': dmax, 'dmaxer': d_res[1], 'rg': rg, 'rger': rg_res[1], 'i0': i0, 'i0er': i0_res[1], 'chisq': c, 'chisq_er': c_res[1], 'alpha': alpha, 'alpha_er': a_res[1], 'evidence': evidence, 'evidence_er': ev_res[1], 'qmin': q[0], 'qmax': q[-1], 'algorithm': 'BIFT', 'filename': os.path.splitext(filename)[0]+'.ift'}
            final_results = {'pr': pr, 'r': r, 'pr_err': err_interp, 'fit': fit, 'q': q, 'i': i, 'err': err, 'details': results, 'fit_extrap': fit_extrap, 'q_extrap': q_extrap}
            return final_results
    return None

class PythonBiftThread(QThread):
    """Worker thread for running the internal python 'bift' command."""
    bift_finished = Signal(dict)
    bift_error = Signal(str)

    def __init__(self, q, i, err, dmax, q_min=None, q_max=None):
        super().__init__()
        self.q, self.i, self.err = q, i, err
        self.dmax_user = dmax
        self.q_min_user = q_min
        self.q_max_user = q_max

    def run(self):
        try:
            q, i, err = self.q, self.i, self.err

            mask = np.ones_like(q, dtype=bool)
            if self.q_min_user:
                mask &= (q >= self.q_min_user)
            if self.q_max_user:
                mask &= (q <= self.q_max_user)
            q, i, err = q[mask], i[mask], err[mask]

            dmax_min = self.dmax_user * 0.5 if self.dmax_user else 10
            dmax_max = self.dmax_user * 1.5 if self.dmax_user else 500
            
            iftm = doBift(
                q, i, err, filename="temp", npts=51,
                alpha_min=0.001, alpha_max=100, alpha_n=20,
                dmax_min=dmax_min, dmax_max=dmax_max, dmax_n=20,
                mc_runs=100, single_proc=True
            )
            if iftm:
                self.bift_finished.emit(iftm)
            else:
                self.bift_error.emit("BIFT calculation failed or was canceled.")
        except Exception as e:
            self.bift_error.emit(f"An unexpected error occurred during BIFT execution: {str(e)}")

#
# Main Window Class
#
class MainWindow(QMainWindow):
    """The main window of the SEXIER application."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SEXIER - SAXS Data Processing v7.2")
        self.setGeometry(100, 100, 1800, 900)
        
        self.data = None; self.bift_data = None; self.ax_integral = None
        self.output_dir = None; self.source_file_prefix = None
        self.auto_fit_score = None 

        main_layout = QHBoxLayout()

        self.figure4panel, axes_flat = plt.subplots(2, 2, figsize=((8, 8)))
        
        self.axA = axes_flat[0, 0]
        self.axC = axes_flat[1, 0]
        self.axD = axes_flat[1, 1]
        
        gs = axes_flat[0, 1].get_gridspec()
        axes_flat[0, 1].remove()
        gs_nested = gs[0, 1].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        self.axB = self.figure4panel.add_subplot(gs_nested[0]) # Top: Guinier
        self.axB_res = self.figure4panel.add_subplot(gs_nested[1], sharex=self.axB) # Bottom: Residual
        
        self.all_4panel_axes = [self.axA, self.axB, self.axB_res, self.axC, self.axD]
        
        self.canvas4panel = FigureCanvas(self.figure4panel)

        control_widget = QFrame()
        control_widget.setFrameStyle(QFrame.StyledPanel | QFrame.Raised)
        control_widget.setFixedWidth(300)
        control_layout = QVBoxLayout(control_widget)
        self._setup_control_panel(control_layout)
        
        self.figure_bift, self.ax_bift = plt.subplots(2, 1, figsize=(8, 8))
        self.canvas_bift = FigureCanvas(self.figure_bift)

        main_layout.addWidget(self.canvas4panel, 2)
        main_layout.addWidget(control_widget)
        main_layout.addWidget(self.canvas_bift, 1)

        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)

        self.canvas4panel.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.canvas_bift.mpl_connect('motion_notify_event', self.on_plot_hover)

    def _setup_control_panel(self, layout):
        # File Selection
        file_group = QtWidgets.QGroupBox("1. Load Data")
        file_layout = QVBoxLayout()
        self.data_file_entry = QLineEdit()
        self.browse_button = QPushButton("Browse File")
        self.browse_button.clicked.connect(self._browse_file)
        file_layout.addWidget(self.data_file_entry)
        file_layout.addWidget(self.browse_button)
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        # Guinier Region
        guinier_group = QtWidgets.QGroupBox("2. Guinier Region")
        guinier_layout = QGridLayout()
        self.qmin_offset_entry = QLineEdit()
        self.qmax_offset_entry = QLineEdit()
        self.process_button = QPushButton("Re-Process Data")
        self.process_button.clicked.connect(self._process_data)
        guinier_layout.addWidget(QLabel("qmin Index:"), 0, 0)
        guinier_layout.addWidget(self.qmin_offset_entry, 0, 1)
        guinier_layout.addWidget(QLabel("qmax Index:"), 1, 0)
        guinier_layout.addWidget(self.qmax_offset_entry, 1, 1)
        guinier_layout.addWidget(self.process_button, 2, 0, 1, 2)
        guinier_group.setLayout(guinier_layout)
        layout.addWidget(guinier_group)

        # Guinier Results
        results_group = QtWidgets.QGroupBox("3. Guinier Results")
        results_layout = QVBoxLayout()
        self.rg_label = QLabel("Rg:")
        self.i0_label = QLabel("I(0):")
        self.q_range_label = QLabel("q range:")
        self.q_rg_label = QLabel("q·Rg:")
        
        for label in [self.rg_label, self.i0_label, self.q_range_label, self.q_rg_label]:
            results_layout.addWidget(label)
        results_group.setLayout(results_layout)
        layout.addWidget(results_group)

        # MW Analysis
        mw_group = QtWidgets.QGroupBox("4. Molecular Weight (Vc)") 
        mw_layout = QGridLayout()
        
        self.mw_qlimit_entry = QLineEdit() 
        self.mw_qlimit_entry.editingFinished.connect(self._process_data) 
        self.mw_vc_result_label = QLabel("MW (Vc):") 
        self.theoretical_mw_entry = QLineEdit()
        self.theoretical_mw_entry.editingFinished.connect(self._process_data)
        self.oligomer_state_label = QLabel("Oligomer State:")
        
        mw_layout.addWidget(QLabel("q_limit for Vc MW (Å⁻¹):"), 0, 0) 
        mw_layout.addWidget(self.mw_qlimit_entry, 0, 1) 
        mw_layout.addWidget(self.mw_vc_result_label, 1, 0, 1, 2) 
        mw_layout.addWidget(QLabel("MW from sequence (kDa):"), 2, 0) 
        mw_layout.addWidget(self.theoretical_mw_entry, 2, 1)
        mw_layout.addWidget(self.oligomer_state_label, 3, 0, 1, 2) 
        
        mw_group.setLayout(mw_layout)
        layout.addWidget(mw_group)

        # P(r) Analysis
        bift_group = QtWidgets.QGroupBox("5. P(r) Analysis (BIFT)")
        bift_layout = QGridLayout()
        self.bift_dmax_entry = QLineEdit() 
        self.bift_qmin_entry = QLineEdit()
        self.bift_qmax_entry = QLineEdit()
        self.run_bift_button = QPushButton("Run BIFT")
        self.run_bift_button.clicked.connect(self._run_bift)
        self.bift_rg_label = QLabel("Rg from P(r):")
        self.bift_dmax_label = QLabel("Dmax from P(r):")
        bift_layout.addWidget(QLabel("Dmax (Å):"), 0, 0) 
        bift_layout.addWidget(self.bift_dmax_entry, 0, 1) 
        bift_layout.addWidget(QLabel("q_min (Å⁻¹):"), 1, 0) 
        bift_layout.addWidget(self.bift_qmin_entry, 1, 1) 
        bift_layout.addWidget(QLabel("q_max (Å⁻¹):"), 2, 0) 
        bift_layout.addWidget(self.bift_qmax_entry, 2, 1) 
        bift_layout.addWidget(self.run_bift_button, 3, 0, 1, 2) 
        bift_layout.addWidget(self.bift_rg_label, 4, 0, 1, 2) 
        bift_layout.addWidget(self.bift_dmax_label, 5, 0, 1, 2) 
        bift_group.setLayout(bift_layout)
        layout.addWidget(bift_group)
        self.run_bift_button.setEnabled(False) 

        layout.addStretch()

        # Plot Options
        options_group = QtWidgets.QGroupBox("Plot Options")
        options_layout = QVBoxLayout()
        self.log_q_checkbox = QCheckBox("Form Factor: log(I) vs log(q)")
        self.log_q_checkbox.toggled.connect(self._process_data)
        options_layout.addWidget(self.log_q_checkbox)
        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # General Controls
        self.reset_button = QPushButton("Reset All")
        self.reset_button.clicked.connect(self._reset)
        layout.addWidget(self.reset_button)
        self.quit_button = QPushButton("Quit")
        self.quit_button.clicked.connect(self.close)
        layout.addWidget(self.quit_button)
        
        layout.addWidget(QLabel("JMB-Scripts Sexier-v7.2", alignment=Qt.AlignCenter))
        self.coord_label = QLabel("Cursor: (x=, y=)", alignment=Qt.AlignCenter)
        layout.addWidget(self.coord_label)

    def on_plot_hover(self, event):
        if event.inaxes: self.coord_label.setText(f"Cursor: (x={event.xdata:.3e}, y={event.ydata:.3e})")

    def _browse_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Select Data File", "", "Data Files (*.dat);;All Files (*.*)")
        if filename:
            try:
                self._reset()
                self.data_file_entry.setText(filename)
                data = load_data(filename)
                q = data[:, 0]
                if np.max(q) > 2.0:
                    QMessageBox.information(self, "Unit Detection", "q-values appear to be in nm⁻¹. Converted to Å⁻¹.")
                    data[:, 0] = q / 10.0
                self.data = data
                self.source_file_prefix = os.path.splitext(os.path.basename(filename))[0]
                self.output_dir = os.path.join(os.path.dirname(filename), f'{self.source_file_prefix}_Sexier_out')
                os.makedirs(self.output_dir, exist_ok=True)
                self.run_bift_button.setEnabled(True)
                self._estimate_guinier_region()
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e)); self.data = None

    def _estimate_guinier_region(self):
        if self.data is None: return
        progress_dialog = QProgressDialog("Finding Guinier region...", None, 0, 100, self)
        progress_dialog.setWindowModality(Qt.WindowModality.WindowModal); progress_dialog.setMinimumDuration(0); progress_dialog.show()
        try:
            q, I = self.data[:, 0], self.data[:, 1]
            positive_mask = I > 0
            if not np.any(positive_mask):
                QMessageBox.warning(self, "Data Error", "No positive intensity values found.")
                return
            q_safe, I_safe = q[positive_mask], I[positive_mask]
            with np.errstate(invalid='ignore'):
                q_squared_safe, ln_I_safe = q_safe**2, np.log(I_safe)
            
            min_pts, max_pts, max_start = 5, 160, 40
            best_score, best_qmin_idx, best_qmax_idx = -np.inf, 0, min_pts - 1
            total_iter = min(max_pts, len(q_squared_safe)) - min_pts

            for size in range(min_pts, min(max_pts, len(q_squared_safe))):
                for start in range(0, min(max_start, len(q_squared_safe) - size)):
                    end = start + size
                    q_win, lnI_win = q_squared_safe[start:end], ln_I_safe[start:end]
                    
                    if np.any(np.isinf(lnI_win)) or np.any(np.isnan(lnI_win)): continue

                    slope, intercept, r, _, _ = linregress(q_win, lnI_win)
                    if slope >= 0: continue
                    rg = np.sqrt(-3 * slope)
                    qmin_rg, qmax_rg = q_safe[start] * rg, q_safe[end - 1] * rg
                    if not (1.1 <= qmax_rg <= 1.5 and qmin_rg <= 0.8): continue
                    score = (r**8) + (1 - abs(qmax_rg - 1.3))
                    if score > best_score:
                        best_score, best_qmin_idx, best_qmax_idx = score, start, end - 1
                if total_iter > 0: progress_dialog.setValue(int((size / total_iter) * 100))

            if best_score > -np.inf:
                original_qmin_val = q_safe[best_qmin_idx]
                original_qmax_val = q_safe[best_qmax_idx]
                original_qmin_offset = np.where(q == original_qmin_val)[0][0]
                original_qmax_offset = np.where(q == original_qmax_val)[0][0]
                self.qmin_offset_entry.setText(str(original_qmin_offset))
                self.qmax_offset_entry.setText(str(original_qmax_offset))
                self.auto_fit_score = best_score 
                self._process_data()
            else: QMessageBox.warning(self, "Estimation Failed", "Could not find a satisfactory Guinier region.")
        finally: progress_dialog.close()

    def _process_data(self):
        if self.data is None: return
        try:
            q, I, E = self.data[:, 0], self.data[:, 1], self.data[:, 2]
            qmin_idx, qmax_idx = int(self.qmin_offset_entry.text()), int(self.qmax_offset_entry.text())
            
            q_range_full, I_range_full, E_range_full = q[qmin_idx:qmax_idx+1], I[qmin_idx:qmax_idx+1], E[qmin_idx:qmax_idx+1]
            valid_mask = I_range_full > 0
            if np.count_nonzero(valid_mask) < 2:
                QMessageBox.warning(self, "Analysis Error", "Not enough valid data points (I > 0) in selected range.")
                return

            q_range_safe, I_range_safe, E_range_safe = q_range_full[valid_mask], I_range_full[valid_mask], E_range_full[valid_mask]
            with np.errstate(invalid='ignore'):
                slope, intercept, r_val, _, std_err = linregress(q_range_safe**2, np.log(I_range_safe))

            if slope >= 0:
                QMessageBox.warning(self, "Fit Error", "Guinier fit failed: slope is not negative. Try a different range.")
                return

            Rg, I0 = np.sqrt(-3 * slope), np.exp(intercept)
            rg_err = np.sqrt((-3 * std_err / (2 * slope))**2) if slope != 0 else 0
            qmin_Rg, qmax_Rg = q[qmin_idx] * Rg, q[qmax_idx] * Rg
            with np.errstate(invalid='ignore'):
                residuals = (np.log(I_range_safe) - (slope*q_range_safe**2 + intercept)) / (E_range_safe / np.sqrt(I_range_safe))
            
            # --- MW Calculations ---
            q_lt_03 = q <= 0.3; Intgr_03 = simpson(I[q_lt_03]*q[q_lt_03], x=q[q_lt_03])
            MW1 = (I0 / Intgr_03)**2 / Rg / 0.1231 if Intgr_03 != 0 else 0
            
            q_vc_limit = 8 / Rg; q_lt_8Rg = q <= q_vc_limit
            Intgr_8Rg = simpson(I[q_lt_8Rg]*q[q_lt_8Rg], x=q[q_lt_8Rg])
            MW_alt = (I0 / Intgr_8Rg)**2 / Rg / 0.1231 if Intgr_8Rg != 0 else 0
            
            q_sq_I = q[q_lt_8Rg]**2 * I[q_lt_8Rg]; Q_p = simpson(q_sq_I, x=q[q_lt_8Rg])
            V_p = (2*np.pi**2 * I0) / Q_p if Q_p != 0 else 0; MV_P = V_p / 1.2
            
            mw_qlimit_val = None
            MW_qlimit = 0
            mw_qlimit_str = self.mw_qlimit_entry.text()
            if mw_qlimit_str:
                try:
                    mw_qlimit_val = float(mw_qlimit_str)
                    if mw_qlimit_val > 0:
                        q_lt_qlimit = q <= mw_qlimit_val
                        Intgr_qlimit = simpson(I[q_lt_qlimit]*q[q_lt_qlimit], x=q[q_lt_qlimit])
                        MW_qlimit = (I0 / Intgr_qlimit)**2 / Rg / 0.1231 if Intgr_qlimit != 0 else 0
                except ValueError:
                    pass 
            
            mw_for_oligomer = MW_qlimit if MW_qlimit > 0 else MW1
            
            # Oligomer state logic 
            theoretical_mw_str = self.theoretical_mw_entry.text()
            if theoretical_mw_str:
                try:
                    theoretical_mw_kda = float(theoretical_mw_str)
                    ratio = mw_for_oligomer / (theoretical_mw_kda * 1000) 
                    oligomer_n = round(ratio)
                    if abs(ratio - oligomer_n) < 0.2: 
                        state = "Monomer" if oligomer_n == 1 else f"Oligomer n≈{oligomer_n}"
                        self.oligomer_state_label.setText(f"Oligomer State: <font color='green'>{state}</font>")
                    else:
                        self.oligomer_state_label.setText("Oligomer State: <font color='red'>MW is suspect</font>")
                except ValueError:
                    self.oligomer_state_label.setText("Oligomer State: <font color='orange'>Invalid MW</font>")
            else:
                 self.oligomer_state_label.setText("Oligomer State:")

            # Update GUI labels
            self.rg_label.setText(f"Rg: {Rg:.2f} ± {rg_err:.2f} Å")
            self.i0_label.setText(f"I(0): {I0:.2e}")
            self.mw_vc_result_label.setText(f"MW (Vc): {mw_for_oligomer:,.0f} Da") 
            self.q_range_label.setText(f"q range: {q[qmin_idx]:.4f} - {q[qmax_idx]:.4f}")
            self.q_rg_label.setText(f"q·Rg: {qmin_Rg:.2f} – {qmax_Rg:.2f}")

            # Automatic saving of Guinier-related data
            self._save_guinier_data(q, I, q_range_safe, I_range_safe, residuals, slope, intercept, Rg, I0, MW1, MW_alt, V_p, MV_P, qmin_Rg, qmax_Rg, r_val**2, mw_qlimit_val, MW_qlimit)
            
            self._plot_4_panel(q, I, q_range_full, I_range_full, E_range_full, q_range_safe, residuals, slope, intercept, Rg, I0, MW1, MW_alt, MV_P, mw_qlimit_val, MW_qlimit)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Processing failed: {e}")

    def _plot_4_panel(self, q, I, q_range_full, I_range_full, E_range_full, q_range_safe, residuals, slope, intercept, Rg, I0, MW1, MW_alt, MV_P, mw_qlimit_val, MW_qlimit):
        for ax in self.all_4panel_axes: ax.clear()
        
        if self.ax_integral: self.ax_integral.remove(); self.ax_integral = None
        
        axA, axB, axB_res, axC, axD = self.axA, self.axB, self.axB_res, self.axC, self.axD
        c1,c2,c3,c4,c5 = '#0D92F4','#77CDFF','#F95454','#C62E2E','#F3F3E0'
        
        # Form Factor (axA)
        axA.plot(q,I,c=c1); axA.set_yscale('log'); axA.grid(True)
        axA.set_title(f'{self.source_file_prefix} Form Factor', fontsize='small') 
        if self.log_q_checkbox.isChecked(): 
            axA.set_xscale('log'); axA.set_xlabel('log(q) (Å⁻¹)')
        else: 
            axA.set_xscale('linear'); axA.set_xlabel('q (Å⁻¹)')
            # --- ADDED ---
            axA.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
        # Guinier Plot (axB)
        valid_plot_mask = I_range_full > 0
        with np.errstate(invalid='ignore'):
            axB.errorbar(q_range_full[valid_plot_mask]**2, np.log(I_range_full[valid_plot_mask]), yerr=E_range_full[valid_plot_mask]/I_range_full[valid_plot_mask], fmt='o', ms=4, c=c1, label='exp.')
        axB.plot(q_range_safe**2, slope*q_range_safe**2+intercept, c=c4, label='Fit')
        
        axB.set(ylabel='ln(I(q))'); axB.grid(False) 
        
        axB.set_title('Guinier Plot', fontsize='small')
        axB.legend(loc='upper center'); 
        
        # --- MODIFIED: Text moved to bottom-left ---
        axB.text(0.05, 0.05, f'Rg = {Rg:.2f} Å\nI0 = {I0:.2e}\nq·Rg: {q_range_safe[0]*Rg:.2f}–{q_range_safe[-1]*Rg:.2f}', 
                 transform=axB.transAxes, ha='left', va='bottom', 
                 bbox=dict(boxstyle='round', fc=c5, alpha=0.7))
        
        plt.setp(axB.get_xticklabels(), visible=False)
        # --- ADDED ---
        axB.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        # Residuals Plot (axB_res)
        axB_res.plot(q_range_safe**2, residuals, 'o', ms=2, c=c1)
        axB_res.axhline(0,c=c4,lw=1)
        axB_res.set_xlabel('q² (Å⁻²)', fontsize='small') 
        axB_res.set_ylabel('Res/σ', fontsize='small')
        axB_res.tick_params(axis='both', which='major', labelsize='x-small')

        axB_res.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        # --- ADDED ---
        axB_res.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
        axB.xaxis.get_offset_text().set_visible(False)
        if axB_res.xaxis.get_offset_text(): # Check if text exists
            axB_res.xaxis.get_offset_text().set_fontsize('x-small')
        
        # Kratky Plot (axC)
        axC.plot(q*Rg, (q*Rg)**2*I/I0, c=c3); axC.grid(True)
        axC.set(xlabel='q·Rg',ylabel='(q·Rg)²·I(q)/I(0)')
        axC.set_title('Kratky', fontsize='small') 
        axC.set_ylim(0,max(2.0, min(np.max((q*Rg)**2*I/I0)*1.1, 5.0))); axC.axvline(np.sqrt(3),c='grey',ls='--'); axC.axhline(1.1,c='grey',ls='--')
        # --- ADDED ---
        axC.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        
        # VC Plot (axD)
        self.ax_integral=axD.twinx()
        line1, = axD.plot(q,I*q,c=c1, label='q·I(q)')
        line2, = self.ax_integral.plot(q,[simpson(I[:i]*q[:i],x=q[:i]) for i in range(1,len(q)+1)],c=c3, label='Integral')
        axD.set(xlabel='q (Å⁻¹)',ylabel='q·I(q)'); axD.grid(True)
        axD.set_title('Volume of Corr.', fontsize='small') 
        axD.tick_params(axis='y',labelcolor=c1); self.ax_integral.tick_params(axis='y',labelcolor=c3)
        axD.axvline(0.3,c=c1,ls='--',alpha=0.8); axD.axvline(8/Rg,c=c2,ls='--',alpha=0.8)
        
        if mw_qlimit_val is not None:
            axD.axvline(mw_qlimit_val, c='purple', ls=':', alpha=0.8)
        mw_qlimit_line = f'MW(q<{mw_qlimit_val:.2f}): {MW_qlimit:,.0f} Da\n' if mw_qlimit_val is not None and MW_qlimit > 0 else ''
        
        self.ax_integral.legend(handles=[line1, line2], loc='upper left')
        
        # This was already in the bottom-left, so no change needed to position
        axD.text(0.05, 0.05, f'{mw_qlimit_line}MW(q<0.3): {MW1:,.0f} Da\nMW(q<8/Rg): {MW_alt:,.0f} Da\nMW(Porod): {MV_P:,.0f} Da', 
                 transform=axD.transAxes, ha='left', va='bottom', 
                 bbox=dict(boxstyle='round', fc=c5, alpha=0.7), zorder=10)
        
        axD.set_ylim(bottom=0)
        self.ax_integral.set_ylim(bottom=0)
        # --- ADDED ---
        axD.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        self.ax_integral.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        self.figure4panel.tight_layout(); self.canvas4panel.draw()
        # Auto-save plots
        self.figure4panel.savefig(os.path.join(self.output_dir, f"{self.source_file_prefix}_Guinier_Graphs.png"), dpi=300)
        self.figure4panel.savefig(os.path.join(self.output_dir, f"{self.source_file_prefix}_Guinier_Graphs.svg"))

        
    def _save_guinier_data(self, q, I, q_range_safe, I_range_safe, residuals, slope, intercept, Rg, I0, MW1, MW_alt, V_p, MV_P, qmin_Rg, qmax_Rg, r2, mw_qlimit_val, MW_qlimit):
        # Save all the text files associated with the Guinier analysis
        np.savetxt(os.path.join(self.output_dir, f'{self.source_file_prefix}_01_Rg.txt'), 
                   np.c_[q_range_safe**2, np.log(I_range_safe), slope * q_range_safe**2 + intercept, residuals],
                   header="q^2\tln(I_exp)\tln(I_theo)\tnormalized_residuals")
        np.savetxt(os.path.join(self.output_dir, f'{self.source_file_prefix}_02_Norm-Kratky.txt'), 
                   np.c_[q*Rg, (q*Rg)**2*I/I0], header="q*Rg\t(q*Rg)^2*(I(q)/I(0))")
        np.savetxt(os.path.join(self.output_dir, f'{self.source_file_prefix}_03_VC.txt'), 
                   np.c_[q, I * q], header="q\tI(q)*q")
        np.savetxt(os.path.join(self.output_dir, f'{self.source_file_prefix}_04_VC_integral.txt'),
                   np.c_[q, [simpson(I[:i]*q[:i], x=q[:i]) for i in range(1,len(q)+1)]], header="q\tCumulative_Integral")
        with open(os.path.join(self.output_dir, f'{self.source_file_prefix}_05_Summary.txt'), 'w') as f:
            f.write(f"Rg\t{Rg:.4e}\nI0\t{I0:.4e}\n")
            f.write(f"qmin*Rg\t{qmin_Rg:.4e}\nqmax*Rg\t{qmax_Rg:.4e}\n")
            f.write(f"R2_value\t{r2:.4f}\n") 
            if hasattr(self, 'auto_fit_score') and self.auto_fit_score is not None:
                f.write(f"Auto_fit_Score\t{self.auto_fit_score:.2f}\n") 
            f.write(f"MW_Vc(q<0.3)\t{MW1:.4e}\n")
            if mw_qlimit_val is not None and MW_qlimit > 0: 
                f.write(f"MW_Vc(q<{mw_qlimit_val:.4f})\t{MW_qlimit:.4e}\n")
            f.write(f"MW_Vc(q<8/Rg)\t{MW_alt:.4e}\n")
            f.write(f"Vp_Porod\t{V_p:.4e}\nMW_Porod\t{MV_P:.4e}\n")
            f.write(f"Oligomeric State\t{self.oligomer_state_label.text().split(': ')[-1]}\n")
        print(f"✅ Guinier analysis results saved to {self.output_dir}")

    def _run_bift(self):
        if self.data is None: return
        self.run_bift_button.setEnabled(False); self.run_bift_button.setText("Running BIFT...")
        q, i, err = self.data[:, 0], self.data[:, 1], self.data[:, 2]
        
        try:
            dmax_user = float(self.bift_dmax_entry.text()) if self.bift_dmax_entry.text() else None 
            qmin_user = float(self.bift_qmin_entry.text()) if self.bift_qmin_entry.text() else None
            qmax_user = float(self.bift_qmax_entry.text()) if self.bift_qmax_entry.text() else None
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Dmax, q_min, and q_max must be numbers.") 
            self.run_bift_button.setEnabled(True); self.run_bift_button.setText("Run BIFT")
            return

        if dmax_user is None:
            try: 
                rg_val = float(self.rg_label.text().split(' ')[1])
                dmax_user = rg_val * 3
                print(f"Info: Dmax not provided, estimating from Guinier Rg: {dmax_user:.1f} Å")
            except: 
                dmax_user = 150.0 # Default fallback
                print("Warning: Could not get Rg from label, defaulting Dmax to 150.0")
        else:
             print(f"Info: Using user-provided Dmax: {dmax_user:.1f} Å")


        self.bift_worker = PythonBiftThread(q, i, err, dmax_user, qmin_user, qmax_user)
        self.bift_worker.bift_finished.connect(self._on_bift_finished)
        self.bift_worker.bift_error.connect(self._on_bift_error)
        self.bift_worker.start()

    @Slot(dict)
    def _on_bift_finished(self, results):
        self.run_bift_button.setEnabled(True); self.run_bift_button.setText("Run BIFT")
        self.bift_data = results
        print("✅ BIFT analysis successful.")
        self._plot_bift()
        self._save_bift_data()
        self._save_bift_plot()
        self._save_bift_fit_data()
        self._update_summary_with_bift()
        
    @Slot(str)
    def _on_bift_error(self, error_message):
        self.run_bift_button.setEnabled(True); self.run_bift_button.setText("Run BIFT")
        QMessageBox.critical(self, "BIFT Error", error_message)

    def _plot_bift(self):
        if self.bift_data is None: return
        try:
            r, Pr, err = self.bift_data['r'], self.bift_data['pr'], self.bift_data['pr_err']
            q_fit, I_fit = self.bift_data['q'], self.bift_data['fit']
            q_exp, I_exp, E_exp = self.bift_data['q'], self.bift_data['i'], self.bift_data['err']
            details = self.bift_data['details']
            rg, dmax, chisq = details.get('rg'), details.get('dmax'), details.get('chisq')
            
            rg_guinier_str = self.rg_label.text().split(' ')[1]
            rg_guinier = float(rg_guinier_str)
            
            color = 'green' if abs(rg - rg_guinier) / rg_guinier < 0.05 else 'black'
            
            if rg is not None: self.bift_rg_label.setText(f"Rg from P(r): <font color='{color}'>{rg:.2f} Å</font>")
            if dmax is not None: self.bift_dmax_label.setText(f"Dmax from P(r): {dmax:.2f} Å")
            
            self.ax_bift[0].clear(); self.ax_bift[1].clear()
            
            pr_legend = f'Dmax={dmax:.1f} Å, Rg={rg:.2f} Å'
            self.ax_bift[0].errorbar(r, Pr, yerr=err, fmt='-o', color='#0D92F4', ms=3, label=pr_legend)
            self.ax_bift[0].axhline(0, color='grey', ls='--')
            self.ax_bift[0].set(xlabel='r (Å)', ylabel='P(r)')
            self.ax_bift[0].set_title('Pair-Distance Distribution Function', fontsize='small') 
            self.ax_bift[0].grid(True); self.ax_bift[0].legend()

            fit_legend = f'BIFT Fit (χ²={chisq:.2f})'
            # Swapped colors. Data is blue, fit is red.
            self.ax_bift[1].errorbar(q_exp, I_exp, yerr=E_exp, fmt='.', color='#0D92F4', label='Experimental Data', alpha=0.5)
            self.ax_bift[1].plot(q_fit, I_fit, color='#C62E2E', label=fit_legend)
            self.ax_bift[1].set_yscale('log')
            self.ax_bift[1].set_xscale('linear')
            self.ax_bift[1].set(xlabel='q (Å⁻¹)', ylabel='log(I(q))')
            self.ax_bift[1].set_title('BIFT Fit to Data', fontsize='small') 
            self.ax_bift[1].grid(True); self.ax_bift[1].legend()

            self.figure_bift.tight_layout(); self.canvas_bift.draw()
        except Exception as e:
            QMessageBox.critical(self, "Plotting Error", f"Failed to plot BIFT results.\n{e}")

    def _save_bift_plot(self):
        if self.bift_data is None: return
        try:
            filename_png = os.path.join(self.output_dir, f"{self.source_file_prefix}_BIFT_Graphs.png")
            filename_svg = os.path.join(self.output_dir, f"{self.source_file_prefix}_BIFT_Graphs.svg")
            self.figure_bift.savefig(filename_png, dpi=300)
            self.figure_bift.savefig(filename_svg)
            print(f"✅ BIFT plot saved to {self.output_dir}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save BIFT plot.\n{e}")

    def _save_bift_data(self):
        if self.bift_data is None: return
        try:
            filename = os.path.join(self.output_dir, f"{self.source_file_prefix}_pr.txt")
            data_to_save = np.c_[self.bift_data['r'], self.bift_data['pr'], self.bift_data['pr_err']]
            np.savetxt(filename, data_to_save, header="r\tP(r)\terror", fmt="%.6e")
            print(f"✅ P(r) data saved to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save P(r) data.\n{e}")

    def _save_bift_fit_data(self):
        if self.bift_data is None: return
        filename = os.path.join(self.output_dir, f"{self.source_file_prefix}_06_BIFT_fit.txt")
        try:
            data_to_save = np.c_[
                self.bift_data['q'], 
                self.bift_data['i'], 
                self.bift_data['err'], 
                self.bift_data['fit']
            ]
            np.savetxt(filename, data_to_save, header="q\tI_exp\tErr_exp\tI_fit", fmt="%.6e")
            print(f"✅ BIFT fit data saved to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save BIFT fit data.\n{e}")

    def _update_summary_with_bift(self):
        if self.bift_data is None: return
        filename = os.path.join(self.output_dir, f'{self.source_file_prefix}_05_Summary.txt')
        try:
            details = self.bift_data['details']
            dmax, dmaxer = details.get('dmax'), details.get('dmaxer', 0)
            rg, rger = details.get('rg'), details.get('rger', 0)
            with open(filename, 'a') as f:
                f.write("\n# BIFT Results\n")
                f.write(f"Rg_BIFT\t{rg:.4e}\t{rger:.4e}\n")
                f.write(f"Dmax_BIFT\t{dmax:.4e}\t{dmaxer:.4e}\n")
            print(f"✅ Summary file updated with BIFT results.")
        except Exception as e:
            QMessageBox.critical(self, "File Error", f"Failed to update summary file.\n{e}")

    def _reset(self):
        self.data, self.bift_data, self.output_dir, self.source_file_prefix = None, None, None, None
        self.auto_fit_score = None 
        
        for entry in [self.data_file_entry, self.qmin_offset_entry, self.qmax_offset_entry, 
                      self.theoretical_mw_entry, self.mw_qlimit_entry, self.bift_dmax_entry, 
                      self.bift_qmin_entry, self.bift_qmax_entry]: 
            entry.clear()
            
        for label in [self.rg_label, self.i0_label, self.q_range_label, self.q_rg_label, 
                      self.oligomer_state_label, self.bift_rg_label, 
                      self.bift_dmax_label, self.mw_vc_result_label]: 
            label.setText(label.text().split(':')[0] + ':')
            
        self.log_q_checkbox.setChecked(False)
        
        for ax in self.all_4panel_axes: ax.clear()
        
        if self.ax_integral: self.ax_integral.remove(); self.ax_integral = None
        self.canvas4panel.draw()
        for ax in self.ax_bift: ax.clear()
        self.canvas_bift.draw()
        self.run_bift_button.setEnabled(False)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())