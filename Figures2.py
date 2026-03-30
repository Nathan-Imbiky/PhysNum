import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math
from scipy.signal import find_peaks

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"/home/imbiky/Desktop/MyFiles/PhysNum/Exercise2_student/rotatingpendulum/problème/Scan_r_pendulum_kappa_0.01_r_1_Omega_7" #n0 outputs ?

plot_layout = {
    "theta_time": False,
    "phase_space": False,
    "energy": False,
    "real_space": False,
    "power": True,
    "energy_balance": False,
    "dt":True,
    "delta_theta":True,
    "Poincare":False
}

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures") #n1 pathjoin
os.makedirs(fig_dir, exist_ok=True) #n2 exist_ok

# ============================================================
# Scan files
# ============================================================
files = sorted(glob.glob(os.path.join(folder, "*.txt")))

datasets = []
param_values = []
param_name = None #n3 comment ça none

for f in files:

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".txt"

    parts = name.split("_")

    param_name = parts[-2]          # scanned parameter
    value = float(parts[-1])        # parameter value

    data = np.loadtxt(f)

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Helper: colored line
# ============================================================

def colored_line(x, y, t, ax, vmin=None, vmax=None):

    points = np.array([x, y]).T.reshape(-1,1,2) #n4 comment ça -1 eft jcomprends pas la question
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='viridis')
    lc.set_array(t)

    if vmin is not None and vmax is not None:
        lc.set_clim(vmin, vmax)

    lc.set_linewidth(2)

    ax.add_collection(lc)
    ax.autoscale()

    return lc

# ============================================================
# Shared time range for color maps
# ============================================================

tmin = min(data[:,0].min() for data in datasets)
tmax = max(data[:,0].max() for data in datasets)

# ============================================================
# Axis layout helper
# ============================================================

def get_axes(plot_key, title): #n5 rien compris

    if plot_layout[plot_key]:

        fig, ax = plt.subplots()
        axes = [ax]*len(datasets)

    else:

        n = len(datasets)
        ncols = min(3, n)
        nrows = math.ceil(n/3)

        fig, axarr = plt.subplots(nrows, ncols,
                                  figsize=(5*ncols,4*nrows))

        axes = np.array(axarr).reshape(-1)

        for j in range(n, len(axes)):
            fig.delaxes(axes[j])

        axes = axes[:n]

    fig.suptitle(title)

    return fig, axes


cmap = plt.get_cmap("tab10") #n6 c qui la colored map ptn
plt.rcParams['font.size'] = 12
plt.rcParams['legend.fontsize'] = 1

# ============================================================
# Plot 9 : Attracteurs
# ============================================================

fig, axes = get_axes("Poincare", "Poincaré section")

lc_ref = None

for i,data in enumerate(datasets):
    
    t = data[:,0]
    #theta = data[:,1]
    theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi
    thetadot = data[:,2]
    
    if(param_name!="Omega"):
        T = 2*np.pi*np.sqrt(9.81/0.2)
        
    else :
        T = 2*np.pi/param_values[i]
    N = int(t[-1]/T)
   
    
    lex = T*np.linspace(0, N, N)
    
    thetaNeo = np.interp(lex, t, theta)
    #, period=T)
    thetaNeo = (thetaNeo + np.pi)%(2*np.pi) -np.pi
    
    thetadotNeo = np.interp(lex, t, thetadot)
    #, period=T)
    
    for k in range(len(lex)) :
        axes[i].plot(thetaNeo[k], thetadotNeo[k], markersize=2, marker='x', color='red')

    #lc = colored_line(thetaNeo, thetadotNeo, t, axes[i], tmin, tmax)
    
    #axes[i].plot(t, theta, color='red', linestyle=None, marker="x")

    axes[i].set_xlabel("theta")
    axes[i].set_ylabel("thetadot")
    axes[i].grid()

    if not plot_layout["Poincare"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")
        marsh =0

    #if lc_ref is None:
        #lc_ref = lc

#cbar = fig.colorbar(lc_ref, ax=axes)
#cbar.set_label("time")

fig.savefig(os.path.join(fig_dir,"poincare.png"), dpi=300)


plt.show()
