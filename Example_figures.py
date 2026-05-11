import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math


R=0.1
# ============================================================
# USER SETTINGS
# ============================================================

folder = r"Scan_N1_electrostatics_b_0.05_R_0.1_trivial_true"

plot_layout = {
    "phi_r": True,
    "Er_r": True,
}

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.out")))

datasets = []
datasetsphi = []
datasetsErDr = []
datasetsdivD = []
param_values = []
param_name = None

for f in files:

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".txt"

    parts = name.split("_")
    
    data = np.loadtxt(f)

     
    if(parts[-1]=="rho") :
        value = float(parts[-3])  # parameter value 
        datasetsdivD.append(data)
        param_values.append(value)      # scanned parameter
        param_name = parts[-4] 
    
    elif(parts[-1]=="phi") :
       datasetsphi.append(data)      

    else :
       datasetsErDr.append(data)  
       
    datasets.append(data)


param_values = np.sort(param_values)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]

# ============================================================
# Helper: colored line
# ============================================================

def colored_line(x, y, t, ax, vmin=None, vmax=None):

    points = np.array([x, y]).T.reshape(-1,1,2)
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

tmin = min([data[:,0].min() for data in datasetsphi])
tmax = max(data[:,0].max() for data in datasetsphi)

# ============================================================
# Axis layout helper
# ============================================================

def get_axes(plot_key, title):

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


cmap = plt.get_cmap("tab10")

#//////////////////////hmin et vmax, pmax et accmax, dTL p et em, x vs y

plt.rcParams.update({'font.size': 16})
plt.rcParams['figure.figsize'] = [9,9]


# ============================================================
# Plot 1 : phi vs r
# ============================================================

fig, axes = get_axes("phi_r", "phi(r)")
phi0 =  []
print(param_values)

for i,data in enumerate(datasetsphi):

    phi = data[:,1] 
    r = data[:,0]
    phi0.append(phi[0])


    
   	 
    color = cmap(i % 10)

    axes[i].plot(r, phi, color=color, label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("r")
    axes[i].set_ylabel("phi")
    axes[i].grid()

    if not plot_layout["phi_r"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["phi_r"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"phi_r.png"), dpi=300, bbox_inches='tight')

# ============================================================
# Plot 2 : Er vs r
# ============================================================

fig, axes = get_axes("Er_r", "Er(r)")


for i,data in enumerate(datasetsErDr):

    Er = data[:,1] 
    r = data[:,0]

    
   	 
    color = cmap(i % 10)

    axes[i].plot(r, Er, color=color, label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("r")
    axes[i].set_ylabel("Er")
    axes[i].grid()

    if not plot_layout["Er_r"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["Er_r"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"Er_r.png"), dpi=300, bbox_inches='tight')

# ============================================================
# Plot 3 : phi0 vs N
# ============================================================


lc_ref = None

phi0_ex = R*R/4

plt.figure()
plt.plot(param_values, phi0, color='r', linestyle='-', marker = 'x', label="phi0 numérique")
plt.axhline(phi0_ex, color='k', linestyle='--', label="phi0 exact")




plt.ylabel("phi0")
plt.xlabel(f"N")
plt.grid()


plt.savefig(os.path.join(fig_dir,"phi0_N.png"), dpi=300, bbox_inches='tight')

phi0_err = np.abs(1 - np.array(phi0) / phi0_ex)

plt.figure()
plt.plot(param_values**(-1), phi0_err, 'r+-', label="phi0 numérique")
plt.loglog(param_values**(-1), param_values**(-1), 'k--', label=f"O(N^{-1})")
plt.loglog(param_values**(-1), param_values**(-2), 'k-.', label=f"O(N^{-2})")
plt.xlabel(f"N")
plt.ylabel(r"phi0")
plt.xscale('log')
#plt.ylim(0, tf/10)  # Set y-limits to focus on the relevant range
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,"phi0_conv.png"), dpi=300, bbox_inches='tight')

plt.tight_layout(pad=2.0)
plt.show()





