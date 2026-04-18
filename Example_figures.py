import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import glob
import re
import math

G = 6.674*1e-11
mA = 8500
d = 5.02
r0 = 314159*1e3
v0 = 1.2*1e3
h = 10000
mT = 5.972e24
mL = 7.3477e22*1e-14
rho_0 = 1.2
R_T =6378.1e3
lambd = 7238.2
Cx = 0.3
dTL = 384748e3
R_L = 1737.4e3
    
vmax = np.sqrt(v0**2 + 2*G*mT*((1/(h+R_T))-(1/r0)))
vt = (vmax*(h+R_T))/r0
vr = np.sqrt(v0**2-vt**2)

#*outputFile << t << " " << y[ix(0)] << " " << y[iy(0)] << " " << y[ivx(0)] << " " << y[ivy(0)] <<" " << emec << " " << pnc <<" "<<quantite_mvmt[0]<<" "<<quantite_mvmt[1]<<endl;

# ============================================================
# USER SETTINGS
# ============================================================

folder = r"Scan_nsteps_pendulum_x1_3.1e+08_y1_0_vx1_-1.2e+03_vy1_2.3e+02"

plot_layout = {
    "theta_time": True,
    "x_y": True,
    "energy": True,
    "real_space": False,
    "power": True,
    "vmax" : True,
    "energy_balance": True
}

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder, "*.txt")))

datasets = []
param_values = []
param_name = None

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

tmin = min(data[:,0].min() for data in datasets)
tmax = max(data[:,0].max() for data in datasets)

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



# ============================================================
# Plot 1 : hmin vs dt
# ============================================================

fig, axes = get_axes("x_y", "hmin")

lc_ref = None

hmin = []
vmaxx = []

for i,data in enumerate(datasets):

    x = data[:,1] 
    y = data[:,2]
    vx = data[:,3] 
    vy = data[:,4]
    
    h = 10000
    hneu = np.sqrt(x**2+y**2) -R_T
    vneu = np.sqrt(vx**2+vy**2)
    
    hmin.append(np.min(hneu))
    vmaxx.append(np.max(vneu))

axes[0].plot(param_values, hmin, color='r', linestyle='-', marker = 'x', label="hmin numérique")
plt.axhline(h, color='k', linestyle='--', label="hmin exact")




axes[0].set_ylabel("hmin")
axes[0].set_xlabel(f"{param_name}")
axes[0].grid()



fig.savefig(os.path.join(fig_dir,"hmin.png"), dpi=300)

h_err = np.abs(1 - np.array(hmin) / h)

plt.figure()
plt.plot(param_values**(-1), h_err, 'r+-', label="hmin numérique")
plt.loglog(param_values**(-1), param_values**(-1), 'k--', label="O(dt)")
plt.loglog(param_values**(-1), param_values**(-2), 'k-.', label="O(dt^2)")
plt.xlabel(r"dt")
plt.ylabel(r"hmin")
plt.xscale('log')
#plt.ylim(0, tf/10)  # Set y-limits to focus on the relevant range
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,"hmin_conv.png"), dpi=300)


# ============================================================
# Plot 2 : vmax vs dt
# ============================================================

plt.figure()

plt.plot(param_values, vmaxx, color='r', linestyle='-', marker = 'x', label="vmax numérique")
plt.axhline(vmax, color='k', linestyle='--', label="hmin exact")




plt.ylabel("vmax")
plt.xlabel(f"{param_name}")
plt.grid()



plt.savefig(os.path.join(fig_dir,"vmax.png"), dpi=300)

v_err = np.abs(1 - np.array(vmaxx) / vmax)

plt.figure()
plt.plot(param_values**(-1), v_err, 'r+-', label="vmax numérique")
plt.loglog(param_values**(-1), param_values**(-1), 'k--', label="O(dt)")
plt.loglog(param_values**(-1), param_values**(-2), 'k-.', label="O(dt^2)")
plt.xlabel(r"dt")
plt.ylabel(r"vmax")
plt.xscale('log')
#plt.ylim(0, tf/10)  # Set y-limits to focus on the relevant range
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,"vmax_conv.png"), dpi=300)
# ============================================================
# Plot 3 : energy
# ============================================================

fig, axes = get_axes("energy", "Energy evolution")

for i,data in enumerate(datasets):

    t = data[:,0]
    E = data[:,3]

    color = cmap(i % 10)

    axes[i].plot(t, E, color=color,
                 label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Energy")
    axes[i].grid()

    if not plot_layout["energy"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["energy"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"energy_all.png"), dpi=300)


# ============================================================
# Plot 4 : real space trajectory
# ============================================================

L = 0.2

fig, axes = get_axes("real_space", "Real space trajectory")

lc_ref = None

for i,data in enumerate(datasets):

    t = data[:,0]
    x = data[:,1]
    y = data[:, 2]

    lc = colored_line(x, y, t, axes[i], tmin, tmax)

    axes[i].set_xlabel("x")
    axes[i].set_ylabel("y")

    if not plot_layout["real_space"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

    if lc_ref is None:
        lc_ref = lc

cbar = fig.colorbar(lc_ref, ax=axes)
cbar.set_label("time")

fig.savefig(os.path.join(fig_dir,"real_space_all.png"), dpi=300)


# ============================================================
# Plot 5 : power
# ============================================================

fig, axes = get_axes("power", "Non-conservative power")

for i,data in enumerate(datasets):

    t = data[:,0]
    Pnc = data[:,4]

    color = cmap(i % 10)

    axes[i].plot(t, Pnc, color=color,
                 label=f"{param_name}={param_values[i]}")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Pnc")
    axes[i].grid()

    if not plot_layout["power"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")

if plot_layout["power"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"power_all.png"), dpi=300)


# ============================================================
# Plot 6 : energy balance
# ============================================================

fig, axes = get_axes("energy_balance", "Energy balance")

for i,data in enumerate(datasets):

    color = cmap(i % 10)

    t = data[:,0]
    E = data[:,3]
    Pnc = data[:,4]

    dEdt = np.gradient(E,t)

    axes[i].plot(t, dEdt, linestyle='-', color=color, label="dE/dt")
    axes[i].plot(t, Pnc, linestyle='--', color=color, label="Pnc")

    axes[i].set_xlabel("t")
    axes[i].set_ylabel("Power")
    axes[i].grid()

    if not plot_layout["energy_balance"]:
        axes[i].set_title(f"{param_name} = {param_values[i]}")
        axes[i].legend()

if plot_layout["energy_balance"]:
    axes[0].legend()

fig.savefig(os.path.join(fig_dir,"energy_balance_all.png"), dpi=300)


plt.show()
