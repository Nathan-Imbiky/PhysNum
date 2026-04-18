import numpy as np
import subprocess
import os

# Parameters
repertoire = ''
executable = './engineexo4'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp

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


input_parameters = {
    'tf': 2, # t final (overwritten if N >0)
    'nsteps': 10000, # number of time steps per period (if N>0), number of timesteps total if N=0
    'G' : 6.674*1e-11, 
    'mA' :8500, 
    'd' :5.02, 
    'r0' :314159*1e3, 
    'v0' :1.2*1e3, 
    'h' :10000, 
    'mT' :5.972e24, 
    'mL' :7.3477e22*1e-14, 
    'rho_0' :1.2, 
    'R_T' :6378.1e3, 
    'lambda' :7238.2, 
    'Cx' :0.3,
    'dTL':384748e3,
    'R_L' :1737.4e3,
    
    'x1': r0,
    'x2':0.0,
    'x3':dTL,
    'y1': 0,
    'y2':0.0,
    'y3':0,
    'vx1': vr,
    'vx2':0.0,
    'vx3':0.0,
    'vy1': vt,
    'vy2':0.0,
    'vy3':0.0,
    'sampling': 1
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'nsteps' # The parameter to scan, must be one of the keys in input_parameters

variable_array = 2**np.arange(3, 15)  # Example values for the parameter scan


outstr = f"pendulum_x1_{input_parameters['x1']:.2g}_y1_{input_parameters['y1']:.2g}_vx1_{input_parameters['vx1']:.2g}_vy1_{input_parameters['vy1']:.2g}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)


for i in range(len(variable_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = variable_array[i]

    output_file = f"{outstr}_{paramstr}_{variable_array[i]}.txt"
    output_path = os.path.join(outdir, output_file)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")

