import numpy as np
import subprocess
import os

# Parameters
repertoire = ''
executable = './engineexo2'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'tf': 2, # t final (overwritten if N >0)
    'N': 7000, # number of excitation periods
    'nsteps': 256, # number of time steps per period (if N>0), number of timesteps total if N=0
    'r': 1,
    'kappa': 0.01,
    'm': 0.1,
    'L': 0.2,
    'g': 9.81,
    'Omega': np.sqrt(9.81/0.2),
    'theta0': 0,
    'thetadot0': 0.,
    'sampling': 1
}

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'r' # The parameter to scan, must be one of the keys in input_parameters

delta = 1e-10

variable_array = np.array([0.01, 0.08, 0.1, 0.12, 0.5, 1, 2, 5, 10])
#0*np.array([1])
#np.linspace(1e-8, np.pi-1e-2, 6)
#np.array([1, 2, 5, 7, 10])
#2**np.arange(15, 20)
#2**np.arange(3, 15)  # Example values for the parameter scan


outstr = f"pendulum_kappa_{input_parameters['kappa']:.2g}_r_{input_parameters['r']:.2g}_Omega_{input_parameters['Omega']:.2g}"

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

