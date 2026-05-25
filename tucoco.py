import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
 
# -----------------------------------------------------------------------
# Script d'analyse – Exercice 5.3 : Tsunami sur la grande barrière
# -----------------------------------------------------------------------
 
repertoire     = ''
executable     = './enginecl'
input_filename = 'configuration.in'
 
# Paramètres physiques
g      = 9.81
L      = 800.0e3   # 800 km
hocean = 8000.0    # m
hrecif = 20.0      # m
xa     = 200.0e3
xb     = 370.0e3
xc     = 430.0e3
xd     = 600.0e3
A      = 1.0       # amplitude 1 m
T      = 15 * 60   # période 15 minutes en secondes
om     = 2*np.pi/T # fréquence angulaire
 
input_parameters = {
    'L'             : L,
    'nx'            : 1000,
    'CFL'           : 0.9,
    'tfin'          : 3 * L / np.sqrt(g * hocean),  # 3 transits océan
    'h00'           : hocean,
    'v_uniform'     : 'false',
    'hL'            : hocean,
    'hR'            : hrecif,
    'xa'            : xa,
    'xb'            : xb,
    'xc'            : xc,
    'xd'            : xd,
    'cb_gauche'     : 'harmonique',
    'cb_droite'     : 'sortie',
    'A'             : A,
    'om'            : om,
    'ecrire_f'      : 'true',
    'n_stride'      : 50,
    'impose_nsteps' : 'false',
    'nsteps'        : 100,
}
 
def lancer(output, params):
    param_string = " ".join(f"{k}={v}" for k, v in params.items())
    cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output}"
    subprocess.run(cmd, shell=True, capture_output=True)
 
def supprimer(output):
    for ext in ['_x', '_v', '_f', '_en']:
        fichier = output + ext
        if os.path.exists(fichier):
            os.remove(fichier)
 
 
# =======================================================================
# (a) VÉRIFICATION VITESSE ET AMPLITUDE
# =======================================================================
print("=== (a) Vitesse et amplitude ===")
 
for eq in ['A', 'B', 'C']:
    params = input_parameters.copy()
    params['equation_type'] = eq
    lancer(f"tsunami_{eq}", params)
 
# Lire les données
x       = np.loadtxt("tsunami_A_x")
vel2    = np.loadtxt("tsunami_A_v")
u_x     = np.sqrt(vel2)          # vitesse u(x) = sqrt(g*h(x))
h_x     = vel2 / g               # profondeur h(x)
 
# -----------------------------------------------------------------------
# Figure 1 : profil de profondeur et vitesse
# -----------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 6))
axes[0].plot(x/1e3, h_x)
axes[0].set_xlabel("x (km)")
axes[0].set_ylabel("h(x) (m)")
axes[0].set_title("Profil de profondeur")
axes[0].grid(True, ls=':')
 
axes[1].plot(x/1e3, u_x/1e3)
axes[1].set_xlabel("x (km)")
axes[1].set_ylabel("u(x) (km/s)")
axes[1].set_title("Vitesse de propagation u(x) = √(g·h(x))")
axes[1].grid(True, ls=':')
plt.tight_layout()
plt.show()
 
# -----------------------------------------------------------------------
# Vitesse de propagation : trouver le temps du maximum pour chaque x
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
 
for eq in ['A', 'B', 'C']:
    data  = np.loadtxt(f"tsunami_{eq}_f")
    temps = data[:, 0]
    F     = data[:, 1:]
 
    # Pour chaque point x, trouver le premier temps où f est maximum
    t_crete = np.zeros(len(x))
    for i in range(len(x)):
        t_crete[i] = temps[np.argmax(F[:, i])]
 
    # Vitesse numérique = dx/dt de la crête
    # On calcule la vitesse entre points consécutifs
    dx_arr = np.diff(x)
    dt_arr = np.diff(t_crete)
    # Éviter division par zéro
    mask = dt_arr > 0
    v_num = np.zeros(len(x)-1)
    v_num[mask] = dx_arr[mask] / dt_arr[mask]
 
    ax.plot(x[:-1]/1e3, v_num/1e3, label=f"Équation {eq} (numérique)", lw=1)
 
ax.plot(x/1e3, u_x/1e3, 'k--', label="u(x) = √(g·h(x)) (théorique)", lw=2)
ax.set_xlabel("x (km)")
ax.set_ylabel("Vitesse (km/s)")
ax.set_title("Vitesse de propagation de la crête")
ax.legend()
ax.grid(True, ls=':')
ax.set_ylim(0, 0.35)
plt.tight_layout()
plt.show()
 
# -----------------------------------------------------------------------
# Amplitude : valeur de f à la crête pour chaque x
# -----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
 
exposants = {'A': 1/4, 'B': -1/4, 'C': -3/4}
 
for eq in ['A', 'B', 'C']:
    data  = np.loadtxt(f"tsunami_{eq}_f")
    temps = data[:, 0]
    F     = data[:, 1:]
 
    # Amplitude = valeur max de f pour chaque x
    amplitude = np.max(np.abs(F), axis=0)
 
    # Normaliser par la valeur au premier point
    if amplitude[0] > 0:
        amplitude = amplitude / amplitude[0]
 
    # Prédiction WKB normalisée
    h_norm = h_x / h_x[0]
    wkb    = h_norm ** exposants[eq]
 
    ax.plot(x/1e3, amplitude, label=f"Équation {eq} (numérique)")
    ax.plot(x/1e3, wkb, '--', label=f"WKB : h^({exposants[eq]:.2f})")
 
ax.set_xlabel("x (km)")
ax.set_ylabel("Amplitude normalisée")
ax.set_title("Amplitude de la vague en fonction de x")
ax.legend(fontsize=8)
ax.grid(True, ls=':')
plt.tight_layout()
plt.show()
 
# Nettoyer les fichiers
for eq in ['A', 'B', 'C']:
    supprimer(f"tsunami_{eq}")
 
# =======================================================================
# (b) BORDS DE PLUS EN PLUS RAIDES (cas B)
# =======================================================================
print("\n=== (b) Bords raides – cas B ===")
 
# xa s'approche de xb : bords de plus en plus raides
xa_vals = [200.0e3, 300.0e3, 340.0e3, 360.0e3, 368.0e3]
 
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
 
for xa_i in xa_vals:
    params_b = input_parameters.copy()
    params_b['equation_type'] = 'B'
    params_b['xa']            = xa_i
    label = f"xa = {xa_i/1e3:.0f} km"
    lancer(f"raide_{xa_i:.0f}", params_b)
 
    data_b  = np.loadtxt(f"raide_{xa_i:.0f}_f")
    F_b     = data_b[:, 1:]
    amplitude_b = np.max(np.abs(F_b), axis=0)
    if amplitude_b[0] > 0:
        amplitude_b = amplitude_b / amplitude_b[0]
 
    axes[0].plot(x/1e3, amplitude_b, label=label)
    supprimer(f"raide_{xa_i:.0f}")
 
# WKB cas B pour référence
h_norm = h_x / h_x[0]
wkb_B  = h_norm ** (-1/4)
axes[0].plot(x/1e3, wkb_B, 'k--', lw=2, label="WKB : h^(-1/4)")
 
axes[0].set_xlabel("x (km)")
axes[0].set_ylabel("Amplitude normalisée")
axes[0].set_title("Amplitude pour différentes raideurs (cas B)")
axes[0].legend(fontsize=8)
axes[0].grid(True, ls=':')
 
# Largeur de la barrière = xb - xa
largeurs = [(xb - xa_i)/1e3 for xa_i in xa_vals]
# Amplitude max sur le récif pour chaque xa
amp_recif = []
for xa_i in xa_vals:
    params_b = input_parameters.copy()
    params_b['equation_type'] = 'B'
    params_b['xa']            = xa_i
    lancer(f"raide2_{xa_i:.0f}", params_b)
    data_b = np.loadtxt(f"raide2_{xa_i:.0f}_f")
    F_b    = data_b[:, 1:]
    amp    = np.max(np.abs(F_b), axis=0)
    # Prendre la valeur sur le récif (entre xb et xc)
    mask_recif = (x >= xb) & (x <= xc)
    amp_recif.append(np.mean(amp[mask_recif]) / amp[0])
    supprimer(f"raide2_{xa_i:.0f}")
 
axes[1].plot(largeurs, amp_recif, 'o-')
axes[1].set_xlabel("Largeur de la montée xb - xa (km)")
axes[1].set_ylabel("Amplitude moyenne sur le récif")
axes[1].set_title("Effet de la raideur sur l'amplitude")
axes[1].grid(True, ls=':')
 
plt.suptitle("Bords de plus en plus raides – cas B", fontsize=12)
plt.tight_layout()
plt.show()
 
print("\nDone !")
