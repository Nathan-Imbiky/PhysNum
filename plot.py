"""
Plots pour 4.3(b)(c)(d) — FEM elements finis
Prerequis : avoir compile et execute engine avec trivial.in et nontrivial.in
  ./engine trivial.in
  ./engine nontrivial.in
"""
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# ---------------------------------------------------------------
# Parametres physiques
# ---------------------------------------------------------------
R  = 0.1
b_trivial    = 0.05
b_nontrivial = 0.02
a0 = 1e4
eps0 = 8.854187817e-12

# ---------------------------------------------------------------
# Solutions analytiques pour le cas uniforme (trivial)
# eps_r=1, rho_lib=eps0 (rho/eps0=1), V0=0
# phi(r) = (R^2 - r^2)/4
# E_r(r) = r/2
# ---------------------------------------------------------------
def phi_uniform(r):
    return (R**2 - r**2) / 4.0

def Er_uniform(r):
    return r / 2.0

# ---------------------------------------------------------------
# Lecture des fichiers de sortie
# ---------------------------------------------------------------
def load(prefix):
    phi_d   = np.loadtxt(prefix + '_phi.out')
    erdr_d  = np.loadtxt(prefix + '_ErDr.out')
    divd_d  = np.loadtxt(prefix + '_divD_rho.out')
    return phi_d, erdr_d, divd_d

# ---------------------------------------------------------------
# 4.3(b) - cas uniforme : comparaison et convergence
# ---------------------------------------------------------------
def plot_trivial():
    phi_d, erdr_d, divd_d = load('trivial')
    r   = phi_d[:,0];   phi = phi_d[:,1]
    rm  = erdr_d[:,0];  Er  = erdr_d[:,1]; Dr = erdr_d[:,2]
    rmm = divd_d[:,0];  dD  = divd_d[:,1]; ro = divd_d[:,2]
    r_f = np.linspace(0, R, 1000)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Cas uniforme (trivial=true, N=5)', fontsize=13)

    axes[0].plot(r*100, phi, 'bo-', ms=5, label='FEM')
    axes[0].plot(r_f*100, phi_uniform(r_f), 'r--', lw=2, label='Exact')
    axes[0].set_xlabel('r [cm]'); axes[0].set_ylabel('φ [V]')
    axes[0].set_title('Potentiel'); axes[0].legend(); axes[0].grid(True)

    axes[1].plot(rm*100, Er, 'bo-', ms=5, label='FEM $E_r$')
    axes[1].plot(rm*100, Dr, 'g^-', ms=5, label='FEM $D_r/\\varepsilon_0$')
    axes[1].plot(r_f*100, Er_uniform(r_f), 'r--', lw=2, label='Exact $E_r$')
    axes[1].set_xlabel('r [cm]'); axes[1].set_ylabel('[V/m]')
    axes[1].set_title('Champ et déplacement'); axes[1].legend(); axes[1].grid(True)

    axes[2].plot(rmm*100, dD, 'bs-', ms=5, label='$\\nabla\\cdot D/\\varepsilon_0$')
    axes[2].plot(rmm*100, ro, 'r--', lw=2, label='$\\rho_{lib}/\\varepsilon_0$')
    axes[2].set_xlabel('r [cm]')
    axes[2].set_title('Vérification $\\nabla\\cdot D = \\rho_{lib}$')
    axes[2].legend(); axes[2].grid(True)

    plt.tight_layout()
    plt.savefig('trivial_solution.png', dpi=150)
    plt.show()

# ---------------------------------------------------------------
# 4.3(b)(ii) - etude de convergence de phi(0)
# Relance le solver avec differents N
# ---------------------------------------------------------------
def convergence_study(engine_path='./engine', config='trivial.in'):
    Ns   = [2, 4, 8, 16, 32, 64, 128]
    phi0 = []
    phi_exact_0 = phi_uniform(0.0)

    for N in Ns:
        cmd = [engine_path, config, f'N1={N}', f'N2={N}', 'output=conv_tmp']
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode != 0:
            print(f"Erreur pour N={N}: {result.stderr.decode()}")
            continue
        data = np.loadtxt('conv_tmp_phi.out')
        phi0.append(data[0, 1])    # phi at r=0 (first node)

    phi0 = np.array(phi0)
    err  = np.abs(phi0 - phi_exact_0)
    hs   = [b_trivial/N for N in Ns[:len(phi0)]]

    print(f"\nConvergence de phi(0) — valeur exacte = {phi_exact_0:.6f} V")
    for N, e, h in zip(Ns, err, hs):
        print(f"  N={N:4d}  h={h:.4f}  phi(0)={phi0[Ns.index(N)]:.8f}  err={e:.3e}")

    # Ordre de convergence
    if len(err) > 2:
        orders = [np.log(err[i]/err[i+1])/np.log(hs[i]/hs[i+1])
                  for i in range(len(err)-1)]
        print(f"  Ordres locaux de convergence: {[f'{o:.2f}' for o in orders]}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(hs, err, 'bo-', ms=6, label='Erreur $|\\phi(0)|$')
    ref2 = err[2] * (np.array(hs)/hs[2])**2
    ax.loglog(hs, ref2, 'k--', label='Pente 2')
    ax.set_xlabel('Pas h [m]'); ax.set_ylabel('Erreur $|\\phi(0)|$ [V]')
    ax.set_title('Convergence FEM — cas uniforme'); ax.legend(); ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig('convergence_phi0.png', dpi=150)
    plt.show()

    # Nettoyer les fichiers temporaires
    for ext in ['_phi.out', '_ErDr.out', '_divD_rho.out']:
        try: os.remove('conv_tmp' + ext)
        except: pass

# ---------------------------------------------------------------
# 4.3(c) - cas non-trivial
# ---------------------------------------------------------------
def plot_nontrivial():
    phi_d, erdr_d, divd_d = load('nontrivial')
    r   = phi_d[:,0];   phi = phi_d[:,1]
    rm  = erdr_d[:,0];  Er  = erdr_d[:,1]; Dr = erdr_d[:,2]
    rmm = divd_d[:,0];  dD  = divd_d[:,1]; ro = divd_d[:,2]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Cas non-trivial (b=2 cm, R=10 cm, N1=N2=20)', fontsize=13)

    axes[0,0].plot(r*100, phi, 'b-o', ms=3)
    axes[0,0].axvline(b_nontrivial*100, color='gray', ls='--', alpha=0.5, label='r=b')
    axes[0,0].set_xlabel('r [cm]'); axes[0,0].set_ylabel('φ [V]')
    axes[0,0].set_title('Potentiel φ(r)'); axes[0,0].legend(); axes[0,0].grid(True)

    axes[0,1].plot(rm*100, Er, 'b-o', ms=3, label='$E_r$')
    axes[0,1].axvline(b_nontrivial*100, color='gray', ls='--', alpha=0.5, label='r=b')
    axes[0,1].set_xlabel('r [cm]'); axes[0,1].set_ylabel('[V/m]')
    axes[0,1].set_title('Champ $E_r(r)$'); axes[0,1].legend(); axes[0,1].grid(True)

    axes[1,0].plot(rm*100, Dr, 'g-^', ms=3, label='$D_r/\\varepsilon_0$')
    axes[1,0].axvline(b_nontrivial*100, color='gray', ls='--', alpha=0.5, label='r=b')
    axes[1,0].set_xlabel('r [cm]'); axes[1,0].set_ylabel('[V/m]')
    axes[1,0].set_title('Déplacement $D_r/\\varepsilon_0$'); axes[1,0].legend(); axes[1,0].grid(True)

    axes[1,1].plot(rmm*100, dD, 'bs-', ms=3, label='$\\nabla\\cdot D/\\varepsilon_0$')
    axes[1,1].plot(rmm*100, ro, 'r--', lw=2, label='$\\rho_{lib}/\\varepsilon_0$')
    axes[1,1].axvline(b_nontrivial*100, color='gray', ls='--', alpha=0.5, label='r=b')
    axes[1,1].set_xlabel('r [cm]')
    axes[1,1].set_title('Vérification $\\nabla\\cdot D = \\rho_{lib}$')
    axes[1,1].legend(); axes[1,1].grid(True)

    plt.tight_layout()
    plt.savefig('nontrivial_solution.png', dpi=150)
    plt.show()

# ---------------------------------------------------------------
# 4.3(d)(ii) - theoreme de Gauss : charges totales
# ---------------------------------------------------------------
def gauss_check(prefix='nontrivial', b=b_nontrivial):
    erdr_d = np.loadtxt(prefix + '_ErDr.out')
    rm = erdr_d[:,0]; Dr = erdr_d[:,1]; # Dr ici c'est D_r/eps0

    # Charge libre totale : Q_lib/eps0 = int_0^R rho_lib/eps0 * 2pi*r*Lz dr
    # Par Gauss sur toute la section: Q/eps0 = 2*pi*Lz * r*D_r|_{r=R}
    # On prend Lz=1 m (par unite de longueur)
    Dr_at_R  = erdr_d[-1, 2]   # D_r/eps0 a r=R (dernier milieu)
    r_at_R   = erdr_d[-1, 0]
    Q_total_gauss = 2*np.pi * r_at_R * Dr_at_R  # Q_lib/(eps0*Lz)
    print(f"\n--- Theoreme de Gauss ---")
    print(f"  D_r(R)/eps0 (dernier midpoint) = {Dr_at_R:.4f} V/m")
    print(f"  Q_libre/(eps0*Lz) [Gauss surface r=R] = {Q_total_gauss:.4f} V")

    # Charge libre numerique par integration directe
    # Q_lib/eps0 = 2*pi * sum_k rho(r_mid_k)/eps0 * r_mid_k * h_k
    phi_d  = np.loadtxt(prefix + '_phi.out')
    r_grid = phi_d[:,0]
    Npts   = len(r_grid)
    Nint   = Npts - 1
    Q_lib  = 0.0
    for k in range(Nint):
        hk     = r_grid[k+1] - r_grid[k]
        rmid_k = 0.5*(r_grid[k] + r_grid[k+1])
        Q_lib += rho_at(rmid_k, b) * rmid_k * hk
    Q_lib *= 2*np.pi
    print(f"  Q_libre/(eps0*Lz) [integration directe] = {Q_lib:.4f} V")

    # Charge de polarisation a r=b : saut de D_r
    # Trouver l'index du midpoint le plus proche de b de chaque cote
    rm = erdr_d[:,0]; Dr_arr = erdr_d[:,2]
    idx_in  = np.searchsorted(rm, b) - 1   # dernier midpoint dans r < b
    idx_out = idx_in + 1                   # premier midpoint dans r > b
    if idx_in >= 0 and idx_out < len(rm):
        D_in  = rm[idx_in]  * Dr_arr[idx_in]
        D_out = rm[idx_out] * Dr_arr[idx_out]
        # Saut: sigma_pol/eps0 = -(D_out - D_in)/b  (discontinuite normale)
        sigma_pol = -(Dr_arr[idx_out] - Dr_arr[idx_in])
        print(f"  Saut de D_r/eps0 a r=b : Delta(D_r/eps0) = {Dr_arr[idx_out]-Dr_arr[idx_in]:.4f}")
        print(f"  Charge de polarisation sigma_pol/eps0 ~ {sigma_pol:.4f} V/m")

def rho_at(r, b):
    if r <= b:
        return a0 * np.sin(np.pi * r / b)
    return 0.0

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
if __name__ == '__main__':
    import sys

    print("=== 4.3(b) Cas uniforme ===")
    try:
        plot_trivial()
    except FileNotFoundError:
        print("Fichiers trivial_*.out introuvables — executer d'abord ./engine trivial.in")

    print("\n=== 4.3(b)(ii) Convergence ===")
    try:
        convergence_study()
    except Exception as e:
        print(f"Convergence skippee ({e}) — verifier le chemin vers ./engine")

    print("\n=== 4.3(c) Cas non-trivial ===")
    try:
        plot_nontrivial()
    except FileNotFoundError:
        print("Fichiers nontrivial_*.out introuvables — executer d'abord ./engine nontrivial.in")

    print("\n=== 4.3(d) Verification Gauss ===")
    try:
        gauss_check()
    except FileNotFoundError:
        print("Fichiers nontrivial_*.out introuvables")
