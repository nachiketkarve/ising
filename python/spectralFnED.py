import os

# Recommended when using joblib parallelism with scipy.linalg.eigh.
# Otherwise BLAS may oversubscribe your CPU.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import eigh
from scipy.ndimage import gaussian_filter1d
from joblib import Parallel, delayed


# ============================================================
# Bit utilities
# ============================================================

def translate_left(s, N):
    """
    Translate bitstring by one site.

    Site convention does not matter as long as it is consistent.
    """
    mask = (1 << N) - 1
    return ((s << 1) & mask) | (s >> (N - 1))


def diagonal_energy_z(s, N, J, hz):
    """
    Diagonal energy of

        J sum_i sigma_i^z sigma_{i+1}^z + hz sum_i sigma_i^z

    for computational basis state s.

    Convention:
        bit 1 -> sigma^z = +1
        bit 0 -> sigma^z = -1
    """
    e = 0.0

    for i in range(N):
        zi = 1 if ((s >> i) & 1) else -1
        zj = 1 if ((s >> ((i + 1) % N)) & 1) else -1

        e += J * zi * zj
        e += hz * zi

    return e


# ============================================================
# Translation orbit construction
# ============================================================

def build_translation_orbits(N):
    """
    Build all translation orbits of the computational basis.

    Returns
    -------
    reps : list[int]
        Representative of each orbit.

    period : dict[int, int]
        period[rep] = orbit length.

    rep_of : ndarray
        rep_of[s] = representative of the orbit containing s.

    shift_of : ndarray
        shift_of[s] = r such that s = T^r rep_of[s].
    """
    D = 1 << N

    seen = np.zeros(D, dtype=bool)
    rep_of = np.empty(D, dtype=np.int64)
    shift_of = np.empty(D, dtype=np.int16)

    reps = []
    period = {}

    for s in range(D):
        if seen[s]:
            continue

        orbit = []
        x = s

        while not seen[x]:
            seen[x] = True
            orbit.append(x)
            x = translate_left(x, N)

        rep = min(orbit)

        # Rebuild orbit starting from canonical representative
        canonical_orbit = []
        x = rep

        while True:
            canonical_orbit.append(x)
            x = translate_left(x, N)
            if x == rep:
                break

        R = len(canonical_orbit)

        reps.append(rep)
        period[rep] = R

        for r, x in enumerate(canonical_orbit):
            rep_of[x] = rep
            shift_of[x] = r

    return reps, period, rep_of, shift_of


# ============================================================
# Momentum-sector construction
# ============================================================

def allowed_representatives_for_momentum(N, q, reps, period):
    """
    Momentum k = 2 pi q / N.

    A translation orbit of period R contributes to momentum q only if

        exp(-i k R) = 1,

    i.e.

        q R = 0 mod N.
    """
    return [rep for rep in reps if (q * period[rep]) % N == 0]


def build_V_sector(N, q, reps, period, rep_of, shift_of):
    """
    Build V = sum_i sigma_i^x in fixed momentum sector q.

    Basis convention:

        |rep, k> = 1/sqrt(R) sum_{r=0}^{R-1} exp(-i k r) T^r |rep>

    Since V is translationally invariant, it preserves momentum.
    """
    k = 2.0 * np.pi * q / N

    sector_reps = allowed_representatives_for_momentum(N, q, reps, period)
    dim = len(sector_reps)

    index = {rep: a for a, rep in enumerate(sector_reps)}

    V = np.zeros((dim, dim), dtype=np.complex128)

    for a, rep_a in enumerate(sector_reps):
        Ra = period[rep_a]

        for i in range(N):
            y = rep_a ^ (1 << i)

            rep_b = int(rep_of[y])

            b = index.get(rep_b)
            if b is None:
                continue

            Rb = period[rep_b]
            shift = int(shift_of[y])

            # Matrix element:
            # <rep_b,k| V |rep_a,k>
            amp = np.sqrt(Ra / Rb) * np.exp(1j * k * shift)

            V[b, a] += amp

    # Remove tiny numerical phase noise
    V = 0.5 * (V + V.conj().T)
    V[np.abs(V) < 1e-14] = 0.0

    return sector_reps, V


def build_H_and_V_sector(N, q, J, hx, hz, reps, period, rep_of, shift_of):
    """
    Build H and V in momentum sector q.
    """
    sector_reps, V = build_V_sector(N, q, reps, period, rep_of, shift_of)

    dim = len(sector_reps)

    H = np.zeros((dim, dim), dtype=np.complex128)

    for a, rep in enumerate(sector_reps):
        H[a, a] = diagonal_energy_z(rep, N, J, hz)

    H += hx * V

    H = 0.5 * (H + H.conj().T)
    H[np.abs(H) < 1e-14] = 0.0

    return H, V


def diagonalize_sector(N, q, J, hx, hz, reps, period, rep_of, shift_of):
    """
    Worker function for one momentum sector.
    """
    H, V = build_H_and_V_sector(
        N, q, J, hx, hz, reps, period, rep_of, shift_of
    )

    if H.shape[0] == 0:
        return None

    E, U = eigh(H)

    # V in the energy eigenbasis
    VE = U.conj().T @ V @ U
    VE = 0.5 * (VE + VE.conj().T)

    return {
        "q": q,
        "E": E,
        "VE": VE,
        "dim": len(E),
    }


# ============================================================
# Spectral accumulation
# ============================================================

def make_frequency_grid(omega_max, n_omega):
    """
    Symmetric frequency grid and histogram edges.
    """
    omega = np.linspace(-omega_max, omega_max, n_omega)
    d_omega = omega[1] - omega[0]

    edges = np.empty(n_omega + 1)
    edges[1:-1] = 0.5 * (omega[:-1] + omega[1:])
    edges[0] = omega[0] - 0.5 * d_omega
    edges[-1] = omega[-1] + 0.5 * d_omega

    return omega, edges, d_omega


def accumulate_sector_histogram(
    sector,
    omega_edges,
    beta,
    Emin,
    Z,
    Vmean,
    drop_diagonal_delta=True,
    block_size=256,
):
    """
    Accumulate the delta-function spectral weight into a histogram.

    This computes integrated weights per frequency bin. Gaussian broadening
    is applied later.
    """
    E = sector["E"]
    VE = sector["VE"]

    dim = len(E)

    hist = np.zeros(len(omega_edges) - 1, dtype=np.float64)

    p = np.exp(-beta * (E - Emin))

    Vc = VE.copy()

    # Connected subtraction: V -> V - <V>
    diag = np.diag(Vc).copy()
    np.fill_diagonal(Vc, diag - Vmean)

    if drop_diagonal_delta:
        # Removes exact finite-size diagonal omega=0 delta peak.
        np.fill_diagonal(Vc, 0.0)

    absV2 = np.abs(Vc) ** 2

    for m0 in range(0, dim, block_size):
        m1 = min(m0 + block_size, dim)

        Em = E[m0:m1][:, None]
        En = E[None, :]

        dE = En - Em

        pm = p[m0:m1][:, None]
        pn = p[None, :]

        weights = (np.pi / Z) * (pm + pn) * absV2[m0:m1, :]

        dE_flat = dE.ravel()
        w_flat = weights.ravel()

        tol = 1e-10
        keep = (w_flat > 1e-28) & (np.abs(dE_flat) > tol)

        if np.any(keep):
            h, _ = np.histogram(
                dE_flat[keep],
                bins=omega_edges,
                weights=w_flat[keep],
            )
            hist += h

    return hist


# ============================================================
# Main function
# ============================================================

def connected_symmetrized_spectral_function(
    N,
    J,
    hx,
    hz,
    beta=0.0,
    eta=0.05,
    n_omega=2001,
    omega_max=None,
    drop_diagonal_delta=True,
    per_site=True,
    n_jobs=-1,
    block_size=256,
    verbose=10,
):
    """
    Compute the connected, symmetrized spectral function of

        V = sum_i sigma_i^x

    for the mixed-field Ising model

        H = J sum_i sigma_i^z sigma_{i+1}^z
          + hx sum_i sigma_i^x
          + hz sum_i sigma_i^z

    using ED with translation symmetry.

    The spectral convention is

        S_sym(omega)
        =
        pi / Z sum_mn (p_m + p_n)
        |<m|V_c|n>|^2 delta(omega - (E_n - E_m)),

    where

        V_c = V - <V>.

    Parameters
    ----------
    N : int
        System size.

    J, hx, hz : float
        Hamiltonian parameters.

    beta : float
        Inverse temperature. beta=0 gives infinite temperature.

    eta : float
        Gaussian broadening width.

    n_omega : int
        Number of frequency grid points.

    omega_max : float or None
        Frequency cutoff. If None, use the many-body bandwidth.

    drop_diagonal_delta : bool
        If True, removes the exact finite-size omega=0 diagonal contribution.

    per_site : bool
        If True, returns S(omega) / N.

    n_jobs : int
        Number of parallel workers. Use -1 for all available cores.

    block_size : int
        Block size for spectral accumulation.

    verbose : int
        joblib verbosity.

    Returns
    -------
    omega : ndarray
        Frequency grid.

    S : ndarray
        Connected symmetrized spectral function.

    info : dict
        Diagnostics.
    """
    D = 1 << N

    print(f"Building translation orbits for N={N}, D={D}...")
    reps, period, rep_of, shift_of = build_translation_orbits(N)

    print(f"Number of translation orbits: {len(reps)}")

    print("Diagonalizing momentum sectors...")

    sectors = Parallel(n_jobs=n_jobs, backend="loky", verbose=verbose)(
        delayed(diagonalize_sector)(
            N, q, J, hx, hz, reps, period, rep_of, shift_of
        )
        for q in range(N)
    )

    sectors = [s for s in sectors if s is not None]

    all_E = np.concatenate([s["E"] for s in sectors])
    Emin = np.min(all_E)
    Emax = np.max(all_E)
    bandwidth = Emax - Emin

    if omega_max is None:
        omega_max = bandwidth

    omega, omega_edges, d_omega = make_frequency_grid(omega_max, n_omega)

    # Partition function and thermal expectation value <V>
    Z = 0.0
    Vmean_num = 0.0

    for s in sectors:
        E = s["E"]
        VE = s["VE"]

        p = np.exp(-beta * (E - Emin))

        Z += np.sum(p)
        Vmean_num += np.sum(p * np.real(np.diag(VE)))

    Vmean = Vmean_num / Z

    print(f"Total Hilbert dimension check: {sum(s['dim'] for s in sectors)}")
    print(f"Emin = {Emin:.12g}, Emax = {Emax:.12g}, bandwidth = {bandwidth:.12g}")
    print(f"Shifted partition function Z = {Z:.12g}")
    print(f"<V> = {Vmean:.12g}")

    print("Accumulating spectral weights...")

    hist_parts = Parallel(n_jobs=n_jobs, backend="loky", verbose=verbose)(
        delayed(accumulate_sector_histogram)(
            s,
            omega_edges,
            beta,
            Emin,
            Z,
            Vmean,
            drop_diagonal_delta,
            block_size,
        )
        for s in sectors
    )

    hist = np.sum(hist_parts, axis=0)

    # Gaussian broadening.
    #
    # hist contains integrated weight per bin.
    # gaussian_filter1d preserves the sum of hist.
    # Divide by d_omega at the end to get a density.
    sigma_bins = eta / d_omega
    S = gaussian_filter1d(hist, sigma=sigma_bins, mode="constant")
    S = S / d_omega

    if per_site:
        S = S / N

    info = {
        "N": N,
        "Hilbert_dimension": D,
        "num_translation_orbits": len(reps),
        "sector_dimensions": {s["q"]: s["dim"] for s in sectors},
        "J": J,
        "hx": hx,
        "hz": hz,
        "beta": beta,
        "eta": eta,
        "omega_max": omega_max,
        "n_omega": n_omega,
        "d_omega": d_omega,
        "Emin": Emin,
        "Emax": Emax,
        "bandwidth": bandwidth,
        "Z_shifted": Z,
        "Vmean": Vmean,
        "drop_diagonal_delta": drop_diagonal_delta,
        "per_site": per_site,
        "n_jobs": n_jobs,
        "block_size": block_size,
    }

    return omega, S, info


# ============================================================
# Example run
# ============================================================

if __name__ == "__main__":

    N = 14

    J = 1.0
    hx = 1.05
    hz = 0.5

    beta = 0.01

    omega, S, info = connected_symmetrized_spectral_function(
        N=N,
        J=J,
        hx=hx,
        hz=hz,
        beta=beta,
        eta=0.005,
        n_omega=4001,
        omega_max=None,
        drop_diagonal_delta=True,
        per_site=True,
        n_jobs=8,
        block_size=256,
        verbose=10,
    )

    print("\nDiagnostics:")
    for key, val in info.items():
        if key != "sector_dimensions":
            print(f"{key}: {val}")

    # Save data
    out = np.column_stack([omega, S])
    np.savetxt(
        f"./../data/spectralFnED_N{N}_J{J:.6f}_hx{hx:.6f}_hz{hz:.6f}_beta{beta:.6f}.csv",
        out,
        delimiter=",",
        header="omega,S_connected_symmetrized",
        comments="",
    )

    # Plot positive frequencies
    mask = omega > 0

    plt.figure()
    plt.plot(omega[mask], S[mask], lw=1.5)
    plt.xlabel(r"$\omega$")
    plt.ylabel(r"$S_{\rm sym}^{\rm conn}(\omega)/N$")
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.show()