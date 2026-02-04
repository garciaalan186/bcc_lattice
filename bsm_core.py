#!/usr/bin/env python3
"""
BSM Core Module — Generator functions and constants
All BSM calculations centralized here

Updated: February 2026 (CODATA 2022 era)
  - α equation extended to four loops (c₄ = (2ρ/(2n+1))·π³)
  - Proton formula in n-explicit factored form (denominators = n/2 and 2n)
  - Muon formula extended to O(α³) with proton-parallel vertex
  - All results benchmarked against CODATA 2022
"""

import numpy as np

# =============================================================================
# GENERATOR FUNCTIONS
# =============================================================================

def rho(n):
    """Radial generator: ρ(n) = n + 1"""
    return n + 1

def sigma(n):
    """Surface generator: σ(n) = n(2n + 1)"""
    return n * (2 * n + 1)

def tau(n):
    """Topological generator: τ(n) = n(2n + 1) + 1"""
    return n * (2 * n + 1) + 1

def mu(n):
    """Mass generator: μ(n) = (3/2)σ(n)ρ(n)"""
    return 1.5 * sigma(n) * rho(n)

def compute_B(n):
    """
    Compute base value B for fine structure constant (four-loop).

    B = τ + c₁/τ - 1/(2τ²) - 1/((n-1)τ³) + c₄/τ⁴

    Loop coefficients and their spectral origins:
      c₁ = π²/2     Brillouin zone momentum integral ⟨|q|²⟩_BZ
      c₂ = -1/2     Spectral mean: -d/⟨D²⟩
      c₃ = -1/(n-1) Spectral variance
      c₄ = (2ρ/(2n+1))·π³  Spectral skewness × zone integral
    """
    t = tau(n)
    r = rho(n)
    c1 = np.pi**2 / 2
    c4 = (2 * r / (2*n + 1)) * np.pi**3
    return t + c1/t - 1/(2*t**2) - 1/((n-1)*t**3) + c4/t**4

def compute_alpha_inverse(n):
    """
    Compute α⁻¹ using cubic Dyson equation:
    α⁻¹ = B + 1/[(n+2)B²]

    This comes from solving: y³ - By² - 1/(n+2) = 0
    interpreted as G⁻¹ = G₀⁻¹ - Σ(G) for the dressed inverse coupling.
    """
    B = compute_B(n)
    return B + 1 / ((n + 2) * B**2)

def compute_alpha(n):
    """Compute fine structure constant α"""
    return 1.0 / compute_alpha_inverse(n)

def compute_mass_ratio(n):
    """
    Compute proton-electron mass ratio (n-explicit factored form):

    m_p/m_e = μ·(1 + 2(2n+1)·π·α³) + (1 + 1/(2n))·π²·α·(1-α²) - α²

    The denominators are functions of n:
      Vertex:  original /4  = /(n/2)  → prefactor 2/n
      VP:      original /16 = /(2n)   → prefactor (2n+1)/(2n) = 1 + 1/(2n)

    Correction taxonomy:
      μ·2(2n+1)·π·α³        — vertex (tree mass × 2·dim(D(n)) dressing)
      (1+1/(2n))·π²·α(1-α²) — VP (universal + lattice correction)
      -α²                    — self-intersection (universal)
    """
    m = mu(n)
    alpha = compute_alpha(n)

    vertex_dressed = m * (1 + 2*(2*n+1) * np.pi * alpha**3)
    vp = (1 + 1/(2*n)) * np.pi**2 * alpha * (1 - alpha**2)
    self_int = -alpha**2

    return vertex_dressed + vp + self_int

def compute_muon_ratio(n, d=3):
    """
    Compute muon-electron mass ratio (d-explicit form):

    m_μ/m_e = (d/2)(τ + πα(1-α²)) + c₁/4 + (d/2)πσα³ - α²

    The muon's structural parameter is d (spatial dimension),
    in contrast to the proton's n (coordination number):
      Proton vertex prefactor: 2/n  (lattice topology)
      Muon vertex prefactor:   d/2  (spatial dimension)

    Correction taxonomy (parallel to proton):
      (d/2)·π·α·(1-α²)  — VP (binary excitation, dressed)
      (d/2)·π·σ·α³       — vertex (= proton vertex × 2d/μ)
      -α²                — self-intersection (universal)
    """
    t = tau(n)
    s = sigma(n)
    c1 = np.pi**2 / 2
    alpha = compute_alpha(n)

    tree = (d / 2) * t + c1 / 4
    vp = (d / 2) * np.pi * alpha * (1 - alpha**2)
    vertex = (d / 2) * np.pi * s * alpha**3
    self_int = -alpha**2

    return tree + vp + vertex + self_int

# =============================================================================
# CONSTANTS AT n = 8
# =============================================================================

N = 8   # BCC coordination number
d = 3   # Spatial dimensions

# Generator values
TAU_8 = tau(N)      # 137
SIGMA_8 = sigma(N)  # 136
RHO_8 = rho(N)      # 9
MU_8 = mu(N)        # 1836

# Derived constants
B_8 = compute_B(N)
ALPHA_INV_BSM = compute_alpha_inverse(N)
ALPHA_BSM = compute_alpha(N)
MASS_RATIO_BSM = compute_mass_ratio(N)
MUON_RATIO_BSM = compute_muon_ratio(N)

# Measured values (CODATA 2022)
ALPHA_INV_MEASURED = 137.035999177
MASS_RATIO_MEASURED = 1836.152673426
MUON_RATIO_MEASURED = 206.7682827

# =============================================================================
# VERIFICATION
# =============================================================================

if __name__ == '__main__':
    print("=" * 72)
    print("BSM CORE MODULE — VERIFICATION (CODATA 2022)")
    print("=" * 72)
    print(f"\nAt n = {N} (BCC coordination):")
    print(f"  ρ(8) = {RHO_8}")
    print(f"  σ(8) = {SIGMA_8}")
    print(f"  τ(8) = {TAU_8}")
    print(f"  μ(8) = {MU_8:.0f}")

    c1 = np.pi**2 / 2
    c4 = (2*RHO_8/(2*N+1)) * np.pi**3
    alpha = ALPHA_BSM

    print(f"\nFine Structure Constant (four-loop):")
    print(f"  B = τ + c₁/τ - 1/(2τ²) - 1/((n-1)τ³) + c₄/τ⁴")
    print(f"  c₁ = π²/2 = {c1:.6f}  (BZ integral)")
    print(f"  c₂ = -1/2             (spectral mean)")
    print(f"  c₃ = -1/(n-1) = -1/7  (spectral variance)")
    print(f"  c₄ = (2ρ/(2n+1))·π³ = (18/17)·π³ = {c4:.6f}  (spectral skewness)")
    print(f"  B = {B_8:.12f}")
    print(f"  α⁻¹ = B + 1/[(n+2)B²] = {ALPHA_INV_BSM:.12f}")
    print(f"  CODATA 2022            = {ALPHA_INV_MEASURED}(21)")
    err_a = abs(ALPHA_INV_BSM - ALPHA_INV_MEASURED)
    print(f"  Error = {err_a:.3e} ({err_a/ALPHA_INV_MEASURED*1e9:.4f} ppb)")

    print(f"\nProton-Electron Mass Ratio (n-explicit):")
    print(f"  m_p/m_e = μ(1 + 2(2n+1)πα³) + (1 + 1/(2n))π²α(1-α²) - α²")
    print(f"  Vertex dressing: 2(2n+1)πα³ = 2×{2*N+1}×π×α³ = {2*(2*N+1)*np.pi*alpha**3:.12f}")
    print(f"  VP (universal):  π²α(1-α²) = {np.pi**2*alpha*(1-alpha**2):.12f}")
    print(f"  VP (lattice):    1/(2n) = 1/{2*N} correction")
    print(f"  Self-int:        -α² = {-alpha**2:.12f}")
    print(f"  m_p/m_e = {MASS_RATIO_BSM:.9f}")
    print(f"  CODATA 2022 = {MASS_RATIO_MEASURED}(32)")
    err_m = abs(MASS_RATIO_BSM - MASS_RATIO_MEASURED)
    print(f"  Error = {err_m:.3e} ({err_m/MASS_RATIO_MEASURED*1e9:.2f} ppb)")

    print(f"\nMuon-Electron Mass Ratio (d-explicit):")
    tree = (d/2)*TAU_8 + c1/4
    mu_vp = (d/2)*np.pi*alpha*(1-alpha**2)
    mu_vt = (d/2)*np.pi*SIGMA_8*alpha**3
    mu_si = -alpha**2
    print(f"  m_μ/m_e = (d/2)(τ + πα(1-α²)) + c₁/4 + (d/2)πσα³ - α²")
    print(f"  Tree:     (d/2)τ + c₁/4        = {tree:.10f}")
    print(f"  VP:       (d/2)πα(1-α²)         = {mu_vp:.12f}")
    print(f"  Vertex:   (d/2)πσα³              = {mu_vt:.12f}")
    print(f"  Self-int: -α²                    = {mu_si:.12f}")
    print(f"  m_μ/m_e = {MUON_RATIO_BSM:.10f}")
    print(f"  CODATA 2022 = {MUON_RATIO_MEASURED}(46)")
    err_mu = abs(MUON_RATIO_BSM - MUON_RATIO_MEASURED)
    print(f"  Error = {err_mu:.3e} ({err_mu/MUON_RATIO_MEASURED*1e9:.1f} ppb)")

    print(f"\nStructural contrast:")
    print(f"  Proton structural parameter: n = {N} (coordination number)")
    print(f"  Muon structural parameter:   d = {d} (spatial dimension)")
    print(f"  Proton vertex prefactor: 2/n = {2/N}")
    print(f"  Muon vertex prefactor:   d/2 = {d/2}")
    print(f"  Proton VP prefactor: (2n+1)/(2n) = {(2*N+1)/(2*N):.6f}")
    print(f"  Muon VP prefactor:   d/2 = {d/2}")

    print(f"\n{'='*72}")
    print(f"SCORECARD: Three constants, two inputs (n=8, π), zero free parameters")
    print(f"  α⁻¹    → {err_a/ALPHA_INV_MEASURED*1e9:.4f} ppb")
    print(f"  m_p/m_e → {err_m/MASS_RATIO_MEASURED*1e9:.2f} ppb")
    print(f"  m_μ/m_e → {err_mu/MUON_RATIO_MEASURED*1e9:.1f} ppb")
    print(f"{'='*72}")
