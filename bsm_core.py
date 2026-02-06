#!/usr/bin/env python3
"""
BSM Core Module — Generator functions and constants
All BSM calculations centralized here

Updated: February 2026 (Session XIII — Higgs Mass Prediction)
  - Switched from numpy to mpmath for arbitrary-precision arithmetic
  - α equation extended to four loops (c₄ = (2ρ/(2n+1))·π³)
  - Muon formula extended to O(α³) with proton-parallel vertex
  - Tau prediction via dressed Koide with closed-form quadratic solution
  - Higgs mass ratio m_H/m_p = (σ-d)(1 + πα/ρ) predicts 125.108 GeV
  - All results benchmarked against CODATA 2022 / LHC Run 2
"""

from mpmath import mp, mpf, pi, sqrt

# Set default precision (50 decimal places; float64 sufficient for physics,
# but mpmath eliminates any doubt about floating-point artifacts)
mp.dps = 50

# =============================================================================
# GENERATOR FUNCTIONS
# =============================================================================

def rho(n):
    """Radial generator: ρ(n) = n + 1"""
    return mpf(n) + 1

def sigma(n):
    """Surface generator: σ(n) = n(2n + 1)"""
    return mpf(n) * (2 * n + 1)

def tau(n):
    """Topological generator: τ(n) = n(2n + 1) + 1"""
    return mpf(n) * (2 * n + 1) + 1

def mu(n):
    """Mass generator: μ(n) = (3/2)σ(n)ρ(n)"""
    return mpf(3) / 2 * sigma(n) * rho(n)

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
    c1 = pi**2 / 2
    c4 = (2 * r / (2*n + 1)) * pi**3
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
    return 1 / compute_alpha_inverse(n)

def compute_mass_ratio(n):
    """
    Compute proton-electron mass ratio (n-explicit factored form):

    m_p/m_e = μ·(1 + 2(2n+1)·π·α³) + (1 + 1/(2n))·π²·α·(1-α²) - α²

    Correction taxonomy:
      μ·2(2n+1)·π·α³        — vertex (tree mass × 2·dim(D(n)) dressing)
      (1+1/(2n))·π²·α(1-α²) — VP (universal + lattice correction)
      -α²                    — self-intersection (universal)
    """
    m = mu(n)
    alpha = compute_alpha(n)

    vertex_dressed = m * (1 + 2*(2*n+1) * pi * alpha**3)
    vp = (1 + 1/(2*mpf(n))) * pi**2 * alpha * (1 - alpha**2)
    self_int = -alpha**2

    return vertex_dressed + vp + self_int

def compute_muon_ratio(n, d=3):
    """
    Compute muon-electron mass ratio (d-explicit form):

    m_μ/m_e = (d/2)(τ + πα(1-α²)) + c₁/4 + (d/2)πσα³ - α²

    Correction taxonomy (parallel to proton):
      (d/2)·π·α·(1-α²)  — VP (binary excitation, dressed)
      (d/2)·π·σ·α³       — vertex (binary, = proton vertex × 2d/μ)
      -α²                — self-intersection (universal)

    Tree level: (d/2)·τ + c₁/4 (binary defect + BZ vacuum energy)
    """
    d = mpf(d)
    t = tau(n)
    s = sigma(n)
    c1 = pi**2 / 2
    alpha = compute_alpha(n)

    tree = (d / 2) * t + c1 / 4
    vp = (d / 2) * pi * alpha * (1 - alpha**2)
    vertex = (d / 2) * pi * s * alpha**3
    self_int = -alpha**2

    return tree + vp + vertex + self_int

def compute_koide_Q(n, d=3):
    """
    Compute dressed Koide parameter Q.

    Q = (d-1)/d + (d-1)·π²·α² / ((n-1)·σ)

    Tree level: (d-1)/d = 2/3 (transverse fraction of line defect)
    Correction: (d-1)·π²·α²/((n-1)·σ) (transverse normal-bundle dressing)

    Factor decomposition of the correction:
      (d-1)     — transverse directions (codimension)
      π²        — BZ momentum integral
      α²        — two-loop EM coupling
      1/(n-1)   — spectral variance coefficient (three-loop, = 1/7 at n=8)
      1/σ       — pairwise channel normalization

    Session XI: identified by systematic scan, interpreted as
    leading Seeley-DeWitt heat kernel coefficient of the transverse D².
    """
    d = mpf(d)
    s = sigma(n)
    alpha = compute_alpha(n)

    Q_tree = (d - 1) / d
    delta_Q = (d - 1) * pi**2 * alpha**2 / ((n - 1) * s)

    return Q_tree + delta_Q

def compute_tau_ratio(n, d=3):
    """
    Compute tau-electron mass ratio via dressed Koide formula
    using the closed-form quadratic solution (Session XII).

    The Koide equation:
      (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = Q

    with m_e = 1, reduces to a quadratic in x = √m_τ:

      (1 - Q)x² - 2Q(1 + b)x + [(1 + m_μ)(1 - Q) - Q·b²] = 0

    where b = √m_μ. The discriminant factors as:

      D = (1 + m_μ)(2Q - 1) + 2Qb

    and the physical (positive) root is:

      x = [Q(1 + b) + √D] / (1 - Q)
      m_τ = x²

    Q = (d-1)/d + (d-1)π²α²/((n-1)σ)

    Physical origin: Koide relation is a partition function constraint
    on the transverse oscillator of the binary line defect, with
    Q = (d-1)/d being the transverse fraction, dressed by the leading
    normal-bundle heat kernel correction.
    """
    m_mu = compute_muon_ratio(n, d)
    Q = compute_koide_Q(n, d)

    b = sqrt(m_mu)
    D = (1 + m_mu) * (2*Q - 1) + 2*Q*b
    x = (Q * (1 + b) + sqrt(D)) / (1 - Q)

    return x**2

def compute_higgs_ratio(n, d=3):
    """
    Compute Higgs-proton mass ratio m_H/m_p.

    m_H/m_p = (σ - d)(1 + πα/ρ)

    Tree level: σ - d = 133 (scalar pairwise modes minus Goldstone bosons)
      σ = n(2n+1) = 136 pairwise scalar couplings at each lattice site
      d = 3 Goldstone modes eaten by W±, Z (one per broken translational generator)

    Radiative correction: πα/ρ (vacuum polarization per radial mode)
      πα  — one-loop vacuum polarization (universal BSM VP factor)
      1/ρ — inverse radial shell count; breathing mode is outermost radial oscillation

    Physical interpretation: The Higgs is the collective breathing mode
    (uniform oscillation of lattice spacing) of the 133 surviving scalar
    channels after Goldstone subtraction.
    """
    s = sigma(n)
    r = rho(n)
    alpha = compute_alpha(n)
    d = mpf(d)

    tree = s - d
    correction = pi * alpha / r

    return tree * (1 + correction)

# Proton mass in GeV (CODATA 2022) for absolute Higgs mass
M_PROTON_GEV = mpf('0.938272088')

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
KOIDE_Q_BSM = compute_koide_Q(N)
TAU_RATIO_BSM = compute_tau_ratio(N)
HIGGS_RATIO_BSM = compute_higgs_ratio(N)
HIGGS_GEV_BSM = HIGGS_RATIO_BSM * M_PROTON_GEV

# Measured values (CODATA 2022 / LHC Run 2)
ALPHA_INV_MEASURED = mpf('137.035999177')
MASS_RATIO_MEASURED = mpf('1836.152673426')
MUON_RATIO_MEASURED = mpf('206.7682827')
TAU_RATIO_MEASURED = mpf('3477.48')
TAU_RATIO_UNCERTAINTY = mpf('0.57')
HIGGS_ATLAS_GEV = mpf('125.11')
HIGGS_ATLAS_UNC = mpf('0.11')
HIGGS_CMS_GEV = mpf('125.35')
HIGGS_CMS_UNC = mpf('0.15')

# =============================================================================
# VERIFICATION
# =============================================================================

if __name__ == '__main__':
    print("=" * 72)
    print("BSM CORE MODULE — VERIFICATION (CODATA 2022, Session XII)")
    print(f"Precision: {mp.dps} decimal places (mpmath)")
    print("=" * 72)
    print(f"\nAt n = {N} (BCC coordination), d = {d}:")
    print(f"  ρ(8) = {int(RHO_8)}")
    print(f"  σ(8) = {int(SIGMA_8)}")
    print(f"  τ(8) = {int(TAU_8)}")
    print(f"  μ(8) = {int(MU_8)}")

    c1 = pi**2 / 2
    c4 = (2*RHO_8/(2*N+1)) * pi**3
    alpha = ALPHA_BSM

    print(f"\n{'─'*72}")
    print(f"Fine Structure Constant (four-loop Dyson):")
    print(f"  B = τ + c₁/τ - 1/(2τ²) - 1/((n-1)τ³) + c₄/τ⁴")
    print(f"  c₁ = π²/2 = {float(c1):.6f}  (BZ integral)")
    print(f"  c₂ = -1/2             (spectral mean)")
    print(f"  c₃ = -1/(n-1) = -1/7  (spectral variance)")
    print(f"  c₄ = (2ρ/(2n+1))·π³ = (18/17)·π³ = {float(c4):.6f}  (spectral skewness)")
    print(f"  B = {B_8}")
    print(f"  α⁻¹ = B + 1/[(n+2)B²] = {ALPHA_INV_BSM}")
    print(f"  CODATA 2022            = {ALPHA_INV_MEASURED}(21)")
    err_a = abs(ALPHA_INV_BSM - ALPHA_INV_MEASURED)
    print(f"  Error = {float(err_a):.3e} ({float(err_a/ALPHA_INV_MEASURED*1e9):.4f} ppb)")

    print(f"\n{'─'*72}")
    print(f"Proton-Electron Mass Ratio (n-explicit):")
    print(f"  m_p/m_e = μ(1 + 2(2n+1)πα³) + (1 + 1/(2n))π²α(1-α²) - α²")
    print(f"  m_p/m_e = {MASS_RATIO_BSM}")
    print(f"  CODATA 2022 = {MASS_RATIO_MEASURED}(32)")
    err_m = abs(MASS_RATIO_BSM - MASS_RATIO_MEASURED)
    print(f"  Error = {float(err_m):.3e} ({float(err_m/MASS_RATIO_MEASURED*1e9):.2f} ppb)")

    print(f"\n{'─'*72}")
    print(f"Muon-Electron Mass Ratio (d-explicit):")
    print(f"  m_μ/m_e = (d/2)(τ + πα(1-α²)) + c₁/4 + (d/2)πσα³ - α²")
    print(f"  m_μ/m_e = {MUON_RATIO_BSM}")
    print(f"  CODATA 2022 = {MUON_RATIO_MEASURED}(46)")
    err_mu = abs(MUON_RATIO_BSM - MUON_RATIO_MEASURED)
    print(f"  Error = {float(err_mu):.3e} ({float(err_mu/MUON_RATIO_MEASURED*1e9):.1f} ppb)")

    print(f"\n{'─'*72}")
    print(f"Tau-Electron Mass Ratio (dressed Koide, closed-form):")
    Q_tree = mpf(d-1)/d
    delta_Q = mpf(d-1) * pi**2 * alpha**2 / ((N-1) * SIGMA_8)
    print(f"  Q = (d-1)/d + (d-1)π²α²/((n-1)σ)")
    print(f"  Q_tree     = {Q_tree}")
    print(f"  δQ         = {delta_Q}")
    print(f"  Q_dressed  = {KOIDE_Q_BSM}")
    print(f"  m_τ/m_e = {TAU_RATIO_BSM}")
    print(f"  CODATA 2022 = {TAU_RATIO_MEASURED}(57)")
    err_tau = abs(TAU_RATIO_BSM - TAU_RATIO_MEASURED)
    print(f"  Error = {float(err_tau):.4f} ({float(err_tau/TAU_RATIO_MEASURED*1e6):.3f} ppm, {float(err_tau/TAU_RATIO_UNCERTAINTY):.4f}σ)")

    print(f"\n{'─'*72}")
    print(f"Higgs Boson Mass (breathing mode):")
    print(f"  m_H/m_p = (σ - d)(1 + πα/ρ)")
    print(f"  Tree level: σ - d = {int(SIGMA_8)} - {d} = {int(SIGMA_8 - d)}")
    print(f"  Correction: πα/ρ = π × {float(alpha):.10f} / {int(RHO_8)} = {float(pi*alpha/RHO_8):.6f}")
    print(f"  m_H/m_p = {float(HIGGS_RATIO_BSM):.6f}")
    print(f"  m_H     = {float(HIGGS_GEV_BSM):.6f} GeV")
    print(f"  ATLAS   = {HIGGS_ATLAS_GEV} ± {HIGGS_ATLAS_UNC} GeV")
    print(f"  CMS     = {HIGGS_CMS_GEV} ± {HIGGS_CMS_UNC} GeV")
    err_atlas = abs(HIGGS_GEV_BSM - HIGGS_ATLAS_GEV)
    err_cms = abs(HIGGS_GEV_BSM - HIGGS_CMS_GEV)
    print(f"  Deviation (ATLAS) = {float(err_atlas):.3f} GeV ({float(err_atlas/HIGGS_ATLAS_UNC):.2f}σ)")
    print(f"  Deviation (CMS)   = {float(err_cms):.3f} GeV ({float(err_cms/HIGGS_CMS_UNC):.1f}σ)")

    print(f"\n{'='*72}")
    print(f"SCORECARD: Five constants, two inputs (n=8, d=3, π), zero free parameters")
    print(f"  α⁻¹    → {float(err_a/ALPHA_INV_MEASURED*1e9):.4f} ppb")
    print(f"  m_p/m_e → {float(err_m/MASS_RATIO_MEASURED*1e9):.2f} ppb")
    print(f"  m_μ/m_e → {float(err_mu/MUON_RATIO_MEASURED*1e9):.1f} ppb")
    print(f"  m_τ/m_e → {float(err_tau/TAU_RATIO_MEASURED*1e6):.3f} ppm  ({float(err_tau/TAU_RATIO_UNCERTAINTY):.4f}σ)")
    print(f"  m_H/m_p → {float(err_atlas):.3f} GeV from ATLAS ({float(err_atlas/HIGGS_ATLAS_UNC):.2f}σ)")
    print(f"{'='*72}")
