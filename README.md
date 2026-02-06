# Numerical Coincidences from a BCC Lattice Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18478082.svg)](https://doi.org/10.5281/zenodo.18478082)

**Author:** Alan Garcia — Independent Researcher  
**Contact:** alan.javier.garcia@gmail.com

## How This Started

This project began with a simple question about wave–particle duality: what would have to be true for a photon to behave as both a wave and a particle?

I imagined two massive objects in a tight helical orbit at a scale far below anything we can currently measure. Viewed from the side, their orbit traces a wave pattern; viewed head-on, they look like a single point — particle-like. I recalled that probing the Planck length with light would require wavelengths so energetic they'd collapse into a black hole, so I pushed the picture further: a sub-Planckian binary system of singularities.

The question then became one of self-consistency. How small would those singularities have to be? What constraints would the framework need to satisfy? Over about three days of iterative computation, the framework converged on a body-centered cubic (BCC) lattice structure that derives the inverse fine structure constant (α⁻¹), the proton-electron mass ratio (m_p/m_e), the muon-electron mass ratio (m_μ/m_e), the tau-electron mass ratio (m_τ/m_e), and the Higgs boson mass (m_H) from two inputs — the BCC coordination number *n* = 8, the spatial dimension *d* = 3, and the transcendental π — with zero free parameters. The Higgs prediction (125.108 GeV) is genuinely predictive: it matches ATLAS to 0.02σ and will be confirmed or ruled out by the HL-LHC.

| Constant | Prediction | Experimental | Agreement |
|----------|-----------|--------------|-----------|
| α⁻¹ | 137.035 999 177 | 137.035 999 177(21) | < 0.001 ppb |
| m_p/m_e | 1836.152 673 5 | 1836.152 673 426(32) | 0.03 ppb |
| m_μ/m_e | 206.768 282 5 | 206.768 282 7(46) | 1.1 ppb (0.05σ) |
| m_τ/m_e | 3477.4799 | 3477.48 ± 0.57 | 0.027 ppm (0.0002σ) |
| m_H | 125.108 GeV | 125.11 ± 0.11 (ATLAS) | 0.02σ (ATLAS) |

The physical picture underlying the lattice structure is still not clear, but the predictive accuracy is difficult to dismiss as coincidence, and I am seeking input from the physics and mathematics communities.

## What's in This Repo

| File | Description |
|------|-------------|
| `bsm_inquiry.tex` | The paper (LaTeX source) |
| `bsm_inquiry.pdf` | Compiled PDF |
| `bsm_streamlined.py` | Minimal verification script (all five constants) |
| `bsm_core.py` | Full computation module with documentation |

## Quick Verification

The core results can be verified in a few lines of Python (requires `mpmath`):

```python
from mpmath import mp, mpf, pi, sqrt

mp.dps = 50  # 50 decimal places; float64 also sufficient

def tau(n): return n*(2*n+1)+1
def sigma(n): return n*(2*n+1)
def mu(n): return mpf(3)/2*sigma(n)*(n+1)

N = 8; d = 3
TAU, SIGMA, RHO, MU = mpf(tau(N)), mpf(sigma(N)), mpf(N+1), mu(N)
c1 = pi**2 / 2

c4 = (2*RHO/(2*N+1)) * pi**3
B = TAU + c1/TAU - 1/(2*TAU**2) - 1/((N-1)*TAU**3) + c4/TAU**4
alpha_inv = B + 1/((N+2)*B**2)
alpha = 1/alpha_inv

mass = MU*(1 + 2*(2*N+1)*pi*alpha**3) + (1 + 1/(2*N))*pi**2*alpha*(1-alpha**2) - alpha**2

muon = mpf(d)/2*(TAU + pi*alpha*(1-alpha**2)) + c1/4 + mpf(d)/2*pi*SIGMA*alpha**3 - alpha**2

Q = mpf(d-1)/d + mpf(d-1)*pi**2*alpha**2/((N-1)*SIGMA)
b = sqrt(muon)
D = (1 + muon)*(2*Q - 1) + 2*Q*b
tau_mass = ((Q*(1+b) + sqrt(D)) / (1-Q))**2

higgs = (SIGMA - d) * (1 + pi*alpha/RHO)  # m_H/m_p: breathing mode

print(f"alpha^-1  = {alpha_inv}")   # 137.035999177...  (CODATA: 137.035999177(21))
print(f"m_p/m_e   = {mass}")        # 1836.15267348...  (CODATA: 1836.152673426(32))
print(f"m_mu/m_e  = {muon}")        # 206.768282475...  (CODATA: 206.7682827(46))
print(f"m_tau/m_e = {tau_mass}")    # 3477.47990762...  (CODATA: 3477.48(57))
print(f"m_H/m_p   = {higgs}")       # 133.339 → 125.108 GeV (ATLAS: 125.11 ± 0.11)
```

## Citation

If you use or reference this work, please cite:

> Garcia, A. (2026). *Numerical Coincidences from a BCC Lattice Framework*. Zenodo. [https://doi.org/10.5281/zenodo.18478082](https://doi.org/10.5281/zenodo.18478082)

Or use the `CITATION.cff` file included in this repository.

## Status

This is an independent investigation, not a claim of proof. I am actively seeking feedback — including identification of a trivial explanation for the numerical coincidences — from anyone with expertise in lattice field theory, representation theory, or mathematical physics.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
