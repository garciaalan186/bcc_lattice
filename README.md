# Numerical Coincidences from a BCC Lattice Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18478082.svg)](https://doi.org/10.5281/zenodo.18478082)

**Author:** Alan Garcia — Independent Researcher  
**Contact:** alan.javier.garcia@gmail.com

## How This Started

This project began with a simple question about wave–particle duality: what would have to be true for a photon to behave as both a wave and a particle?

I imagined two massive objects in a tight helical orbit at a scale far below anything we can currently measure. Viewed from the side, their orbit traces a wave pattern; viewed head-on, they look like a single point — particle-like. I recalled that probing the Planck length with light would require wavelengths so energetic they'd collapse into a black hole, so I pushed the picture further: a sub-Planckian binary system of singularities.

The question then became one of self-consistency. How small would those singularities have to be? What constraints would the framework need to satisfy? Over about three days of iterative computation, the framework converged on a body-centered cubic (BCC) lattice structure that derives the inverse fine structure constant (α⁻¹), the proton-electron mass ratio (m_p/m_e), and the muon-electron mass ratio (m_μ/m_e) from two inputs — the BCC coordination number *n* = 8 and the transcendental π — with zero free parameters.

| Constant | Prediction | CODATA 2022 | Agreement |
|----------|-----------|-------------|-----------|
| α⁻¹ | 137.035 999 177 | 137.035 999 177(21) | < 0.001 ppb |
| m_p/m_e | 1836.152 673 5 | 1836.152 673 426(32) | < 0.03 ppb |
| m_μ/m_e | 206.768 282 5 | 206.768 282 7(46) | 1.1 ppb (0.05σ) |

The physical picture underlying the lattice structure is still not clear, but the predictive accuracy is difficult to dismiss as coincidence, and I am seeking input from the physics and mathematics communities.

## The Equations

**Fine structure constant** (four-loop Dyson equation):

```
α⁻¹ = B + 1/((n+2)·B²)

B = τ + c₁/τ - 1/(2τ²) - 1/((n-1)τ³) + c₄/τ⁴
```

where c₁ = π²/2 (BZ momentum integral) and c₄ = (2ρ/(2n+1))·π³ (spectral skewness).

**Proton mass ratio** (n-explicit factored form):

```
m_p/m_e = μ·(1 + 2(2n+1)·π·α³) + (1 + 1/(2n))·π²·α·(1-α²) - α²
```

The denominators that appear as 4 and 16 at n = 8 are n/2 and 2n — the coordination number, halved and doubled. The vertex dresses the tree mass μ by 2·dim(D(n)) = 2(2n+1). The VP splits into a universal piece π²α(1−α²) plus a lattice correction 1/(2n) that vanishes as n → ∞.

**Muon mass ratio** (d-explicit form):

```
m_μ/m_e = (d/2)·(τ + π·α·(1-α²)) + c₁/4 + (d/2)·π·σ·α³ - α²
```

The muon follows the same correction taxonomy as the proton — vertex, dressed VP, self-intersection — but its structural parameter is the spatial dimension *d* = 3 rather than the coordination number *n* = 8. The proton knows it lives on a lattice; the muon only knows it lives in space.

## What's in This Repo

| File | Description |
|------|-------------|
| `bsm_inquiry.tex` | The paper (LaTeX source) |
| `bsm_inquiry.pdf` | Compiled PDF |
| `bsm_core.py` | Full module with all generator functions and verification |
| `bsm_streamlined.py` | Minimal 12-line verification script |

## Quick Verification

The core results can be verified in a few lines of Python:

```python
import numpy as np

def tau(n): return n*(2*n+1)+1
def sigma(n): return n*(2*n+1)
def mu(n): return 1.5*sigma(n)*(n+1)

N = 8; d = 3
TAU, SIGMA, RHO, MU = tau(N), sigma(N), N+1, mu(N)
c1 = np.pi**2 / 2

c4 = (2*RHO/(2*N+1)) * np.pi**3
B = TAU + c1/TAU - 1/(2*TAU**2) - 1/((N-1)*TAU**3) + c4/TAU**4
alpha_inv = B + 1/((N+2)*B**2)
alpha = 1/alpha_inv

mass = MU*(1 + 2*(2*N+1)*np.pi*alpha**3) + (1 + 1/(2*N))*np.pi**2*alpha*(1-alpha**2) - alpha**2

muon = (d/2)*(TAU + np.pi*alpha*(1-alpha**2)) + c1/4 + (d/2)*np.pi*SIGMA*alpha**3 - alpha**2

print(f"alpha^-1  = {alpha_inv:.12f}")  # 137.035999177  (CODATA 2022: 137.035999177(21))
print(f"m_p/m_e   = {mass:.9f}")        # 1836.152673485 (CODATA 2022: 1836.152673426(32))
print(f"m_mu/m_e  = {muon:.10f}")       # 206.7682824754 (CODATA 2022: 206.7682827(46))
```

## Citation

If you use or reference this work, please cite:

> Garcia, A. (2026). *Numerical Coincidences from a BCC Lattice Framework*. Zenodo. [https://doi.org/10.5281/zenodo.18478082](https://doi.org/10.5281/zenodo.18478082)

Or use the `CITATION.cff` file included in this repository.

## Status

This is an independent investigation, not a claim of proof. I am actively seeking feedback — including identification of a trivial explanation for the numerical coincidences — from anyone with expertise in lattice field theory, representation theory, or mathematical physics.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
