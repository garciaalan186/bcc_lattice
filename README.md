# Numerical Coincidences from a BCC Lattice Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18478082.svg)](https://doi.org/10.5281/zenodo.18478082)

**Author:** Alan Garcia — Independent Researcher  
**Contact:** alan.javier.garcia@gmail.com

## How This Started

This project began with a simple question about wave–particle duality: what would have to be true for a photon to behave as both a wave and a particle?

I imagined two massive objects in a tight helical orbit at a scale far below anything we can currently measure. Viewed from the side, their orbit traces a wave pattern; viewed head-on, they look like a single point — particle-like. I recalled that probing the Planck length with light would require wavelengths so energetic they'd collapse into a black hole, so I pushed the picture further: a sub-Planckian binary system of singularities.

The question then became one of self-consistency. How small would those singularities have to be? What constraints would the framework need to satisfy? Over about three days of iterative computation, the framework converged on a body-centered cubic (BCC) lattice structure that derives the inverse fine structure constant (α⁻¹), the proton-electron mass ratio (m_p/m_e), and the muon-electron mass ratio (m_μ/m_e) from two inputs — the BCC coordination number *n* = 8 and the transcendental π — with zero free parameters.

| Constant | Prediction | Experimental | Agreement |
|----------|-----------|--------------|-----------|
| α⁻¹ | 137.035 999 084 | 137.035 999 084(21) | < 0.005 ppb |
| m_p/m_e | 1836.152 674 | 1836.152 673 43(11) | < 0.03 ppb |
| m_μ/m_e | 206.768 286 | 206.768 283 0(46) | 12.6 ppb (0.6σ) |

The physical picture underlying the lattice structure is still not clear, but the predictive accuracy is difficult to dismiss as coincidence, and I am seeking input from the physics and mathematics communities.

## What's in This Repo

| File | Description |
|------|-------------|
| `bsm_inquiry.tex` | The paper (LaTeX source, REVTeX 4.2 / APS format) |
| `bsm_inquiry.pdf` | Compiled PDF |

## Quick Verification

The core results can be verified in a few lines of Python:

```python
import numpy as np

def tau(n): return n*(2*n+1)+1
def sigma(n): return n*(2*n+1)
def mu(n): return 1.5*sigma(n)*(n+1)

N = 8; d = 3
TAU, SIGMA, MU = tau(N), sigma(N), mu(N)
c1 = np.pi**2 / 2

B = TAU + c1/TAU - 1/(2*TAU**2) - 1/((N-1)*TAU**3)
alpha_inv = B + 1/((N+2)*B**2)
alpha = 1/alpha_inv

t1 = np.pi*MU*SIGMA*alpha**3/4
t2 = (2*N+1)*np.pi**2*alpha*(1-alpha**2)/16
mass = MU + t1 + t2 - alpha**2

muon = (d/2)*(TAU + np.pi*alpha) + (c1/4)*(1 + d*alpha**2)

print(f"alpha^-1  = {alpha_inv:.9f}")  # 137.035999084  (CODATA: 137.035999084(21))
print(f"m_p/m_e   = {mass:.6f}")       # 1836.152674    (CODATA: 1836.15267343(11))
print(f"m_mu/m_e  = {muon:.6f}")       # 206.768286     (Expt:   206.7682830(46))
```

## Citation

If you use or reference this work, please cite:

> Garcia, A. (2026). *Numerical Coincidences from a BCC Lattice Framework*. Zenodo. [https://doi.org/10.5281/zenodo.18478082](https://doi.org/10.5281/zenodo.18478082)

Or use the `CITATION.cff` file included in this repository.

## Status

This is an independent investigation, not a claim of proof. I am actively seeking feedback — including identification of a trivial explanation for the numerical coincidences — from anyone with expertise in lattice field theory, representation theory, or mathematical physics.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
