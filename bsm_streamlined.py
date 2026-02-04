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
