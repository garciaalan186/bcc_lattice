from mpmath import mp, mpf, pi, sqrt

mp.dps = 50  # 50 decimal places (float64 sufficient, but mpmath is exact)

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

print(f"alpha^-1  = {alpha_inv}")   # 137.035999177  (CODATA 2022: 137.035999177(21))
print(f"m_p/m_e   = {mass}")        # 1836.152673485 (CODATA 2022: 1836.152673426(32))
print(f"m_mu/m_e  = {muon}")        # 206.7682824754 (CODATA 2022: 206.7682827(46))
print(f"m_tau/m_e = {tau_mass}")    # 3477.4799      (CODATA 2022: 3477.48(57))
