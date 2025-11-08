"""
Follows Bret Plasma talks #9.

Solves eq. 11 for alpha
"""

import sympy as sp

x = sp.symbols("x")
alpha, Z, gamma_b, gamma_p, beta = sp.symbols("alpha Z gamma_b gamma_p beta", real=True)
e_xx = 1 - alpha / (x**2 * gamma_b) - 1 / (x**2 * gamma_p)
e_zz = 1 - alpha / (x**2 * gamma_b**3) - 1 / (x**2 * gamma_p**3) \
    - alpha * Z**2 / (x**4 * gamma_b) - alpha**2 * Z**2 / (x**4 * gamma_p)
e_xz = alpha * Z / x**3 * (1 / gamma_p - 1 / gamma_b)

dispersion_eq = sp.Eq(e_xx * (e_zz - Z**2 / x**2 / beta**2), e_xz)
solved_x4 = sp.solve(dispersion_eq, x**4)
assert len(solved_x4) == 1

x4lim = sp.limit(solved_x4[0], Z, sp.oo).doit()

solved_xlim = sp.solve(sp.Eq(x**4, x4lim), x)
sp.pprint(solved_xlim)

beta_from_gamma = lambda gamma: sp.sqrt(1 - 1 / gamma)

params = [(x, sp.I * x), (alpha, 1), (beta, beta_from_gamma(gamma_b)), (gamma_p, gamma_b), (gamma_b, 10)]
sp.plot_implicit(dispersion_eq.subs(params), x_var=(Z, 0, 10), y_var=(x, 0, 1))
