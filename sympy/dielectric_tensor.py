"""
[1] Bret "Beam-plasma dielectric tensor with Mathematica" (2007).
"""


import sympy as sp

kx, kz = sp.symbols("k_x k_z")
E1x, E1y, E1z = sp.symbols("E_{1x} E_{1y} E_{1z}")
vb0z, vp0z = sp.symbols("v_{b0z} v_{p0z}")
m, c, omega_b, q, omega = sp.symbols("m c omega_b q omega")
vb1x, vb1y, vb1z = sp.symbols("v_{b1x} v_{b1y} v_{b1z}")
vp1x, vp1y, vp1z = sp.symbols("v_{p1x} v_{p1y} v_{p1z}")
gamma_b, gamma_p = sp.symbols("gamma_b gamma_p")

k = sp.Matrix([kx, 0, kz])
E1 = sp.Matrix([E1x, E1y, E1z])
Vb0 = sp.Matrix([0, 0, vb0z])
Vp0 = sp.Matrix([0, 0, vp0z])
B0 = sp.Matrix([0, 0, m * c * omega_b / q])
B1 = c * k.cross(E1) / omega
Vb1 = sp.Matrix([vb1x, vb1y, vb1z])
Vp1 = sp.Matrix([vp1x, vp1y, vp1z])

def gamma(v):
    return 1 / (1 - v.dot(v) / c**2)

def eq4(V0, V1, g):
    """
    Eq. 4 of [1].
    Note that [1] defines q < 0.
    """

    lhs = sp.I * m * g * (k.dot(V0) - omega) * (V1 + g**2 / c**2 * V0.dot(V1) * V0)
    rhs = q * (E1 + ((V0 + V1).cross(B0) + V0.cross(B1)) / c)
    return sp.Eq(lhs, rhs)

EqVb = eq4(Vb0, Vb1, gamma_b).simplify()
solved_Vb1_c = sp.solve(EqVb, Vb1)
solved_Vb1 = sp.Matrix([solved_Vb1_c[vb1x], solved_Vb1_c[vb1y], solved_Vb1_c[vb1z]])

EqVp = eq4(Vp0, Vp1, gamma_p).simplify()
solved_Vp1_c = sp.solve(EqVp, Vp1)
solved_Vp1 = sp.Matrix([solved_Vp1_c[vp1x], solved_Vp1_c[vp1y], solved_Vp1_c[vp1z]])

omega_pb, omega_pp = sp.symbols("omega_{pb} omega_{pp}")

nb0 = omega_pb**2 * m / (4 * sp.pi * q**2)
np0 = omega_pp**2 * m / (4 * sp.pi * q**2)

# Eq. 3 of [1].
nb1 = (nb0 * k.dot(solved_Vb1) / (omega - k.dot(Vb0)))
np1 = (np0 * k.dot(solved_Vp1) / (omega - k.dot(Vp0)))

# Eq. 5 of [1].
J = q * (nb0 * solved_Vb1 + nb1 * Vb0 + np0 * solved_Vp1 + np1 * Vp0)

def E1_coeff_matrix(X):
    """
    Given vector X, which is written in terms of E1,
    this functions returns a matrix M, s.t. M * E1 = X
    """
    rows = [[sp.collect(X[i].expand(), E1x).coeff(E1x),
             sp.collect(X[i].expand(), E1y).coeff(E1y),
             sp.collect(X[i].expand(), E1z).coeff(E1z)] for i in range(3)]

    return sp.Matrix(rows)

# k x (k x E1) || see eq. 6 of [1].
kk = E1_coeff_matrix(k.cross(k.cross(E1)))
sp.pprint(kk)

# epsilon := E1 + 4 * i * pi * J / omega
epsilon = E1_coeff_matrix(E1 + 4 * sp.I * sp.pi * J / omega)

# Substitute dimensional variables from eq. 7 of [1].
Zx, Zz, x, alpha = sp.symbols("Zx Zz x alpha")

epsilon_sub = epsilon.subs(vp0z, -alpha * vb0z) \
                     .subs(kz, omega_pp * Zz / vb0z) \
                     .subs(kx, omega_pp * Zx / vb0z) \
                     .subs(omega_pb, sp.sqrt(alpha) * omega_pp) \
                     .subs(omega, x * omega_pp) \
                     .subs(omega_b, 0)

# Used to store wich terms are non-zero without Zz = 0 substitution
m = sp.eye(3)

# This has to be manually simplified to get nice expression.
# But even this does not do everything.
for i in range(3):
    for j in range(3):
        m[i, j] = int(epsilon_sub[i, j] != 0)
        epsilon_sub[i, j] = epsilon_sub[i, j].subs(Zz, 0)
        epsilon_sub[i, j] = sp.Add(*[sp.simplify(sp.cancel(t)) for t in epsilon_sub[i, j].args])

sp.pprint(epsilon_sub)

print("Which terms of epsilon are non-zero without substituing Zz = 0:")
sp.pprint(m)
