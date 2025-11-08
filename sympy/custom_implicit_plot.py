import numpy as np
import scipy

v = lambda gamma: np.sqrt(1 - 1 / gamma**2)
gamma = lambda v: 1.0 / np.sqrt(1 - v**2)

def newton_solve_disperion_eq(
    known_good, # tuple (x0, i, j) which contains known good guess for Zx[i], Zz[j]
    gamma_b: float,
    alpha: float,
    Zx_space,
    Zz_space,
):
    vb = v(gamma_b)
    vp = -alpha * vb
    gamma_p = gamma(vp)

    c = 1

    # define units in which plasma frequency is one
    omega = 1

    def dispersion_eq(x, Zx, Zz):
        k_x = omega * Zx / vb
        k_z = omega * Zz / vb
        return (-k_x**2 + omega**2*x**2*(-Zx**2*Zz*alpha**3/(gamma_p*x**2*(Zz**3*alpha**3 + 3*Zz**2*alpha**2*x + 3*Zz*alpha*x**2 + x**3)) + Zx**2*Zz*alpha/(gamma_b*x**2*(-Zz**3 + 3*Zz**2*x - 3*Zz*x**2 + x**3)) - Zx**2*alpha**2/(gamma_p*x*(Zz**3*alpha**3 + 3*Zz**2*alpha**2*x + 3*Zz*alpha*x**2 + x**3)) + Zx**2*alpha/(gamma_b*x*(Zz**3 - 3*Zz**2*x + 3*Zz*x**2 - x**3)) + Zz*alpha*c**2/(gamma_p*x*(Zz**2*alpha**4*gamma_p**2*vb**2 + Zz**2*alpha**2*c**2 + 2*Zz*alpha**3*gamma_p**2*vb**2*x + 2*Zz*alpha*c**2*x + alpha**2*gamma_p**2*vb**2*x**2 + c**2*x**2)) - Zz*alpha*c**2/(gamma_b*x*(Zz**2*c**2 + Zz**2*gamma_b**2*vb**2 - 2*Zz*c**2*x - 2*Zz*gamma_b**2*vb**2*x + c**2*x**2 + gamma_b**2*vb**2*x**2)) + alpha*c**2/(gamma_b*x*(Zz*c**2 + Zz*gamma_b**2*vb**2 - c**2*x - gamma_b**2*vb**2*x)) - c**2/(gamma_p*x*(Zz*alpha**3*gamma_p**2*vb**2 + Zz*alpha*c**2 + alpha**2*gamma_p**2*vb**2*x + c**2*x)) + 1)/c**2)*(-k_z**2 + omega**2*x**2*(-Zz**2*alpha**2/(gamma_p*x**2*(Zz**2*alpha**2 + 2*Zz*alpha*x + x**2)) - Zz**2*alpha/(gamma_b*x**2*(Zz**2 - 2*Zz*x + x**2)) - 2*Zz*alpha/(gamma_p*x*(Zz**2*alpha**2 + 2*Zz*alpha*x + x**2)) + 2*Zz*alpha/(gamma_b*x*(Zz**2 - 2*Zz*x + x**2)) - alpha/(gamma_b*(Zz**2 - 2*Zz*x + x**2)) + 1 - 1/(gamma_p*(Zz**2*alpha**2 + 2*Zz*alpha*x + x**2)))/c**2) - (k_x*k_z + omega**2*x**2*(Zx*Zz**2*alpha**3/(gamma_p*x**2*(Zz**3*alpha**3 + 3*Zz**2*alpha**2*x + 3*Zz*alpha*x**2 + x**3)) + Zx*Zz**2*alpha/(gamma_b*x**2*(Zz**3 - 3*Zz**2*x + 3*Zz*x**2 - x**3)) + 2*Zx*Zz*alpha**2/(gamma_p*x*(Zz**3*alpha**3 + 3*Zz**2*alpha**2*x + 3*Zz*alpha*x**2 + x**3)) + 2*Zx*Zz*alpha/(gamma_b*x*(-Zz**3 + 3*Zz**2*x - 3*Zz*x**2 + x**3)) + Zx*alpha/(gamma_p*(Zz**3*alpha**3 + 3*Zz**2*alpha**2*x + 3*Zz*alpha*x**2 + x**3)) + Zx*alpha/(gamma_b*(Zz**3 - 3*Zz**2*x + 3*Zz*x**2 - x**3)))/c**2)**2

    x_g, i_g, j_g = known_good

    x_arr = np.zeros((Zx_space.size, Zz_space.size), dtype=complex)

    def opt(x0, i, j):
        try:
            return scipy.optimize.newton(lambda x: dispersion_eq(x, Zx_space[i], Zz_space[j]),
                                         x0=x0,
                                         tol=1e-9,
                                         maxiter=1000)
        except:
            return np.nan

    def solve(dir, i, j):
        di, dj = dir

        xa = x_arr[i - di, j - dj]
        xb = x_arr[i - 2 * di, j - 2 * dj]

        delta_x = xa - xb

        x0 = xa + delta_x
        return opt(x0, i, j)

    x_arr[i_g, j_g] = opt(x_g, i_g, j_g)
    x_arr[i_g + 1, j_g] = opt(x_g, i_g + 1, j_g)
    x_arr[i_g, j_g + 1] = opt(x_g, i_g, j_g + 1)
    x_arr[i_g + 1, j_g + 1] = opt(x_g, i_g + 1, j_g + 1)

    for j in reversed(range(0, j_g)):
        x_arr[i_g, j] = solve((0, -1), i_g, j)
        x_arr[i_g + 1, j] = solve((0, -1), i_g + 1, j)

    for j in range(j_g + 1, Zz_space.size):
        x_arr[i_g, j] = solve((0, 1), i_g, j)
        x_arr[i_g + 1, j] = solve((0, 1), i_g + 1, j)

    for j in range(0, Zz_space.size):
        for i in reversed(range(0, i_g)):
            x_arr[i, j] = solve((-1, 0), i, j)

        for i in range(i_g + 2, Zx_space.size):
            x_arr[i, j] = solve((1, 0), i, j)

    Zx_arr, Zz_arr = np.meshgrid(Zx_space, Zz_space, indexing="ij")
    return x_arr, dispersion_eq(x_arr, Zx_arr, Zz_arr), Zx_arr, Zz_arr



def argfind_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

case = "A"

if case == "A":
    Zx_space = np.linspace(0, 4, 20)
    Zz_space = np.linspace(0, 1.5, 40)

    gamma_b = 3
    alpha = 1e-2

    Zx_good = 1.08
    Zz_good = 1
    Zx_good_idx = argfind_nearest(Zx_space, Zx_good)
    Zz_good_idx = argfind_nearest(Zz_space, Zz_good)
    x_good = 0.972 + 0.094j

    largest_value = np.sqrt(3) / 2**(4 / 3) * (alpha / gamma_b)**(1 / 3)

elif case == "B":
    Zx_space = np.linspace(0, 4, 20)
    Zz_space = np.linspace(0, 1.5, 80)

    gamma_b = 1.1
    alpha = 1

    Zx_good = 3
    Zz_good = 1.1
    Zx_good_idx = argfind_nearest(Zx_space, Zx_good)
    Zz_good_idx = argfind_nearest(Zz_space, Zz_good)
    x_good = 1 + 0.45j

    vb = v(gamma_b)
    largest_value = vb * np.sqrt(2 / gamma_b)


x_plot, _, Zx_plot, Zz_plot = newton_solve_disperion_eq(known_good=(x_good,
                                                                    Zx_good_idx,
                                                                    Zz_good_idx),
                                                        gamma_b=gamma_b,
                                                        alpha=alpha,
                                                        Zx_space=Zx_space,
                                                        Zz_space=Zz_space,
                                                        )

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()

ax = fig.add_subplot(121, projection='3d')
axre = fig.add_subplot(122, projection='3d')

ax.plot_surface(Zz_plot, Zx_plot, np.imag(x_plot))
ax.plot_surface(Zz_plot, Zx_plot, np.ones_like(Zz_plot) * largest_value)
axre.plot_surface(Zz_plot, Zx_plot, np.real(x_plot))

plt.show()
