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
    print(x_arr[i_g, j_g])
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

def solve(Zx_space,
          Zz_space,
          gamma_b,
          alpha,
          Zx_good,
          Zz_good,
          x_good,
):
    Zx_good_idx = argfind_nearest(Zx_space, Zx_good)
    Zz_good_idx = argfind_nearest(Zz_space, Zz_good)
    return newton_solve_disperion_eq(known_good=(x_good, Zx_good_idx, Zz_good_idx),
                                     gamma_b=gamma_b,
                                     alpha=alpha,
                                     Zx_space=Zx_space,
                                     Zz_space=Zz_space,
                                     )




import matplotlib
import matplotlib.pyplot as plt

def setup_figure(columns=1,
                 nrow_fig=1,
                 ncol_fig=1,
                 # control these (in units of [0,1]) to position the figure
                 axleft=0.08,
                 axbottom=0.24,
                 axright=0.96,
                 axtop=0.92,
                 figsize_scale=1,
                 colorbar=True,
):
    """
    Based on: https://natj.github.io/teaching/figs/
    """

    if colorbar:
        axtop *= 0.8

    if columns == 1:
        fig = plt.figure(1, figsize=(figsize_scale * 3.25, figsize_scale * 2.2))
    elif columns == 2:
        fig = plt.figure(1, figsize=(figsize_scale * 7.0, figsize_scale * 4))
    else:
        raise RuntimeError("This only supports single- or two-column figures")

    # add ticks to both sides
    plt.rc('xtick', top   = True)
    plt.rc('ytick', right = True)

    plt.rc('font',  family='serif',)
    plt.rc('text',  usetex=False)

    # make labels slightly smaller
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=8)
    plt.rc('legend', handlelength=4.0)

    gs = plt.GridSpec(nrow_fig, ncol_fig)
    gs.update(wspace = 0.25)
    gs.update(hspace = 0.35)

    axs = np.empty((nrow_fig, ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)

    if colorbar:
        baxs = dict()
        for j in range(ncol_fig):
            for i in range(nrow_fig):
                p = axs[i, j].get_position()
                axwidth  = p.x1 - p.x0
                axheight  = (p.y1 - p.y0) * 0.03
                axpad = 0.02
                baxs[i, j] = fig.add_axes([p.x0, p.y1 + axpad, axwidth, axheight])

        return fig, axs, baxs


    return fig, axs

fig, axs, baxs = setup_figure(columns=2, ncol_fig=2, figsize_scale=0.7, colorbar=True)
cmap = matplotlib.colormaps['jet']

res = 80

x_plot, _, Zx_plot, Zz_plot = solve(
    Zx_space = np.linspace(0, 4, res),
    Zz_space = np.linspace(0, 1.5, res),
    gamma_b = 3,
    alpha = 1,
    Zx_good = 3,
    Zz_good = 0.1,
    x_good = 0 + 0.55j,
)

norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(np.imag(x_plot)))

E = [np.min(Zz_plot), np.max(Zz_plot), np.min(Zx_plot), np.max(Zx_plot)]
axs[0, 0].imshow(np.flip(np.imag(x_plot), axis=0), cmap='jet', extent=E, aspect='auto')
axs[0, 0].set(xlabel='$Z_z$', ylabel='$Z_x$')

cbar = matplotlib.colorbar.ColorbarBase(baxs[0, 0],
                                        cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal',
                                        ticklocation='top')
cbar.set_label('$\\Im(x)$')

x_plot, _, Zx_plot, Zz_plot = solve(
    Zx_space = np.linspace(0, 4, res),
    Zz_space = np.linspace(0, 1.5, res),
    gamma_b = 3,
    alpha = 0.1,
    Zx_good = 2.95,
    Zz_good = 1.27,
    x_good = 0.98 + 0.121j,
)

E = [np.min(Zz_plot), np.max(Zz_plot), np.min(Zx_plot), np.max(Zx_plot)]
axs[0, 1].imshow(np.flip(np.imag(x_plot), axis=0), cmap='jet', extent=E, aspect='auto')
axs[0, 1].set(xlabel='$Z_z$', ylabel='$Z_x$')

norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(np.imag(x_plot)))
cbar = matplotlib.colorbar.ColorbarBase(baxs[0, 1],
                                        cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal',
                                        ticklocation='top')
cbar.set_label('$\\Im(x)$')

plt.show()
fig.savefig('growth-rate.pdf', dpi=300)
