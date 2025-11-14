import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt

def setup_figure(columns=1,
                 nrow_fig=1,
                 ncol_fig=1,
                 # control these (in units of [0,1]) to position the figure
                 axleft=0.12,
                 axbottom=0.14,
                 axright=0.96,
                 axtop=1.07,
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

    axs = np.empty((nrow_fig, ncol_fig), dtype=object)

    for j in range(ncol_fig):
        for i in range(nrow_fig):
            axs[i,j] = plt.subplot(gs[i,j])
            axs[i,j].minorticks_on()

    fig.subplots_adjust(left=axleft, bottom=axbottom, right=axright, top=axtop)


    return fig, axs

def solve(K_arr, CFL):
    return np.acos(CFL**2 * (np.cos(K_arr * np.pi) - 1) + 1) / CFL



fig, axs = setup_figure(columns=1, ncol_fig=1, figsize_scale=2, colorbar=True)

N = 20
cfl_plot_stable = np.linspace(0, 1, N)
cfl_plot_unstable = np.linspace(1, 1.1, N)
cfl_plot = np.concatenate((cfl_plot_stable, cfl_plot_unstable))

norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(cfl_plot), vcenter=1, vmax=np.max(cfl_plot))
norm_stable = matplotlib.colors.Normalize(vmin=np.min(cfl_plot_stable),
                                          vmax=np.max(cfl_plot_stable))
norm_unstable = matplotlib.colors.Normalize(vmin=np.min(cfl_plot_unstable),
                                            vmax=np.max(cfl_plot_unstable))


cmap_stable = matplotlib.colormaps['winter']
cmap_unstable = matplotlib.colormaps['autumn']

all_colors = np.vstack((cmap_stable(norm_stable(cfl_plot_stable)),
                        cmap_unstable(norm_unstable(cfl_plot_unstable))))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom', all_colors)

colors = np.concatenate((cmap_stable(norm_stable(cfl_plot_stable)),
                         cmap_unstable(norm_unstable(cfl_plot_unstable))))


K_plot = np.linspace(0, 2, 200, dtype=complex)
for j, cfl in enumerate(reversed(cfl_plot)):
    i = cfl_plot.size - j - 1

    T_plot = solve(K_plot, cfl) / np.pi
    lw = 1
    for sign in (-1, 1):
        axs[0, 0].plot(K_plot,
                       sign * np.imag(T_plot),
                       color=colors[i],
                       alpha=1,
                       lw=lw,
                       linestyle=(0, (5, 1)),
                       )

        axs[0, 0].plot(K_plot,
                       sign * np.real(T_plot),
                       color=colors[i],
                       alpha=1,
                       lw=lw,
                       linestyle='solid',
                       )

cfl = 1
T_plot = solve(K_plot, cfl) / np.pi
axs[0, 0].plot(K_plot,
               -1 * np.real(T_plot),
               color="black",
               alpha=1,
               lw=2,
               linestyle='solid',
               )
axs[0, 0].plot(K_plot,
               np.real(T_plot),
               color="black",
               alpha=1,
               lw=2,
               linestyle='solid',
               )

axs[0, 0].set(ylabel=r'$\frac{\omega}{c} \frac{\Delta x}{\pi}$',
              xlabel=r'$k \frac{\Delta x}{\pi}$',
              )


p = axs[0, 0].get_position()
axwidth  = p.x1 - p.x0
axheight  = (p.y1 - p.y0) * 0.03
axpad = 0.02
baxs = fig.add_axes([p.x0, p.y1 + axpad, axwidth, axheight])

cbar = matplotlib.colorbar.ColorbarBase(baxs,
                                        cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal',
                                        ticklocation='top')
cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1, np.max(cfl_plot)])

cbar.set_label('$\\frac{c\\Delta t}{\\Delta x}$')

plt.show()
fig.savefig('cfl-plot.pdf', dpi=300)
