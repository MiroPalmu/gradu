import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt

def setup_figure(columns=1,
                 nrow_fig=1,
                 ncol_fig=1,
                 # control these (in units of [0,1]) to position the figure
                 axleft=0.14,
                 axbottom=0.18,
                 axright=0.96,
                 axtop=0.92,
                 figsize_scale=1,
):
    """
    Based on: https://natj.github.io/teaching/figs/
    """

    if columns == 1:
        fig = plt.figure(1, figsize=(figsize_scale * 3.25, figsize_scale * 2.2))
    elif columns == 2:
        fig = plt.figure(1, figsize=(figsize_scale * 7.0, figsize_scale * 2.5))
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

    return fig, axs

def B(x, k, i, t):
   if k == 0:
      return 1.0 if t[i] <= x < t[i+1] else 0.0
   if t[i+k] == t[i]:
      c1 = 0.0
   else:
      c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
   if t[i+k+1] == t[i+1]:
      c2 = 0.0
   else:
      c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
   return c1 + c2

def divide_to_knots(x, k):
    if k < 2:
        raise RuntimeError()

    N = x.size
    M = int(N / (k - 1))

    knots = np.array([x[i] for i in np.arange(0, N, M)])
    if knots[-1] == x[-1]:
        return knots
    else:
        return np.concatenate((knots, np.array([x[-1]])))


fig, axs = setup_figure(columns=1, ncol_fig=1, figsize_scale=1)
cmap = matplotlib.colormaps['jet']

res = 400
k = 4
norm_stable = matplotlib.colors.Normalize(vmin=0,
                                          vmax=k - 1)

for i in range(k):
    x_plot = np.linspace(-(i + 1) / 2., (i + 1) / 2., res)
    knots = divide_to_knots(x_plot[1:], i + 2)
    y_plot = [B(x, i, 0, knots) for x in x_plot]

    axs[0, 0].plot(x_plot, y_plot, color=cmap(i / k), alpha=1, linestyle='solid', lw=1.5)


pos_ticks = np.array([0.5 * i for i in range(1, k + 1)])
axs[0, 0].set(xlabel=r'$\Delta x$',
              ylabel=r'$\Delta x S(x)$',
              xticks=np.concatenate((-pos_ticks, [0], pos_ticks))
              )

plt.show()
fig.savefig('spline.pdf', dpi=300)
