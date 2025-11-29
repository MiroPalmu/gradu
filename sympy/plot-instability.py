"""
usage: <this-script> <filepath> [shift-x [shift-y]]

Plots the output csv file from pic-beam-instabilities/beam.py.

filepath: path to the csv file.
shift-x:  how much growth rate lines are shifted in x-direction.
shift-x:  how much growth rate lines are shifted in y-direction.
"""


import numpy as np
import sys
import argparse

import matplotlib
import matplotlib.pyplot as plt


def setup_figure(columns=1,
                 nrow_fig=1,
                 ncol_fig=1,
                 # control these (in units of [0,1]) to position the figure
                 axleft=0.20,
                 axbottom=0.18,
                 axright=0.96,
                 axtop=0.92,
                 figsize_scale=1,
):
    """
    Based on: https://natj.github.io/teaching/figs/
    """

    if columns == 1:
        fig = plt.figure(1, figsize=(figsize_scale * 3.25, figsize_scale * 2.2), dpi=300)
    elif columns == 2:
        fig = plt.figure(1, figsize=(figsize_scale * 7.0, figsize_scale * 2.5), dpi=300)
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


if __name__ == "__main__":

    fig, axs = setup_figure(figsize_scale=1)

    lw = 1.5
    ymin = 1e-5
    plot_shift_y = ymin

    for i in range(1, len(sys.argv), 2):
        file = sys.argv[i]
        plot_shift_x = float(sys.argv[i + 1])

        data = np.loadtxt(file, delimiter=",")

        t = data[:, 0]
        B2_parallel = data[:, 1]
        B2_perpendicular = data[:, 2]

        filamentation_growth = data[:, 7]


        C = "C" + str(int((i - 1) / 2))
        axs[0, 0].plot(t, B2_parallel + B2_perpendicular, alpha=1, lw=lw, linestyle="solid", color=C)




        axs[0, 0].plot(t + plot_shift_x,
                       plot_shift_y * filamentation_growth,
                       '--',
                       alpha=1,
                       lw=lw,
                       linestyle="dashed",
                       color="black")


    axs[0,0].set(yscale="log",
                 ylim=(ymin, 1),
                 ylabel=r'$\left\langle\frac{\rho_B}{\gamma_b n_0 m_e c^2}\right\rangle$',
                 xlabel=r'$t\omega_p$',
                 )

    plt.show()
    fig.savefig("instability-plot.pdf")
