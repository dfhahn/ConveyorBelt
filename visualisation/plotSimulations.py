import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import figaspect

from ConveyorBelt.src.system import system


def static_sim_plots(sys: system, x_range: tuple = None, title: str = "", out_path: str = None) -> str:
    """
    Plot giving the sampled space, position distribution and forces
    :param sys:
    :param x_range:
    :param title:
    :param out_path:
    :return:
    """
    # gather data
    x = [state.position for state in sys.trajectory]
    y = [state.totPotEnergy for state in sys.trajectory]
    shift = [state.dhdpos for state in sys.trajectory]

    if (x_range == None):
        x_pot = np.arange(min(x) + min(x) * 0.25, max(x) + max(x) * 0.25)
    elif (type(x_range) == range):
        x_pot = x_range
    else:
        x_pot = np.arange(min(x_range), max(x_range)+1)

    ytot_space = sys.potential.ene(x_pot)

    # plot
    w, h = figaspect(0.5)
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=[w, h])

    ax1.scatter(x, y, c="orange", alpha=0.8)
    ax1.scatter(x[0], y[0], c="green", alpha=0.8)
    ax1.scatter(x[-1], y[-1], c="red", alpha=0.8)
    ax1.plot(x_pot, ytot_space)
    ax1.set_ylabel("$V_pot$")
    ax1.set_xlabel("$x$")
    ax1.set_title("Potential Sampling")

    ax2.boxplot(x)
    ax2.set_ylabel("$x$")
    ax2.set_xlabel("$simulation$")
    ax2.set_title("x-Distribution")

    ax3.plot(range(len(x)), shift)
    ax3.set_ylabel("$dhdpos$")
    ax3.set_xlabel("$t$")
    ax3.set_title("Forces/shifts")

    fig.suptitle(title, y=1.08)
    fig.tight_layout()
    if (out_path == None):
        plt.show()
    else:
        fig.savefig(out_path)
        plt.close(fig)
    return out_path