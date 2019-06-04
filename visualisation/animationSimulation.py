import os
import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from IPython.display import HTML

from src.system import system

def animation_trajectory(sys: system, x_range=None, title:str=None, out_path:str=None)-> (animation.Animation, (str or None)):
    # plotting
    x1data = [state.position for state in sys.trajectory]
    y1data = [state.totPotEnergy for state in sys.trajectory]
    shift = [state.dhdpos for state in sys.trajectory]
    x_max = max(x1data)
    x_min = min(x1data)

    if (x_range == None):
        xtot_space = np.arange(x_min + 0.2 * x_min, x_max + 0.2 * x_max + 1)
    else:
        xtot_space = np.arange(min(x_range), max(x_range) + 1, 1)
    ytot_space = sys.potential.ene(xtot_space)

    tmax = len(y1data) - 1
    t0 = 0
    step_size = 1

    xdata, ydata = [], []

    # figures
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xtot_space, ytot_space, label="wholePot")
    # data structures in ani
    line, = ax.plot([], [], c="skyblue", lw=2)
    scatter, = ax.plot([], [], "bo", c="orange", alpha=0.6, ms=6)
    start_p, = ax.plot([], [], "bo", c="g", ms=10)
    end_p, = ax.plot([], [], "bo", c="r", ms=10)
    curr_p, = ax.plot([], [], "bo", c="k", ms=10)

    #Params
    ax.set_xlabel("$r$")
    ax.set_ylabel("$V$")
    if (title != None):
        fig.suptitle(title)


    def init():
        del xdata[:], ydata[:]

        line.set_data(xtot_space, ytot_space)
        start_p.set_data(x1data[0], y1data[0])
        scatter.set_data([], [])
        end_p.set_data([], [])
        curr_p.set_data([], [])

        if (x_range != None):
            ax.set_xlim(x_range)
        return line,

    def data_gen(t=t0):
        while t < tmax:
            t += step_size
            yield x1data[t], y1data[t]

    def run(data):
        # update the data
        x, V = data

        if (x == x1data[-1]):  # last step of traj
            curr_p.set_data([], [])
            end_p.set_data(x1data[-1], y1data[-1])
        else:
            curr_p.set_data([x], [V])
            xdata.append(x)
            ydata.append(V)
            scatter.set_data(xdata, ydata)

        if(min(ax.get_ylim())> V):
            ax.set_ylim(V+V*0.1, max(ax.get_ylim()))
        return line

    ani = animation.FuncAnimation(fig=fig, func=run, frames=data_gen, init_func=init, blit=False,
                                  interval=20, repeat=False, save_count=len(x1data))
    if(out_path != None):
        if( not "mp4" in os.path.splitext(out_path)[1]):
            out_path += ".mp4"
        # Set up formatting for the movie files
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='animationsMD1D'), bitrate=1800)
        ani.save('lines.mp4', writer=writer)

    return ani, out_path

def animation_EDS_trajectory(sys: system, x_range=None, title:str=None, out_path:str=None, s_values:list=[1.0])-> (animation.Animation, (str or None)):
    # plotting
    x1data = [state.position for state in sys.trajectory]
    y1data = [state.totPotEnergy for state in sys.trajectory]
    shift = [state.dhdpos for state in sys.trajectory]
    x_max = max(x1data)
    x_min = min(x1data)

    if (x_range == None):
        xtot_space = np.arange(x_min + 0.2 * x_min, x_max + 0.2 * x_max + 1)
    else:
        xtot_space = np.arange(min(x_range), max(x_range) + 1, 1)

    tmax = len(y1data) - 1
    t0 = 0
    step_size = 1

    xdata, ydata = [], []

    # figures
    fig = plt.figure()
    ax = fig.add_subplot(111)

    from visualisation.plotPotentials import envPot_differentS_overlay_plot

    _, ax = envPot_differentS_overlay_plot(eds_potential=sys.potential, s_values=s_values, positions=xtot_space, axes=ax)

    # data structures in ani
    #line, = ax.plot([], [], c="skyblue", lw=2)
    scatter, = ax.plot([], [], "bo", c="orange", alpha=0.6, ms=6)
    start_p, = ax.plot([], [], "bo", c="g", ms=10)
    end_p, = ax.plot([], [], "bo", c="r", ms=10)
    curr_p, = ax.plot([], [], "bo", c="k", ms=10)

    #Params
    ax.set_xlabel("$r$")
    ax.set_ylabel("$V$")
    if (title != None):
        fig.suptitle(title)


    def init():
        del xdata[:], ydata[:]

        #line.set_data(xtot_space, ytot_space)
        start_p.set_data(x1data[0], y1data[0])
        scatter.set_data([], [])
        end_p.set_data([], [])
        curr_p.set_data([], [])

        if (x_range != None):
            ax.set_xlim(x_range)
        #return line,

    def data_gen(t=t0):
        while t < tmax:
            t += step_size
            yield x1data[t], y1data[t]

    def run(data):
        # update the data
        x, V = data

        if (x == x1data[-1]):  # last step of traj
            curr_p.set_data([], [])
            end_p.set_data(x1data[-1], y1data[-1])
        else:
            curr_p.set_data([x], [V])
            xdata.append(x)
            ydata.append(V)
            scatter.set_data(xdata, ydata)

        if(min(ax.get_ylim())> V):
            ax.set_ylim(V+V*0.1, max(ax.get_ylim()))
        #return line

    ani = animation.FuncAnimation(fig=fig, func=run, frames=data_gen, init_func=init, blit=False,
                                  interval=20, repeat=False, save_count=len(x1data))
    if(out_path != None):
        if( not "mp4" in os.path.splitext(out_path)[1]):
            out_path += ".mp4"
        # Set up formatting for the movie files
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='animationsMD1D'), bitrate=1800)
        ani.save('lines.mp4', writer=writer)

    return ani, out_path