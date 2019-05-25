from src import potential1D as pot
from matplotlib import pyplot as plt


#UTIL FUNCTIONS
def significant_decimals(s:float)->float:
    significant_decimal=2
    if(s % 1 != 0):
        decimals = str(float(s)).split(".")[-1]
        for digit in decimals:
            if(digit == "0"):
                significant_decimal +=1
            else:
                return round(s, significant_decimal)
    else:
        return s

#show feature landscape per s
def envPot_differentS_overlay_min0_plot(eds_potential:pot.envelopedPotential, s_values:list, positions:list,
                                        y_range:tuple=None, hide_legend:bool=False, title:str=None, out_path:str=None):
    #generate energy values
    ys = []
    scale = 1 # 0.1
    for s in s_values:
        eds_potential.s=s
        enes = eds_potential.ene(positions)
        y_min =min(enes)
        y=list(map(lambda z: (z-y_min)*scale, enes))
        ys.append(y)

    #plotting
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20,10))
    for s, y in reversed(list(zip(s_values, ys))):
        axes.plot(positions, y, label="s_"+str(significant_decimals(s)))

    if (y_range != None):
        axes.set_ylim(y_range)
    axes.set_xlim(min(positions),max(positions))

    #styling
    axes.set_ylabel("Vr/[kJ]")
    axes.set_xlabel("r")
    axes.set_title("different Vrs aligned at 0 with different s-values overlayed ")

    ##optionals
    if(not hide_legend): axes.legend()
    if(title):    fig.suptitle(title)
    if(out_path): fig.savefig(out_path)
    fig.show()

    return fig, axes

#show feature landscape per s
def envPot_differentS_overlay_plot(eds_potential:pot.envelopedPotential, s_values:list, positions:list,
                                   y_range:tuple=None, hide_legend:bool=False, title:str=None, out_path:str=None):
    #generate energy values
    ys = []
    for s in s_values:
        eds_potential.s=s
        enes = eds_potential.ene(positions)
        ys.append(enes)

    #plotting
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20,10))
    for s, y in reversed(list(zip(s_values, ys))):
        axes.plot(positions, y, label="s_"+str(significant_decimals(s)))

    #styling
    axes.set_xlim(min(positions),max(positions))
    axes.set_ylabel("Vr/[kJ]")
    axes.set_xlabel("r")
    axes.set_title("different Vrs with different s-values overlayed ")

    ##optionals
    if (y_range != None): axes.set_ylim(y_range)
    if(not hide_legend): axes.legend()
    if(title):    fig.suptitle(title)
    if(out_path): fig.savefig(out_path)
    fig.show()

    return fig, axes

def envPot_diffS_compare(eds_potential:pot.envelopedPotential, s_values:list, positions:list,
                         y_range:tuple=None,title:str=None, out_path:str=None):
    ##row/column ratio
    per_row =4
    n_rows = (len(s_values)//per_row)+1 if ((len(s_values)%per_row)>0) else (len(s_values)//per_row)

    ##plot
    fig, axes = plt.subplots(nrows=n_rows, ncols=per_row, figsize=(20,10))
    axes = [ax for ax_row in axes for ax in ax_row]

    for ax, s in zip( axes, s_values):
        eds_potential.s=s
        y=eds_potential.ene(positions)
        ax.plot(positions, y)

        #styling
        ax.set_xlim(min(positions), max(positions))
        ax.set_title("s_"+str(significant_decimals(s)))
        ax.set_ylabel("Vr/[kJ]")
        ax.set_xlabel("r")
        if (y_range != None): ax.set_ylim(y_range)

    ##optionals
    if(title):    fig.suptitle(title)
    if(out_path): fig.savefig(out_path)
    fig.show()
    return fig, axes

def plot_envelopedPotential_system(eds_potential:pot.envelopedPotential, positions:list, s_value:float=None, Eoffi:list=None,
                                   y_range:tuple=None,title:str=None, out_path:str=None):
    if(s_value!=None):
        eds_potential.s = s_value       #set new s
    if(Eoffi!=None):
        if(len(Eoffi) == len(eds_potential.V_is)):
            eds_potential.Eoff_i = Eoffi
        else:
            raise IOError("There are "+str(len(eds_potential.V_is))+" states and "+str(Eoffi)+", but the numbers have to be equal!")

    ##calc energies
    energy_Vr = eds_potential.ene(positions)
    energy_Vis = [state.ene(positions) for state in eds_potential.V_is]
    num_states = len(eds_potential.V_is)

    ##plot nicely
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    axes = [ax for ax_row in axes for ax in ax_row]
    y_values = energy_Vis + [energy_Vr]
    labels = ["state_"+str(ind) for ind in range(1,len(energy_Vis)+1)]+["refState"]

    for ax, y, label in zip(axes, y_values, labels):
        ax.plot(positions, y)
        ax.set_xlim(min(positions), max(positions))
        ax.set_ylim(y_range)
        ax.set_title(label)
        ax.set_ylabel("Vr/[kJ]")
        ax.set_xlabel("r_"+label)

    ##optionals
    if(title):    fig.suptitle(title)
    if(out_path): fig.savefig(out_path)
    fig.show()
    return fig, axes

def plot_envelopedPotential_2State_System(eds_potential:pot.envelopedPotential, positions:list, s_value:float=None, Eoffi:list=None,
                                          title:str=None, out_path:str=None, V_max:float=600, V_min:float=None):
    if(len(eds_potential.V_is)>2):
        raise IOError(__name__+" can only be used with two states in the potential!")

    if(s_value!=None):
        eds_potential.s = s_value
    if (Eoffi != None):
        if (len(Eoffi) == len(eds_potential.V_is)):
            eds_potential.Eoff_i = Eoffi
        else:
            raise IOError("There are " + str(len(eds_potential.V_is)) + " states and " + str(
                Eoffi) + ", but the numbers have to be equal!")

    #Calculate energies
    energy_Vr = eds_potential.ene(positions)
    energy_Vis = [state.ene(positions) for state in eds_potential.V_is]
    energy_map = []
    min_e = 0
    for x in positions:
        row = eds_potential.ene([[x for i in positions], positions])
        row_cut = list(map(lambda x:  V_max if(V_max != None and x > V_max) else x, row))
        energy_map.append(row_cut)
        if(min(row)< min_e):
            min_e=min(row)
    if(V_min==None):
        V_min=min_e

    ##plot nicely
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    axes = [ax for ax_row in axes for ax in ax_row]
    y_values = energy_Vis + [energy_Vr]
    labels = ["State_" + str(ind) for ind in range(1, len(energy_Vis) + 1)] + ["State_R"]

    #plot the line potentials
    for ax, y, label in zip(axes, y_values, labels):
        ax.plot(positions, y)
        ax.set_xlim(min(positions), max(positions))
        ax.set_ylim([V_min, V_max])
        ax.set_title("Potential $"+label+"$")
        ax.set_ylabel("$V/[kJ]$")
        ax.set_xlabel("$r_{ " + label+"} $")

    #plot phase space surface
    ax = axes[-1]
    surf = ax.imshow(energy_map, cmap="viridis", interpolation="nearest",
                     origin='center', extent=[ min(positions), max(positions), min(positions), max(positions)],
                     vmax=V_max, vmin=V_min)
    ax.set_xlabel("$r_{"+labels[0]+"}$")
    ax.set_ylabel("$r_{"+labels[1]+"}$")
    ax.set_title("complete phaseSpace of $state_R$")
    fig.colorbar(surf, aspect=5, label='Energy/kJ')

    ##optionals
    if(title):    fig.suptitle(title)
    if(out_path): fig.savefig(out_path)
    fig.show()
    return fig, axes

if __name__ == "__main__":
    pass