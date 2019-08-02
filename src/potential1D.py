"""
Module: Potential
    This module shall be used to implement subclasses of Potential. This module contains all available potentials."""

import math
import numpy as np
import typing as t
import numbers
import scipy.constants as const
from collections.abc import Iterable

class _potentialCls:
    '''
    potential base class
    '''

    name:str = "Unknown"
    nDim:int =1
    def __init__(self, *kargs):
        return

    def __name__(self)->str:
        return str("Unknown")

    def _check_positions_type(self, positions:t.Union[t.Iterable[numbers.Number], numbers.Number, t.Iterable[t.Iterable[numbers.Number]]])->t.List[float]:
        if (isinstance(positions, Iterable)):
            if(all([isinstance(x, numbers.Number) for x in positions])):
                positions = [positions]
                pass
                #positions = list(map(float, list(positions)))
            elif(len(positions) == self.nDim and all([isinstance(x, Iterable) for x in positions])):
                pass
                #positions = list(map(lambda x: list(map(float, x)), positions))[:self.nDim]
            else:
                raise Exception("list dimensionality does not fit to potential dimensionality! len(list)="+str(len(list))+" potDims "+str(self.nDim))
        else:
            positions = [[float(positions)]]
        return positions

    def _calculate_energies(self, positions:t.List[float], *kargs):
        raise NotImplementedError("Function " + __name__ + " was not implemented for class " + str(__class__) + "")

    def _calculate_dhdpos(self, positions:t.List[float], *kargs):
        raise NotImplementedError("Function " + __name__ + " was not implemented for class " + str(__class__) + "")

    def ene(self, positions:(t.List[float] or float), *kargs) -> (t.List[float] or float):
        '''
        calculates energy of particle
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: energy
        '''
        positions = self._check_positions_type(positions)
        return self._calculate_energies(positions)

    def dhdpos(self, positions:(t.List[float] or float), *kargs) -> (t.List[float] or float):
        '''
        calculates derivative with respect to position
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: derivative dh/dpos
        '''
        positions = self._check_positions_type(positions)
        return self._calculate_dhdpos(positions)


class flat_well(_potentialCls):
    '''
        .. autoclass:: flat well potential
    '''
    name:str = "Flat Well"
    x_min: float = None
    x_max: float = None
    y_max:float = None
    y_min:float = None

    def __init__(self, x_range: list = [0, 1], y_max: float = 1000, y_min:float=0):
        '''
        initializes flat well potential class

        '''
        super(flat_well, self).__init__()
        self.x_min = min(x_range)
        self.x_max = max(x_range)
        self.y_max = y_max
        self.y_min = y_min

    def _calculate_energies(self, positions, *kargs):
        return [self.y_min if (pos >= self.x_min and pos <= self.x_max) else self.y_max for pos in positions[0]]

    def _calculate_dhdpos(self, positions: (t.List[float] or float), *kargs) -> (t.List[float] or float):
        return [0 for pos in positions[0]]

"""
    SIMPLE POTENTIALS
"""

class harmonicOsc1D(_potentialCls):
    '''
        .. autoclass:: harmonic oscillator potential
    '''
    name:str = "harmonicOscilator"
    x_shift = None
    fc = None

    def __init__(self, fc:float=1.0, x_shift:float=0.0, y_shift:float=0.0):
        '''
        initializes harmonicOsc1D class
        :param fc: force constant
        :param x_shift: minimum position of harmonic oscillator
        '''
        super(harmonicOsc1D, self).__init__()
        self.fc = fc
        self.x_shift = x_shift
        self.y_shift = y_shift

    def _calculate_energies(self, positions:t.List[float], *kargs) -> (t.List[float]):
        return [list(map(lambda pos: 0.5 * self.fc * (pos - self.x_shift) ** 2 - self.y_shift, dimPos)) for dimPos in positions]

    def _calculate_dhdpos(self, positions:t.List[float], *kargs) -> (t.List[float]):
        return [list(map(lambda pos: self.fc * (pos - self.x_shift), dimPos)) for dimPos in positions]

class wavePotential(_potentialCls):
    '''
        .. autoclass:: Wave Potential
    '''
    name:str = "Wave Potential"

    phase_shift:float = 0.0
    amplitude:float = 1.0
    multiplicity:float = 1.0
    radians:bool = False

    def __init__(self,  phase_shift:float=0.0, multiplicity:float=1.0, amplitude:float=1.0, radians:bool=False):
        '''
        initializes wavePotential potential class
        '''
        super().__init__()
        self.amplitude = amplitude
        self.multiplicity = multiplicity
        self.set_radians(radians)
        if(radians):
            self.phase_shift = phase_shift
        else:
            self.phase_shift = np.deg2rad(phase_shift)

    def set_degrees(self, degrees:bool=True):
        self.set_radians(radians=not degrees)
    def set_radians(self, radians:bool=True):
        self.radians=radians
        if(radians):
            self._calculate_energies = lambda positions: list(map(lambda x: self.amplitude*math.cos(self.multiplicity*(x + self.phase_shift)), positions))
            self._calculate_dhdpos = lambda positions: list(map(lambda x: self.amplitude*math.sin(self.multiplicity*(x + self.phase_shift)), positions))
        else:
            self._calculate_energies = lambda positions: list(
                map(lambda x: self.amplitude * math.cos(self.multiplicity * (x + self.phase_shift)), np.deg2rad(positions)))
            self._calculate_dhdpos = lambda positions: list(
                map(lambda x: self.amplitude * math.sin(self.multiplicity * (x + self.phase_shift)), np.deg2rad(positions)))

class torsionPotential(_potentialCls):
    '''
        .. autoclass:: Torsion Potential
    '''
    name:str = "Torsion Potential"

    phase:float=1.0
    wave_potentials:t.List[wavePotential]=[]

    def __init__(self, wave_potentials:t.List[wavePotential]):
        '''
        initializes torsions Potential
        '''
        super().__init__()
        if(type(wave_potentials) == list):
            self.wave_potentials=wave_potentials
        else:
            self.wave_potentials=[wave_potentials]

        if(len(self.wave_potentials)>1):
            self._calculate_energies = lambda positions: np.add(*map(lambda x: np.array(x.ene(positions)), self.wave_potentials))
            self._caclulcate_dhdpos = lambda positions: np.add(*map(lambda x: np.array(x.dhdpos(positions)), self.wave_potentials))
        elif(len(self.wave_potentials)==1):
            self._calculate_energies = lambda positions: np.array(self.wave_potentials[0].ene(positions))
            self._caclulcate_dhdpos = lambda positions: np.array(self.wave_potentials[0].dhdpos(positions))

class coulombPotential(_potentialCls):

    epsilon:float

    coulombLaw = lambda q1,q2,r,epsilon: (q1*q2/(r*epsilon*4*math.pi))
    def __init__(self, q1=1, q2=1, epsilon=1):
        self.q1=q1
        self.q2=q2
        self.epsilon=epsilon
        pass

    def _calculate_energies(self, distances:t.List[float], *kargs):
        coulombLaw_curry = lambda x: self.coulombLaw(self.q1, self.q2, x, self.epsilon)
        return list(map(coulombLaw_curry, distances))
    def _calculate_dhdpos(self, distances:t.List[float], *kargs):
        raise NotImplemented("deviation not implemented yet")

class lennardJonesPotential(_potentialCls):
    '''
        .. autoclass:: Lennard Jones Potential
    '''
    name:str = "Lennard Jones Potential"

    c6: float = None
    c12: float = None
    x_shift: float = None
    y_shift: float = None

    def __init__(self, c6: float=0.2, c12: float=0.0001, x_shift: float = 0, y_shift=0):
        '''
        initializes flat well potential class
        '''
        super().__init__()
        self.c6 = c6
        self.c12 = c12
        x_shift = x_shift
        y_shift = y_shift

    def _calculate_energies(self, positions, *kargs) -> t.List[float]:
        return [(self.c12 / pos ** 12) - (self.c6 / pos ** 6) for pos in positions]

    def _calculate_dhdpos(self, positions: (t.List[float] or float), *kargs) -> t.List[float]:
        return [6 * ((2 * self.c12) - (pos ** 6 * self.c6)) / pos ** 13 for pos in positions]


class doubleWellPot1D(_potentialCls):
    '''
        .. autoclass:: unperturbed double well potential
    '''
    name:str = "Double Well Potential"

    a = None
    b = None
    Vmax = None

    def __init__(self, Vmax=100.0, a=0.0, b=0.5):
        '''
        initializes double well potential
        :param Vmax:
        :param a:
        :param b:
        '''
        super(_potentialCls).__init__()
        self.Vmax = Vmax
        self.a = a
        self.b = b


    def _calculate_energies(self, positions:t.List[float], *kargs) -> (t.List[float]):
        return [self.Vmax / self.b ** 4 * ((pos - self.a / 2) ** 2 - self.b ** 2) ** 2 for pos in positions]

    def _calculate_dhdpos(self, positions:t.List[float], *kargs) -> (t.List[float]):
        return [4 * self.Vmax / self.b ** 4 * ((pos - self.a / 2) ** 2 - self.b ** 2) * (pos - self.a / 2) for pos in positions]

"""
    PERTURBED POTENTIALS
"""

class _perturbedPotentialCls(_potentialCls):
    """
        .. autoclass:: perturbedPotentialCls
    """
    lam:float

    def __init__(self, lam:float=0.0):
        '''
        Initializes a potential of the form V = 0.5 * (1 + alpha * lam) * fc * (pos - gamma * lam) ** 2
        :param fc: force constant
        :param alpha: perturbation parameter for width of harmonic oscillator
        :param gamma: perturbation parameter for position of harmonic oscillator
        '''
        super(_potentialCls).__init__()
        self.lam = lam

    def _calculate_energies(self, positions:t.List[float], lam:float=1.0):
        raise NotImplementedError("Function " + __name__ + " was not implemented for class " + str(__class__) + "")

    def _calculate_dhdpos(self, positions:t.List[float], lam:float=1.0):
        raise NotImplementedError("Function " + __name__ + " was not implemented for class " + str(__class__) + "")

    def _calculate_dhdlam(self, positions:t.List[float], lam:float=1.0):
        raise NotImplementedError("Function " + __name__ + " was not implemented for class " + str(__class__) + "")

    def set_lam(self, lam:float):
        self.lam = lam

    def ene(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        '''
        calculates energy of particle
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: energy
        '''
        positions = self._check_positions_type(positions)
        enes =  self._calculate_energies(positions)
        return enes


    def dhdpos(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        '''
        calculates derivative with respect to position
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: derivative dh/dpos
        '''
        positions = self._check_positions_type(positions)
        return self._calculate_dhdpos(positions)

    def dhdlam(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        '''
        calculates derivative with respect to lambda value
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: derivative dh/dlan
        '''
        positions = self._check_positions_type(positions)
        return self._calculate_dhdlam(positions)

class linCoupledHosc(_perturbedPotentialCls):
    def __init__(self, ha=harmonicOsc1D(fc=1.0, x_shift=0.0), hb=harmonicOsc1D(fc=11.0, x_shift=0.0), lam=0):
        super(_potentialCls).__init__()

        self.ha = ha
        self.hb = hb
        self.lam = lam
        self.couple_H = lambda Vab_t: (1.0 - self.lam) * Vab_t[0] + self.lam * Vab_t[1]


    def _calculate_energies(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        return list(map(self.couple_H,  zip(self.ha.ene(positions), self.hb.ene(positions))))

    def _calculate_dhdpos(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        return list(map(self.couple_H,  zip(self.ha.dhdpos(self.lam, positions), self.hb.dhdpos(self.lam, positions))))

    def _calculate_dhdlam(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        return list(map(lambda Va_t, Vb_t: Vb_t-Va_t, self.ha.dhdpos(self.lam, positions), self.hb.dhdpos(self.lam, positions)))



class expCoupledHosc(_perturbedPotentialCls):
    def __init__(self, ha=harmonicOsc1D(fc=1.0, x_shift=0.0), hb=harmonicOsc1D(fc=11.0, x_shift=0.0), s=1.0, temp=300.0, lam:float = 0.0):
        super(_potentialCls).__init__()

        self.ha = ha
        self.hb = hb
        self.s = s
        self.beta = const.gas_constant/1000.0 * temp
        self.lam = lam

        self.couple_H = lambda V_t: -1.0 / (self.beta * self.s) * np.log(
            self.lam * np.exp(-self.beta * self.s * V_t[0]) + (1.0 - self.lam) * np.exp(-self.beta * self.s * V_t[1]))
        self.couple_H_dhdpos = lambda V_t: self.lam * V_t[1] * np.exp(-self.beta * self.s * V_t[1]) + (1.0 - self.lam) * V_t[0] * np.exp(-self.beta * self.s *V_t[0])
        self.couple_H_dhdlam = lambda V_t: -1.0 / (self.beta * self.s) * 1.0 / (self.lam * np.exp(-self.beta * self.s * V_t[1]) + (1.0 - self.lam) * np.exp(-self.beta * self.s * V_t[0])) * (np.exp(-self.beta * self.s * V_t[1]) - np.exp(-self.beta * self.s *V_t[0]))

    def _calculate_energies(self, positions:(t.List[float] or float)) -> (t.List[float] or float):
        return list(map(self.couple_H, zip(self.ha.ene(positions),self.hb.ene(positions))))

    def _calculate_dhdpos(self, positions:(t.List[float] or float), lam:float=1.0) -> (t.List[float] or float):
        return list(map(self.couple_H_dhdpos, zip(self.ha.ene(positions),self.hb.ene(positions))))

    def _calculate_dhdlam(self, positions:(t.List[float] or float), lam:float=1.0) -> (t.List[float] or float):
        return list(map(self.couple_H_dhdlam, zip(self.ha.ene(positions),self.hb.ene(positions))))

class pertHarmonicOsc1D(_perturbedPotentialCls):
    """
        .. autoclass:: pertHarmonixsOsc1D
    """
    name = "perturbed Harmonic Oscilator"
    def __init__(self, fc=1.0, alpha=10.0, gamma=0.0, lam:float=0.0):
        '''
        Initializes a potential of the form V = 0.5 * (1 + alpha * lam) * fc * (pos - gamma * lam) ** 2
        :param fc: force constant
        :param alpha: perturbation parameter for width of harmonic oscillator
        :param gamma: perturbation parameter for position of harmonic oscillator
        '''
        super(_perturbedPotentialCls).__init__()

        self.fc = fc
        self.alpha = alpha
        self.gamma = gamma
        self.lam = lam

    def _calculate_energies(self, positions:(t.List[float] or float)) -> t.List[float]:
        return [0.5 * (1.0 + self.alpha * self.lam) * self.fc * (pos - self.gamma * self.lam) ** 2 for pos in positions]

    def _calculate_dhdpos(self, positions:(t.List[float] or float)) -> t.List[float]:
        return [(1.0 + self.alpha * self.lam) * self.fc * (pos - self.gamma * self.lam) for pos in positions]

    def _calculate_dhdlam(self, positions:(t.List[float] or float))  -> t.List[float]:
        return [0.5 * self.alpha * self.fc * (pos - self.gamma * self.lam) ** 2 - (
                1.0 + self.alpha * self.lam) * self.fc * self.gamma * (pos - self.gamma * self.lam) for pos in positions]


"""
    ENVELOPED POTENTIALS
"""
class envelopedPotential(_potentialCls):
    """
    .. autoclass:: envelopedPotential
    """
    V_is:t.List[_potentialCls]=None
    E_is:t.List[float]=None
    numStates:int=None
    s:float=None

    def __init__(self, V_is: t.List[_potentialCls], s: float = 1.0, Eoff_i: t.List[float] = None):
        """

        :param V_is:
        :param s:
        :param Eoff_i:
        """
        super(_potentialCls).__init__()
        self.numStates = len(V_is)
        if (self.numStates < 2):
            raise IOError("It does not make sense enveloping less than two potentials!")
        if (Eoff_i == None):
            Eoff_i = [0.0 for state in range(len(V_is))]
        elif (len(Eoff_i) != self.numStates):
            raise IOError(
                "Energy offset Vector and state potentials don't have the same length!\n states in Eoff " + str(
                    len(Eoff_i)) + "\t states in Vi" + str(len(V_is)))

        #Todo: think about n states with each m dims.
        self.nDim = V_is[0].nDim
        if(any([V.nDim != self.nDim for V in V_is])):
            raise Exception("Not all endstates have the same dimensionality! This is not imnplemented.\n Dims: "+str([V.nDim != nDim for V in V_is]))

        self.V_is = V_is
        self.s = s
        self.Eoff_i = Eoff_i

    #each state gets a position list
    def _check_positions_type(self, positions:t.List[float])->t.List[float]:
        if (type(positions) in [float, int, str]):
            if(len(positions) != self.numStates):
                positions = [positions for state in range(self.numStates + 1)]
            else:
                positions = [float(positions) for state in range(self.numStates+1)]
        elif(isinstance(positions, Iterable)):
            if(len(positions)!= self.numStates):
                if (isinstance(positions[0], Iterable)):    # Ndimensional case
                    positions = [list(map(lambda dimlist: np.array(list(map(float, dimlist))), positions)) for state in range(self.numStates)]
                else:   #onedimensional
                    positions = [list(map(float, positions)) for state in range(self.numStates)]
            else:   #TODO: insert check here! for fitting numstates
                return positions
        else:
            raise Exception("This is an unknown type of Data structure: "+str(type(positions))+"\n"+str(positions))
        return positions

    def _calculate_energies(self, positions:(t.List[float] or float), *kargs) -> list:
        partA = [-self.s*(Vit-self.Eoff_i[0]) for Vit in self.V_is[0].ene(positions[0])]
        partB = [-self.s*(Vit-self.Eoff_i[1]) for Vit in self.V_is[1].ene(positions[1])]
        sum_prefactors = [max(A_t,B_t)+math.log(1+math.exp(min(A_t,B_t)-max(A_t,B_t))) for A_t, B_t in zip(partA, partB)]

        #more than two states!
        for state in range(2, self.numStates):
            partN = [-self.s * (Vit - self.Eoff_i[state]) for Vit in self.V_is[state].ene(positions[state])]
            sum_prefactors = [max(sum_prefactors_t, N_t) + math.log(1 + math.exp(min(sum_prefactors_t, N_t) - max(sum_prefactors_t, N_t))) for sum_prefactors_t, N_t in zip(sum_prefactors, partN)]

        Vr = [-1/float(self.s)*partitionF for partitionF in sum_prefactors]
        return Vr

    
    def _calculate_dhdpos(self,  positions:(t.List[float] or float), *kargs):
        ###CHECK!THIS FUNC!!! not correct
        V_R_ene = self.ene(positions)
        V_Is_ene = [statePot.ene(state_pos) for statePot, state_pos in zip(self.V_is, positions)]
        V_Is_dhdpos = [statePot.dhdpos(state_pos) for statePot, state_pos in zip(self.V_is, positions)]
        dhdpos=[]

        for position in range(len(positions[0])):
            dhdpos_R = 0
            dhdpos_state=[]
            for V_state_ene, V_state_dhdpos in zip(V_Is_ene, V_Is_dhdpos):
                #den = sum([math.exp(-const.k *Vstate[position]) for Vstate in V_Is_ene])
                #prefactor = (math.exp(-const.k *V_state_ene[position]))/den if (den!=0) else 0
                if(V_state_ene[position]==0): 
                    prefactor = 0 
                else:
                    prefactor = 1-(V_state_ene[position]/(sum([ Vstate[position] for Vstate in V_Is_ene]))) if (sum([ Vstate[position] for Vstate in V_Is_ene])!=0) else 0
                #print(round(positions[0][position],2),"\t",round(prefactor,2),"\t" , round(V_state_dhdpos[position]), "\t", round(V_R_ene[position]))
                dhdpos_state.append(prefactor*V_state_dhdpos[position])
                dhdpos_R += prefactor*V_state_dhdpos[position]
                dlist = [dhdpos_R]
                dlist.extend(dhdpos_state)
            dhdpos.append(dlist)   
        return dhdpos

class envelopedDoubleWellPotential(envelopedPotential):
    def __init__(self, y_shifts: list = None, x_shifts=None,
                 smoothing: float = 1, fcs=None):
        if (y_shifts == None):
            y_shifts = [0, 0]
        if (x_shifts == None):
            x_shifts = [-1, 1]
        if (fcs == None):
            fcs = [1, 1]

        V_is = [harmonicOsc1D(x_shift=x_shift, y_shift=y_shift, fc=fc)
                for y_shift, x_shift, fc in zip(y_shifts, x_shifts, fcs)]
        super().__init__(V_is=V_is, s=smoothing)


class envelopedPotentialMultiS(envelopedPotential):
    """
    .. autoclass:: envelopedPotential
    """
    V_is:t.List[_potentialCls] = None
    E_is:t.List[float] = None
    numStates:int = None
    s:t.List[float] = None

    def __init__(self, V_is: t.List[_potentialCls], s: t.List[float], Eoff_i: t.List[float] = None):
        """

        :param V_is:
        :param s:
        :param Eoff_i:
        """
        super(_potentialCls).__init__()
        self.numStates = len(V_is)
        if (self.numStates < 2):
            raise IOError("It does not make sense enveloping less than two potentials!")
        if (Eoff_i == None):
            Eoff_i = [0.0 for state in range(len(V_is))]

        elif (len(Eoff_i) != self.numStates):
            raise IOError(
                "Energy offset Vector and state potentials don't have the same length!\n states in Eoff " + str(
                    len(Eoff_i)) + "\t states in Vi" + str(len(V_is)))

        # Todo: think about n states with each m dims.
        self.nDim = V_is[0].nDim
        if (any([V.nDim != self.nDim for V in V_is])):
            raise Exception("Not all endstates have the same dimensionality! This is not imnplemented.\n Dims: " + str(
                [V.nDim != nDim for V in V_is]))

        self.V_is = V_is
        self.s = s
        self.Eoff_i = Eoff_i

    def _calculate_energies(self, positions: (t.List[float] or float), *kargs) -> list:
        partA = [-self.s[0] * (Vit - self.Eoff_i[0]) for Vit in self.V_is[0].ene(positions[0])]
        partB = [-self.s[1] * (Vit - self.Eoff_i[1]) for Vit in self.V_is[1].ene(positions[1])]
        sum_prefactors = [max(A_t, B_t) + math.log(1 + math.exp(min(A_t, B_t) - max(A_t, B_t))) for A_t, B_t in
                          zip(partA, partB)]

        # more than two states!
        for state in range(2, self.numStates):
            partN = [-self.s[state] * (Vit - self.Eoff_i[state]) for Vit in self.V_is[state].ene(positions[state])]
            sum_prefactors = [max(sum_prefactors_t, N_t) + math.log(
                1 + math.exp(min(sum_prefactors_t, N_t) - max(sum_prefactors_t, N_t))) for sum_prefactors_t, N_t in
                              zip(sum_prefactors, partN)]

        Vr = [- partitionF for state,partitionF in enumerate(sum_prefactors)]    #1/ float(self.s[0]) *
        return Vr

    def _calculate_dhdpos(self, positions: (t.List[float] or float), *kargs):
        ###CHECK!THIS FUNC!!! not correct
        V_R_ene = self.ene(positions)
        V_Is_ene = [statePot.ene(state_pos) for statePot, state_pos in zip(self.V_is, positions)]
        V_Is_dhdpos = [statePot.dhdpos(state_pos) for statePot, state_pos in zip(self.V_is, positions)]
        dhdpos = []

        for position in range(len(positions[0])):
            dhdpos_R = 0
            dhdpos_state = []
            for V_state_ene, V_state_dhdpos in zip(V_Is_ene, V_Is_dhdpos):
                # den = sum([math.exp(-const.k *Vstate[position]) for Vstate in V_Is_ene])
                # prefactor = (math.exp(-const.k *V_state_ene[position]))/den if (den!=0) else 0
                if (V_state_ene[position] == 0):
                    prefactor = 0
                else:
                    prefactor = 1 - (V_state_ene[position] / (sum([Vstate[position] for Vstate in V_Is_ene]))) if (
                                sum([Vstate[position] for Vstate in V_Is_ene]) != 0) else 0
                # print(round(positions[0][position],2),"\t",round(prefactor,2),"\t" , round(V_state_dhdpos[position]), "\t", round(V_R_ene[position]))
                dhdpos_state.append(prefactor * V_state_dhdpos[position])
                dhdpos_R += prefactor * V_state_dhdpos[position]
                dlist = [dhdpos_R]
                dlist.extend(dhdpos_state)
            dhdpos.append(dlist)
        return dhdpos


if __name__ == "__main__":
    pass

