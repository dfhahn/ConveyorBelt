"""
Module: integrator
    This module shall be used to implement subclasses of integrator. The integrators are use for propagating simulatoins.
    Think about how to link conditions to integrator???
"""
import numpy as np
from typing import Tuple
import scipy.constants as const

#from src.system import system
system = None

class integrator:
    """
    autoclass: integrator
        This class is the parent class to all other classes.
    """
    #general:
    verbose:bool = False
    #Params
    maxStepSize:float = None
    minStepSize:float = None
    spaceRange:tuple = None

    #calculation
    posShift:float = 0

    #Limits:
    _critInSpaceRange = lambda self,pos: self.spaceRange == None or (self.spaceRange != None and pos >= min(self.spaceRange) and pos <= max(self.spaceRange))

    def __init__(self):
        raise NotImplementedError("This "+__class__+" class is not implemented")
    
    def step(self, system:system):
        """
        ..autofunction: step
            This is the parent function that is the interface for all integrator step functions.

        :param system: This is a system, that should be integrated.
        :type system: src.system.system
        :return: (new Position, new velocity, position Shift/ force)
        :rtype: (float, float, float)
        """
        raise NotImplementedError("The function step in "+__class__+" class is not implemented")

    def integrate(self, system:system, steps:int):
        for step in range(steps):
            (newPosition, newVelocity, newForces) = self.step(system=system)
            system.append_state(newPosition=newPosition, newVelocity=newVelocity, newForces=newForces)

    def setVerbose(self, verbose:bool=True):
        self.verbose = verbose

"""
Stochastic Integrators
"""
class monteCarloIntegrator(integrator):
    """
    ..autoclass: monteCarloIntegrator
        This class implements the classic monte carlo integrator.
        It choses its moves purely randomly.
    """
    resolution:float = 0.01   #increase the ammount of different possible values = between 0 and 10 there are 10/0.01 different positions.

    def __init__(self, maxStepSize:float=None, minStepSize:float=None, spaceRange:tuple=None):
        self.maxStepSize = maxStepSize
        self.minStepSize = minStepSize
        self.spaceRange = spaceRange
        pass
    
    def step(self, system)-> Tuple[float, None, float]:
        """
        ..autofunction: step
            This function is performing an integration step in MonteCarlo fashion.
        :param system: This is a system, that should be integrated.
        :type system: src.system.system
        :return: (new Position, None, position Shift)
        :rtype: (float, None, float)
        """
        # integrate
        # while no value in spaceRange was found, terminates in first run if no spaceRange
        while(True):
            current_state = system.currentState
            self.oldpos = current_state.position
            self.randomShift()
            self.newPos = self.oldpos + self.posShift
            
            #only get positions in certain range or accept if no range
            if(self._critInSpaceRange(self.newPos)):
                break
            else:
                self.newPos = self.oldpos           #reject step outside of range
        return self.newPos, None, self.posShift
    
    def randomShift(self, nDim)->float:
        """
        ..autofunction: randomShift
            This function calculates the shift for the current position.

        :return: position shift
        :rtype: float
        """
        #which sign will the shift have?
        sign = [-1 if(np.random.randint(low=0, high=100) <50) else 1 for x in range(nDim)]
        #Check if there is a space restriction? - converges faster
        if(self.spaceRange!=None):
            shift = [abs(np.random.randint(low=self.spaceRange[0]/self.resolution, high=self.spaceRange[1]/self.resolution)*self.resolution) for x in range(nDim)]
        else:
            shift = [abs(np.random.rand()) for x in range(nDim)]

        #Is the step shift in the allowed area?
        if(self.maxStepSize != None and shift > self.maxStepSize):#is there a maximal step size?
            self.posShift = np.multiply(sign, self.maxStepSize)
        elif(self.minStepSize != None and shift < self.minStepSize):
            self.posShift = np.multiply(sign, self.minStepSize)
        else:
            self.posShift = np.multiply(sign, shift)

        return self.posShift

class metropolisMonteCarloIntegrator(monteCarloIntegrator):
    """
    ..autoclass: metropolisMonteCarloInegrator
        This class is implementing a metropolis monte carlo Integrator.
        In opposite to the Monte Carlo Integrator, that is completley random, this integrator has limitations to the randomness.
        Theis limitation is expressed in the Metropolis Criterion.

        There is a standard Metropolis Criterion implemented, but it can also be exchanged with a different one.

        Default Metropolis Criterion:
            $ decision =  (E_{t} < E_{t-1}) ||  ( rand <= e^{(-1/(R/T*1000))*(E_t-E_{t-1})}$
            with:
                - $R$ as universal gas constant

        The original Metropolis Criterion (Nicholas Metropolis et al.; J. Chem. Phys.; 1953 ;doi: https://doi.org/10.1063/1.1699114):

            $ p_A(E_{t}, E_{t-1}, T) = min(1, e^{-1/(k_b*T) * (E_{t} - E_{t-1})})
            $ decision:  True if( 0.5 < p_A(E_{t}, E_{t-1}, T)) else False
            with:
                - $k_b$ as Boltzmann Constant
    """
    #
    #Parameters:
    metropolisCriterion=None    #use a different Criterion
    randomnessIncreaseFactor:float = 1  #tune randomness of your results
    maxIterationTillAccept:float = 100  #how often shall the integrator iterate till it accepts a step forcefully

    #METROPOLIS CRITERION
    ##random part of Metropolis Criterion:
    _defaultRandomness = lambda self, ene_new, currentState: ((1/self.randomnessIncreaseFactor)*np.random.rand() <= np.exp(-1.0 / (const.gas_constant / 1000.0 * currentState.temperature) * (ene_new - currentState.totPotEnergy)))
    ##default Metropolis Criterion
    _defaultMetropolisCriterion = lambda self, ene_new, currentState: (ene_new < currentState.totEnergy or self._defaultRandomness(ene_new, currentState))
    ## original criterion not useful causes overflows:
    #_defaultMetropolisCriterion = lambda self, ene_new, currentState: True if(0.5 > min(1, np.e**(-1/(const.k * currentState.temperature)*(ene_new-currentState.totPotEnergy)))) else False

    def __init__(self, minStepSize:float=None, maxStepSize:float=None, spaceRange:tuple=None, metropolisCriterion=None, randomnessIncreaseFactor=1, maxIterationTillAccept:int=100):
        self.maxStepSize = maxStepSize
        self.minStepSize = minStepSize
        self.spaceRange = spaceRange
        self.randomnessIncreaseFactor = randomnessIncreaseFactor
        self.maxIterationTillAccept = maxIterationTillAccept
        if(metropolisCriterion == None):
            self.metropolisCriterion = self._defaultMetropolisCriterion
        else:
            self.metropolisCriterion = metropolisCriterion

    def step(self, system):
        """
        ..autofunction: step
            This function is performing an integration step in MetropolisMonteCarlo fashion.
        :param system: This is a system, that should be integrated.
        :type system: src.system.system
        :return: (new Position, None, position Shift)
        :rtype: (float, None, float)
        """

        iterstep = 0
        current_state = system.currentState
        self.oldpos = current_state.position
        nDim = system.nDim
        # integrate position
        while(True):    #while no value in spaceRange was found, terminates in first run if no spaceRange
            self.randomShift(nDim)
            #eval new Energy
            system._currentPosition = np.add(self.oldpos, self.posShift)
            ene = system.totPot()

            #MetropolisCriterion
            if ((self._critInSpaceRange(system._currentPosition) and self.metropolisCriterion(ene, current_state)) or iterstep==self.maxIterationTillAccept):
                break
            else:   #not accepted
                iterstep += 1
                continue
        return system._currentPosition , None, self.posShift


 #OLD
"""
Newtonian Integrators
"""
class newtonianIntegrator(integrator):
    currentPosition:float
    currentVelocity:float
    currentForces:float

    dt:float

class verlocityVerletIntegrator(newtonianIntegrator):
    """
        .. autoclass:: Verlocity Verlet Integrator,
        is not implemented yet
    """

    pass

class positionVerletIntegrator(newtonianIntegrator):

    def __init__(self, dt=0.0005):
        self.dt = dt

    def step(self, system):
        #init
        currentPosition = system._currentPosition
        currentVelocity = system._currentVelocities

        #calculation:
        newForces = system.potential.dhdpos(currentPosition)[0]    #Todo: make multi particles possible - use current forces!
        velocity_new = currentVelocity - (newForces / system.mass * self.dt)
        new_velocity = velocity_new #0.5 * (currentVelocity + velocity_new)  #update velo
        new_position = currentPosition + new_velocity * self.dt


        if(self.verbose):
            print("INTEGRATOR: current forces\t ", newForces)
            print("INTEGRATOR: current Velocities\t ", currentVelocity)
            print("INTEGRATOR: current_position\t ", currentPosition)

            print("INTEGRATOR: newVel\t ", new_velocity)
            print("INTEGRATOR: newPosition\t ", new_position)
            print("\n")
        return new_position, new_velocity, newForces


class leapFrogIntegrator(newtonianIntegrator):
    def __init__(self, dt=0.0005):
        self.dt = dt

    def step(self, system):
        pass

"""
OLD Integrators:
class newtonIntegrator(integrator): #LEAPFROG
    def step(self, sys):
        sys.pos = sys.newpos  # t
        sys.force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
        sys.oldvel = sys.new_velocity  # t - 0.5 Dt
        sys.new_velocity += sys.force / sys.mu * self.dt  # t+0.5Dt
        sys.vel = 0.5 * (sys.oldvel + sys.new_velocity)
        sys.newpos += self.dt * sys.new_velocity  # t+Dt
        sys.veltemp = sys.mu / const.gas_constant / 1000.0 * sys.vel ** 2  # t
        sys.updateEne()
        return [sys.pos, sys.vel, sys.veltemp, sys.totkin, sys.totpot, sys.totene, sys.lam, sys.dhdlam]

    def __init__(self, dt=1e-1):
        self.dt = dt
        raise NotImplementedError("This "+__class__+" class is not implemented")

class nhIntegrator(integrator): Nosehover-leapfrog
    def scaleVel(self, sys):
        freetemp = 2.0 / const.gas_constant / 1000.0 * sys.mu * sys.new_velocity ** 2  # t+0.5Dt
        self.oldxi = self.xi  # t-0.5Dt
        self.xi += self.dt / (self.tau * self.tau) * (freetemp / sys.temp - 1.0)  # t+0.5t
        scale = 1.0 - self.xi * self.dt
        return scale

    def step(self, sys):
        sys.pos = sys.newpos  # t
        sys.force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
        sys.oldvel = sys.new_velocity  # t - 0.5 Dt
        sys.new_velocity += sys.force / sys.mu * self.dt  # t+0.5t
        sys.new_velocity *= self.scaleVel(sys)
        sys.vel = 0.5 * (sys.oldvel + sys.new_velocity)
        sys.newpos += self.dt * sys.new_velocity  # t+Dt
        sys.veltemp = sys.mu / const.gas_constant/1000.0 * sys.vel ** 2  # t
        sys.updateEne()
        return [sys.pos, sys.vel, sys.veltemp, sys.totkin, sys.totpot, sys.totene, sys.lam, sys.dhdlam]

    def __init__(self, xi=0.0, tau=0.1, dt=1e-1):
        self.xi = xi
        self.oldxi = xi
        self.tau = tau
        self.dt = dt
        raise NotImplementedError("This "+__class__+" class is not implemented")

class hmcIntegrator(integrator):
    def step(self, sys):
        accept = 0
        oldene = sys.totene
        oldpos = sys.pos  # t
        oldvel = sys.vel
        sys.initVel()  # t-0.5Dt
        for i in range(self.steps):
            sys.pos += self.dt * sys.new_velocity  # t
            force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
            sys.vel = (oldvel + sys.new_velocity) / 2.0
            sys.veltemp = sys.mu / const.gas_constant / 1000.0 * sys.vel ** 2  # t
            sys.updateEne()
            oldvel = sys.new_velocity
            sys.new_velocity += force / sys.mu * self.dt  # t+0.5t
        accept = 0
        if sys.totene < oldene:
            accept = 1
        else:
            if np.random.rand() <= np.exp(-1 / (const.gas_constant / 1000.0* sys.temp) * (sys.totene - oldene)):
                accept = 1
        if accept == 0:
            sys.pos = oldpos
            sys.vel = oldvel
            sys.veltemp = sys.mu / const.gas_constant * sys.vel ** 2  # t
            sys.updateEne()
        return [sys.pos, sys.vel, sys.veltemp, sys.totkin, sys.totpot, sys.totene, sys.lam, sys.dhdlam]

    def __init__(self, steps=5, dt=1e-1):
        self.steps = steps
        self.dt = dt
        raise NotImplementedError("This "+__class__+" class is not implemented")
"""
