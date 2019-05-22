"""
Module: integrator
    This module shall be used to implement subclasses of integrator. The integrators are use for propagating simulatoins.
"""
#import newSystem1D as sys
import numpy as np
import scipy.constants as const

class integrator:

    def __init__(self):
        raise NotImplementedError("This "+__class__+" class is not implemented")
    
    def step(self, system:None):
        raise NotImplementedError("The function step in "+__class__+" class is not implemented")
    
class monteCarloIntegrator(integrator):
    posShift:float = 0
    maxStepSize:float = None
    spaceRange:tuple = None
    resolution:float = 0.0001
    _critInSpaceRange = lambda self,pos: self.spaceRange != None and pos >= min(self.spaceRange) and pos <= max(self.spaceRange)
    
    def __init__(self, maxStepSize:float=None, spaceRange:tuple=None):
        self.maxStepSize = maxStepSize
        self.spaceRange = spaceRange
        pass
    
    def step(self, system):
        while(True):    #while no value in spaceRange was found, terminates in first run if no spaceRange
            #integrate
            current_state = system.currentState
            self.oldpos = current_state.position
            self.randomShift()
            self.newPos = self.oldpos + self.posShift
            
            #only get positions in certain range
            if(self._critInSpaceRange(self.newPos)):
                break
            else:
                self.newPos = self.oldpos           #reject step outside of range
        return self.newPos, None, None
    
    def randomShift(self):
        if(self.spaceRange!=None):   #check for space range, that posShift is not bigger than space Range > converges faster in 
            self.posShift = np.random.randint(low=self.spaceRange[0]/self.resolution, high=self.spaceRange[1]/self.resolution)*self.resolution
        else:
            self.posShift = np.random.rand()
            
        if(self.maxStepSize != None):#is there a maximal step size?
            self.posShift = self.posShift/100*self.maxStepSize

        return self.posShift

    
class metropolisMonteCarloIntegrator(monteCarloIntegrator):
    posShift:float = 0
    maxStepSize:float = None
    spaceRange:tuple = None
    resolution:float = 0.0001
    metropolisCriterion=None
    
    _defaultMetropolisCriterion = lambda self, ene_new, currentState: (ene_new < currentState.totEnergy or 
                                                                       np.random.rand() <= np.exp(-1.0 / (const.gas_constant / 1000.0 * currentState.temperature) * (ene_new - currentState.totPotEnergy)))
    
    def __init__(self, maxStepSize:float=None, spaceRange:tuple=None, metropolisCriterion=None):
        self.maxStepSize = maxStepSize
        self.spaceRange = spaceRange
        
        if(metropolisCriterion == None):
            self.metropolisCriterion = self._defaultMetropolisCriterion
        else:
            self.metropolisCriterion = metropolisCriterion
                
                
    def step(self, system):
        while(True):    #while no value in spaceRange was found, terminates in first run if no spaceRange
            #integrate position
            current_state = system.currentState
            self.oldpos = current_state.position
            self.randomShift()
            self.newPos = self.oldpos + self.posShift

            #MetropolisCriterion
            ene = system.pot()
            if (self._critInSpaceRange(self.newPos) and self.metropolisCriterion(ene, current_state)):
                break
            else:
                self.newPos = self.oldpos  # not accepted
        return self.newPos, None, None
 
    
class newtonIntegrator(integrator):
    def step(self, sys):
        sys.pos = sys.newpos  # t
        sys.force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
        sys.oldvel = sys.newvel  # t - 0.5 Dt
        sys.newvel += sys.force / sys.mu * self.dt  # t+0.5Dt
        sys.vel = 0.5 * (sys.oldvel + sys.newvel)
        sys.newpos += self.dt * sys.newvel  # t+Dt
        sys.veltemp = sys.mu / const.gas_constant / 1000.0 * sys.vel ** 2  # t
        sys.updateEne()
        return [sys.pos, sys.vel, sys.veltemp, sys.totkin, sys.totpot, sys.totene, sys.lam, sys.dhdlam]

    def __init__(self, dt=1e-1):
        self.dt = dt
        raise NotImplementedError("This "+__class__+" class is not implemented")

class nhIntegrator(integrator):
    def scaleVel(self, sys):
        freetemp = 2.0 / const.gas_constant / 1000.0 * sys.mu * sys.newvel ** 2  # t+0.5Dt
        self.oldxi = self.xi  # t-0.5Dt
        self.xi += self.dt / (self.tau * self.tau) * (freetemp / sys.temp - 1.0)  # t+0.5t
        scale = 1.0 - self.xi * self.dt
        return scale

    def step(self, sys):
        sys.pos = sys.newpos  # t
        sys.force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
        sys.oldvel = sys.newvel  # t - 0.5 Dt
        sys.newvel += sys.force / sys.mu * self.dt  # t+0.5t
        sys.newvel *= self.scaleVel(sys)
        sys.vel = 0.5 * (sys.oldvel + sys.newvel)
        sys.newpos += self.dt * sys.newvel  # t+Dt
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
            sys.pos += self.dt * sys.newvel  # t
            force = -sys.pot.dhdpos(sys.lam, sys.pos)  # t
            sys.vel = (oldvel + sys.newvel) / 2.0
            sys.veltemp = sys.mu / const.gas_constant / 1000.0 * sys.vel ** 2  # t
            sys.updateEne()
            oldvel = sys.newvel
            sys.newvel += force / sys.mu * self.dt  # t+0.5t
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
