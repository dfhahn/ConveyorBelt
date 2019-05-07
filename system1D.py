import numpy as np
import scipy.constants as const
import potential1D

class system1D:
    '''
    Class of a (possibly perturbed) 1D system on a
    potential energy surface defined by a potential1D.potential object.
    Different integrators/propagators can be chosen to create an ensemble of states.
    '''
    def kin(self):
        return 0.5 * self.mu * self.vel ** 2

    @staticmethod
    def randomPos():
        pos = np.random.rand() * 20.0 - 10.0
        return pos

    def randomShift(self):
        posShift = (np.random.rand() * 2.0 - 1.0) * self.posShift
        return posShift

    def initVel(self):
        self.vel = np.sqrt(const.gas_constant / 1000.0 * self.temp / self.mu) * np.random.normal()
        self.newvel = self.vel
        self.veltemp = self.mu / const.gas_constant / 1000.0 * self.vel ** 2  # t
        return self.vel

    def update(self, lam):
        self.lam = lam
        self.omega = np.sqrt((1.0 + self.alpha * self.lam) * self.fc / self.mu)
        self.updateEne()

    def updateTemp(self, temp):
        self.temp = temp
        self.posShift = np.sqrt(0.5 * (1.0 + self.alpha * 0.5) * self.fc / (const.gas_constant / 1000.0 * self.temp))
        self.updateEne()

    def updateEne(self):
        self.totpot = self.pot.ene(self.lam, self.pos)
        self.dhdlam = self.pot.dhdlam(self.lam, self.pos)
        self.totkin = self.kin()
        self.totene = self.totpot + self.totkin
        self.redene = self.totpot / (const.gas_constant / 1000.0 * self.temp)

    def getState(self):
        return [self.pos, self.temp, 0, 0, self.totpot, self.totene, self.lam, self.dhdlam]

    def propagate(self):
        return self.integrator.step(self)

    def revert(self):
        self.pos = self.oldpos
        self.updateEne()
        return

    class sdIntegrator:
        @staticmethod
        def step(sys):
            sys.oldpos = sys.pos
            sys.pos += sys.randomShift()
            ene = sys.pot.ene(sys.lam, sys.pos)
            if ene < sys.totpot:
                sys.updateEne()
            elif np.random.rand() <= np.exp(-1.0 / (const.gas_constant / 1000.0 * sys.temp) * (ene - sys.totpot)):
                sys.updateEne()
            else:
                sys.pos = sys.oldpos  # not accepted
            return sys.getState()

    class newtonIntegrator:
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

    class nhIntegrator:
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

    class hmcIntegrator:
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

    def __init__(self, temp=300.0, fc=1.0, lam=0.0, alpha=10.0, mass=None, potential=potential1D.potential(),
                 integrator='sd'):
        self.temp = temp
        self.fc = fc
        self.lam = lam
        self.alpha = alpha
        self.posShift = np.sqrt(0.5 * (1.0 + self.alpha * 0.5) * self.fc / (const.gas_constant /1000.0 * self.temp))
        if mass is None:
            mass = [1.0, 1.0]
        self.mu = mass[0] * mass[1] / (mass[0] + mass[1])
        self.omega = np.sqrt(
            (1.0 + self.alpha * self.lam) * self.fc / const.Avogadro * 1e21 / (self.mu * const.atomic_mass))
        self.newpos = 0.0
        self.oldpos = 0.0
        self.pos = self.newpos
        self.vel = None
        self.newvel = None
        self.veltemp = None
        self.totpot = None
        self.dhdlam = None
        self.totkin = None
        self.totene = None
        self.redene = None
        self.oldxi = None
        self.initVel()
        self.pot = potential
        self.updateEne()
        if integrator == 'sd':
            self.vel = 0.0
            self.integrator = self.sdIntegrator()
        elif integrator == 'nose-hoover':
            self.integrator = self.nhIntegrator()
            self.propagate()
        elif integrator == 'newton':
            self.integrator = self.newtonIntegrator()
            self.propagate()
        elif integrator == 'hmc':
            self.integrator = self.hmcIntegrator()
            self.propagate()

