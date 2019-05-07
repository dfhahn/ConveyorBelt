import numpy as np
import scipy.constants as const

class potential:
    '''
    potential base class
    '''

    def __init__(self):
        return

    def ene(self, lam, pos) -> float:
        '''
        calculates energy of particle
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: energy
        '''
        return 0.0

    def dhdpos(self, lam, pos) -> float:
        '''
        calculates derivative with respect to position
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: derivative dh/dpos
        '''
        return 0.0

    def dhdlam(self, lam, pos) -> float:
        '''
        calculates derivative with respect to lambda value
        :param lam: alchemical parameter lambda
        :param pos: position on 1D potential energy surface
        :return: derivative dh/dlan
        '''
        return 0.0


class harmonicOsc1D(potential):
    '''
    unperturbed harmonic oscillator potential
    '''
    shift = None
    fc = None

    def ene(self, lam, pos):
        return 0.5 * self.fc * (pos - self.shift) ** 2

    def dhdpos(self, lam, pos):
        return self.fc * (pos - self.shift)

    @staticmethod
    def dhdlam(lam, pos):
        return 0.0

    def __init__(self, fc=1.0, shift=0.0):
        '''
        initializes harmonicOsc1D class
        :param fc: force constant
        :param shift: minimum position of harmonic oscillator
        '''
        super(harmonicOsc1D, self).__init__()
        self.fc = fc
        self.shift = shift


class doubleWellPot1D(potential):
    '''
    unperturbed double well potential
    '''
    a = None
    b = None
    Vmax = None

    def ene(self, lam, pos):
        return self.Vmax / self.b ** 4 * ((pos - self.a / 2) ** 2 - self.b ** 2) ** 2

    def dhdq(self, lam, pos):
        return 4 * self.Vmax / self.b ** 4 * ((pos - self.a / 2) ** 2 - self.b ** 2) * (pos - self.a / 2)

    @staticmethod
    def dhdlam(lam, pos):
        return 0.0

    def __init__(self, Vmax=100.0, a=0.0, b=0.5):
        '''
        initializes double well potential
        :param Vmax:
        :param a:
        :param b:
        '''
        self.Vmax = Vmax
        self.a = a
        self.b = b


class pertHarmonicOsc1D(potential):
    def ene(self, lam, pos):
        return 0.5 * (1.0 + self.alpha * lam) * self.fc * (pos - self.gamma * lam) ** 2

    def dhdpos(self, lam, pos):
        return (1.0 + self.alpha * lam) * self.fc * (pos - self.gamma * lam)

    def dhdlam(self, lam, pos):
        return 0.5 * self.alpha * self.fc * (pos - self.gamma * lam) ** 2 - (
                1.0 + self.alpha * lam) * self.fc * self.gamma * (pos - self.gamma * lam)

    def __init__(self, fc=1.0, alpha=10.0, gamma=0.0):
        '''
        Initializes a potential of the form V = 0.5 * (1 + alpha * lam) * fc * (pos - gamma * lam) ** 2
        :param fc: force constant
        :param alpha: perturbation parameter for width of harmonic oscillator
        :param gamma: perturbation parameter for position of harmonic oscillator
        '''
        self.fc = fc
        self.alpha = alpha
        self.gamma = gamma


class linCoupledHosc(potential):
    def ene(self, lam, pos):
        return (1.0 - lam) * self.ha.ene(lam, pos) + lam * self.hb.ene(lam, pos)

    def dhdpos(self, lam, pos):
        return (1.0 - lam) * self.ha.dhdpos(lam, pos) + lam * self.hb.dhdpos(lam, pos)

    def dhdlam(self, lam, pos):
        return self.hb.ene(lam, pos) - self.ha.ene(lam, pos)

    def __init__(self, ha=harmonicOsc1D(fc=1.0, shift=0.0), hb=harmonicOsc1D(fc=11.0, shift=0.0)):
        self.ha = ha
        self.hb = hb


class expCoupledHosc(potential):
    def ene(self, lam, pos):
        return -1.0 / (self.beta * self.s) * np.log(
            lam * np.exp(-self.beta * self.s * self.hb.ene(lam, pos)) + (1.0 - lam) * np.exp(
                -self.beta * self.s * self.ha.ene(lam, pos)))

    def dhdpos(self, lam, pos):
        return 1.0 / (lam * np.exp(-self.beta * self.s * self.hb.ene(lam, pos)) + (1.0 - lam) *
                      np.exp(-self.beta * self.s * self.ha.ene(lam, pos))) * \
               (lam * self.hb.dhdpos(lam, pos) * np.exp(-self.beta * self.s * self.hb.ene(lam, pos)) + (1.0 - lam) *
                self.ha.dhdpos(lam, pos) * np.exp(-self.beta * self.s * self.ha.ene(lam, pos)))

    def dhdlam(self, lam, pos):
        return -1.0 / (self.beta * self.s) * 1.0 / (
                lam * np.exp(-self.beta * self.s * self.hb.ene(lam, pos)) + (1.0 - lam) * np.exp(
            -self.beta * self.s * self.ha.ene(lam, pos))) * (
                       np.exp(-self.beta * self.s * self.hb.ene(lam, pos)) - np.exp(
                   -self.beta * self.s * self.ha.ene(lam, pos)))

    def __init__(self, ha=harmonicOsc1D(fc=1.0, shift=0.0), hb=harmonicOsc1D(fc=11.0, shift=0.0), s=1.0, temp=300.0):
        self.ha = ha
        self.hb = hb
        self.s = s
        self.beta = const.gas_constant/1000.0 * temp