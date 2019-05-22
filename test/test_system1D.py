import unittest
from src import system, potential1D as pot


class testHarmonicOsc1D(unittest.TestCase):
    def testFc(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0
        for fc in [1.0, 1.5, 10.0, 20.0]:
            with self.subTest(fc=fc):
                sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                      potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha))
                sys.pos = 1.0
                sys.updateEne()
                self.assertAlmostEqual(sys.totpot, fc)
                self.assertAlmostEqual(sys.pot.dhdpos(sys.lam, sys.pos), 2.0 * fc)
                self.assertAlmostEqual(sys.pot.dhdl(sys.lam, sys.pos), 0.5 * fc)

    def testAlpha(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0

        for alpha in [1.0, 10.0, 100.0]:
            with self.subTest(alpha=alpha):
                sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                      potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha))
                sys.pos = 1.0
                sys.updateEne()
                self.assertAlmostEqual(sys.totpot, 0.5 + 0.5 * alpha)
                self.assertAlmostEqual(sys.pot.dhdpos(sys.lam, sys.pos), 1.0 + alpha)
                self.assertAlmostEqual(sys.pot.dhdl(sys.lam, sys.pos), 0.5 * alpha)

    def testLam(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0

        for lam in [0.0, 0.5, 1.0]:
            with self.subTest(lam=lam):
                sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                      potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha))
                sys.pos = 1.0
                sys.updateEne()
                self.assertAlmostEqual(sys.totpot, 0.5 + 0.5 * lam)
                self.assertAlmostEqual(sys.pot.dhdpos(sys.lam, sys.pos), 1.0 + lam)
                self.assertAlmostEqual(sys.pot.dhdl(sys.lam, sys.pos), 0.5)

    def testKin(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0

        for temp in [250, 300, 350]:
            with self.subTest(temp=temp):
                for sd in [0, 1]:
                    for nose in [0, 1]:
                        for mass in [[1.0, 1.0], [1.0, 2.0], [2.0, 2.0]]:
                            sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                                  potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha))

    def testNoseHoover(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0
        for vel in [1.0, 2.0, 3.0]:
            with self.subTest(vel=vel):
                sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                      potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha))
                sys.vel = vel
                sys.updateEne()
                self.assertAlmostEqual(sys.totkin, vel ** 2 / 4.0)
                sys.propagate()

    def testNewton(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0
        for vel in [1.0, 2.0, 3.0]:
            with self.subTest(vel=vel):
                sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                      potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha), integrator='newton')
                sys.pos = 1.0
                sys.newpos = 1.0
                sys.newvel = vel
                sys.vel = vel
                sys.updateEne()
                self.assertAlmostEqual(sys.totkin, vel ** 2 / 4.0)
                sys.propagate()
                self.assertAlmostEqual(sys.newvel, vel - 0.4)
                self.assertAlmostEqual(sys.vel, vel - 0.2)
                self.assertAlmostEqual(sys.newpos, 1.0 + (vel - 0.4) * 0.1)
        vel = 1.0
        sys = sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                                    potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha), integrator='newton')
        sys.pos = 1.0
        sys.newpos = 1.0
        sys.vel = 1.0
        sys.newvel = 1.0
        for i, j in zip(
                [1.0, 0.8, 3.848716665874072e-05, 0.16000000000000003, 1.0, 1.1600000000000001, 1.0, 0.0],
                sys.propagate()):
            self.assertAlmostEqual(i, j)

    def testHMC(self):
        temp = 300.0
        mass = [1.0, 1.0]
        alpha = 1.0
        lam = 1.0
        fc = 1.0
        sys = system.system1D(temp=temp, fc=fc, lam=lam, alpha=alpha, mass=mass,
                              potential=pot.pertHarmonicOsc1D(fc=fc, alpha=alpha), integrator='hmc')


if __name__ == '__main__':
    unittest.main()
