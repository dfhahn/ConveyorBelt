import unittest
import numpy as np
from ConveyorBelt import potential1D as pot

class test_Integrators(unittest.TestCase):
    def test_Monte_Carlo_Integrator(self):
        expected_norm_dist = np

        s = 1
        Eoffs = [0, 0]
        V_is = [pot.harmonicOsc1D(x_shift=-10), pot.harmonicOsc1D(x_shift=10)]
        eds_pot = pot.envelopedPotential(V_is=V_is, s=s, Eoff_i=Eoffs)

        positions = list(range(-100,100))
        energies = eds_pot.ene(positions)

        self.assert_(any([ex==ene for ex, ene in zip(expected, energies)]))

if __name__ == '__main__':
    unittest.main()