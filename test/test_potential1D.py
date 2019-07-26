import os,sys
import unittest
import numpy as np
import numbers
from collections.abc import Iterable
sys.path.append(os.path.dirname(__file__+"/../.."))

from ConveyorBelt.src import potential1D as pot


class potentialCls(unittest.TestCase):

    """
    TEST for Potential inputs
    """
    def test_check_positions__float_type(self):
        #check single Float
        position = 1.0
        potentialCls = pot.flat_well()
        checked_pos = potentialCls._check_positions_type(positions=position)

        if(not isinstance(checked_pos, Iterable)):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Float] - no list")
        elif(any([not isinstance(x, numbers.Number) for dim in checked_pos for x in dim])):
            print(type(checked_pos), type(checked_pos[0]))
            print(checked_pos)
            raise Exception("The return Type has to be Iterable[Float] - not all list elements are float")

    def test_check_positions_list_type(self):
        #check LIST[Float]
        position = [1.0, 2.0, 3.0]
        potentialCls = pot.flat_well()
        checked_pos = potentialCls._check_positions_type(positions=position)

        if(not isinstance(checked_pos, Iterable)):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Float] - no list")
        elif(any([not isinstance(x, numbers.Number)  for dim in checked_pos for x in dim])):
            print(type(checked_pos), type(checked_pos[0]))
            print(checked_pos)
            raise Exception("The return Type has to be Iterable[Float] - not all list elements are float")

    def test_check_positions_nDlist_type(self):
        position = [[1.0, 2.0, 3.0]]
        potentialCls = pot.flat_well()
        checked_pos = potentialCls._check_positions_type(positions=position)

        if(not isinstance(checked_pos, Iterable)):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Float] - no list")
        if(not any([isinstance(x, Iterable) for x in checked_pos])):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Iterable[Float]] - no list")
        elif(any([not isinstance(x, numbers.Number)  for dim in checked_pos for x in dim])):
            print(type(checked_pos), type(checked_pos[0]))
            print(checked_pos)
            raise Exception("The return Type has to be Iterable[Float] - not all list elements are float")

    def test_check_positions_npArray_type(self):
        #check nparray
        position = np.arange(1,10)
        potentialCls = pot.flat_well()
        checked_pos = potentialCls._check_positions_type(positions=position)

        if(not isinstance(checked_pos, Iterable)):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Float] - no list")
        elif(any([not isinstance(x, numbers.Number)for dim in checked_pos for x in dim])):
            print(type(checked_pos), type(checked_pos[0]))
            raise Exception("The return Type has to be Iterable[Float] - not all list elements are float")


    """
        Check all constructor defaults
    """
    def test_check_init_potentials(self):
        """
        Very coarse test for all potentials!
        :return:
        """
        import inspect

        is_class_member = lambda member: inspect.isclass(member) and member.__module__ == pot.__name__
        for member in inspect.getmembers(pot, is_class_member):
            if(not str(member[0]).startswith("_") and not any([x in str(member[0]) for x in ["envelopedPotential", "torsionPotential"]])):
                print(member)
                test = member[1]()

        #Test The Rest
        poteinf = pot.wavePotential()
        pot.envelopedPotential(V_is=[poteinf, poteinf])
        pot.envelopedPotentialMultiS(V_is=[poteinf, poteinf],s=[1.0])
        pot.torsionPotential(wave_potentials=[poteinf, poteinf])
        print("initiated all potential 1D classes")

"""
TEST for Potentials 1D
"""

class potentialCls_flatwell(unittest.TestCase):
    def test_energies(self):
        x_range = [0,1]
        y_max = 10
        y_min = 0
        positions = [0,2,1,0.5]
        expected_result = [0, 10, 0, 0]

        potential = pot.flat_well(x_range=x_range, y_max=y_max, y_min=y_min)

        energies = potential.ene(positions)

        self.assertEquals(first=expected_result, second=energies, msg="The results were not correct!")

    def test_dHdpos(self):
        x_range = [0,1]
        y_max = 10
        y_min = 0
        positions = [0,2,1,0.5]
        potential = pot.flat_well(x_range=x_range, y_max=y_max, y_min=y_min)
        expected_result = [0, 0, 0, 0]

        energies = potential.dhdpos(positions)

        self.assertEquals(first=expected_result, second=energies, msg="The results were not correct!")

class potentialCls_harmonicOsc1D(unittest.TestCase):

    def test_energies(self):
        fc= 1.0
        x_shift = 0.0
        y_shift = 0.0
        positions = [0,2,1,0.5]
        expected_result = [0, 10, 0, 0]

        potential = pot.harmonicOsc1D(fc=fc, x_shift=x_shift, y_shift=y_shift)

        energies = potential.ene(positions)

        self.assertEquals(first=expected_result, second=energies, msg="The results were not correct!")

    def test_dHdpos(self):
        fc: float = 1.0
        x_shift:float = 0.0
        y_shift:float = 0.0
        positions = [0,2,1,0.5]
        expected_result = [0, 10, 0, 0]

        potential = pot.harmonicOsc1D(fc=fc, x_shift=x_shift, y_shift=y_shift)

        energies = potential.dhdpos(positions)

        self.assertEquals(first=expected_result, second=energies, msg="The results were not correct!")

class potentialCls_wavePotential(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")

    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_torsionPotential(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_coulombPotential(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_lennardJonesPotential(unittest.TestCase):

    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_doubleWellPot1D(unittest.TestCase):

    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_perturbedLinCoupledHosc(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_perturbedExpCoupledHosc(unittest.TestCase):

    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_perturbedHarmonicOsc1D(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_envelopedPotential(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")

class potentialCls_envelopedPotentialMultiS(unittest.TestCase):
    def test_energies(self):
        raise NotImplementedError("Implement this test maaaaan!")
    def test_dHdpos(self):
        raise NotImplementedError("Implement this test maaaaan!")


if __name__ == '__main__':
    unittest.main()