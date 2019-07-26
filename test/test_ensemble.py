import unittest
import os,sys
sys.path.append(os.path.dirname(__file__+"/../.."))


from ConveyorBelt.src import ensemble
from ConveyorBelt.src import system, integrator, potential1D as pot


class test_ReplicaExchangeCls(unittest.TestCase):
    integrator = integrator.monteCarloIntegrator()
    potential = pot.harmonicOsc1D()
    sys = system.system(potential=potential, integrator=integrator)

    def test_init_1DREnsemble(self):
        ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))

    def test_init_2DREnsemble(self):
        ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature", "mass"], parameter_Ranges=[range(288, 310), range(1,10)])

    def test_run_1DREnsemble(self):
        group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        group.run()

    def test_getTraj_1DREnsemble(self):
        replicas =22
        nsteps = 100
        group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        group.nsteps = nsteps
        group.run()
        trajectories = group.get_trajectories()


        print(len(trajectories))
        print([len(trajectories[t]) for t in trajectories])

        self.assertEqual(len(trajectories), 22, msg="not enough trajectories were retrieved!")
        self.assertEquals([len(trajectories[t]) for t in trajectories], second=[nsteps for x in range(replicas)], msg="traj lengths are not correct!")

    def test_getTotPot_1DREnsemble(self):
        replicas =22
        nsteps = 100
        group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        group.nsteps = nsteps
        group.run()
        totPots = group.get_TotPot_Energy()


        print(len(totPots))
        print(totPots)
        self.assertEqual(len(totPots), replicas, msg="not enough trajectories were retrieved!")

class test_TemperatureReplicaExchangeCls(unittest.TestCase):
    integrator = integrator.monteCarloIntegrator()
    potential = pot.harmonicOsc1D()
    sys = system.system(potential=potential, integrator=integrator)

    def test_exchange(self):
        replicas =22
        nsteps = 100
        group = ensemble.TemperatureReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        group.nsteps = nsteps
        group.run()
        totPots = group.get_TotPot_Energy()
        group.exchange()



        pass
if __name__ == '__main__':
    unittest.main()