import unittest
import os ,sys
sys.path.append(os.path.dirname(__file__+"/../.."))

from ConveyorBelt.src import system, integrator as integrators, potential1D as pot
from ConveyorBelt.src import ensemble


class test_ReplicaExchangeCls(unittest.TestCase):

    RE = ensemble.ReplicaExchange
    integrator = integrators.monteCarloIntegrator()
    potential = pot.harmonicOsc1D()
    sys = system.system(potential=potential, integrator=integrator)
    group:ensemble.ReplicaExchange = None

    def setUp(self) -> None:
        self.group = None

    def tearDown(self) -> None:
        self.RE.replicas = {}

    def test_init_1DREnsemble(self):
        from ConveyorBelt.src import ensemble
        ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))

    def test_init_2DREnsemble(self):
        from ConveyorBelt.src import ensemble
        ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature", "mass"], parameter_Ranges=[range(288, 310), range(1,10)])

    def test_run_1DREnsemble(self):
        from ConveyorBelt.src import ensemble
        group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        group.run()

    def test_getTraj_1DREnsemble(self):
        replicas =22
        nsteps = 100
        group = None
        from ConveyorBelt.src import ensemble
        self.group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        self.group.nSteps_between_trials = nsteps
        self.group.run()
        trajectories = self.group.get_trajectories()


        print(len(trajectories))
        print([len(trajectories[t]) for t in trajectories])

        self.assertEqual(len(trajectories), 22, msg="not enough trajectories were retrieved!")
        self.assertEquals([len(trajectories[t]) for t in trajectories], second=[nsteps for x in range(replicas)], msg="traj lengths are not correct!")

    def test_getTotPot_1DREnsemble(self):
        replicas =22
        nsteps = 100
        from ConveyorBelt.src import ensemble
        self.group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=range(288, 310))
        self.group.nSteps_between_trials = nsteps
        self.group.run()
        totPots = self.group.get_Total_Energy()


        print(len(totPots))
        print(totPots)
        self.assertEqual(len(totPots), replicas, msg="not enough trajectories were retrieved!")

    def test_setPositionsList_1DREnsemble(self):
        param_range = range(288, 310)
        replicas =len(param_range)
        expected_pos= range(replicas)
        from ConveyorBelt.src import ensemble
        self.group = ensemble.ReplicaExchange(system=self.sys, parameter_Names=["temperature"], parameter_Ranges=param_range)

        initial_positions = sorted([self.group.replicas[replica]._currentPosition for replica in self.group.replicas])
        self.group.setReplicasPositions(expected_pos)
        setted_pos = sorted([self.group.replicas[replica]._currentPosition for replica in self.group.replicas])

        self.assertEqual(len(self.group.replicas), 22, msg="not enough trajectories were retrieved!")
        self.assertNotEqual(initial_positions, setted_pos, msg="Setted positions are the same as before!")
        self.assertEqual(setted_pos, list(expected_pos), msg="The positions were not set correctly!")

class test_TemperatureReplicaExchangeCls(unittest.TestCase):
    TRE = ensemble.TemperatureReplicaExchange
    group:ensemble.TemperatureReplicaExchange

    def tearDown(self) -> None:
        setattr(self, "group", None)


    def setUp(self) -> None:
        self.group = None
        self.TRE.replicas = {}


    def test_init(self):
        integrator = integrators.monteCarloIntegrator()
        potential = pot.harmonicOsc1D()
        sys = system.system(potential=potential, integrator=integrator)

        replicas =22
        nsteps = 100
        T_range=range(288, 310)
        setattr(self, "group", None)
        self.group = ensemble.TemperatureReplicaExchange(system=sys, temperature_Range=T_range)

    def test_run(self):
        integrator = integrators.monteCarloIntegrator()
        potential = pot.harmonicOsc1D()
        sys = system.system(potential=potential, integrator=integrator)

        replicas =22
        nsteps = 100
        T_range=range(288, 310)
        self.group = ensemble.TemperatureReplicaExchange(system=sys, temperature_Range=T_range)
        self.group.nSteps_between_trials = nsteps
        self.group.run()

    def test_exchange_all(self):
        integrator = integrators.monteCarloIntegrator()
        potential = pot.harmonicOsc1D()
        sys = system.system(potential=potential, integrator=integrator)

        T_range=range(1, 10)
        nReplicas = len(T_range)
        positions = list([float(1) for x in range(nReplicas)])
        velocities = list([float(0) for x in range(nReplicas)])


        self.group = ensemble.TemperatureReplicaExchange(system=sys, temperature_Range=T_range)
        self.group.setReplicasPositions(positions)
        self.group.setReplicasVelocities(velocities)
        self.group._defaultRandomness= lambda x,y: False

        self.group.exchange()
        all_exchanges = self.group._current_exchanges
        finpositions = list(self.group.getReplicasPositions().values())
        finvelocities = list(self.group.getReplicasVelocities().values())

        #Checking:
        ##constant params?
        self.assertEqual(len(self.group.replicas), nReplicas, msg="not enough trajectories were retrieved!")
        self.assertListEqual(finpositions, positions, msg="Positions should not change during exchange!")
        self.assertListEqual(finvelocities, velocities, msg="Velocities should not change during exchange!")
        ##exchange process
        self.assertEqual(nReplicas // 2, len(all_exchanges), msg="length of all exchanges is not correct!")
        self.assertTrue(all(list(all_exchanges.values())), msg="not all exchanges are True!!")

    def test_exchange_none(self):
        integrator = integrators.positionVerletIntegrator()
        potential = pot.harmonicOsc1D()
        sys = system.system(potential=potential, integrator=integrator)

        T_range = range(1, 500, 100)
        nReplicas = len(T_range)
        positions = [float(x)*100 for x in range(nReplicas)]
        velocities = list([float(1) for x in range(nReplicas)])

        self.group = ensemble.TemperatureReplicaExchange(system=sys, temperature_Range=T_range)
        #remove Randomness!
        self.group._defaultRandomness= lambda x,y: False
        self.group.setReplicasPositions(positions)
        self.group.setReplicasVelocities(velocities)

        #first round
        self.group.exchange()

        all_exchanges = self.group._current_exchanges
        finpositions = list(self.group.getReplicasPositions().values())
        finvelocities = list(self.group.getReplicasVelocities().values())

        #Checking:
        ##constant params?
        self.assertEqual(len(self.group.replicas), nReplicas, msg="not enough trajectories were retrieved!")
        self.assertListEqual(finpositions, positions, msg="Positions should not change during exchange!")
        self.assertListEqual(finvelocities, velocities, msg="Velocities should not change during exchange!")
        ##exchange process
        self.assertEqual(nReplicas//2, len(all_exchanges), msg="length of all exchanges is not correct!")
        self.assertFalse(all(list(all_exchanges.values())), msg="length of all exchanges is not correct!")


if __name__ == '__main__':
    unittest.main()