"""
.. automodule: ensemble
    This module shall be used to implement subclasses of ensemble.
    It is a class, that is using multiple system. It can be used for RE or Conveyor belt
"""
from collections import Iterable
import numpy as np
import scipy.constants as const
from typing import List, Dict, Tuple
import copy
import itertools as it


class ReplicaExchange:
    parameters:dict
    dimensions:int
    replicas:dict={}
    exchange_information:list = []
    nsteps = 100


    #METROPOLIS CRITERION
    ##random part of Metropolis Criterion:
    randomnessIncreaseFactor = 1
    _defaultRandomness = lambda self, partner1, partner2: ((1 / self.randomnessIncreaseFactor) * np.random.rand() <= np.exp(-1.0 / (const.gas_constant / 1000.0 * partner2.temperature) * (partner1 - partner2.totPotEnergy)))
    ##default Metropolis Criterion
    _defaultMetropolisCriterion = lambda self, partner1, partner2: (partner1 < partner2 or self._defaultRandomness(partner1, partner2))
    exchange_criterium = _defaultMetropolisCriterion
    exchange_offset=0

    def __init__(self, system, parameter_Names, parameter_Ranges, exchange_criterium=None, steps_between_trials=None):
        self.dimensions = len(parameter_Names) #get dimensionality
        #TODO do some fancy parsing
        if(len(parameter_Names) > 1):
            self.parameters = {parameter_Name: parameter_Range for parameter_Name, parameter_Range in zip(parameter_Names, parameter_Ranges)}
        elif(len(parameter_Names)==1):
            self.parameters = {parameter_Names[0]: parameter_Ranges}


        self.system = system
        self.initialise()

        #exchange Criterium
        if(exchange_criterium != None):
            self.exchange_criterium = exchange_criterium

    def exchange(self):
        raise NotImplementedError("This method was not implemented for "+str(__name__)+" do so!")
        pass

    def run(self):
        for replica_coords, replica in self.replicas.items():
            replica.simulate(steps=self.nsteps)
        pass

    def initialise(self):
        # BUILD replicas
        self.initialise_replica_graph()

    def initialise_replica_graph(self):
        coord_dims = list(sorted(self.parameters))

        if(len(self.parameters) > 1):
            coord_it=it.product(*[list(self.parameters[r]) for r in sorted(self.parameters)])
        elif(len(self.parameters)==1):
            coord_it = list(map(lambda x: (x), self.parameters[coord_dims[0]]))
        else:
            raise Exception("Could not find parameters to exchange")

        for coords in coord_it:
            replica = copy.deepcopy(self.system)
            replica.trajectory = [] #is not deepcopied!!!
            for ind, parameter_Name in enumerate(coord_dims):
                if (hasattr(replica, parameter_Name)):
                    if (isinstance(coords, Iterable)):
                        setattr(replica, parameter_Name, coords[ind])
                    else:
                        setattr(replica, parameter_Name, coords)
                else:
                    raise Exception()
            self.replicas.update({coords: replica})
        self.nReplicas = len(self.replicas)

    def simulate(self, ntrials:int):
        for trial in range(ntrials-1):
            self.run()
            self.exchange()

    def get_trajectories(self)->Dict[Tuple, List]:
        return {coord:replica.trajectory for coord, replica in self.replicas.items()}

    def get_TotPot_Energy(self)->Dict[Tuple, float]:
        return {coord:replica._currentTotPot for coord, replica in self.replicas.items()}
    pass

class TemperatureReplicaExchange(ReplicaExchange):
    parameters:list = ["T"]
    dimensions:int = 1

    def exchange(self):

        totPots = self.get_TotPot_Energy()
        replica_keys = list(totPots.keys())
        print("TotPots",totPots)
        check_exchange = [self.exchange_criterium(totPots[partner1], totPots[partner2]) for partner1, partner2 in zip(replica_keys[self.exchange_offset::2], replica_keys[1+self.exchange_offset::2])]

        print("Exchange: ", check_exchange)
        for partner1, partner2 in zip(self.replicas[self.exchange_offset::2], self.replicas[1+self.exchange_offset::2]):
            pass

class HamiltonianReplicaExchange(ReplicaExchange):
    pass

class ReplicaExchangeEnvelopingDistributionSampling(ReplicaExchange):
    pass
