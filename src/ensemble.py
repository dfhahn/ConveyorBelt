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
from ConveyorBelt.src import system

class ReplicaExchange:
    ##Parameters
    parameters:Dict
    parameter_names:List
    coordinate_dimensions:int

    ##Replicas
    replicas:dict={}
    nReplicas:int
    replica_graph_dimensions:int

    ##Exchange params/io
    exchange_information:list = []

    ##simulation Params
    nSteps_between_trials = 100


    #METROPOLIS CRITERION
    ##random part of Metropolis Criterion:
    randomnessIncreaseFactor = 0.01
    _temperature_exchange:float= 298
    _defaultRandomness = lambda self, originalParams, swappedParams: ((1 / self.randomnessIncreaseFactor) * np.random.rand() <= np.exp(-1.0 / (const.gas_constant / 1000.0 * self._temperature_exchange) * (originalParams - swappedParams)))

    ##default Metropolis Criterion
    _defaultMetropolisCriterion = lambda self, originalParams, swappedParams: (np.greater_equal(originalParams, swappedParams) or self._defaultRandomness(originalParams, swappedParams))
    exchange_criterium = _defaultMetropolisCriterion
    exchange_offset=0

    def __init__(self, system:system.system, parameter_Names, parameter_Ranges, exchange_criterium=None, steps_between_trials=None):

        #TODO do some fancy parsing
        #SET PARAMETER FIELDS
        if(len(parameter_Names) > 1):
            self.parameters = {parameter_Name: parameter_Range for parameter_Name, parameter_Range in zip(parameter_Names, parameter_Ranges)}
        elif(len(parameter_Names)==1):
            self.parameters = {parameter_Names[0]: parameter_Ranges}

        self.coordinate_dimensions = len(parameter_Names)  # get dimensionality
        self.parameter_names = list(self.parameters.keys())

        #SET SYSTEM
        self.system = system

        if(steps_between_trials != None):
            self.nSteps_between_trials = steps_between_trials

        #initialize the replica graphs
        self.initialise()

        #exchange Criterium
        if(exchange_criterium != None):
            self.exchange_criterium = exchange_criterium

        #steps_between_trials:

    def exchange(self):
        raise NotImplementedError("This method was not implemented for "+str(__name__)+" do so!")
        pass

    def run(self):
        for replica_coords, replica in self.replicas.items():
            replica.simulate(steps=self.nSteps_between_trials)
        pass

    def initialise(self):
        # BUILD replicas
        self.initialise_replica_graph()

    def initialise_replica_graph(self, verbose:bool=False):
        coord_dims = list(sorted(self.parameters))

        #generate all parameter combinations
        if(len(self.parameters) > 1):
            coord_it=it.product(*[list(self.parameters[r]) for r in sorted(self.parameters)])
        elif(len(self.parameters)==1):
            coord_it = list(map(lambda x: (x), self.parameters[coord_dims[0]]))
        else:
            raise Exception("Could not find parameters to exchange")

        #set all parameters
        self.nReplicas = len(list(coord_it))
        if(verbose):
            print("Coord_prod", list(coord_it))
            print("Coord Dim", self.nReplicas)
        replicas = [copy.deepcopy(self.system) for x in range(self.nReplicas)]

        #final copying and field updating (coords later)
        for replica in replicas:
            replica.trajectory = [] #fields are not deepcopied!!!
            # set steps between trials
            replica.nsteps = self.nSteps_between_trials

        for coords, replica in zip(coord_it, replicas):
            for ind, parameter_Name in enumerate(coord_dims):
                #set parameter set
                if (hasattr(replica, parameter_Name)):
                    if (isinstance(coords, Iterable)):
                        setattr(replica, parameter_Name, coords[ind])
                    else:
                        setattr(replica, parameter_Name, coords)
                else:
                    raise Exception("REPLICA INIT FAILDE: Replica does not have a field: "+parameter_Name+"\n")
            replica.initVel()
            self.replicas.update({coords: replica})
        self.nReplicas = len(self.replicas)

    def getReplicasPositions(self)->Dict:
        """
        .. autofunction:: getReplicaPositions
            If a list is passed to this function, the replicas will be aligned by sorted and then sequentially filled up with the new positions.
            Else a Dictionary can be passed and the directed position is transferred to the replica at the replica coordinate(key)
        :param positions:
        :return:
        """
        vals_dict = {}
        for replicaName, replica in self.replicas.items():
            vals_dict.update({replicaName: getattr(replica, "_currentPosition")})
        return vals_dict

    def getReplicasVelocities(self) -> Dict:
        """
        .. autofunction:: getReplicaPositions
            If a list is passed to this function, the replicas will be aligned by sorted and then sequentially filled up with the new positions.
            Else a Dictionary can be passed and the directed position is transferred to the replica at the replica coordinate(key)
        :param positions:
        :return:
        """
        vals_dict = {}
        for replicaName, replica in self.replicas.items():
            vals_dict.update({replicaName: getattr(replica, "_currentVelocities")})
        return vals_dict

    def setReplicasPositions(self, positions:(List or Dict)):
        """
        .. autofunction:: setReplicaPositions
            If a list is passed to this function, the replicas will be aligned by sorted and then sequentially filled up with the new positions.
            Else a Dictionary can be passed and the directed position is transferred to the replica at the replica coordinate(key)
        :param positions:
        :return:
        """
        if(len(positions)==self.nReplicas):
            if (type(positions) == dict):
                for replicaName, position in positions.items():
                    setattr(self.replicas[replicaName], "_currentPosition", position)
                    self.replicas[replicaName].totPot()
            elif(isinstance(positions, Iterable)):
                for replicaName, position in zip(sorted(self.replicas), positions):
                    setattr(self.replicas[replicaName], "_currentPosition", position)
                    self.replicas[replicaName].totPot()
            else:
                raise Exception("Did not understand the the type of the new positions "+str(type(positions)))
        else:
            raise ValueError("Not enough positions got passed to setReplicapositions\n replicas: "+str(self.nReplicas)+"\n positions: "+str(len(positions))+"\n"+str(positions))

    def setReplicasVelocities(self, velocities:(List or Dict)):
        """
        .. autofunction:: setReplicasVelocities
            If a list is passed to this function, the replicas will be aligned by sorted and then sequentially filled up with the new velocity.
            Else a Dictionary can be passed and the directed position is transferred to the replica at the replica coordinate(key)
        :param positions:
        :return:
        """
        if(len(velocities)==self.nReplicas):
            if (type(velocities) == dict):
                for replicaName, velocity in velocities.items():
                    setattr(self.replicas[replicaName], "_currentVelocities", velocity)
                    self.replicas[replicaName].totKin()
            elif(isinstance(velocities, Iterable)):
                for replicaName, velocity in zip(sorted(self.replicas), velocities):
                    setattr(self.replicas[replicaName], "_currentVelocities", velocity)
                    self.replicas[replicaName].totKin()
            else:
                raise Exception("Did not understand the the type of the new positions "+str(type(velocities)))
        else:
            raise ValueError("Not enough positions got passed to setReplicapositions\n replicas: "+str(self.nReplicas)+"\n positions: "+str(len(velocities))+"\n"+str(velocities))

    def setParameterSet(self, coordinates:List, replicas:List):
        """
            ..autofunction:: set ParameterSet
            This function is setting new coordinates to the replicas in the replica lists.
            The coordinates will be assigned sequentially in the same order to the replicas List.

        :warning: This function is Overwritting old coordinates!
        :param coordinates:
        :return:
        """

        if(self.coordinate_dimensions>1):
            self.replicas = {}
            for coords, replica in zip(coordinates, replicas):
                for ind, parameter_Name in enumerate(self.parameters):
                    #set parameter set
                    if (hasattr(replica, parameter_Name)):
                        setattr(replica, parameter_Name, coords[ind])
                    else:
                        raise Exception("REPLICA INIT FAILDE: Replica does not have a field: "+parameter_Name+"\n")
                self.replicas.update({coords: replica})
        else:
            self.replicas = {}
            for coords, replica in zip(coordinates, replicas):
                #set parameter set
                if (hasattr(replica, self.parameter_names[0])):
                    if (isinstance(coords, Iterable)):
                        setattr(replica, self.parameter_names[0], coords[0])
                    else:
                        setattr(replica, self.parameter_names[0], coords)

                else:
                    raise Exception("REPLICA INIT FAILDE: Replica does not have a field: "+self.parameters[0]+"\n")
                self.replicas.update({coords: replica})

    def simulate(self, ntrials:int):
        for trial in range(ntrials-1):
            self.run()
            self.exchange()

    def get_trajectories(self)->Dict[Tuple, List]:
        return {coord:replica.trajectory for coord, replica in self.replicas.items()}

    def get_Total_Energy(self)->Dict[Tuple, float]:
        return {coord:replica.getTotEnergy() for coord, replica in self.replicas.items()}
    pass

class TemperatureReplicaExchange(ReplicaExchange):

    parameter_names:list = ["temperature"]
    coordinate_dimensions:int = 1
    replica_graph_dimensions:int = 1

    def __init__(self, system, temperature_Range:Iterable=range(298,320), exchange_criterium=None, steps_between_trials=None, exchange_trajs:bool=True):
        super().__init__(system=system, parameter_Names=self.parameter_names, parameter_Ranges=temperature_Range, exchange_criterium=exchange_criterium, steps_between_trials=steps_between_trials)

        if(exchange_trajs):
            self.exchange_param = "trajectory"
        else:
            self.exchange_param = self.parameter_names[0]

    def exchange(self, verbose:bool=True):
        """
        .. autofunction:: Exchange the Trajectory of T-replicas in pairwise fashion
        :param verbose:
        :return:
        """
        original_totPots = self.get_Total_Energy()
        original_T = list(sorted(original_totPots.keys()))
        replica_values = list([self.replicas[key] for key in original_T])

        #SWAP temperature params pairwise
        ##take care of offset situations and border replicas
        swapped_T = [] if self.exchange_offset==0 else [original_T[0]]
        ##generate sequence with swapped params
        for partner1, partner2 in zip(original_T[self.exchange_offset::2], original_T[1+self.exchange_offset::2]):
            swapped_T.extend([partner2, partner1])
        ##last replica on the border?
        if(self.exchange_offset == 0):
            swapped_T.append(original_T[-1])

        if(verbose):
            print("original T ", original_T)
            print("swapped T ", swapped_T)

        #SWAP Params to calculate energies in swapped case
        self.setParameterSet(coordinates=swapped_T, replicas=replica_values)    #swap parameters

        #scaleVel
        if( not any([getattr(self.replicas[replica], "_currentVelocities") == None for replica in self.replicas])):
            [setattr(self.replicas[replica], "_currentVelocities", np.multiply(getattr(self.replicas[replica], "_currentVelocities"), np.divide(swapped_T[i], original_T[i]))) for i, replica in enumerate(self.replicas)]

        swapped_totPots = self.get_Total_Energy()  #calc swapped parameter Energies
        #scale Vel
        if (not any([getattr(self.replicas[replica], "_currentVelocities") == None for replica in self.replicas])):
            [setattr(self.replicas[replica], "_currentVelocities", np.multiply(self.replicas[replica]._currentVelocities, np.divide(original_T[i], swapped_T[i]))) for i, replica in enumerate(self.replicas)]

        self.setParameterSet(coordinates=original_T, replicas=replica_values)        #swap back parameters

        if(verbose):
            print("origTotE ", [original_totPots[key] for key in sorted(original_T)])
            print("SWPTotE ", [swapped_totPots[key] for key in sorted(swapped_T)])

        #compare along 1D axis
        exchanges_to_make={}
        for partner1, partner2 in zip(original_T[self.exchange_offset::2], original_T[1+self.exchange_offset::2]):
            originalEnergies = np.add(original_totPots.get(partner1), original_totPots.get(partner2))
            swapEnergies = np.add(swapped_totPots.get(partner1), swapped_totPots.get(partner2))
            print("partners: "+str(partner1)+"/"+str(partner2)+" \t originalEnergies "+str(originalEnergies)+" / Swap "+str(swapEnergies))
            print("randomness part: "+str(self._defaultRandomness(originalEnergies, swapEnergies)))

            exchanges_to_make.update({(partner1, partner2): self.exchange_criterium(originalEnergies, swapEnergies)})

        #Acutal Exchange of params (actually trajs here
        if(verbose):
            print("Exchange: ", exchanges_to_make)
            print("exchaning param: ", self.exchange_param)
        for (partner1ID, partner2ID), exchange in exchanges_to_make.items():
            if(exchange):
                if(verbose): print("Exchanging: "+str(partner1ID)+"\t"+str(partner2ID)+"\t"+str(exchange)+"\n")
                partner1 = self.replicas[partner1ID]
                partner2 = self.replicas[partner2ID]

                #T=self.parameter_names[0]
                exchange_param = "trajectory"
                param = getattr(partner1, self.exchange_param)
                setattr(partner1, self.exchange_param,  getattr(partner2, self.exchange_param))
                setattr(partner2, self.exchange_param, param)
            else:
                if (verbose): print("not Exchanging: "+str(partner1ID)+" / "+str(partner2ID)+" \n")
        self._current_exchanges = exchanges_to_make

        #update the offset
        self.exchange_offset = (self.exchange_offset+1)%2

class HamiltonianReplicaExchange(ReplicaExchange):
    pass

class ReplicaExchangeEnvelopingDistributionSampling(ReplicaExchange):
    pass
