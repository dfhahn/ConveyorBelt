import numpy as np
import typing as t
from collections import Iterable
from ConveyorBelt.src.conditions.conditions import Condition
from ConveyorBelt.src.system import system as sys

class periodicBoundaryCondition(Condition):
    """
        ..autoclass:: periodicBoundaryCondition
            This class allows to enable sampling in mirror images and projects the coordinates to the restricted space.
    """

    lowerbounds:Iterable
    higherbounds:Iterable
    _tau:int = 1    #each _tau steps apply

    def __init__(self, boundary:Iterable, system:sys=None):
        self._parse_boundary(boundary)
        if(system != None):
            self.system = system

    def apply(self):
        new_current_position = []
        for dim_pos, dimlBound, dimhBound in zip(self.system._currentPosition[0], self.lowerbounds, self.higherbounds):
            if(dim_pos < dimlBound):
                new_current_position.append(dimhBound - (dimlBound - dim_pos))
                #new_current_position.append(np.subtract(dimhBound, np.subtract(dimlBound, dim_pos%dimhBound)))
            elif(dim_pos > dimhBound):
                new_current_position.append(dimlBound + (dim_pos- dimhBound))
                #new_current_position.append(np.add(dimlBound, np.subtract(dim_pos%dimlBound, dimhBound)))
            else:
                new_current_position.append(dim_pos)
        self.system._currentPosition = [new_current_position]

    def _parse_boundary(self, boundary):  #Todo: really needed?
        if(isinstance(boundary, Iterable)):
            if (isinstance(boundary[0], Iterable)):
                self.higherbounds = np.array(list(map(max, boundary )))
                self.lowerbounds = np.array(list(map(min, boundary )))

            else:
                self.higherbounds = max(boundary)
                self.lowerbounds = min(boundary)
        return True