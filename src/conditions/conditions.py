"""
Module: Conditions
    This module shall be used to implement subclasses of conditions like, thermostat or distance restraints
"""

#from ConveyorBelt.src.system import system as sys

class condition:
    _tau:float  #tau = apply every tau steps

    def __init__(self , sys):   #system):
        raise NotImplementedError("This " + __class__ + " class is not implemented")

    def apply(self):#, system:system):
        raise NotImplementedError("The function step in " + __class__ + " class is not implemented")

    def coupleSystem(self, system): #sys):
        self.system = system