"""
Module: Conditions
    This module shall be used to implement subclasses of conditions like, thermostat or distance restraints
"""

from src.system import system1D

class condition:
    def __init__(self):
        raise NotImplementedError("This " + __class__ + " class is not implemented")

    def apply(self, system:system1D):
        raise NotImplementedError("The function step in " + __class__ + " class is not implemented")

