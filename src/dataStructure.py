"""
Module: dataStructure
    This module shall be used to implement all needed data Structures for the project.
"""
from collections import namedtuple


#states:
basicState = namedtuple("State", ["position", "temperature",
                                  "totEnergy", "totPotEnergy", "totKinEnergy",
                                  "dhdpos", "velocity"])
lambdaState = namedtuple("LState", ["position", "temperature",
                                   "totEnergy", "totPotEnergy", "totKinEnergy",
                                   "dhdpos",  "velocity",
                                   "lamb", "dhdlam"])
envelopedPStstate = namedtuple("EState", ["position", "temperature",
                                         "totEnergy", "totPotEnergy", "totKinEnergy",
                                         "dhdpos",  "velocity",
                                         "s", "Eoff"])
