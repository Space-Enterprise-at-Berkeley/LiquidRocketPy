from unittest.mock import patch
import os

import numpy as np
import matplotlib.pyplot as plt

from rocketpy.motors.LiquidMotor import MassBasedTank
import Fluid


# @PBales1
def test_mass_based_motor():
    lox = Fluid(name = "LOx", density = 1141, quality = 1.0) #Placeholder quality value
    propane = Fluid(name = "Propane", density = 493, quality = 1.0) #Placeholder quality value
    n2 = Fluid(name = "Nitrogen Gas", density = 51.75, quality = 1.0) #Placeholder quality value; density value may be estimate
    
    example_motor = MassBasedTank("Example Tank", 0.1540, 0.66, 0.7, "Placeholder", "Placeholder", lox, n2) 
    #Need docs to be pushed + tank dimension values
