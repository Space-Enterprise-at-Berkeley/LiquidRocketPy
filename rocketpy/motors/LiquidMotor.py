# -*- coding: utf-8 -*-

__author__ = "Giovani Hidalgo Ceotto, Oscar Mauricio Prada Ramirez, João Lemes Gribel Soares, Mateus Stano and Pedro Henrique Marinho Bressan"
__copyright__ = "Copyright 20XX, RocketPy Team"
__license__ = "MIT"

from abc import ABC, abstractmethod

import numpy as np
from scipy import integrate

from rocketpy.Function import Function
from rocketpy.motors import Motor
from rocketpy.supplement import Disk, Cylinder, Hemisphere
from rocketpy.motors import Fluid

# @Stano
class LiquidMotor(Motor):
    def __init__(
        self,
        thrustSource,
        burnOut,
        nozzleRadius,
        throatRadius,
        reshapeThrustCurve=False,
        interpolationMethod="linear",
    ):

        super.__init__()
        self.tanks = []
        pass

    def evaluateMassFlowRate(self):
        pass

    def evaluateCenterOfMass(self):
        pass

    def evaluateInertiaTensor(self):
        pass

    def addTank(self, tank, position):
        self.tanks.append({"tank": tank, "position": position})


class Tank(ABC):
    def __init__(self, name: str, diameter, height, gas: Fluid, liquid: Fluid = 0, endcap: str = "flat"):
        self.name = name
        self.diameter = diameter
        self.height = height
        self.gas = gas
        self.liquid = liquid

        self.capMap = {
            "flat": Disk(diameter / 2),
            "spherical": Hemisphere(diameter / 2),
        }
        self.cap = self.capMap.get(endcap)

        self.tankMap = {
            "flat": Cylinder(diameter/2, height),
            "spherical": SphericalCylinder(diameter/2, height)
        }
        self.tank = self.tankMap.get(endcap)

    @abstractmethod
    def mass(self, t):
        """Returns the total mass of liquid and gases inside the tank as a
        function of time.

        Parameters
        ----------
        time : float
            Time in seconds.

        Returns
        -------
        Function
            Mass of the tank as a function of time. Units in kg.
        """
        pass

    @abstractmethod
    def netMassFlowRate(self, t):
        """Returns the net mass flow rate of the tank as a function of time.
        Net mass flow rate is the mass flow rate exiting the tank minus the
        mass flow rate entering the tank, including liquids and gases.

        Parameters
        ----------
        time : float
            Time in seconds.

        Returns
        -------
        Function
            Net mass flow rate of the tank as a function of time.
        """
        pass

    @abstractmethod
    def liquidVolume(self, t):
        """Returns the volume of liquid inside the tank as a function
        of time.

        Parameters
        ----------
        time : float
            Time in seconds.

        Returns
        -------
        Function
            Tank's liquid volume as a function of time.
        """
        pass

    def centerOfMass(self, t):
        """Returns the center of mass of the tank's fluids as a function of
        time.

        Parameters
        ----------
        time : float
            Time in seconds.

        Returns
        -------
        Function
            Center of mass of the tank's fluids as a function of time.
        """
        liquid_volume = self.liquidVolume(t)
        if liquid_volume < self.cap.volume:
            self.cap.filled_volume = liquid_volume
            return self.cap.filled_centroid
        else:
            self.cylinder.filled_volume = liquid_volume - self.cap.volume

            cylinder_mass = self.cylinder.filled_volume * self.liquid.density
            cap_mass = self.cap.volume * self.liquid.density

            return (
                self.cap.centroid * cap_mass + self.cylinder.centroid * cylinder_mass
            ) / (cap_mass + cylinder_mass)

    def inertiaTensor(self, t):
        """Returns the inertia tensor of the tank's fluids as a function of
        time.

        Parameters
        ----------
        time : float
            Time in seconds.

        Returns
        -------
        Function
            Inertia tensor of the tank's fluids as a function of time.
        """
        ...

# @MrGribel
class MassFlowRateBasedTank(Tank):
    def __init__(
        self,
        name,
        diameter,
        height,
        endcap,
        initial_liquid_mass,
        initial_gas_mass,
        liquid_mass_flow_rate_in,
        gas_mass_flow_rate_in,
        liquid_mass_flow_rate_out,
        gas_mass_flow_rate_out,
        liquid,
        gas,
    ):
        super().__init__(name, diameter, height, endcap, gas, liquid)


# @phmbressan
class UllageBasedTank(Tank):
    def __init__(
        self,
        name,
        diameter,
        height,
        endcap: str,
        liquid: Fluid,
        gas: Fluid,
        ullage,
    ):
        super().__init__(name, diameter, height, endcap, gas, liquid)

        # ullage v. time
        self.ullage = Function(ullage)
        
        # We are not accounting for pressurized fluids

    def mass(self, t):
        liquid_volume = self.liquidVolume(t)

        liquid_mass = liquid_volume * self.liquid.density
        gas_mass = (self.tank.volume() - liquid_volume) * self.gas.density

        return gas_mass + liquid_mass

    def netMassFlowRate(self, t):
        current_height = self.full_height - self.ullage.getValue(t)
        delta_volume = Function(self.tank.filled_volume).differentiate(current_height)

        return delta_volume * self.liquid.density + delta_volume * self.gas.density

    def liquidVolume(self, t):
        return self.tank.filled_volume(self.height-self.ullage.getValue(t))


# @ompro07
class MassBasedTank(Tank):
    def __init__(
        self,
        name,
        diameter,
        height,
        endcap,
        liquid_mass,
        gas_mass,
        liquid,
        gas,
    ):
        super().__init__(name, diameter, height, endcap, gas, liquid)
        pass
