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
    def __init__(self, name, diameter, height, gas, liquid=0, endcap="flat"):
        self.name = name
        self.diameter = diameter
        self.height = height
        self.gas = gas
        self.liquid = liquid

        self.capMap = {
            "flat": Disk(diameter / 2),
            "spherical": Hemisphere(diameter / 2),
        }
        self.cylinder = Cylinder(diameter / 2, height)
        self.cap = self.capMap.get(endcap)
        pass

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
        endcap,
        liquid,
        gas,
        ullage,
    ):
        super().__init__(name, diameter, height, endcap, gas, liquid)
    
        self.name = name
        self.diameter = diameter
        self.height = height
        self.endcap = endcap
        self.liquid = liquid
        self.gas = gas
        self.ullage = ullage

    def mass(self, t):
        # self.height is height of cylinder
        # assuming self.liquid/self.gas is density
        # self.cap should either be a Disk or Hemisphere class

        # template code. Still buggy
        def return_mass(time):
            volume_liquid = 0
            if self.cap is Disk:
                cross_sectional_area = self.cap.area()
                volume_liquid = (self.height-self.ullage) * cross_sectional_area
                volume_gas = self.cylinder.volume() - volume_liquid
                return volume_gas*self.gas + volume_liquid*self.liquid
        return return_mass

    def netMassFlowRate(self, t):
        # gas going in, mass leaving
        pass

    def liquidVolume(self, t):
        def return_liquid_volume(time):
            cross_sectional_area = self.cap.area()
            volume_liquid = (self.height-self.ullage) * cross_sectional_area
            return volume_liquid
        return return_liquid_volume



# @ompro07
# @PBales1
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
        super().__init__(
            name, 
            diameter, 
            height, 
            endcap, 
            gas, 
            liquid)

        #Define tank dimensions, fluid properties, and current fluid state
        self.name = name
        self.diameter = diameter
        self.height = height
        self.endcap = endcap
        self.liquid_mass = liquid_mass
        self.gas_mass = gas_mass
        self.liquid = liquid
        self.gas = gas

    def mass(self, t):
        current_mass = Function.getValue(self.liquid_mass, t) + Function.getValue(self.gas_mass, t)
        return Function(current_mass, t)

    def netMassFlowRate(self, t):
        tank_volume = 1/4*np.pi*self.height*(self.diameter ** 2) + 2*self.endcap

        liquid_initial_mass = self.liquid.denisty * tank_volume
        liquid_current_mass = Function.getValue(self.liquid_mass, t)
        liquid_mfr = (liquid_initial_mass - liquid_current_mass) / t

        gas_current_mass = Function.getValue(self.gas_mass, t)
        gas_mfr = (0 - gas_current_mass) / t #initial gas mass is 0

        return Function(liquid_mfr + gas_mfr, t)

    def liquidVolume(self, t):
        # density = m/V -> V_current = m_current/density
        current_liquid_mass = Function.getValue(self.liquid_mass, t)
        return Function(current_liquid_mass / self.liquid.density, t)
