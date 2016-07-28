# -*- coding: utf-8 -*-

"""Routine for calculating DM collision rates with a generic direct detection experiment."""

__author__ = 'Grace Lawrence'
__email__ = 'gracie2084@gmail.com'
__version__ = '0.0.1'

import astropy as ap
from astropy import constants as const
import scipy 
import numpy as np
import math as math
import matplotlib.pyplot as plt

detvar=input("Which Detector Would You Like to Simulate?")
dmvar=input("Which Particle Candidate/Model would you like to use: 'CDM' or 'Other'?")

class Detector:
        def __init__(self, name, atomic_mass, mass_of_detector):
            """Class for variables of Detector properties"""
            self.name=name
            self.atomic_mass=atomic_mass * ap.units.g / ap.units.mol
            self.mass_of_detector=mass_of_detector * ap.units.kg

class DM:
        def __init__(self, density, velocity, mass, cross_section):
            """Class for variables of DM model"""
            self.density=density * ap.units.GeV / ap.units.cm**3 
            self.velocity=velocity * ap.units.cm / ap.units.s
            self.mass=mass * ap.units.GeV #/ const.c**2
            self.cross_section=cross_section * ap.units.cm**2

def run(detvar, dmvar='CDM'):
    """Accepts input and pre-set values for DM and Detector classes"""
    if isinstance(detvar, dict):
        det = input((Detector(detvar['name'],detvar['atomic_mass'],detvar['mass_of_detector'])))
    elif detvar == 'SABRE': 
        det = Detector('SABRE', 149.89, 50)
        print("Detector: ",det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    elif detvar == "Xenon10":
    	det = Detector('Xenon10', 131.29, 15)
    	print("Detector: ",det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    elif detvar == "DAMA":
    	det = Detector('DAMA', 149.89, 87.3)
    	print("Detector: ",det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    else:
    	detvar= ""
    while detvar != "SABRE" and detvar != "Xenon10" and detvar != "DAMA":
    	detvar=input("Please Pass a Valid Detector: ")
    

    if isinstance(dmvar, dict):
         dm_p = DM(dmvar['density'], dmvar['velocity'], dmvar['mass'], dmvar['cross_section'])
    elif dmvar == 'CDM':
        dm_p = DM(0.3, 2.3e7, 100, 1e-36)
        print("DM Density: ", dm_p.density,"; ", "DM Velocity: ", dm_p.velocity,"; ", "Particle Mass: ",dm_p.mass,"; ", "Cross Section",dm_p.cross_section)
    elif dmvar == 'Other':
        dm_p = DM(float(input("DM Density?(GeV/cm**2): ")),(input("DM Velocity?(cm/s): ")), (input("Particle Mass?(GeV): ")), (input("Cross Section?(cm**2): "))  )
        print("DM Density: ", dm_p.density,"; ", "DM Velocity: ", dm_p.velocity,"; ", "Particle Mass: ",dm_p.mass,"; ", "Cross Section",dm_p.cross_section)
    else:
    	dmvar=""
    while dmvar != "CDM" and dmvar != "Other":
    	dmvar=input("Please specify a valid candidate/model: ")
    
    return (det, dm_p)
    #return crash statement.         

def calcrate(det, dm_p):
    """Calculates expected Reaction Rates and Number Densities"""
    NaInumd = ((ap.constants.N_A/det.atomic_mass))
    DMnumd = ((dm_p.density/dm_p.mass))
    print ("Number density of NaI atoms per kg: ",NaInumd), (ap.constants.N_A / det.atomic_mass).cgs
    print ("Number density of DM particles: ", DMnumd) #(dm_p.density / dm_p.mass).cgs
    rate = (dm_p.cross_section*dm_p.velocity*((dm_p.density/dm_p.mass)*(ap.constants.N_A/det.atomic_mass)))*((1000*86400))
    print("Expected Count Rate: ", rate.value, "counts/kg/day")
    print ("Rate per day:", rate.value*det.mass_of_detector.value, "counts/day")
    print("Rate per year: ", (rate.value*det.mass_of_detector.value)*365, "counts/year")
    
    return dm_p, det
    
def calcrate_1(dm_p):
    """Accepts time period, calculates time-dependent velocity, graphs expected modulation"""
    time=int(input("Please Enter the Observing Period (in days): "))
    t=np.arange(0, time)
    v_esc=232+15*np.cos((2*np.pi)*(((t)-152.5)/365.25))
    rate_1 = dm_p.cross_section*v_esc*(dm_p.density/dm_p.mass)*(ap.constants.N_A/det.atomic_mass)*det.mass_of_detector
    plt.plot(t, rate_1*det.mass_of_detector)
    plt.ylabel('Counts (per day)')
    plt.xlabel('Time (days)')
    plt.title('Predicted Modulation')
    plt.grid(True)
    plt.show()

    
det,dm_p =run(detvar, dmvar)
calcrate(det, dm_p)
calcrate_1(dm_p)


