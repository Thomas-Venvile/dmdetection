# -*- coding: utf-8 -*-

"""Routine for calculating DM collision rates with a generic direct detection experiment."""

__author__ = 'Grace Lawrence'
__email__ = 'gracie2084@gmail.com'
__version__ = '0.0.1'


import astropy as ap
from astropy import constants as const
from astropy import units as units
import scipy 
import numpy as np
from random import randrange, uniform, randint
import math as math
import matplotlib.pyplot as plt
import emcee 



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
        det= input((Detector(detvar['name'],detvar['atomic_mass'],detvar['mass_of_detector'])))
    elif detvar == 'SABRE': 
        det= Detector('SABRE', 149.89, 50)
        print('Detector: ',det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    elif detvar == "Xenon10":
        det= Detector('Xenon10', 131.29, 15)
        print("Detector: ",det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    elif detvar == "DAMA":
        det= Detector('DAMA', 149.89, 87.3)
        print("Detector: ",det.name,"; ",  "Atomic Mass: ", det.atomic_mass,"; ", "Detector Mass: ", det.mass_of_detector)
    else:
        detvar= ""
    while detvar != "SABRE" and detvar != "Xenon10" and detvar != "DAMA":
        detvar=raw_input("Please Pass a Valid Detector: ")
    
    if isinstance(dmvar, dict):
        dm_p = DM(dmvar['density'], dmvar['velocity'], dmvar['mass'], dmvar['cross_section'])
    elif dmvar == 'CDM':
        dm_p = DM(0.3, 2.3e7, 100, 1e-36)
        print("DM Density: ", dm_p.density,"; ", "DM Velocity: ", dm_p.velocity,"; ", "Particle Mass: ",dm_p.mass,"; ", "Cross Section",dm_p.cross_section)
    elif dmvar == 'Other':
        dm_p = DM(float(input("DM Density?(GeV/cm**2): ")),(input("DM Velocity?(cm/s): ")), (input("Particle Mass?(GeV): ")), (input("Cross Section?(cm**2): "))  )
        print("DM Density: ", dm_p.density,"; ", "DM Velocity: ", dm_p.velocity,"; ", "Particle Mass: ",dm_p.mass,"; ", "Cross Section",dm_p.cross_section)
    else:
        dmvar= " "
    while dmvar != "CDM" and dmvar != "Other":
        dmvar=raw_input("The particle was not recognised. Please specify a valid candidate/model: ")

    return det,dm_p

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
    plt.clf()

#create a maxwellian curve that randomly distributes a velocity between 0 and like, 232 (vmax).  
#Assuming the detector is on the solar sytem reference frame (stationary), then the recoil energy taken from Lewin and Smith ;
def recoil_energy(dm_p,det,N_particles):
    step_size = 232/int(N_particles)
    vel = np.arange(1,232, step_size) #This sample of velocities will be replaced with MCMC sample
    M_d = 107.35
    M_T =  det.atomic_mass.value
    #Define the energy
    E = 0.5*dm_p.mass.value*(vel)**2
    #Define the Kinematic Factor
    r = (4*M_d*M_T)/(M_d+M_T)**2
    #print("r", r)

    #Define the velocity distribution (Maxwellian)
    k_T = 0.5*(dm_p.mass.value)*(vel)**2
    pi = np.pi
    maxwell_dist = np.sqrt((dm_p.mass.value/(2*pi*k_T))**3)*4*pi*vel**2*np.exp(4)
    #plt.plot(maxwell_dist)
    #plt.title("Maxwell Distribution, single velocity")
    #plt.show()
    
    #Sample cos(theta) between zero and one, for x=100,10000,1e6 particles.
    theta_range = int(N_particles)
    i = np.arange(0,theta_range)
    angle_range = np.cos(i)
    #print angle_range
    
    angles = np.random.choice(angle_range, theta_range, replace=True)
    print("The particle angles are: ", angles)
    
    #Calculate Recoil Energy
    E_r = (E*r)*(1-angles)/(2)
    print("Recoil Energy", E_r)
    plt.plot(E_r)
    plt.title('Recoil_energies')
    plt.show()
    return E_r

def helm_form(dm_p, det, E_r): #Constants taken from Shield's thesis:
    M_d = dm_p.mass.value
    M_T =  det.atomic_mass.value 
    
    
    #Momentum Transfer
    q = np.sqrt(2*M_T*E_r)
    print("q values", q)
    #r_n is the effective nuclear radius
    c=((1.23)*(det.atomic_mass.value**(1/3))-0.6)*10**-13 #this conversts femtometres, 10^-15m, to cm. 
    a = 0.52*10**-13
    s = 0.9*10**-13#Nuclear skin thickness
    r_n = (c**2)+((7/3)*np.pi**2*a**2)-(5*s**2)
    print("R_N", r_n)

    #Calculate the Helm Holtz Form Factor
    Helm_form =3*np.exp((-q**2*s**2)/(2))*(np.degrees(np.sin(q*r_n))-(q*r_n*(np.degrees(np.cos(q*r_n))))/(q*r_n))
    print("Helm form factor; ", Helm_form)
    return Helm_form


#The helm factor is then used to calculate the effective collision cross section:
def Effective_cross_section(Helm_form):
    effective_cs = (10**(-38)) * (Helm_form)**2
    print ("effective cross section", effective_cs)

#Calculate the Quenching Factor, Q.
def quenching_fact(E_r):
    p_0= 2.126
    p_1=5.632
    Q = p_0 * np.sqrt(E_r) + p_1
    Electron_Equiv_Energy = Q * E_r
    print ("electron equivalent energies for the sample are: ", Electron_Equiv_Energy)
    plt.plot(Electron_Equiv_Energy) #Electron Equivalent energies are just the energy of the scintillation event. 
    plt.title("Electron Equivalent Energies") #then integrate through in certain energy levels. 
    plt.show()
    



N_particles = input("How many particles would you like to simulate? (Enter 231)")
detvar=raw_input("Which Detector Would You Like to Simulate: SABRE, Xenon10, DAMA?")
dmvar=raw_input("Which Particle Candidate/Model would you like to use: 'CDM' or 'Other'?")  
det,dm_p=run(detvar, dmvar)
#calcrate(det, dm_p)
#calcrate_1(dm_p)
E_r = recoil_energy(dm_p,det,N_particles)
Helm_form = helm_form(dm_p, det,E_r)
Effective_cross_section(Helm_form)
quenching_fact(E_r)

