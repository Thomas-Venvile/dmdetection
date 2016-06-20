# -*- coding: utf-8 -*-

"""Routine for calculating DM collision rates with a generic direct detection experiment."""

__author__ = 'Alan Duffy'
__email__ = 'mail@alanrduffy.com'
__version__ = '0.1.0'

import astropy as ap
import numpy as np

class Detector:
        def __init__(self, name, atomic_mass, mass_of_detector):
            """Class for variables of Detector properties"""
            self.name=name
            self.atomic_mass=atomic_mass * ap.units.g / ap.units.mol
            self.mass_of_detector=mass_of_detector * ap.units.kg

class DM:
        def __init__(self, density, velocity, mass, cross_section):
            """Class for variables of DM particle"""
            self.density=density * ap.units.GeV / ap.units.cm**3 #/ ap.constants.c**2 #ap.units.msolMass / ap.units.mpc**3
            self.velocity=velocity * ap.units.km / ap.units.s
            self.mass=mass * ap.units.GeV #/ ap.constants.c**2
            self.cross_section=cross_section * ap.units.pbarn

def calcrate(det,dm_p):
    """ Calculate rates per second """
    rate = dm_p.cross_section*dm_p.velocity*(dm_p.density/dm_p.mass)*(ap.constants.N_A/det.atomic_mass)*det.mass_of_detector

    return rate

def run(detvar, dmvar='CDM'):
    """ Calculate rates and plot """

    if isinstance(detvar, dict):
        ## user is providing own parameters, i.e. det = ['name' : 'SABRE', 'atomic_mass' : 150., 'mass_of_detector' : 50.]
        det = Detector(detvar['name'],detvar['atomic_mass'],detvar['mass_of_detector'])
    elif detvar == 'SABRE':
        det = Detector('SABRE', 149.89, 50.)
    else:
        print "What detector do you want? You said ", detvar
        ## Add a return crash statement

    if isinstance(dmvar, dict):
        ## user is providing own parameters, i.e. particle = ['density' : 0.3 'GeV/cm^3', 'velocity' : 230., 'mass' : 100., 'cross_section' : 1e-5]
        dm_p = DM(dmvar['density'], dmvar['velocity'], dmvar['mass'], dmvar['cross_section'])
    elif dmvar == 'CDM':
        dm_p = DM(0.3, 230., 10., 1e-5)
    else:
        print "What particle is this? You passed ", dmvar

    ## Check numvers
    print "Number density of NaI atoms per kg ", (ap.constants.N_A / det.atomic_mass)
    print "Number density of DM particles ", (dm_p.density / dm_p.mass).cgs

    rate = calcrate(det, dm_p)
    print "Rate per s ",rate.cgs
    print "Rate per year ",(rate * ap.units.yr).cgs

    return rate


