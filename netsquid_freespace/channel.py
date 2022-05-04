#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 20:01:15 2021

@author: schiavon

This script implements the actual moving satellite channel to be used with the
FixedSatelliteLossModel NetSQUID class. 
"""

import orekit

from orekit.pyhelpers import download_orekit_data_curdir, setup_orekit_curdir, datetime_to_absolutedate

import os

import lowtran
import numpy as np

from math import radians

# setup orekit
vm = orekit.initVM()

if not os.path.isfile('orekit-data.zip'):
    download_orekit_data_curdir()

setup_orekit_curdir()

# Orekit imports
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.orekit.utils import IERSConventions, Constants


# Prepare global variables - time and space reference systems
ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                         Constants.WGS84_EARTH_FLATTENING, ITRF)
inertialFrame = FramesFactory.getEME2000()


# utc = TimeScalesFactory.getUTC()


# Define classes for the different objects
class orbitModelError(Exception):
    pass

class Satellite:
    """ Satellite orbit
    
    This class implements a real satellite orbit.
    
    The class defines the satellite using different simulation orbit techniques:
        * 'tle': use the two-line elements (TLE) of an existing satellite
        * 'polOrbPass': use the model of polar orbit passage described in
            [Moll et al., PRA 99, 053830 (2019)]
            
    Paramters
    ---------
    incAngle: float
        Inclination angle of the satellite w.r.t. the ground station [deg].
    satAlt: float
        Altitude of the satellite [km].    
    """
    
    simTypeAllowed = ('tle','polOrbPass','')
    
    def __init__(self,tle='',simType='tle',incAngle=0,satAlt=0):
        self.simType = simType
        
        if simType == 'tle':
            """ Initialization of the satellite orbit from a TLE list """
            self.tleList = tle
            self.tleObject = TLE(*self.tleList)
            self.propagator = TLEPropagator.selectExtrapolator(self.tleObject)
        elif simType == 'polOrbPass':
            self.incAngle = radians(incAngle)
            self.satAlt = satAlt*1e3
        else:
            self.simType = ''
                   
    def setSimTLE(self,tle):
        self.simType = 'tle'
        self.tleList = tle
        self.tleObject = TLE(*self.tleList)
        self.propagator = TLEPropagator.selectExtrapolator(self.tleObject)
        
    def setSimPolOrbPass(self,incAngle,satAlt):
        self.simType = 'polOrbPass'
        self.incAngle = radians(incAngle)
        self.satAlt = satAlt*1e3
        
    def isPolOrbPass(self):
        return (self.simType == 'polOrbPass')
    
    def isTLE(self):
        return (self.simType == 'tle')
    


class GroundStation:
    """ Ground station
    
    This class implements a ground station.
    
    Parameters
    ----------
    latitude: float
        Latitude of the ground station [degrees]. Positive towards North.
    longitude: float
        Longitude of the ground station [degrees]. Positive towards East.
    altitude: float
        Altitude of the ground station [m].
    name: str
        Name of the station
    """
    
    def __init__(self,lat,long,alt,name):
        self.latitude = radians(lat)
        self.longitude = radians(long)
        self.altitude = float(alt)
        self.name = name
        
        # define the station Topocentric Frame
        self.geodeticPoint = GeodeticPoint(self.latitude, self.longitude, self.altitude)
        self.frame = TopocentricFrame(earth, self.geodeticPoint, self.name)

# atmospheric models
atmModel = {'TROPICAL': 1,
            'MIDLAT_SUMMER': 2,
            'MIDLAT_WINTER': 3,
            'SUBARCTIC_SUMMER': 4,
            'SUBARCTIC_WINTER': 5,
            'US_STANDARD': 6
            }

# aerosol models
aeroModel = {'NO_AEROSOLS': 0,
             'RURAL_23KM': 1,
             'RURAL_5KM': 2,
             'NAVY': 3,
             'MARITIME': 4,
             'URBAN_5KM': 5,
             'TROPOSPHERIC': 6,
             'USER_DEFINED': 7,
             'FOG1': 8,
             'FOG2': 9}

# mode of execution IEMSCT
execMode = {'TRANSMITTANCE_MODE': 0,
            'THERMAL_RADIANCE_MODE': 1,
            'RADIANCE_MODE': 2,
            'TRANSMITTED_SOLAR_IRRADIANCE': 3}


class AtmosphereTransmittanceModel:
    """ Atmospheric transmittance model
    
    This class contains the atmospheric transmittance model used for calculating
    Tatm. It uses the LOWTRAN7 model.
    
    Parameters
    ----------
    wavelength: float
        Wavelength of the radiation [m].
    altitude: float
        Altitude of the ground station [m].
    model: string
        Atmospheric model [MODEL in LOWTRAN].
    aerosolModel: string
        Aerosol model [IHAZE in LOWTRAN].
    visibility: float
        Surface meteorological range [VIS in LOWTRAN] [km].
    """   
    
    def __init__(self, wavelength, altitude, model='MIDLAT_SUMMER', aerosolModel='NO_AEROSOLS', visibility=0. ):
        
        self.wavelength = wavelength
        self.altitude = altitude

        self.c1 = {'h1': float(self.altitude/1000),
                  'wlshort': float(self.wavelength*1e9),
                  'wllong': float(self.wavelength*1e9),
                  'wlstep': 10,
                  'vis': float(visibility),
                  'gndAlt': float(self.altitude/1000)
                  }

        # set the atmospheric model
        if model in atmModel:
            self.c1['model'] = atmModel[model]
            self.model = model
        else:
            print('[WARNING] Model '+model+' does not exist: use default MIDLAT_SUMMER model.')
            self.c1['model'] = atmModel['MIDLAT_SUMMER']
            self.model = 'MIDLAT_SUMMER'
            
        # set the aerosol model
        if aerosolModel in aeroModel:
            self.c1['ihaze'] = aeroModel[aerosolModel]
            self.aerosolModel = aerosolModel
        else:
            print('[WARNING] Aerosol model '+aerosolModel+' does not exist: use default NO_AEROSOLS model.')
            self.c1['ihaze'] = atmModel['NO_AEROSOLS']
            self.aerosolModel = 'NO_AEROSOLS'


    def calculateTransmittance(self,elevation):
        """ Calculate the atmospheric transmittance
        
        Parameters
        ----------
        elevation: float
            Elevation of the satellite [degrees].
            
        Returns
        -------
        float
            The atmospheric transmittance of the elevation is larger than 0, 
            otherwise it returns 0.
        """
        
        self.c1['itype'] = 3
        self.c1['iemsct'] = 0
        
        if elevation > 0:
            self.c1['angle'] = 90 - elevation
            return float(lowtran.golowtran(self.c1).transmission.values)
        else:
            return 0.
    
    def calculateTransmittanceHorizontal(self,distance):
        """ Calculate the atmospheric transmittance for an horizontal path
        
        The height of the two stations is defined by the altitude parameter.
        
        Parameters
        ----------
        distance: float
            Distance of the two stations [m].
            
        Returns
        -------
        float
            The atmospheric transmittance for the path.
        """
        
        self.c1['itype'] = 1
        self.c1['iemsct'] = 0
        
        self.c1['range_km'] = float(distance/1000)
        
        return float(lowtran.golowtran(self.c1).transmission.values)
    
    def calculateTransmittanceSlant(self,h2,distance):
        """ Calculate the atmospheric transmittance for a slant path
        
        The height of the first station is defined by the altitude parameter 
        while the second station is at an altitude h2.
        
        Parameters
        ----------
        h2: float
            Altitude of the second station [m].
        distance: float
            Distance of the two stations [m].
        
        Returns
        -------
        float
            The atmospheric transmittance for the required path.
        """
        
        self.c1['itype'] = 2
        self.c1['iemsct'] = 0
        
        self.c1['h2'] = float(h2/1000)
        self.c1['range_km'] = float(distance/1000)
        
        return float(lowtran.golowtran(self.c1).transmission.values)
        
        
        
class SimpleDownlinkChannel:
    """ Downlink channel
    
    This class describes a simple downlink channel.
    
    Parameters
    ----------
    satellite: Satellite
        Instance of the Satellite class.
    groundStation: GroundStation
        Instance of the GroundStation class.
    wavelength: float
        Wavelength of the radiation [m].
    atmModel: AtmosphereTransmittanceModel
        Atmospheric transmittance model used in the channel.
    """
    
    def __init__(self,sat,gs,wl,atmModel=None):
        self.satellite = sat
        self.groundStation = gs
        self.wavelength = float(wl)
        
        # initialize the atmospheric transmittance model
        if atmModel is None:
            self.atm = AtmosphereTransmittanceModel(self.wavelength, self.groundStation.altitude)
        else:
            self.atm = atmModel
        
    def calculateChannelParameters(self, timeList):
        """ Calculate channel paramters
        
        This function calculates the parameters of the channels for the times
        given in timeList.
        
        Parameters
        ----------
        timeList : list of :obj:`~datetime.datetime`
            List of times at which to calculate satellite parameters.
            
        Returns
        -------
        tuple (np.array, np.array, np.array)
            The elements of the tuple are the parameters of the channel for the
            different times in timeList:
            - length [m]
            - Tatm
            - elevation [degrees]
        """
        
        if self.satellite.isPolOrbPass():
            # calculate the orbit parameters using the [Moll et al.] model.
            psi = self.groundStation.latitude
            
            deltaI = self.satellite.incAngle
            hs = self.satellite.satAlt
            
            deltaMin = np.arccos( np.cos(psi)*np.cos(deltaI) / np.sqrt(1 - ( np.cos(psi)*np.sin(deltaI) )**2) )
            tMin = timeList[ int(len(timeList)/2) ]
            
            # Earth parameters (from Daniele's code)
            Rt = 6.37e6     # Earth radius
            M = 5.97e24     # Earth mass
            G = 6.67e-11    # Gravitational constant
            
            Omega = np.sqrt( G*M / (Rt + hs)**3)
            
            relTime = np.array([(timeList[i] - tMin).total_seconds() for i in range(len(timeList))])
            delta = Omega * relTime + deltaMin
            
            Zc = np.arccos( np.sin(psi)*np.sin(delta) + np.cos(psi)*np.cos(delta)*np.cos(deltaI) )
            Z = np.arcsin( (Rt+hs)*np.sin(Zc) / np.sqrt( Rt**2 + (Rt+hs)**2 - 2*Rt*(Rt+hs)*np.cos(Zc)) )
            elevation = 90 - np.rad2deg(Z)
            
            channelLength = -Rt*np.cos(Z) + np.sqrt( (Rt*np.cos(Z))**2 + 2*Rt*hs + hs**2)
            
            atmTrans = np.array([self.atm.calculateTransmittance(elevation[i]) for i in range(len(timeList))])
            
            
        elif self.satellite.isTLE():
            channelLength = np.zeros( (len(timeList),) )
            atmTrans = np.zeros( (len(timeList),) )
            elevation = np.zeros( (len(timeList),) )
        
            # calculate the orbit parameters using the TLE
            absDateList = [datetime_to_absolutedate(i) for i in timeList]
            
            for i in range(len(absDateList)):
                pv = self.satellite.propagator.getPVCoordinates(absDateList[i], inertialFrame)
                pos_tmp = pv.getPosition()
    
                frameTrans = inertialFrame.getTransformTo(self.groundStation.frame, absDateList[i])
                
                channelLength[i] = frameTrans.transformPosition(pos_tmp).getNorm()
                
                elevation[i] = np.rad2deg(self.groundStation.frame.getElevation(pv.getPosition(),inertialFrame,absDateList[i]))
                
                atmTrans[i] = self.atm.calculateTransmittance(elevation[i])
            
        return (channelLength, atmTrans, elevation)