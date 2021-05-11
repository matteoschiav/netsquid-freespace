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

class Satellite:
    """ Satellite orbit
    
    This class implements a real satellite orbit.
    """
    
    def __init__(self,tle):
        """ Initialization of the satellite orbit from a TLE list """
        self.tleList = tle
        self.tleObject = TLE(*self.tleList)
        self.propagator = TLEPropagator.selectExtrapolator(self.tleObject)



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


class AtmosphereTransmittanceModel:
    """ Atmospheric transmittance model
    
    This class contains the atmospheric transmittance model used for calculating
    Tatm.
    
    Parameters
    ----------
    wavelength: float
        Wavelength of the radiation [m].
    altitude: float
        Altitude of the ground station [m].
    """
    
    def __init__(self, wavelength, altitude):
        MIDLAT_SUMMER = 2
        
        self.wavelength = wavelength
        self.altitude = altitude

        self.c1 = {'model': MIDLAT_SUMMER,
                  'h1': float(self.altitude/1000),
                  'angle': 60,
                  'wlshort': float(self.wavelength*1e9),
                  'wllong': float(self.wavelength*1e9),
                  'wlstep': 10,
                  'itype': 3,
                  'iemsct': 0}

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
        
        if elevation > 0:
            self.c1['angle'] = 90 - elevation
            return float(lowtran.golowtran(self.c1).transmission.values)
        else:
            return 0.
    
        
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
    """
    
    def __init__(self,sat,gs,wl):
        self.satellite = sat
        self.groundStation = gs
        self.wavelength = float(wl)
        
        # initialize the atmospheric transmittance model
        self.atm = AtmosphereTransmittanceModel(self.wavelength, self.groundStation.altitude)
        
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
        
        absDateList = [datetime_to_absolutedate(i) for i in timeList]
        
        channelLength = np.zeros( (len(absDateList),) )
        atmTrans = np.zeros( (len(absDateList),) )
        elevation = np.zeros( (len(absDateList),) )
        
        for i in range(len(absDateList)):
            pv = self.satellite.propagator.getPVCoordinates(absDateList[i], inertialFrame)
            pos_tmp = pv.getPosition()

            frameTrans = inertialFrame.getTransformTo(self.groundStation.frame, absDateList[i])
            
            channelLength[i] = frameTrans.transformPosition(pos_tmp).getNorm()
            
            elevation[i] = np.rad2deg(self.groundStation.frame.getElevation(pv.getPosition(),inertialFrame,absDateList[i]))
            
            atmTrans[i] = self.atm.calculateTransmittance(elevation[i])
            
        return (channelLength, atmTrans, elevation)