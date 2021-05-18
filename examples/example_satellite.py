#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 21:07:34 2021

@author: schiavon

This examples calculates the orbital parameters for the passage of a real satellite
(in this case QSS-Micius), using the two-line elements. It then calculates the
elevation, the channel length and the atmospheric losses with respect to two ground
stations, placed in Paris and Delft. These parameters can be input to the
FixedSatelliteLossModel to implement the corresponding channel on netsquid.

"""

from netsquid_freespace import channel

from datetime import datetime, timedelta

import numpy as np
from matplotlib import pyplot as plt

#%% Initialize channel paramters
wavelength = 1550e-9

#%% Initialize the satellite

# Micius (QSS)
tleMicius = ['1 41731U 16051A   21117.42314584  .00000696  00000-0  30260-4 0  9998',
             '2 41731  97.3499  30.8507 0012844 347.0485 124.2616 15.25507799261429']

satMicius = channel.Satellite(tleMicius)

#%% Initialize the ground stations

# Paris
latParis = 48.857
longParis = 2.352
altParis = 80.

staParis = channel.GroundStation(latParis, longParis, altParis, 'Paris')

# Delft
latDelft = 52.012
longDelft = 4.357
altDelft = 0.

staDelft = channel.GroundStation(latDelft, longDelft, altDelft, 'Delft')

#%% Initialize the downlink channels

downSatParis = channel.SimpleDownlinkChannel(satMicius, staParis, wavelength)

downSatDelft = channel.SimpleDownlinkChannel(satMicius, staDelft, wavelength)

#%% Initialize the time array

dt = startTime = datetime(2021, 5, 15, 0, 0, 0)
endTime = datetime(2021, 5, 15, 0, 20, 0)
timeStep = timedelta(seconds = 10.)

timeList = []
while dt < endTime:
    timeList.append(dt)
    dt += timeStep
    
#%% Calculate the orbital parameters for the two channels

lenSatParis, tSatParis, elSatParis = downSatParis.calculateChannelParameters(timeList)

lenSatDelft, tSatDelft, elSatDelft = downSatDelft.calculateChannelParameters(timeList)

#%% Plot data

times = np.array([ (timeList[i]-timeList[0]).seconds  for i in range(len(timeList)) ])

plt.figure(figsize=(18,6))

plt.subplot(131)
plt.plot(times/60,elSatParis,'b')
plt.plot(times/60,elSatDelft,'r')
plt.ylim([0,90])
plt.ylabel('Elevation [degrees]')
plt.xlabel('Passage time [minutes]')
plt.legend(['Paris','Delft'])

plt.subplot(132)
plt.plot(times/60,lenSatParis/1000,'b')
plt.plot(times/60,lenSatDelft/1000,'r')
plt.ylabel('Channel length [km]')
plt.xlabel('Passage time [minutes]')
plt.legend(['Paris','Delft'])

plt.subplot(133)
plt.plot(times/60,tSatParis,'b')
plt.plot(times/60,tSatDelft,'r')
plt.ylabel('Tatm')
plt.xlabel('Passage time [minutes]')
plt.legend(['Paris','Delft'])

plt.suptitle('Micius satellite - startDate 15/05/2021 00:00')