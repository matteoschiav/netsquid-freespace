#!/usr/bin/env python
# coding: utf-8
"""
Created on Tue May 11 21:07:34 2021

@author: schiavon

This example implements a simple quantum network consisting on a transmitter, placed
on a fixed satellite (or ground telescope) and a receiver on a ground station. The 
channel is implemented using the FixedSatelliteLossModel (or FreeSpaceLossModel), 
whose parameters can be changed by commenting/uncommenting the respective lines.
The two stations perform several runs of the BB84 protocol, calculating the number of
received bits for each run and, from that, the average transmission of the channel.

"""


import netsquid as ns

import netsquid.components.instructions as instr
from scipy.stats import bernoulli
import numpy as np

from netsquid.components import QuantumChannel
from netsquid.protocols import NodeProtocol
from netsquid.components.qprocessor import QuantumProcessor
from netsquid.nodes.network import Network
from netsquid.components.qprocessor import PhysicalInstruction

from netsquid_freespace.lossmodel import FreeSpaceLossModel, FixedSatelliteLossModel


#### Free space channel parameter
W0 = 0.1
rx_aperture_freespace = 1
Tatm = 1
# Atmospheric turbulence
Cn2_freespace = 0
# Cn2_freespace = 1e-15


### Satellite to Ground channel parameters
distsat = 450 #distance from the satellite to the ground station in km
txDiv = 10e-6 #beam divergence
rx_aperture_sat = 1 #aperture of the receiving telescope
# Atmospheric turbulence
Cn2_sat = 0 
# Cn2_sat = 1e-15
# Pointing error
sigmaPoint = 0.5e-6 
# sigmaPoint = 0


#Light parameter
wavelength = 1550e-9 # telecom wavelength
# wavelength = 850e-9 # Micius satellite
c = 299792.458 #speed of light in km/s

#Node parameters
init_time = 100 #time to create a qubit in ns 


#physical instructions
physical_instructions = [
    PhysicalInstruction(instr.INSTR_INIT, duration=init_time),
    PhysicalInstruction(instr.INSTR_H, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_X, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_Z, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_S, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_I, duration=1, parallel=True, topology=[0]),
    PhysicalInstruction(instr.INSTR_CNOT, duration=4, parallel=True),
    PhysicalInstruction(instr.INSTR_MEASURE, duration=1, parallel=True,
                        topology=[0]),
    PhysicalInstruction(instr.INSTR_MEASURE_BELL, duration = 1, parallel=True),
    PhysicalInstruction(instr.INSTR_SWAP, duration = 1, parallel=True)
]

#Protocol definitions
class SendProtocol(NodeProtocol):
    
    #Protocol performed by a node to send a BB84 qubit. It assumes the network is already created and 
    # the node are connected
    
    #Parameter:
    # node: sending node
    # othernode: receiving node
    
    def __init__(self,othernode, node):
        super().__init__(node=node)
        self._othernode = othernode
    
    def run(self):
            if self.node.name[0:9]== 'Satellite':
                mem = self.node.subcomponents["SatelliteMemory"]
            else:
                mem = self.node.subcomponents["StationMemory"]
        
            while True:
                mem.reset()
                mem.execute_instruction(instr.INSTR_INIT,[0])
                yield self.await_program(mem,await_done=True,await_fail=True)
            
                base = bernoulli.rvs(0.5) #random choice of a basis
                if base <0.5:
                    mem.execute_instruction(instr.INSTR_H,[0])
                    base = "plusmoins"
                else:
                    mem.execute_instruction(instr.INSTR_I,[0])
                    base = "zeroun"
                
                yield self.await_program(mem,await_done=True,await_fail=True)

                bit = bernoulli.rvs(0.5) #random choice of a bit
                if bit < 0.5:
                    mem.execute_instruction(instr.INSTR_I, [0], physical=False)
                else:
                    if base == "zeroun":
                        mem.execute_instruction(instr.INSTR_X, [0], physical=False)
                    elif base == "plusmoins":
                        mem.execute_instruction(instr.INSTR_Z, [0], physical=False)
                
                qubit, = mem.pop([0])

class ReceiveProtocol(NodeProtocol):
    
    # Protocol performed by a node to receive a state a measure it
    
    #Parameters:
    # othernode: node from which a qubit is expected
    # port : receiving port
    
        def __init__(self, othernode, port, node):
            super().__init__(node=node)
            self._othernode = othernode
            self._port = port

        def run(self):
            if self.node.name[0:9]== 'Satellite':
                mem = self.node.subcomponents["SatelliteMemory"]
            else:
                mem = self.node.subcomponents["StationMemory"]
            while True:
                yield self.await_port_input(self._port)

                base = bernoulli.rvs(0.5) #choose a random basis
                if base < 0.5:
                    mem.execute_instruction(instr.INSTR_H, [0], physical = False)
                    base = "plusmoins"
                else:
                    mem.execute_instruction(instr.INSTR_I, [0],physical = False)
                    base = "zeroun"
                        
                m,_,_ = mem.execute_instruction(instr.INSTR_MEASURE,[0],output_key="M1")
                yield self.await_program(mem,await_done=True,await_fail=True)
                        
                Key.append(m['M1'][0])
                mem.reset()
                        
#Network and node initialisation
ns.sim_reset()
network = Network("TestNetwork")
Satellite, GroundStation = network.add_nodes(["Satellite", "GroundStation"])
qmem1 = QuantumProcessor("SatelliteMemory", num_positions=4,
                                phys_instructions=physical_instructions)
Satellite.add_subcomponent(qmem1)

qmem2 = QuantumProcessor("StationMemory", num_position=4, 
                         phys_instructions = physical_instructions)
GroundStation.add_subcomponent(qmem2)

Key = [] 


#Channel initialisation with the satellite loss model
qchannel = QuantumChannel("Downlink",length=distsat, delay=1, models={"quantum_loss_model": FixedSatelliteLossModel(txDiv, sigmaPoint,
                                                                            rx_aperture_sat, Cn2_sat, wavelength,Tatm)})


#Channel initialisation with the free space loss model (uncomment to apply, and comment previous channel initialisation part)
# qchannel = QuantumChannel("FreespaceLink",length=distsat, delay=1, models={"quantum_loss_model": FreeSpaceLossModel(W0, rx_aperture_freespace,
#       Cn2_freespace, wavelength, Tatm)})

#port management
Sat_send, Station_receive = network.add_connection(
                    Satellite, GroundStation, channel_to=qchannel, label="Downlink")

GroundStation.ports[Station_receive].forward_input(qmem2.ports["qin"]) 
qmem1.ports["qout"].forward_output(Satellite.ports[Sat_send])

#Initialisation of the protocol
SatProtocol = SendProtocol(GroundStation, Satellite)
SatProtocol.start()

StationProtocol = ReceiveProtocol(Satellite,GroundStation.ports[Station_receive], GroundStation)
StationProtocol.start()

#%% Run the simulation several times

numSimulations = 100
siftBits = np.zeros((numSimulations,))

simDuration = 10000 # simulation time [ns]

print('Transmitter beam radius:',qchannel.models['quantum_loss_model'].W0)

sentBits = int(simDuration/init_time)
print('Qubits sent:',sentBits)

for i in range(numSimulations):
    # a list to store the measurement outputs of qubit arriving from the satellite
    Key = []
    
    #Run simulation
    stat =ns.sim_run(duration=simDuration)
    
    #Show resulting key
    print('Sifted bits:',len(Key),'->',Key)
    
    siftBits[i] = len(Key)

print('Mean sifted bit ratio:',np.mean(siftBits)/sentBits)