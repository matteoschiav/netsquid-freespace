# This is an example module file with an associated unit test.
# You should replace these files with your own content.
#
# We suggest you use Numpy style docstrings in case your work could be integrated
# into NetSquid in the future:
# https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
"""Atmospheric loss model for the free space channel

This file contains the different loss models implemented for the free space channel.

"""
import numpy as np
from scipy.special import i0, i1

import netsquid as ns
from netsquid.components.models.qerrormodels import QuantumErrorModel
import netsquid.util.simtools as simtools
from netsquid.util.simlog import warn_deprecated



class FreeSpaceLossModel(QuantumErrorModel):
    """Model for photon loss on a free space channel
    
    Uses beam-wandering PDTC from [Vasylyev et al., PRL 108, 220501 (2012)] to
    sample the loss probability of the photon.
    
    Parameters
    ----------
    rng : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    
    """
    def __init__(self, length, W0, rx_aperture, Cn2, wavelength, rng=None):
        super().__init__()
        self.rng = rng if rng else simtools.get_random_state()
        
        # TODO: substitute kwargs with parameters of the model
        # -> length could stay in kwargs, the other ones must be in the parameter of init
        
    def PDTC(self,T,**kwargs):
        # TODO: put the initialization of these paramters in the __init__ function
        z = kwargs['length']*1e3
        W = kwargs['W0']*np.sqrt(1 + (z*kwargs['wavelength']/(np.pi*kwargs['W0']**2))**2)
        X = (kwargs['rx_aperture']/W)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        sigma = np.sqrt(1.919 * kwargs['Cn2'] * kwargs['length']**3 * (2*kwargs['W0'])**(-1./3.))
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X)) / np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)))
        R = kwargs['rx_aperture'] * np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)) )**(-1./l)
        
        return 2*R**2/(sigma**2*l*T) * (2*np.log(T0/T))**((2./l)-1) * np.exp( -R**2 * (2*np.log(T0/T))**(2./l) / (2*sigma**2) )
    

class FibreLossModel(QuantumErrorModel):
    """Model for exponential photon loss on fibre optic channels.

    Uses length of transmitting channel to sample an
    exponential loss probability.

    Parameters
    ----------
    p_loss_init : float, optional
        Initial probability of losing a photon once it enters a channel.
        e.g. due to frequency conversion.
    p_loss_length : float, optional
        Photon survival probability per channel length [dB/km].
    rng : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.

    """
    def __init__(self, p_loss_init=0.2, p_loss_length=0.25, rng=None):
        super().__init__()
        self.p_loss_init = p_loss_init
        self.p_loss_length = p_loss_length
        self.rng = rng if rng else simtools.get_random_state()
        self.required_properties = ["length"]

    @property
    def rng(self):
        """ :obj:`~numpy.random.RandomState`: Random number generator."""
        return self.properties['rng']

    @rng.setter
    def rng(self, value):
        if not isinstance(value, np.random.RandomState):
            raise TypeError("{} is not a valid numpy RandomState".format(value))
        self.properties['rng'] = value

    @property
    def p_loss_init(self):
        """float: initial probability of losing a photon when it enters channel."""
        return self.properties['p_loss_init']

    @p_loss_init.setter
    def p_loss_init(self, value):
        if not 0 <= value <= 1:
            raise ValueError
        self.properties['p_loss_init'] = value

    @property
    def p_loss_length(self):
        """float: photon survival probability per channel length [dB/km]."""
        return self.properties['p_loss_length']

    @p_loss_length.setter
    def p_loss_length(self, value):
        if value < 0:
            raise ValueError
        self.properties['p_loss_length'] = value

    def error_operation(self, qubits, delta_time=0, **kwargs):
        """Error operation to apply to qubits.

        Parameters
        ----------
        qubits : tuple of :obj:`~netsquid.qubits.qubit.Qubit`
            Qubits to apply noise to.
        delta_time : float, optional
            Time qubits have spent on a component [ns].

        """
        if 'channel' in kwargs:
            warn_deprecated("channel parameter is deprecated. "
                            "Pass length parameter directly instead.",
                            key="FibreLossModel.compute_model.channel")
            kwargs['length'] = kwargs['channel'].properties["length"]
            del kwargs['channel']
        #self.apply_loss(qubits, delta_time, **kwargs)
        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            prob_loss = 1 - (1 - self.p_loss_init) * np.power(10, - kwargs['length'] * self.p_loss_length / 10)
            self.lose_qubit(qubits, idx, prob_loss, rng=self.properties['rng'])

    def prob_item_lost(self, item, delta_time=0, **kwargs):
        # DEPRECATED
        return 1 - (1 - self.p_loss_init) * np.power(10, - kwargs['length'] * self.p_loss_length / 10)
