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

from netsquid.components.models.qerrormodels import QuantumErrorModel
import netsquid.util.simtools as simtools
from netsquid.util.simlog import warn_deprecated

def randomSample(pdf,xmin=0,xmax=1,rng=None):
    """ Sample a random number from the probability distribution pdf
    
    This function returns a random value in the interval [xmin,xmax] distributed
    according to the distribution pdf(x).
    """
    
    if rng is None:
        rng = simtools.get_random_state()
        
    x = xmin + rng.random_sample() * (xmax-xmin)
    y = rng.random_sample()
    
    while y > pdf(x):
        x = xmin + rng.random_sample() * (xmax-xmin)
        y = rng.random_sample()
    
    return x
    

class FreeSpaceLossModel(QuantumErrorModel):
    """Model for photon loss on a free space channel
    
    Uses beam-wandering PDTC from [Vasylyev et al., PRL 108, 220501 (2012)] to
    sample the loss probability of the photon.
    
    Parameters
    ----------
    W0 : float
        Waist of the beam at the transmitter [m].
    rx_aperture : float
        Radius of the receiving telescope [m].
    Cn2 : float
        Index of refraction structure constant [m**(-2/3)].
    wavelength : float
        Wavelength of the radiation [m].
    rng : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    
    """
    def __init__(self, W0, rx_aperture, Cn2, wavelength, rng=None):
        super().__init__()
        self.rng = rng if rng else simtools.get_random_state()
        self.W0 = W0
        self.rx_aperture = rx_aperture
        self.Cn2 = Cn2
        self.wavelength = wavelength
        self.required_properties = ['length']

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
    def W0(self):
        """float: beam waist at the transmitter [m]."""
        return self.properties['W0']

    @W0.setter
    def W0(self, value):
        if value < 0:
            raise ValueError
        self.properties['W0'] = value

    @property
    def rx_aperture(self):
        """float: radius of the receiving telescope [m]."""
        return self.properties['rx_aperture']

    @rx_aperture.setter
    def rx_aperture(self, value):
        if value < 0:
            raise ValueError
        self.properties['rx_aperture'] = value

    @property
    def Cn2(self):
        """float: index of refraction structure constant [m**(-2/3)]."""
        return self.properties['Cn2']

    @Cn2.setter
    def Cn2(self, value):
        if value < 0:
            raise ValueError
        self.properties['Cn2'] = value

    @property
    def wavelength(self):
        """float: wavelength of the radiation [m]."""
        return self.properties['wavelength']

    @wavelength.setter
    def wavelength(self, value):
        if value < 0:
            raise ValueError
        self.properties['wavelength'] = value
            
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
                            key="FreeSpaceLossModel.compute_model.channel")
            kwargs['length'] = kwargs['channel'].properties["length"]
            del kwargs['channel']
            
        # calculate the parameters used in the PDTC
        z = kwargs['length']*1e3
        W = self.W0*np.sqrt(1 + (z*self.wavelength/(np.pi*self.W0**2))**2)
        X = (self.rx_aperture/W)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        sigma = np.sqrt(1.919 * self.Cn2 * kwargs['length']**3 * (2*self.W0)**(-1./3.))
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X)) / np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)))
        R = self.rx_aperture * np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)) )**(-1./l)
        # define the PDTC
        PDTC = lambda T: 2*R**2/(sigma**2*l*T) * (2*np.log(T0/T))**((2./l)-1) * np.exp( -R**2 * (2*np.log(T0/T))**(2./l) / (2*sigma**2) )

        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            # extract the value of the transmission coefficient
            T = randomSample(PDTC,0,T0,self.rng)
            # calculate the probability of losing the qubit
            prob_loss = 1 - T**2
            self.lose_qubit(qubits, idx, prob_loss, rng=self.properties['rng'])
