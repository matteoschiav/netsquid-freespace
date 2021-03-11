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
from numpy.random import weibull

from netsquid.components.models.qerrormodels import QuantumErrorModel
import netsquid.util.simtools as simtools
from netsquid.util.simlog import warn_deprecated

__all__ = [
    'FreeSpaceLossModel',
    'FixedSatelliteLossModel'
]


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
    Tatm : float
        Atmospheric transmittance (square of the transmission coefficient).
    rng : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """
    def __init__(self, W0, rx_aperture, Cn2, wavelength, Tatm=1, rng=None):
        super().__init__()
        self.rng = rng if rng else simtools.get_random_state()
        self.W0 = W0
        self.rx_aperture = rx_aperture
        self.Cn2 = Cn2
        self.wavelength = wavelength
        self.Tatm = Tatm
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
    def Tatm(self):
        """ :float: atmosphere transmittance. """
        return self.properties['Tatm']

    @Tatm.setter
    def Tatm(self, value):
        if (value < 0) or (value > 1):
            raise ValueError
        self.properties['Tatm'] = value

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

    def _compute_weibull_loss_model_parameters(self, length):
        """
        Parameters
        ----------
        length: float
            Length of the channel.

        Returns
        -------
        tuple (float, float, float)
            The elements of the tuple are properties of the
            Weibull distribution. From left to right:
            - the 'shape' parameter
            - the 'scale' parameter
            - 'T0'
        """
        # TODO improve explanation in the docstring of this method
        # TODO explain the model below in the main docstring of this class
        # calculate the parameters used in the PDTC
        
        z = length*1e3
        W = self.W0*np.sqrt(1 + (z*self.wavelength/(np.pi*self.W0**2))**2)
        X = (self.rx_aperture/W)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        sigma = np.sqrt(1.919 * self.Cn2 * z**3 * (2*self.W0)**(-1./3.))
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X)) / np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)))
        R = self.rx_aperture * np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)) )**(-1./l)

        # define the parameters of the Weibull distribution
        a = 2/l
        scaleL = (2*(sigma/R)**2)**(l/2)

        return (a, scaleL, T0)

    def _sample_loss_probability(self, length):
        """
        Parameters
        ----------
        length: float
            Length of the channel.

        Returns
        -------
        float
        """
        a, scaleL, T0 = self._compute_weibull_loss_model_parameters(length=length)

        # extract the value of the transmission coefficient
        # print('Transmission coefficient extraction')
        x = weibull(a, 1)
        scaleX = scaleL * x
        T = T0*np.exp(-scaleX/2)
        # print('T =',T)
        # calculate the probability of losing the qubit
        prob_loss = 1 - self.Tatm * T**2
        return prob_loss

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

#        # test the distribution
#        x = weibull(a,10000)
#        scaleX = scaleL * x
#        T = T0*np.exp(-scaleX/2)
#        import matplotlib.pyplot as plt
#        plt.hist(T,1000)

        for idx, qubit in enumerate(qubits):
            if qubit is None:
                continue
            prob_loss = self._sample_loss_probability(length=kwargs['length'])
#            print(prob_loss)
            self.lose_qubit(qubits, idx, prob_loss, rng=self.properties['rng'])


class FixedSatelliteLossModel(FreeSpaceLossModel):
    """Model for photon loss on a satellite-to-ground static channel

    Uses beam-wandering PDTC from [Vasylyev et al., PRL 108, 220501 (2012)] to
    sample the loss probability of the photon.

    Parameters
    ----------
    txDiv : float
        Divergence of the beam sent from the satellite [rad].
    sigmaPoint :
        Pointing error of the satellite, standard deviation [rad].
    rx_aperture : float
        Radius of the receiving telescope [m].
    Cn2 : float
        Index of refraction structure constant [m**(-2/3)].
    wavelength : float
        Wavelength of the radiation [m].
    Tatm : float
        Atmospheric transmittance (square of the transmission coefficient).
    rng : :obj:`~numpy.random.RandomState` or None, optional
        Random number generator to use. If ``None`` then
        :obj:`~netsquid.util.simtools.get_random_state` is used.
    """
    def __init__(self, txDiv, sigmaPoint, rx_aperture, Cn2, wavelength, Tatm=1, rng=None):
        super().__init__(0.21*wavelength/txDiv,rx_aperture,Cn2,wavelength,Tatm,rng)
        self.txDiv = txDiv
        self.sigmaPoint = sigmaPoint
        self.required_properties = ['length']

    @property
    def txDiv(self):
        """float: divergence of the beam at the transmitter (satellite) [m]."""
        return self.properties['txDiv']

    @txDiv.setter
    def txDiv(self, value):
        if value < 0:
            raise ValueError
        self.properties['txDiv'] = value

    @property
    def sigmaPoint(self):
        """float: pointing error at the transmitter (satellite) [m]."""
        return self.properties['sigmaPoint']

    @sigmaPoint.setter
    def sigmaPoint(self, value):
        if value < 0:
            raise ValueError
        self.properties['sigmaPoint'] = value
        
    def _compute_weibull_loss_model_parameters(self, length):
        """
        Parameters
        ----------
        length: float
            Length of the channel.

        Returns
        -------
        tuple (float, float, float)
            The elements of the tuple are properties of the
            Weibull distribution. From left to right:
            - the 'shape' parameter
            - the 'scale' parameter
            - 'T0'
        """
        # TODO improve explanation in the docstring of this method
        # TODO explain the model below in the main docstring of this class
        # calculate the parameters used in the PDTC
        
        # this function cannot be used for range values lower than 10 km
        if length <= 10:
            raise ValueError
        
        z = length*1e3
        W = self.txDiv * z
        X = (self.rx_aperture/W)**2
        T0 = np.sqrt(1 - np.exp(-2*X))
        sigmaTurb = np.sqrt(1.919 * self.Cn2 * 10e3**3 * (2*self.txDiv*(z-10e3))**(-1./3.))
        sigma = np.sqrt( (self.sigmaPoint*z)**2 + sigmaTurb**2 )
        l = 8 * X * np.exp(-4*X) * i1(4*X) / (1 - np.exp(-4*X)*i0(4*X)) / np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)))
        R = self.rx_aperture * np.log( 2*T0**2/(1-np.exp(-4*X)*i0(4*X)) )**(-1./l)

        # define the parameters of the Weibull distribution
        a = 2/l
        scaleL = (2*(sigma/R)**2)**(l/2)

        return (a, scaleL, T0)
