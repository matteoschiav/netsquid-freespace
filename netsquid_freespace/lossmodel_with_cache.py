from netsquid_freespace.lossmodel import FreeSpaceLossModel


def get_freespacelossmodel_class_with_cache(freespacelossmodelcls):
    """Returns a subclass of :cls:`~netsquid_freespace.lossmodel.FreeSpaceLossModel` with
    the same functionality but potentially faster runtime.
    This is achieved by caching the length of the channel, so that the Weibull distribution
    needs not be computed every time the same length is called.

    **Important note:** every time the input parameters to the class `spacelossmodelcls`
    (for example for :cls:`~netsquid_freespace.FixedSatelliteLossModel`: `txDiv`, `sigmaPoint`, etc.)
    set to a new value, the cache should be manually reset by called the method `reset`.

    Parameters
    ----------
    freespacelossmodelcls : subclass of :cls:`~netsquid_freespace.lossmodel.FreeSpaceLossModel`

    Returns
    -------
    subclass of :cls:`~netsquid_freespace.lossmodel.FreeSpaceLossModel`

    Example usage
    -------------

    >>> # Define an example set of parameters
    >>> txDiv = 10e-6
    >>> sigmaPoint = 0
    >>> rx_aperture = 0.75
    >>> Cn2 = 0
    >>> wavelength = 1550e-9
    >>> length = 1000

    >>> # define a new loss model class, with cache
    >>> FixedSatelliteLossModelWithCache = \
    >>>     get_freespacelossmodel_class_with_cache(freespacelossmodelcls=FixedSatelliteLossModel)

    >>> # initialize an object of the new loss model with cache
    >>> loss_model = FixedSatelliteLossModelWithCache(txDiv=txDiv,
    >>>                                               sigmaPoint=sigmaPoint,
    >>>                                               rx_aperture=rx_aperture,
    >>>                                               Cn2=Cn2,
    >>>                                               wavelength=wavelength)
    >>> # now the new loss model can be used as one could use any
    >>> # FreeSpaceLossModel object
    """

    if not issubclass(freespacelossmodelcls, FreeSpaceLossModel):
        raise TypeError("{} is not a subclass of FreeSpaceLossModel".format(freespacelossmodelcls))

    class SpaceLossModelWithCache(freespacelossmodelcls):

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self._cached_length_to_weibull_params = {}
            self._reset_cache()

        def reset(self):
            super().reset()
            self._reset_cache()

        def _reset_cache(self):
            self._cached_length_to_weibull_params = {}

        def _compute_weibull_loss_model_parameters(self, length):
            if length in self._cached_length_to_weibull_params:
                return self._cached_length_to_weibull_params[length]
            else:
                weibull_params = super()._compute_weibull_loss_model_parameters(length=length)
                self._cached_length_to_weibull_params[length] = weibull_params
                return weibull_params

    return SpaceLossModelWithCache
