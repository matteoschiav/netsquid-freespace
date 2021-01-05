import unittest
import netsquid.qubits.qubitapi as qapi
from netsquid_freespace.lossmodel import FreeSpaceLossModel


class TestFreeSpaceLossModel(unittest.TestCase):

    def test_weibull_loss_model_computation(self):

        # Define an example set of parameters
        # TODO replace by values for which you know the outcome
        W0 = 1
        rx_aperture = 1
        Cn2 = 0.5
        wavelength = 1500e-9
        length = 100

        # expected parameters of the distribution
        # TODO replace by values we really expect
        expected_shape = 0
        expected_scale = 1
        expected_T0 = 1

        # computed parameters
        loss_model = FreeSpaceLossModel(W0=W0,
                                        rx_aperture=rx_aperture,
                                        Cn2=Cn2,
                                        wavelength=wavelength)
        shape, scale, T0 = loss_model._compute_weibull_loss_model_parameters(length=length)

        # check that the expected parameters are the same as
        # the computed ones
        self.assertAlmostEqual(shape, expected_shape)
        self.assertAlmostEqual(scale, expected_scale)
        self.assertAlmostEqual(T0, expected_T0)

    def test_loss_is_applied_to_qubits(self):
        """We override the `_sample_loss_probability` method
        of FreeSpaceLossModel to investigate whether FreeSpaceLossModel
        works correctly in the following edge cases:

        - when qubits are always lost
        - when qubits are never lost
        """
        self._test_case_qubits_always_lost()
        self._test_case_qubits_never_lost()

    def _test_case_qubits_always_lost(self):
        # Edge case: qubit is always lost

        # Dummy values
        W0 = 1
        rx_aperture = 1
        Cn2 = 0.5
        wavelength = 1500e-9
        length = 100

        class AlwaysLossFreeSpaceLossModel(FreeSpaceLossModel):
            """
            A loss model which always loses qubits.
            """

            def _sample_loss_probability(self, length):
                return 1

        loss_model = AlwaysLossFreeSpaceLossModel(W0=W0,
                                                  rx_aperture=rx_aperture,
                                                  Cn2=Cn2,
                                                  wavelength=wavelength)

        qubits = qapi.create_qubits(10)
        loss_model.error_operation(qubits=qubits, length=length)
        for qubit in qubits:
            self.assertTrue(qubit is None)

    def _test_case_qubits_never_lost(self):
        # Edge case: qubit is never lost

        # Dummy values
        W0 = 1
        rx_aperture = 1
        Cn2 = 0.5
        wavelength = 1500e-9
        length = 100

        class NeverLossFreeSpaceLossModel(FreeSpaceLossModel):
            """
            A loss model which always loses qubits.
            """

            def _sample_loss_probability(self, length):
                return 0

        loss_model = NeverLossFreeSpaceLossModel(W0=W0,
                                                 rx_aperture=rx_aperture,
                                                 Cn2=Cn2,
                                                 wavelength=wavelength)

        qubits = qapi.create_qubits(10)
        loss_model.error_operation(qubits=qubits, length=length)
        for qubit in qubits:
            self.assertFalse(qubit is None)


if __name__ == "__main__":
    unittest.main()
