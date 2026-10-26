"""Bad samples and unbounded work must fail before optional TLS can run."""
import unittest
from unittest.mock import patch

import numpy as np

from discovery.adapters.tls import search


class TestTLSInputContract(unittest.TestCase):
    def test_invalid_samples_are_never_silently_cleaned(self):
        time = np.linspace(100, 110, 100)
        flux = np.ones(100)
        error = np.full(100, 0.001)
        invalid = [
            (np.where(np.arange(100) == 20, np.nan, time), flux, error),
            (time, flux[:-1], error), (time[::-1], flux, error),
            (time, flux, np.zeros(100)), (time, flux * 2, error),
        ]
        parameters = {"use_threads": 1, "period_min": 1.1, "period_max": 8.0}
        with patch("discovery.adapters.tls.installed_provenance", side_effect=AssertionError("Detector must not run")):
            for t, y, e in invalid:
                with self.subTest(kind=(len(t), len(y), len(e))), self.assertRaises(ValueError):
                    search(t, y, e, parameters)

    def test_unvalidated_baseline_and_compute_expansion_are_rejected(self):
        time = np.linspace(100, 110, 100)
        flux = np.ones(100)
        error = np.full(100, 0.001)
        parameters = {"use_threads": 1, "period_min": 1.1, "period_max": 8.0}
        with patch("discovery.adapters.tls.installed_provenance", side_effect=AssertionError("Detector must not run")):
            with self.assertRaises(ValueError):
                search(np.linspace(100, 500, 100), flux, error, parameters)
            for changes in ({"use_threads": 4}, {"use_threads": True}, {"period_max": 100}, {"period_min": 9}):
                with self.subTest(changes=changes), self.assertRaises(ValueError):
                    search(time, flux, error, {**parameters, **changes})
