import os
import pytest
import numpy as np
from cwest_polymer.fmr_classes import FractionalMRDataset
from piblin.data import Measurement, MeasurementSet


@pytest.fixture()
def test_dir():
    """Returns path to the test_data directory located in the same directory as the test file."""
    test_dir = os.path.join(os.path.dirname(__file__), 'test_data')

    if not os.path.isdir(test_dir):
        raise ValueError(f"Test data directory not found at {test_dir}")

    return test_dir


@pytest.fixture()
def mass():
    """Returns a numpy array of mass values for testing."""

    mass = np.array([50, 100, 150, 200, 250, 63, 113, 163, 213, 263])
    return mass


@pytest.fixture()
def rt(mass):
    rt = np.arange(1, len(mass) + 1)
    return rt


@pytest.fixture()
def abundance(mass):
    abundance = np.ones_like(mass) * 100
    return abundance


@pytest.fixture()
def fmr_values(mass):
    """Returns a numpy array of fractional mass ratio values for testing."""
    fmr_values = np.mod(mass/50, 1)
    return fmr_values


@pytest.fixture()
def clusters(mass):
    """Returns a numpy array of clusters for testing."""
    cluster = np.array([0 if m % 50 == 0 else 1 for m in mass])
    return cluster


@pytest.fixture()
def repeat_units():
    """Returns a list of repeat unit masses for testing."""
    repeat_units = [50, 13]

    return dict(zip(['RU1', 'RU2'], repeat_units))


@pytest.fixture()
def data(mass, rt, abundance):
    """Returns a dictionary representing a dataset for testing."""
    dataset = FractionalMRDataset(
        mass=mass,
        rt=rt,
        abundance=abundance
    )

    measurement = Measurement(
        datasets=[dataset],
        details={},
        conditions={'file_name': 'test_data'}
    )

    return MeasurementSet(
        measurements=[measurement],
        merge_redundant=False
    )
