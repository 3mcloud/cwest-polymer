from cwest_polymer.fmr_filereaders.fmr_mass_spreadsheet_reader import MassSpreadsheetReader as FileReader
from cwest_polymer.fmr_classes.fmr_datasets import FractionalMRDataset
import os


def test_mass_spreadsheet_reader(test_dir):
    data_dir = os.path.join(test_dir, "mass_list.csv")
    data = FileReader().data_from_filepath(filepath=data_dir)
    assert data.num_measurements == 1
    assert isinstance(data.measurements[0].datasets[0], FractionalMRDataset)
