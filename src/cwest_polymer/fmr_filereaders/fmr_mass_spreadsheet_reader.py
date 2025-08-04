"""A file reader module for importing mass spectrometry data from spreadsheet files.

This module provides functionality to read mass spectrometry data from CSV and Excel files
into FractionalMRDataset objects through piblin FileReader class. It supports custom column mapping for mass, m/z,
retention time, spatial position, and abundance values.

Default supported column headers (case-insensitive):
    - mass: 'mass'
    - m/z: 'mz', 'm/z'
    - abundance: 'abundance', 'intensity', 'area'
    - retention time: 'rt', 'time'
    - x position: 'x_pos'
    - y position: 'y_pos'

Parameters are referenced from fmr_parameters.py
Custom column mappings can be provided via read_kwargs using the CUSTOM_COLUMNS parameter.
"""

import pandas as pd
from typing import Dict
from pathlib import Path

from piblin.dataio.fileio.read.file_reader import FileReader
from piblin.data.data_collections.measurement import Measurement
from piblin.data.data_collections.measurement_set import MeasurementSet

from ..fmr_classes.fmr_datasets import FractionalMRDataset
from .. import fmr_parameters as p


class MassSpreadsheetReader(FileReader):
    """
    A piblin file reader for mass spectrometry data in spreadsheet format (csv, xlsx). read_kwargs enable custom column
    reading.
    """
    supported_extensions = {'csv', 'xlsx'}

    def _read_file_contents(self, filepath: str, **read_kwargs):
        """Read the contents of a file into a single object.

        Parameters
        ----------
        filepath : str
            The path to the file to be read.
        read_kwargs : dict
            used for custom column names using 'mass','mz','rt','x','y','abundance' as keys (check fmr_parameters if
            values change)
        """

        column_update_dict = read_kwargs.get(p.CUSTOM_COLUMNS, None)

        if column_update_dict is None:
            column_update_dict = p.COLUMN_UPDATE_DICT
        else:
            column_update_dict = {x.lower(): y for x, y in column_update_dict.items()}

        file_contents = {}
        extension = Path(filepath).suffix.lower()
        if extension == '.csv':
            df = pd.read_csv(filepath)
        elif extension == '.xlsx':
            df = pd.read_excel(filepath)
        else:
            raise ValueError("file extension incorrect")

        for column in df.columns:
            if column.lower() in [x.lower() for x in column_update_dict.keys()]:
                data_key = column_update_dict[column.lower()]
                data_value = df[column].to_numpy()

                file_contents[data_key] = data_value

        return file_contents

    @property
    def default_mode(self) -> str:
        return ''

    @classmethod
    def _data_from_file_contents(cls, file_contents: Dict, file_location=None, file_name=None,
                                 **read_kwargs) -> MeasurementSet:

        mass = file_contents.get(p.MASS_LABEL, None)
        mz = file_contents.get(p.MZ_LABEL, None)
        rt = file_contents.get(p.RT_LABEL, None)
        abundance = file_contents.get(p.ABUNDANCE_LABEL, None)

        result = FractionalMRDataset(
            mass=mass,
            mz=mz,
            rt=rt,
            abundance=abundance
        )

        # removes extension so the file_name can be used in subsequent steps
        file_name = Path(file_name).with_suffix('')

        measurements = [Measurement(datasets=[result], conditions={'file_name': file_name})]
        return MeasurementSet(measurements=measurements, merge_redundant=False)
