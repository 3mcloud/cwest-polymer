"""MALDI mass list reader module.

This module provides functionality to read and parse mass spectrometry data from MALDI Excel exports through the pibin
FileReader class. It handles multi-sheet xlsx files containing mass lists with specific column headers for m/z values,
intensities, and other mass spectrometry parameters found in COLUMN_HEADERS below.
"""

from typing import Dict
import os
import pandas as pd

from piblin.dataio import FileReader
from piblin.data import Measurement, MeasurementSet

from ..fmr_classes.fmr_datasets import FractionalMRDataset
from .. import fmr_parameters as p

COLUMN_HEADERS = ['m/z', 'time', 'Intens.', 'SN', 'Quality Fac.', 'Res.', 'Area', 'Rel. Intens.', 'FWHM', 'Chi^2',
                  'Bk. Peak']


class MaldiMassListReader(FileReader):
    """
    Reads .xlsx files from maldi export. Reader for .xlsx document where the sheets are named based on the sample
    information.
    """
    supported_extensions = {'xlsx'}

    def _read_file_contents(self, filepath: str, **read_kwargs):
        """Read the contents of a file into a single object.

        Parameters
        ----------
        filepath : str
            The path to the file to be read.
        read_kwargs : dict
            Read parameters from piblin file reader class. Also enables custom columns not recommended for MALDI
            results, but can be used to import different values.
        """

        column_update_dict = read_kwargs.get(p.CUSTOM_COLUMNS, None)

        all_sheets = pd.read_excel(filepath, sheet_name=None)

        file_contents = {}
        for sheet_name, df in all_sheets.items():
            if len(df) == 0:
                continue
            if column_update_dict is None:
                column_update_dict = p.COLUMN_UPDATE_DICT

            skip = 0
            while not all([x in df.columns for x in COLUMN_HEADERS]):
                skip += 1
                df = pd.read_excel(filepath, sheet_name=sheet_name, header=skip)

            sheet_contents = {}
            for column in df.columns:
                if column.lower() in column_update_dict.keys():
                    data_key = column_update_dict[column.lower()]
                    data_value = df[column].to_numpy()

                    sheet_contents[data_key] = data_value
            file_contents[sheet_name] = sheet_contents

        return file_contents

    @property
    def default_mode(self) -> str:
        return ''

    @classmethod
    def _data_from_file_contents(cls, file_contents: Dict, file_location=None, file_name=None,
                                 **read_kwargs) -> MeasurementSet:
        measurements = []
        for sheet_name, contents in file_contents.items():
            mass = contents.get(p.MASS_LABEL, None)
            mz = contents.get(p.MZ_LABEL, None)
            rt = contents.get(p.RT_LABEL, None)
            abundance = contents.get(p.ABUNDANCE_LABEL, None)

            result = FractionalMRDataset(
                mass=mass,
                mz=mz,
                rt=rt,
                abundance=abundance
            )

            file_path = os.path.join(file_location, file_name)

            # xlsx sheet name is in place of the file, maintains uniqueness from single file with multiple sheets
            measurements.append(
                Measurement(datasets=[result], conditions={'file_name': sheet_name}, details={'file_path': file_path})
            )

        return MeasurementSet(measurements=measurements, merge_redundant=False)
