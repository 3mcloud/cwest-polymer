"""This module contains constants and default values used in the FMR (Fractional Molecular Repeat) analysis.

This module defines default repeat units, data array labels, units, and other parameters used throughout
the FMR analysis pipeline.

The parameters include:
- Default polymer repeat units
- Column header mappings
- Piblin parameters for data and transforms
"""

# default package repeat units
DEFAULT_REPEAT_UNITS = {
    "PEG": "C2 H4 O",
    "PPG": "C3 H6 O",
    "PTHF": "C4 H8 O",
    "PET": "C10 H8 O4",
    "PE": "C2 H4",
    "PP": "C3 H6",
    "Perfluoro": "C F2",
    "PDMS": "C2 H6 Si O",
    "BPA": "C18 H20 O3",
    "Acrylamide": "C3 H5 N O",
    "Acrylic acid": "C3 H4 O2",
    "Nylon 6 6": "C12 H22 N2 O2",
}

# standardized labels for data arrays
MASS_LABEL = 'mass'
MZ_LABEL = 'mz'
RT_LABEL = 'rt'
X_LABEL = 'x'
Y_LABEL = 'y'
ABUNDANCE_LABEL = 'abundance'
MASS_LIST_LABEL = 'mass_list'
FILTER_LABEL = 'filtered'
FMR_LABEL = 'fmr_list'
CLUSTER_LABEL = 'clusters'

CUSTOM_COLUMNS = 'custom_columns'

# unit names for data arrays
MASS_UNIT = 'Da'
RT_UNIT = 'min'
ABUNDANCE_UNIT = 'AU'
CLUSTER_UNIT = 'count'

# FILEREADER: used to check if columns headers match the following and replaces them with standardized names for data
# arrays
COLUMN_UPDATE_DICT = {
    'mass': MASS_LABEL,
    'mz': MZ_LABEL,
    'm/z': MZ_LABEL,
    'rt': RT_LABEL,
    'retention time': RT_LABEL,
    'abundance': ABUNDANCE_LABEL,
    'intensity': ABUNDANCE_LABEL,
    'Area': ABUNDANCE_LABEL,
    'area': ABUNDANCE_LABEL,
    'x_pos': X_LABEL,
    'y_pos': Y_LABEL
}

ACCEPTED_COLUMN_HEADERS = list(COLUMN_UPDATE_DICT.keys())
RESULT_COLUMN_HEADERS = list(COLUMN_UPDATE_DICT.values())

# DATASET: used to create standardized piblin dataset
PARAMETER_LABELS = [MASS_LIST_LABEL, MASS_LABEL, MZ_LABEL, RT_LABEL, X_LABEL, Y_LABEL, ABUNDANCE_LABEL, FILTER_LABEL,
                    CLUSTER_LABEL, FMR_LABEL]

PARAMETER_UNITS = [MASS_UNIT, MASS_UNIT, MASS_UNIT, RT_UNIT, '', '', ABUNDANCE_UNIT, '', CLUSTER_UNIT, MASS_LABEL]
UNITS_DICT = dict(zip(PARAMETER_LABELS, PARAMETER_UNITS))
SOURCE = 'Spreadsheet containing mass spectral data was read for fractional MR polymer analysis'

# TRANSFORMS: labels for measurements details and conditions
CONDITION_RU_LABEL = 'repeat_unit_value'
DETAIL_RU_LABEL = 'repeat_unit_information'

# default name for repeat units if float is given instead of string
RU_LABEL = "RU"
