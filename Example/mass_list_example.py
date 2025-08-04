from cwest_polymer import MassSpreadsheetReader, transforms
from cwest_polymer import fmr_parameters as p
from piblin.data import MeasurementSet
import pandas as pd
import numpy as np
import os

folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# directory to .csv files
path = os.path.join(folder_path, "Data", "mass_list.csv")
path1 = os.path.join(folder_path, "Data", "custom_columns.csv")

result_path = os.path.join(folder_path, "Results", "spreadsheet_results")

ppm_tol = 10
mz_tol = 0.001
min_samples = 5

custom_columns = {
    'a': 'mass',
    'b': 'rt',
    'c': 'abundance'
}

read_kwargs = {
    p.CUSTOM_COLUMNS: custom_columns
}

# read files using default and custom columns for the mass spec data
data = MassSpreadsheetReader().data_from_filepath(filepath=path)
data1 = MassSpreadsheetReader().data_from_filepath(filepath=path1, **read_kwargs)

# merge into measurementset
data = MeasurementSet.from_measurement_sets(measurement_sets=[data, data1], merge_redundant=False)

# create transforms and pipeline
calc_fmr = transforms.FractionalMRTransform.create(repeat_units='C2 H4 O', fractional_values=1, default_list=False)
calc_fmr_clusters = transforms.ClusterTransform.create(mz_tol=mz_tol, ppm_tol=ppm_tol, min_samples=min_samples)
pipeline = calc_fmr + calc_fmr_clusters

# run the pipeline
results = pipeline(data)

fmr_results = results.split_by_condition_name('file_name')

# analyze by file
for result in fmr_results:
    for measurement in result.measurements:
        if measurement.datasets[0].number_of_points() == 0:
            continue
        name = measurement.conditions['file_name']

        ru_name = measurement.details['repeat_unit_information'][0]
        df = pd.DataFrame(np.array(measurement.datasets[0].data_arrays).T,
                          columns=measurement.datasets[0].data_array_names)
        df.to_csv(os.path.join(result_path, f'result_{name}_{ru_name}.csv'))
        fig, _ = measurement.visualize()
        fig.savefig(os.path.join(result_path, f'plot_{name}_{ru_name}.png'), dpi=1000, bbox_inches='tight')
