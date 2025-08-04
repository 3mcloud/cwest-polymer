from cwest_polymer import MassSpreadsheetReader, transforms
import pandas as pd
import numpy as np
import os

folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# fatty acids calculated based on molecular weights
path = os.path.join(folder_path, "Data", "lipid_data.csv")
result_path = os.path.join(folder_path, "Results", "lipid_results")

ppm_tol = 10
mz_tol = 0.001
min_samples = 2

repeat_units = {
    'alkanes': 'C H2',
    'alkenes': 'C2 H2',
    'dehydrogenate': 'H2',
}

data = MassSpreadsheetReader().data_from_filepath(filepath=path)

calc_fmr = transforms.FractionalMRTransform.create(repeat_units=repeat_units, fractional_values=1, default_list=False)
cluster_data = transforms.ClusterTransform.create(mz_tol=mz_tol, ppm_tol=ppm_tol, min_samples=min_samples)
cluster_filter = transforms.FilterByClusterSize.create(min_samples=min_samples)

pipeline = calc_fmr + cluster_data + cluster_filter

fmr_clusters = pipeline(data)
fmr_results = fmr_clusters.split_by_condition_name('file_name')

for result in fmr_results:
    for measurement in result.measurements:
        if measurement.datasets[0].number_of_points() == 0:
            continue
        name = os.path.basename(measurement.conditions['file_name'])
        ru_name = measurement.details['repeat_unit_information'][0]
        df = pd.DataFrame(np.array(measurement.datasets[0].data_arrays).T,
                          columns=measurement.datasets[0].data_array_names)
        df.to_csv(os.path.join(result_path, f'result_{name}_{ru_name}.csv'))
        fig, _ = measurement.visualize()
        fig.savefig(os.path.join(result_path, f'plot_{name}_{ru_name}.png'), dpi=1000, bbox_inches='tight')
