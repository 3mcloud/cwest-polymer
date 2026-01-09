from cwest_polymer.fmr_transforms import fmr_transforms as transforms
from cwest_polymer import fmr_parameters as p
import numpy as np
import pandas as pd


def test_convert_fromula_to_mass():
    formula = "C2H4O"

    mass = transforms.convert_formula_to_mass(formula)
    assert int(mass) == 44


def test_ppm_metric():
    x1 = np.array([75.1, 0.1])
    x2 = np.array([124.9, 0.9])
    ru_value = 50
    mz_tolerance: float = 0.09
    ppm_tolerance: float = 100

    eps = transforms.ppm_metric(
        x1=x1,
        x2=x2,
        repeat_unit=ru_value,
        mz_tolerance=mz_tolerance,
        ppm_tolerance=ppm_tolerance
    )

    assert round(eps) == 100


def test_sort_cluster_by(clusters, abundance):
    # increase cluster 1 abundances to sort them first
    abundance[5:] = abundance[5:]*10
    result = transforms.sort_cluster_by(clusters=clusters, abundance=abundance)
    assert np.all(result == np.array([1]*5+[0]*5))


def test_update_clusters(clusters, mass, repeat_units):
    ru_value = repeat_units['RU1']
    # clustered when size above min_size
    min_size = 3
    result = transforms.update_clusters(clusters=clusters, masses=mass, ru=ru_value, min_size=min_size)
    assert np.all(result == np.array([0]*5+[1]*5))

    # unclustered when size above min_size
    min_size = 6
    result = transforms.update_clusters(clusters=clusters, masses=mass, ru=ru_value, min_size=min_size)
    assert np.all(result == np.array([-1]*10))


def test_fractional_mr_transform(data, repeat_units):
    transform = transforms.FractionalMRTransform.create(repeat_units=repeat_units, default_list=False)
    assert len(transform.data_independent_parameters) == 1

    result = transform(data)
    assert result.measurements[0].details[p.DETAIL_RU_LABEL][1] == 50

    data_arrays = result.measurements[0].datasets[0].data_arrays
    data_array_names = result.measurements[0].datasets[0].data_array_names
    assert np.all(data_arrays[data_array_names.index(p.FILTER_LABEL)] == np.zeros(len(data_arrays[0])))
    assert np.all(np.isnan(data_arrays[data_array_names.index(p.CLUSTER_LABEL)]))


def test_cluster_transform(data, repeat_units):
    transform = transforms.FractionalMRTransform.create(repeat_units=repeat_units, default_list=False)
    transform += transforms.ClusterTransform.create(mz_tol=0.001, ppm_tol=10, min_samples=3)
    assert len(transform.transforms) == 2

    result = transform(data)
    data_arrays = result.measurements[0].datasets[0].data_arrays
    data_array_names = result.measurements[0].datasets[0].data_array_names
    assert np.all(data_arrays[data_array_names.index(p.CLUSTER_LABEL)] == np.array([1]*5+[0]*5))
    assert np.all(data_arrays[data_array_names.index(p.FILTER_LABEL)] == np.zeros(len(data_arrays[0])))

    data_arrays = result.measurements[1].datasets[0].data_arrays
    data_array_names = result.measurements[1].datasets[0].data_array_names
    assert np.all(data_arrays[data_array_names.index(p.CLUSTER_LABEL)] == np.array(p.UNCLUSTERED_LABELS*10))
    assert np.all(data_arrays[data_array_names.index(p.FILTER_LABEL)] == np.ones(len(data_arrays[0])))


def test_filter_by_cluster_size(data, repeat_units):
    transform = transforms.FractionalMRTransform.create(repeat_units=repeat_units, default_list=False)
    transform += transforms.ClusterTransform.create(mz_tol=0.001, ppm_tol=10, min_samples=3)
    transform += transforms.FilterByClusterSize.create(min_samples=4, remove_list=[0])
    assert len(transform.transforms) == 3

    result = transform(data)
    data_arrays = result.measurements[0].datasets[0].data_arrays
    data_array_names = result.measurements[0].datasets[0].data_array_names
    assert np.all(data_arrays[data_array_names.index(p.CLUSTER_LABEL)] == np.array([1]*5+[0]*5))
    assert np.all(data_arrays[data_array_names.index(p.FILTER_LABEL)] == np.array([1]+[0]*(len(data_arrays[0])-1)))

    data_arrays = result.measurements[1].datasets[0].data_arrays
    data_array_names = result.measurements[1].datasets[0].data_array_names
    assert np.all(data_arrays[data_array_names.index(p.CLUSTER_LABEL)] == -np.ones(len(data_arrays[0])))
    assert np.all(data_arrays[data_array_names.index(p.FILTER_LABEL)] == np.ones(len(data_arrays[0])))


def test_calculate_polymer_groups(data, repeat_units):
    transform = transforms.FractionalMRTransform.create(repeat_units=repeat_units, default_list=False)
    transform += transforms.ClusterTransform.create(mz_tol=0.001, ppm_tol=10, min_samples=3)
    transform += transforms.FilterByClusterSize.create(min_samples=4, remove_list=[0])
    transform += transforms.CalculatePolymerGroups.create()
    assert len(transform.transforms) == 4

    result = transform(data)
    df = pd.DataFrame.from_dict(result, orient='index')
    assert df.loc[('test_data', 50, 0), 'end_group'] == 13
    assert df.shape == (2, 15)
