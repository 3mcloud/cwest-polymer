from cwest_polymer.fmr_classes.fmr_datasets import FractionalMRDataset
from cwest_polymer.fmr_classes.fmr_datasets import FractionalMRMeasurement
from cwest_polymer import fmr_parameters as p


def test_fractional_mr_dataset(mass, rt, abundance):
    fmr_data = FractionalMRDataset(
        mass=mass,
        rt=rt,
        abundance=abundance
    )
    assert fmr_data.number_of_points() == len(mass)


def test_fractional_mr_measurement(mass, rt, abundance, clusters):
    dataset = FractionalMRDataset(
        mass=mass,
        rt=rt,
        abundance=abundance
    )

    dataset.data_arrays[dataset.data_array_names.index(p.CLUSTER_LABEL)] = clusters
    dataset.data_arrays[dataset.data_array_names.index(p.FILTER_LABEL)] = [0]*len(mass)

    measurement = FractionalMRMeasurement(
                dataset=dataset,
                repeat_unit_information=('RU1', 50),
                details={p.DETAIL_RU_LABEL: ('RU1', 50)},
                conditions={'file_name': 'test_data'}
    )
    assert measurement.details[p.DETAIL_RU_LABEL][1] == 50
    assert measurement.datasets[0].number_of_points() == len(mass)
