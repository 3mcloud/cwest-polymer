"""FMR data transformation module.

This module provides transformers for processing Fractional Mass Remainder (FMR) data.

Utility and Distance functions:
    convert_formula_to_mass:
        Convert chemical formulas to masses using molmass package.

    ppm_metric:
        Calculate PPM-based metrics for mass comparisons. This approach uses the average mass of two values to calculate
        the error based on PPM and m/z tolerances.

    sort_cluster_by:
        Sort clusters by various criteria

    update_clusters:
        Update cluster order based on mass values and minimum size criteria

Piblin Transformations:
    FractionalMRTransform:
        Transform measurement data based on FMR values and repeat unit calculations.

    ClusterTransform:
        Cluster mass spectrometry data based on FMR values. DBsCAN is used for clustering, calling the custom PPM
        metric.

    CalculatePolymerGroups:
        Calculate polymer group statistics for each cluster in FMR datasets.

    FilterByClusterSize:
        Filter clusters based on cluster size and abundance
"""

from sklearn.cluster import DBSCAN
from typing import List, Union, Dict
from molmass import Formula
import numpy as np
from fnmatch import fnmatch
from piblin.data import Measurement, MeasurementSet, Dataset
from piblin.transform import MeasurementSetTransform, DatasetTransform
from .. import fmr_parameters as p
from ..fmr_parameters import DEFAULT_REPEAT_UNITS
from ..fmr_classes.fmr_datasets import FractionalMRDataset as FmrDataset


### Utility and Distance functions ###


def convert_formula_to_mass(formula: str):
    """
    Function used by checks valid formula based on molmass package and returns monoisotopic value
    """
    try:
        return float(Formula(formula).monoisotopic_mass)
    except Exception as e:
        raise ValueError(f"Formula {formula} could not be converted to mass: {e}")


def ppm_metric(
        x1: np.array,
        x2: np.array,
        repeat_unit: float,
        mz_tolerance: float = 0,
        ppm_tolerance: float = 5,
        circular: bool = True,
):
    """
    calculates distance between two np.arrays of vectors (mass, FMR) points based on PPM and m/z tolerances.

    Parameters
    ----------
    x1 : np.array
        First vector (mass, FMR)
    x2 : np.array
        Second vector (mass, FMR)
    repeat_unit : float
        The repeat unit value used for calculating the error.
    mz_tolerance : float
        m/z tolerance value in Da (default is 0 mz units)
    ppm_tolerance : float
        PPM tolerance value (default is 5 ppm)
    circular : bool
        If True, the FMR values are treated as a circular variable (default is True).
    """
    c1 = ppm_tolerance * 10 ** -6
    c2 = mz_tolerance
    avg_mass = (x1[0] + x2[0]) / 2
    error = (c1 * avg_mass + c2) / repeat_unit

    if circular:
        result = min(abs(x1[1] - x2[1]), abs(1 - abs(x1[1] - x2[1]))) / error
    else:
        result = abs(x1[1] - x2[1]) / error
    return result


def sort_cluster_by(
        clusters: np.array,
        abundance: np.array = None,
):
    """
    Sort cluster labels based on abundance values, if provided. Otherwise, sort by the count of features in each cluster.
    """
    if abundance is None:
        abundance = np.ones(len(clusters))

    unique_clusters = [x for x in list(set(clusters)) if x not in p.UNCLUSTERED_LABELS]
    abundance_array = np.array([unique_clusters, np.zeros(len(unique_clusters)), np.zeros(len(unique_clusters))]).T
    for cluster in unique_clusters:
        count = abundance[clusters == cluster].sum()
        abundance_array[abundance_array[:, 0] == cluster, 1] = count
    abundance_array = abundance_array[abundance_array[:, 1].argsort()][::-1]

    abundance_array[:, 2] = [x for x in range(len(unique_clusters))]
    new_clusters = -1 * np.ones(len(clusters))
    for cluster in unique_clusters:
        new_cluster = abundance_array[np.where(abundance_array[:, 0] == cluster), 2]
        new_clusters[np.where(clusters == cluster)] = new_cluster
    return new_clusters


def update_clusters(
        clusters: np.array,
        masses: np.array,
        min_size: int,
        ru: float = None
):
    """
    Update clusters based on number of features in a group.
    """
    if ru is None:
        km_values = masses
    else:
        km_values = masses * round(ru) / ru

    new_clusters = clusters
    unique_clusters = [x for x in list(set(clusters)) if x not in p.UNCLUSTERED_LABELS]
    for cluster in unique_clusters:
        if not np.isnan(cluster):
            cluster_km = km_values[clusters == cluster]
            size1 = len(set(np.round(cluster_km)))
            size2 = len(set(np.floor(cluster_km)))

            if min(size1, size2) < min_size:
                new_clusters[clusters == cluster] = -1

    return new_clusters


class FractionalMRTransform(MeasurementSetTransform):
    """
    Transform based repeat unit to generate fmr calculations within fmr datasets.
    """
    def __init__(self, data_independent_parameters: List[object] = None, *args, **kwargs):
        req = 1
        if len(data_independent_parameters) != req:
            raise ValueError(
                f"Incorrect number of data-independent parameter passed to transform (needs {req}): {len(data_independent_parameters)} given"
            )

        self._repeat_units_values: Dict[str, float] = data_independent_parameters[0]
        super().__init__(data_independent_parameters, *args, **kwargs)

    @staticmethod
    def create(
            repeat_units: Union[str, float, List, Dict[str, Union[str, float]]] = None,
            fractional_values: Union[int, List[int]] = 1,
            default_list: bool = True,
            kmd: bool = False,
    ):
        if repeat_units is None:
            repeat_units = []

        # check repeat unit values
        ru_values = {}
        if isinstance(repeat_units, str) or (isinstance(repeat_units, float)):
            repeat_units = [repeat_units]

        if isinstance(repeat_units, list):
            for ru in repeat_units:
                if isinstance(ru, str):
                    ru_str = ru
                    ru: float = convert_formula_to_mass(ru)
                else:
                    ru_str = f"{p.RU_LABEL}{len(ru_values)}"
                ru_values[ru_str] = ru

        elif isinstance(repeat_units, dict):
            # update values to floats only from formulas
            repeat_units = {x: convert_formula_to_mass(y) if isinstance(y, str) else y for x, y in repeat_units.items()}
            ru_values.update(repeat_units)

        if default_list is True:
            default_values: Dict[str, float] = {x: convert_formula_to_mass(x) for x in DEFAULT_REPEAT_UNITS.values()}
            ru_values.update(default_values)

        # check that fractional values are valid
        if isinstance(fractional_values, int):
            fractional_values = [fractional_values]

        fractional_values = list(set([int(x) for x in fractional_values if x > 0]))  # unique values greater than 0

        if len(fractional_values) == 0:
            fractional_values = [1]

        # adds fractional values of each repeat unit to the final list
        final_rus = {}
        temp_dict = ru_values.copy()
        # add fractional values - default 1
        for f_val in fractional_values:
            final_rus.update({f"{k}_{f_val}": v / int(f_val) for k, v in temp_dict.items()})

        # add kmd values if true
        if kmd is True:
            final_rus.update({f"{k}_k": v / round(v) for k, v in temp_dict.items()})

        return FractionalMRTransform(data_independent_parameters=[final_rus])

    def _apply(self, target: MeasurementSet, **kwargs):
        fmr_measurements = []
        for measurement in target.measurements:
            for dataset in measurement.datasets:
                # check if this is an fMR dataset
                if not isinstance(dataset, FmrDataset):
                    raise ValueError("This transform requires a FractionalMRDataset")

                else:
                    for ru_str, ru_value in self._repeat_units_values.items():
                        properties = dataset.to_dict()
                        properties.pop(p.CMPD_LABEL)
                        properties.pop(p.MASS_LIST_LABEL)
                        properties.pop(p.CLUSTER_LABEL)
                        properties.pop(p.FMR_LABEL)
                        properties.pop(p.FILTER_LABEL)

                        new_dataset = FmrDataset(
                            **properties
                        )

                        masses = dataset.data_arrays[dataset.data_array_names.index(p.MASS_LIST_LABEL)]
                        fmr_values = np.array([(x / ru_value) % 1 for x in masses])

                        new_dataset.data_arrays[dataset.data_array_names.index(p.FMR_LABEL)] = fmr_values

                        details = measurement.details.copy()
                        details[p.DETAIL_RU_LABEL] = (ru_str, ru_value)

                        conditions = measurement.conditions.copy()
                        conditions[p.CONDITION_RU_LABEL] = ru_value

                        measurement = Measurement(
                            datasets=[new_dataset],
                            details=details,
                            conditions=conditions
                        )

                        fmr_measurements.append(measurement)

        return MeasurementSet(measurements=fmr_measurements, merge_redundant=False)


class ClusterTransform(MeasurementSetTransform):
    """
    Transform data into clusters based on DBSCAN cluster algorithm using FMR values and PPM/mz tolerances.
    """
    def __init__(self, data_independent_parameters: List[object] = None, *args, **kwargs):
        req = 4
        if len(data_independent_parameters) != req:
            raise ValueError(
                f"Incorrect number of data-independent parameter passed to transform (needs {req}): {len(data_independent_parameters)} given"
            )
        self._align_params = {
            "mz_tolerance": data_independent_parameters[0],
            "ppm_tolerance": data_independent_parameters[1]
        }
        self._min_samples = data_independent_parameters[2]
        self._eps = data_independent_parameters[3]
        super().__init__(data_independent_parameters, *args, **kwargs)

    @staticmethod
    def create(mz_tol: float, ppm_tol: float, min_samples: int, eps: float = 1):
        return ClusterTransform(data_independent_parameters=[mz_tol, ppm_tol, min_samples, eps])

    def _apply(self, target: MeasurementSet, **kwargs):
        cluster_measurements = []
        for n, measurement in enumerate(target.measurements):
            repeat_unit = measurement.conditions['repeat_unit_value']
            for m, dataset in enumerate(measurement.datasets):
                ru = measurement.conditions[p.CONDITION_RU_LABEL]
                if not isinstance(dataset, FmrDataset):
                    continue

                data_arrays = dataset.data_arrays
                data_array_names = dataset.data_array_names
                mass_values = data_arrays[data_array_names.index(p.MASS_LABEL)]
                fmr_values = data_arrays[data_array_names.index(p.FMR_LABEL)]
                abundance = data_arrays[data_array_names.index(p.ABUNDANCE_LABEL)]

                x = np.array([mass_values, fmr_values]).T

                sorted_index = np.argsort(x[:, 1])
                revert_sorted_index = np.argsort(sorted_index)
                x = x[sorted_index]

                params = self._align_params.copy()
                params['repeat_unit'] = repeat_unit
                dbscan = DBSCAN(eps=self._eps, min_samples=self._min_samples, metric=ppm_metric,
                                metric_params=params)
                dbscan.fit(x)

                labels = dbscan.labels_
                cluster = labels

                # revert sort on cluster and X
                cluster = cluster[revert_sorted_index]
                cluster = update_clusters(clusters=cluster, masses=mass_values, min_size=self._min_samples, ru=ru)
                cluster = sort_cluster_by(cluster, abundance)

                # update clusters and filter unclustered points by default
                dataset.data_arrays[dataset.data_array_names.index(p.CLUSTER_LABEL)] = cluster
                dataset.remove_clusters()

                details = measurement.details
                details['alignment parameters'] = self._align_params

                measurement = Measurement(datasets=[dataset], details=details, conditions=measurement.conditions)
                cluster_measurements.append(measurement)

        return MeasurementSet(measurements=cluster_measurements, merge_redundant=False)


class FilterByClusterSize(DatasetTransform):
    """
    Filter clusters based on minimum cluster size and remove list. Recommended to use cluster size OR specified remove list.
    """
    def __init__(self, data_independent_parameters: List[object] = None, *args, **kwargs):
        req = 2
        if len(data_independent_parameters) != req:
            raise ValueError(
                f"Incorrect number of data-independent parameter passed to transform (needs {req}): {len(data_independent_parameters)} given"
            )
        self._min_samples = data_independent_parameters[0]
        self._remove_list = data_independent_parameters[1]
        if self._remove_list is None:
            self._remove_list = []
        self._remove_list = np.array(self._remove_list)

        super().__init__(data_independent_parameters, *args, **kwargs)

    @staticmethod
    def create(min_samples: int, remove_list: List[int] = None):
        return FilterByClusterSize(data_independent_parameters=[min_samples, remove_list])

    def _apply(self, target: Dataset, **kwargs):
        # check for fmr dataset
        if not isinstance(target, FmrDataset):
            return target

        # get data arrays and unique clusters
        filter_list = target.data_arrays[target.data_array_names.index(p.FILTER_LABEL)]
        clusters = target.data_arrays[target.data_array_names.index(p.CLUSTER_LABEL)]
        unique_clusters = [x for x in np.unique(clusters) if x not in p.UNCLUSTERED_LABELS]

        # remove based on mass indices
        if self._remove_list is not None:
            filter_list[self._remove_list] = 1

        # remove based on cluster size
        for cluster in unique_clusters:
            if len(np.where(clusters == cluster)[0]) < self._min_samples:
                filter_list[np.where(clusters == cluster)[0]] = 1

        target.data_arrays[target.data_array_names.index(p.FILTER_LABEL)] = filter_list
        return target


class CalculatePolymerGroups(MeasurementSetTransform):
    """
    Calculate polymer group statistics for each cluster in FMR datasets, including Mn, Mw, Dispersion, End-group, etc.
    """
    def __init__(self, data_independent_parameters: List[object] = None, *args, **kwargs):
        req = 1
        if len(data_independent_parameters) != req:
            raise ValueError(
                f"Incorrect number of data-independent parameter passed to transform (needs {req}): {len(data_independent_parameters)} given"
            )
        self._ignore_list = data_independent_parameters[0]

        super().__init__(data_independent_parameters, *args, **kwargs)

    @staticmethod
    def create(ignore_list: List[str | float] = None):
        return CalculatePolymerGroups(data_independent_parameters=[ignore_list])

    def _apply(self, target: MeasurementSet, **kwargs) -> Dict[int, Dict[str, Union[str, float]]]:
        # track results and file names
        calculated_values = {}
        file_names = []

        # calculate polymer groups for each repeat unit i.e. measurement
        for measurement in target.measurements:
            temp_values = {}
            temp_values.update(measurement.conditions)

            # get file name and track unique names
            file_name = measurement.conditions.get('file_name', None)
            if file_name is None:
                file_name = measurement.details.get('source_filename', None)
            if file_name is None:
                file_name = f'file_{len(file_names)}'
            file_names.append(file_name)
            file_names = list(set(file_names))

            ru_label, ru_value = measurement.details.get(p.DETAIL_RU_LABEL, (None, None))

            temp_values.update({
                'repeat_unit_label': ru_label,
            })

            if ru_value is None:
                measurement.details.get(p.CONDITION_RU_LABEL, None)

            if ru_value is None:
                print(f"This transform requires a repeat unit {file_name}")
                continue

            # skip if in ignore list
            if self._ignore_list is not None:
                if any(fnmatch(ru_label, f"{x}*") for x in self._ignore_list) or (ru_value in self._ignore_list):
                    print(f"Skipping calculation of polymer groups for {ru_label}: {file_name}")
                    continue

            # find FMR datasets in the measurement - should be 1:1
            for dataset in measurement.datasets:
                # handle incorrect datasets
                if not isinstance(dataset, FmrDataset):
                    print(f"This transform requires a FractionalMRDataset {file_name}")
                    continue
                # Retrieve data arrays
                data_arrays = dataset.data_arrays
                data_array_names = dataset.data_array_names
                mass_list = data_arrays[data_array_names.index(p.MASS_LIST_LABEL)]
                rt_list = data_arrays[data_array_names.index(p.RT_LABEL)]
                abundance = data_arrays[data_array_names.index(p.ABUNDANCE_LABEL)]
                clusters = data_arrays[data_array_names.index(p.CLUSTER_LABEL)]
                # handle none values
                if abundance is None:
                    abundance = np.ones(len(mass_list))
                if rt_list is None:
                    rt_list = np.array([np.nan] * len(mass_list))

                unique_clusters = np.array([x for x in np.unique(clusters) if x not in p.UNCLUSTERED_LABELS])

                # calculate each polymer group
                for cluster in unique_clusters:
                    # select points in the specific group
                    cmpd_idx = np.where(clusters == cluster)[0]
                    cluster_mass = mass_list[cmpd_idx]
                    cluster_rt = rt_list[cmpd_idx]
                    cluster_abundance = abundance[cmpd_idx]

                    # calculate mp, mn, mw...
                    mp = cluster_mass[np.argmax(cluster_abundance)]

                    n = cluster_abundance / sum(cluster_abundance)
                    mn = (cluster_mass * n).sum()
                    mw = (cluster_mass ** 2 * n).sum() / mn
                    mz = (cluster_mass ** 3 * n).sum() / (mn * mw)

                    poly_disp = mw / mn

                    # calculate end-group, mass_min, mass_max, rt_min, rt_max
                    end_groups = np.mod(cluster_mass, ru_value)
                    end_group = end_groups.mean()
                    end_group_sd = end_groups.std()
                    mass_min = min(cluster_mass)
                    mass_max = max(cluster_mass)
                    rt_min = min(cluster_rt)
                    rt_max = max(cluster_rt)

                    # add to temp values
                    temp_values.update({
                        'n': len(cmpd_idx),
                        'mp': mp,
                        'mn': mn,
                        'mw': mw,
                        'mz': mz,
                        'pd': poly_disp,
                        'end_group': end_group,
                        'end_group_sd': end_group_sd,
                        'mass_min': mass_min,
                        'mass_max': mass_max,
                        'rt_min': rt_min,
                        'rt_max': rt_max,
                    })
                    # add to final values
                    calculated_values[(file_name, ru_value, int(cluster))] = temp_values.copy()

        # return dict of values
        return calculated_values
