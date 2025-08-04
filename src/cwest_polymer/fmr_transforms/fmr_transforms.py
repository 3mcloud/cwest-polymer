"""FMR data transformation module.

This module provides transformers for processing Fractional Mass Remainder (FMR) data.
Utility and Distance functions:
    convert_formula_to_mass:
        Convert chemical formulas to masses using molmass package.

    ppm_metric:
        Calculate PPM-based metrics for mass comparisons. This approach uses the average mass of two values to calculate
        the error based on PPM and m/z tolerances.

Calculation and cluster transformations:
    FractionalMRTransform:
        Transform measurement data based on FMR values and repeat unit calculations.

    ClusterTransform:
        Cluster mass spectrometry data based on FMR values. DBsCAN is used for clustering, calling the custom PPM
        metric.

Piblin Utility Transformations:
    FilterByClusterSize:
        Filter clusters based on cluster size and abundance

    sort_cluster_by:
        Sort clusters by various criteria

    update_clusters:
        Update cluster order based on mass values and minimum size criteria
"""

from sklearn.cluster import DBSCAN
from typing import List, Union, Dict
from molmass import Formula
import numpy as np
from piblin.data import Measurement, MeasurementSet, Dataset
from piblin.transform import MeasurementSetTransform, DatasetTransform
from .. import fmr_parameters as p
from ..fmr_parameters import DEFAULT_REPEAT_UNITS
from ..fmr_classes.fmr_datasets import FractionalMRDataset as FmrDataset


def convert_formula_to_mass(formula: str):
    """
    checks valid formula based on molmass package and returns monoisotopic value
    """
    try:
        return float(Formula(formula).monoisotopic_mass)
    except Exception:
        raise ValueError(f"Formula {formula} could not be converted to mass")


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
    Sorts clusters by their abundance or size. If abundance is not provided, it defaults to an array of ones.

    Parameters
    ----------
    clusters : np.array
        Array of cluster labels to be sorted.
    abundance : np.array, optional
        Array of abundance values corresponding to each cluster. Default is None.

    Returns
    -------
    np.array New array of sorted cluster labels.
    """

    if abundance is None:
        abundance = np.ones(len(clusters))

    unique_clusters = list(set(clusters))
    abundance_array = np.transpose([unique_clusters, np.zeros(len(unique_clusters)), np.zeros(len(unique_clusters))])
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
    Updates the cluster labels based on the mass values and minimum size criteria. If the cluster size is less than the
    minimum size, it assigns a new cluster label.

    Parameters
    ----------
    clusters : np.array
        Array of cluster labels to be updated.
    masses : np.array
        Array of mass values corresponding to each cluster.
    min_size : int
        Minimum size threshold for clusters.
    ru : float, optional
        Repeat unit value for mass normalization. Default is None.

    Returns
    -------
    np.array Updated array of cluster labels.
    """
    if ru is None:
        km_values = masses
    else:
        km_values = masses * round(ru) / ru

    new_clusters = clusters
    for cluster in set(clusters):
        if not np.isnan(cluster):
            cluster_km = km_values[clusters == cluster]
            size1 = len(set(np.round(cluster_km)))
            size2 = len(set(np.floor(cluster_km)))

            if min(size1, size2) < min_size:
                new_clusters[clusters == cluster] = -1

    result = np.zeros(shape=clusters.shape, dtype=np.int_)
    c = max(clusters)
    for n, value in enumerate(clusters):
        if value == -1:
            c += 1
            result[n] = c
        else:
            result[n] = value
    return result


class FractionalMRTransform(MeasurementSetTransform):
    """A transform class for calculating Fractional Mass Remainder (FMR) values.

    This transform calculates FMR values for mass spectrometry data by dividing masses by repeat unit values.
    FMR values are calculated as (mass / repeat_unit) % 1, resulting in values between 0 and 1.

    The transform supports:
    - Multiple repeat units specified as formulas or masses
    - Fractional repeat units which are applied to all repeat units (e.g., RU/2, RU/3)
    - Option to add Kendrick Mass Defect (KMD) calculations
    - Default repeat units from package parameters

    """
    def __init__(self, data_independent_parameters: List[object] = None, *args, **kwargs):
        req = 1
        if len(data_independent_parameters) != req:
            raise ValueError(
                f"Incorrect number of data-independent parameter passed to transform (needs {req}): "
                f"{len(data_independent_parameters)} given"
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
        """ Method used to create the piblin FractionalMRTransform, which can be applied FMR MeasurementSet.

        This method effectively builds a comprehensive list of repeat units including fractions and optional Kendrick
        Mass Defect calculations.

        Overview:
        1. Creates an empty repeat unit list if no repeat units provided, accepts input as formula (string), repeat
        unit values (float), a list, or dictionary thereof.
        2. Converts repeat units to mass values (using formulas if given), processes default units if enabled.
        3. Validates fractional values (positive integers only) and creates fractional versions of each repeat unit.
        4. Returns new FractionalMRTransform with processed repeat unit dictionary, including optional KMD calculations.

        Parameters
        ----------
        repeat_units : Union[str, float, List, Dict[str, Union[str, float]]], optional
            Chemical formulas or masses to use as repeat units.
        fractional_values : Union[int, List[int]], default=1
            List of fractional values for fractional repeat units.
        default_list : bool, default=True
            Whether to include default repeat units from package parameters.
        kmd : bool, default=False
            Whether to calculate Kendrick Mass Defect values.
        """
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
            final_rus.update({f"{k}_k": v/round(v) for k, v in temp_dict.items()})

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
                        properties = dict(zip(dataset.data_array_names,dataset.data_arrays))
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
    """A transform class for clustering mass spectrometry data based on FMR values.

    This transform performs density-based clustering (DBSCAN) on mass/FMR pairs using a custom PPM-based metric.
    Clusters are formed based on mass accuracy (PPM tolerance), m/z tolerance, and minimum cluster size criteria.

    The clustering process:
    1. Sorts data points by FMR values
    2. Applies DBSCAN with custom PPM metric
    3. Updates cluster labels based on size criteria
    4. Sorts clusters by abundance

    Input datasets must be a measurementset composed of FractionalMRDataset instances with FMR values calculated.
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
        """ Method used to create the piblin ClusterTransform, which can be applied to FMR MeasurementSet.

        Parameters
        ----------
        mz_tol : float
            m/z tolerance value in Da.
        ppm_tol : float
            PPM tolerance value.
        min_samples : int
            Minimum number of samples required to form a cluster.
        eps : float, default=1
            The maximum relative distance between two samples for one to be considered as in the neighborhood of the
            other.

        """
        return ClusterTransform(data_independent_parameters=[mz_tol, ppm_tol, min_samples, eps])

    def _apply(self, target: MeasurementSet, **kwargs):
        cluster_measurements = []
        for n, measurement in enumerate(target.measurements):
            repeat_unit = measurement.conditions['repeat_unit_value']
            for m, dataset in enumerate(measurement.datasets):
                ru = measurement.conditions[p.CONDITION_RU_LABEL]
                if not isinstance(dataset, FmrDataset):
                    continue
                mass_values = dataset.data_arrays[dataset.data_array_names.index(p.MASS_LIST_LABEL)]
                fmr_values = dataset.data_arrays[dataset.data_array_names.index(p.FMR_LABEL)]
                if p.ABUNDANCE_LABEL in dataset.data_array_names:
                    abundance = dataset.data_arrays[dataset.data_array_names.index(p.ABUNDANCE_LABEL)]
                else:
                    abundance = None

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

                dataset.data_arrays[dataset.data_array_names.index(p.CLUSTER_LABEL)] = cluster

                details = measurement.details
                details['alignment parameters'] = self._align_params

                measurement = Measurement(datasets=[dataset], details=details, conditions=measurement.conditions)
                cluster_measurements.append(measurement)

        return MeasurementSet(measurements=cluster_measurements, merge_redundant=False)


class FilterByClusterSize(DatasetTransform):
    """A transform class for filtering clusters based on size and abundance criteria.

    This transform removes clusters that don't meet minimum size requirements or are in a specified removal list.
    It updates the filter status of points in removed clusters while preserving the original cluster assignments.

    The filtering process:
    1. Checks each cluster's size against minimum sample threshold
    3. Updates filter flags for specified cluster numbers in the FMR datasets

    Input datasets must be FractionalMRDataset instances with cluster assignments.
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
        if not isinstance(target, FmrDataset):
            return target

        filter_list = target.data_arrays[target.data_array_names.index(p.FILTER_LABEL)]
        clusters = target.data_arrays[target.data_array_names.index(p.CLUSTER_LABEL)]

        unique_clusters = np.unique(clusters)

        remove_list = self._remove_list

        for cluster in unique_clusters:
            if (len(np.where(clusters == cluster)[0]) < self._min_samples) and (cluster not in remove_list):
                remove_list = np.concat([remove_list, np.array([cluster])])
                filter_list[np.where(clusters == cluster)[0]] = 1

        target.data_arrays[target.data_array_names.index(p.FILTER_LABEL)] = filter_list
        return target
