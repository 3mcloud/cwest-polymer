"""FMR Datasets and Mesaurements within piblin framework.

This module provides classes for handling Fractional Mass Remainder (FMR) data within the piblin framework.
It includes FractionalMRDataset for storing mass spectrometry data with FMR calculations, and
FractionalMRMeasurement for managing measurements with specific repeat unit information.


__Classes and Methods__

FractionalMRDataset:
    - remove_clusters:
        Remove specific clusters from the dataset
    - to_dict:
        Convert the dataset to a dictionary format
    - _plot_on_axes:
        Base plotting method inherited from OneDimensionalCompositeDataset. Called through piblin visualize methods

FractionalMRMeasurement:
    - plot_mass_rt:
        Visualize mass vs retention time plots for each identified cluster

Methods are called directly or through piblin Measurement or Dataset objects to analyze and visualize FMR data.
"""

import numpy as np
from typing import Union, List, Dict
import matplotlib
import matplotlib.pyplot as plt
from piblin.data import OneDimensionalCompositeDataset, Measurement
from .. import fmr_parameters as p


class FractionalMRDataset(OneDimensionalCompositeDataset):
    """
    A piblin dataset for storing Fractional Mass Remainder (FMR) data.

    Parameters
    ----------
    mass : List[float] | np.array, optional
        List of mass values.
    mz : List[float] | np.array, optional
        List of mass-to-charge ratio values.
    rt : List[float] | np.array, optional
        List of retention time values.
    x : List[float] | np.array, optional
        List of x-coordinates.
    y : List[float] | np.array, optional
        List of y-coordinates.
    abundance : List[float] | np.array, optional
        List of abundance values.
    mass_priority : bool, default True
        Whether to mass_list is based on mass or mz values.
    """
    def __init__(
            self,
            mass: Union[List[float], np.array] = None,
            mz: Union[List[float], np.array] = None,
            rt: Union[List[float], np.array] = None,
            x: Union[List[float], np.array] = None,
            y: Union[List[float], np.array] = None,
            abundance: Union[List[float], np.array] = None,
            mass_priority: bool = True,
    ):
        if mass is None:
            mass = []
        if mz is None:
            mz = []
        if rt is None:
            rt = []
        if x is None:
            x = []
        if y is None:
            y = []
        if abundance is None:
            abundance = []

        # check that mass or mz values exist
        if max(len(mass), len(mz)) == 0:
            raise ValueError("Missing any mass data")

        if (len(mass) > 0) and (mass_priority is True):
            mass_list = mass
        else:
            mass_list = mz

        if len(mass_list) == 0:
            raise ValueError("if mass values are not priority, then provide mz values")

        n_pts = len(mass_list)

        data_arrays = []
        data_array_names = []
        data_array_units = []

        data_list = [mass_list, mass, mz, rt, x, y, abundance, [0] * n_pts, [np.nan] * n_pts, [np.nan] * n_pts]
        data_dict = dict(zip(p.PARAMETER_LABELS, data_list))

        for name, data in data_dict.items():
            # checks if data is equal in length
            if len(data) == n_pts:
                if isinstance(data, list):
                    data = np.array(data)

                data_arrays.append(data)
                data_array_names.append(name)
                data_array_units.append(p.UNITS_DICT[name])

        super().__init__(
            data_arrays=data_arrays,
            data_array_names=data_array_names,
            data_array_units=data_array_units,
            default_independent_name=p.PARAMETER_LABELS[0],
            default_dependent_name=p.PARAMETER_LABELS[-1],
            source=p.SOURCE
        )

    def remove_clusters(
            self,
            clusters: List[int]
    ):
        """
        Remove specific clusters from the dataset by setting their filter values to 1.
        """
        filter_list = self.__getattribute__(p.FILTER_LABEL)
        cluster_list = self.__getattribute__(p.CLUSTER_LABEL)
        for cluster in clusters:
            filter_list[np.where(cluster_list == cluster)[0]] = 1
        self.__setattr__(p.FILTER_LABEL, filter_list)

    def to_dict(self) -> Dict:
        """
        Convert the dataset to a dictionary format with relevant parameters.
        """
        result = {}
        for prop in p.PARAMETER_LABELS:
            value = getattr(self, prop, None)
            if value is not None:
                result[prop] = value
        return result

    def _plot_on_axes(
                self,
                axes: matplotlib.axes.Axes,
                **axes_plotting_kwargs
    ) -> None:
        """
        Updates plotting method for piblin visualize function.
        """
        # get variables
        filter_list = self.data_arrays[self.data_array_names.index(p.FILTER_LABEL)]
        mass_list = self.data_arrays[self.data_array_names.index(p.MASS_LIST_LABEL)]
        fmr_list = self.data_arrays[self.data_array_names.index(p.FMR_LABEL)]
        if p.ABUNDANCE_LABEL in self.data_array_names:
            abundance = self.data_arrays[self.data_array_names.index(p.ABUNDANCE_LABEL)]
        else:
            abundance = np.ones(len(mass_list))  # if no abundance, use ones

        group = self.data_arrays[self.data_array_names.index(p.CLUSTER_LABEL)]

        # Set up color mapping for the clusters
        cmap_value = axes_plotting_kwargs.get('colormap', 'gist_rainbow')
        cmap = plt.get_cmap(cmap_value)
        unique_groups = np.unique(group[np.where(filter_list == 0)[0]])
        # get size of pts
        group_numeric = [np.where(unique_groups == g)[0] for g in group]
        group_numeric = np.array([x[0] if len(x) > 0 else -1 for x in group_numeric])
        n_size = axes_plotting_kwargs.get('cmap_split', 256)
        group_numeric = np.round(group_numeric /n_size) + np.mod(group_numeric, n_size) * n_size
        # get max point size
        size = axes_plotting_kwargs.get('point_size', 7)

        if abundance is not None:
            # normalized to max abundance of unfiltered peaks
            if any([x == 0 for x in filter_list]):
                max_abund = max(abundance[np.where(filter_list == 0)[0]])
            else:
                max_abund = max(abundance)
            normalized_size = 2 * abundance / max_abund * size
        else:
            normalized_size = size

        filter_labels = axes_plotting_kwargs.get('filter_labels', [0, 1])
        # plot points
        for value in filter_labels:
            # Add scatter plot separately
            mass_ = mass_list[np.where(filter_list == value)[0]]
            fmr_ = fmr_list[np.where(filter_list == value)[0]]
            group_ = group_numeric[np.where(filter_list == value)[0]]
            size_ = normalized_size[np.where(filter_list == value)[0]]
            if value == 0:
                axes.scatter(x=mass_, y=fmr_, s=size_, c=group_, cmap=cmap)
            else:
                axes.scatter(x=mass_, y=fmr_, s=size_, facecolors='none', edgecolors='black', alpha=0.5, linewidth=0.25)

        axes.set_ylim(0, 1)


class FractionalMRMeasurement(Measurement):
    """
    A piblin measurement class for Fractional Mass Remainder (FMR) data.
    """
    def __init__(
            self,
            dataset: FractionalMRDataset,
            repeat_unit_information: tuple[str, float],
            details: dict = None,
            conditions: dict = None,
    ):
        if details is None:
            details = {}
            details[p.DETAIL_RU_LABEL] = repeat_unit_information
        if conditions is None:
            conditions = {}
            conditions[p.CONDITION_RU_LABEL] = repeat_unit_information[1]

        super().__init__(
            datasets=[dataset],
            conditions=conditions,
            details=details,
        )

    def plot_mass_rt(self) -> List[plt.Figure]:
        """Plot mass vs retention time for each group in the dataset.

        Returns
        -------
        List[matplotlib.figure.Figure]
            A list of figure objects, one for each identified cluster in the dataset.
        """
        result = []
        dataset = self.datasets[0]
        ru_label = [f'{x:.04f}' if isinstance(x, float) else x for x in dataset.details[p.DETAIL_RU_LABEL]]
        if p.RT_LABEL not in dataset.data_array_names:
            return result
        groups = dataset.data_arrays[dataset.data_array_names.index(p.CLUSTER_LABEL)]
        masses = dataset.data_arrays[dataset.data_array_names.index(p.MASS_LIST_LABEL)]
        rt = dataset.data_arrays[dataset.data_array_names.index(p.RT_LABEL)]

        for group in np.unique(groups):
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(rt[groups == group], masses[groups == group], s=5)
            ax.set_xlabel('Mass (Da)')
            ax.set_ylabel('RT (min)')
            ax.set_title(f'{ru_label} Group {group}')
            result.append(fig)

        return result
