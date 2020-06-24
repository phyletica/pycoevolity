#! /usr/bin/env python

import sys
import os
import logging

from munkres import Munkres

from pycoevolity import errors 

_LOG = logging.getLogger(__name__)


def standardize_partition(elements):
    el_map = {}
    next_idx = 0
    partition = []
    values = {}
    for i, el in enumerate(elements):
        if not el in el_map:
            el_map[el] = next_idx
            values[next_idx] = [el]
            next_idx += 1
        partition.append(el_map[el])
    return tuple(partition), values


class Subset(object):
    """
    A set of indices and a corresponding parameter value (float).
    """
    def __init__(self, set_of_indices = set(), value = None):
        self._indices = set()
        for i in set_of_indices:
            self.add_index(i)
        if value is None:
            self._value = value
        else:
            self.value = value
    
    def add_index(self, i):
        if i in self._indices:
            raise Exception("Subset already has index {0}".format(i))
        self._indices.add(i)

    def _get_value(self):
        return self._value

    def _set_value(self, float_value):
        try:
            self._value = float(float_value)
        except ValueError:
            _LOG.error("Invalid subset value: {0!r}".format(float_rate))
            raise

    value = property(_get_value, _set_value)

    def _get_set_of_indices(self):
        return self._indices
    set_of_indices = property(_get_set_of_indices)

    def _get_list_of_indices(self):
        return sorted(self._indices)
    list_of_indices = property(_get_list_of_indices)

    def _get_number_of_elements(self):
        return len(self._indices) 
    number_of_elements = property(_get_number_of_elements)


class SetPartition(object):
    """
    A list of Subset objects.
    """
    def __init__(self, list_of_subset_objects = []):
        self._indices = set()
        self._subsets = []
        for ss in list_of_subset_objects:
            self.add_subset(ss)
        if self._subsets:
            self.check_validity()

    @classmethod
    def get_from_indices(self, partition_indices):
        indices_list, values = standardize_partition(partition_indices)
        indices = sorted(set(indices_list))
        assert indices == list(range(len(indices)))
        subsets = [Subset() for i in indices]
        for el_index, subset_index in enumerate(indices_list):
            subsets[subset_index].add_index(el_index)
        return SetPartition(subsets)

    def add_subset(self, subset_object):
        assert isinstance(subset_object, Subset)
        self._subsets.append(subset_object)
        for i in subset_object.set_of_indices:
            if i in self._indices:
                raise errors.InvalidSetPartitionError("Duplicate element: {}".format(i))
            self._indices.add(i)

    def _get_number_of_subsets(self):
        return len(self._subsets)
    number_of_subsets = property(_get_number_of_subsets)

    def _get_number_of_elements(self):
        n = 0
        for ss in self.subsets:
            n += ss.number_of_elements
        return n
    number_of_elements = property(_get_number_of_elements)

    def _get_subsets(self):
        return self._subsets
    subsets = property(_get_subsets)

    def distance(self, other):
        assert self.number_of_elements == other.number_of_elements
        mat = []
        dim = max(self.number_of_subsets, other.number_of_subsets)
        for x_subset in self.subsets:
            row = [0] * dim
            for i, y_subset in enumerate(other.subsets):
                intersection = x_subset.set_of_indices & y_subset.set_of_indices
                row[i] = len(intersection)
            mat.append(row)
        n_to_add = dim - self.number_of_subsets
        for i in range(n_to_add):
            mat.append([0] * dim)

        cost_matrix = []
        for row in mat:
            cost_row = [sys.maxsize - col for col in row]
            cost_matrix.append(cost_row)

        indexes = Munkres().compute(cost_matrix)

        total = 0
        for row, column in indexes:
            value = mat[row][column]
            total += value
        return self.number_of_elements - total

    def check_validity(self):
        all_indices = []
        for subset in self.subsets:
            for i in subset.set_of_indices:
                all_indices.append(i)
        duplicates = []
        missing = []
        for i in range(0, (max(all_indices) + 1)):
            n = all_indices.count(i)
            if n > 1:
                duplicates.append(i)
            if n < 1:
                missing.append(i)
        if duplicates or missing:
            msg = "ERROR: Invalid SetPartition\n"
            if missing:
                msg += "Missing indices: {0}\n".format(
                        ", ".join(str(x) for x in missing))
            if duplicates:
                msg += "Duplicate indices: {0}\n".format(
                        ", ".join(str(x) for x in duplicates))
            raise errors.InvalidSetPartitionError(msg)


class SetPartitionCollection(object):
    def __init__(self, set_partitions = []):
        self._set_partitons = []
        for p in set_partitions:
            self.add_set_partition(p)

    @classmethod
    def get_from_indices(self, collection_of_partitions):
        spc = SetPartitionCollection()
        for indices in collection_of_partitions:
            spc.add_set_partition(SetPartition.get_from_indices(indices))
        return spc

    def add_set_partition(self, p):
        assert isinstance(p, SetPartition)
        self._set_partitons.append(p)

    def _get_partitions(self):
        return self._set_partitions
    set_partitions = property(_get_partitions)

    def _get_number_of_partitions(self):
        return len(self._set_partitions)
    number_of_partitions = property(_get_number_of_partitions)

    def get_distance_matrix(self):
        lower_triangle = []
        for n, partition in enumerate(self._set_partitions):
            distances = tuple(
                    partition.distance(j) for j in self._set_partitions[:n+1])
            lower_triangle.append(distances)
        return lower_triangle

    def median_sampled_partitions(self):
        distance_matrix = self.get_distance_matrix()
        median_partitions = []
        dim = len(distance_matrix)
        sum_dist = [0] * dim
        for i in range(dim):
            for j in range(i):
                element = distance_matrix[i][j]
                sum_dist[i] += element
                sum_dist[j] += element
        min_dist = min(sum_dist)
        for n, dist in enumerate(sum_dist):
            if dist == min_dist:
                median_partitions.append(self._set_partitions[n])
        return median_partitions, min_dist

    def distances_from(self, set_partition):
        for p in self.set_partitions:
            yield set_partition.distance(p)
