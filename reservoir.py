import math
import random
from collections import namedtuple
import numpy as np


StructEntry = namedtuple('StructEntry', 'index weight')


class Reservoir(object):
    def __init__(self, size):
        self._size = size
        self._contents = []

    def insert(self, struct_entry):
        self._contents.append(struct_entry)
        self._resample()

    def pop_random(self):
        index = random.randrange(self._size)
        struct_entry = self._contents.pop(index)
        self._resample()
        return struct_entry

    def _resample(self):
        weights = [x.weight for x in self._contents]
        total_weight = sum(weights)
        new_weight = total_weight / float(self._size)

        sampled_indices = _systematic_resample(weights, self._size)

        self._contents = [StructEntry(self._contents[i].index, new_weight) for i in sampled_indices]


class HighMarkReservoir(object):
    def __init__(self, target_size, low_mark, high_mark, max_ratio):
        # warning
        # this might not be correct
        #
        raise RuntimeError('Not sure if this works yet.)

        self._target_size = int(target_size)
        self._low_mark = int(low_mark)
        self._high_mark = int(high_mark)
        self._max_ratio = float(max_ratio)

        if not self._target_size > 0:
            raise ValueError('target_size must be > 0')
        if self._low_mark <= 0 or self._low_mark > self._target_size:
            raise ValueError('low_mark must be >0 and <= target size')
        if self._high_mark < self._target_size:
            raise ValueError('low_mark must be >= target_size')
        if not self._max_ratio > 0:
            raise ValueError('max_ration must be > 0')

        self._contents = []
        self._average_weight = 0.0
        self._count = 0

    def insert(self, struct_entry):
        self._update_average_weight(struct_entry.weight)
        n_copies = self._get_n_copies(struct_entry.weight)

        for i in range(n_copies):
            self._contents.append(
                StructEntry(index=struct_entry.index, weight=self._average_weight))

        self._resample_if_needed()

    def pop_random(self):
        index = random.randrange(len(self._contents))
        struct_entry = self._contents.pop(index)
        self._resample_if_needed()
        return struct_entry

    def _update_average_weight(self, weight):
        self._count += 1
        p = float(self._count - 1) / self._count
        self._average_weight = p * self._average_weight + (1.0 - p) * weight

    def _get_n_copies(self, weight):
        ratio = weight / self._average_weight
        lower = math.floor(ratio)
        p = ratio - lower
        if random.random() < p:
            return int(lower)
        else:
            return int(lower + 1)

    def _resample_if_needed(self):
        resample_needed = False
        if len(self._contents) <= self._low_mark:
            resample_needed = True
        if len(self._contents) >= self._high_mark:
            resample_needed = True

        # this could be made more efficient by updating
        # min_ and max_ as we insert and delete
        min_ = min(x.weight for x in self._contents)
        max_ = max(x.weight for x in self._contents)
        if max_ / min_ > self._max_ratio:
            resample_needed = True

        if resample_needed:
            self._resample()

    def _resample(self):
        weights = [x.weight for x in self._contents]
        total_weight = sum(weights)
        new_weight = total_weight / float(self._target_size)

        sampled_indices = _systematic_resample(weights, self._target_size)

        self._contents = [StructEntry(self._contents[i].index, new_weight) for i in sampled_indices]


def _systematic_resample(weights, n):
    weights = np.array(weights)
    probs = weights / np.sum(weights)
    output = []

    u = random.random() / n
    j = 0
    sumQ = probs[j]

    for i in range(n):
        while sumQ < u:
            j = j + 1
            sumQ = sumQ + probs[j]
        output.append(j)
        u = u + 1.0 / n

    return output
