import random
from collections import namedtuple
import numpy as np


StructEntry = namedtuple('StructEntry', 'index weight')


class Reservoir(object):
    def __init__(self, size):
        self.size = size
        self.contents = []

    def insert(self, struct_entry):
        self.contents.append(struct_entry)
        self.resample()

    def pop_random(self):
        index = random.randrange(self.size)
        struct_entry = self.contents.pop(index)
        self.resample()
        return struct_entry

    def resample(self):
        current_size = len(self.contents)
        weights = [x.weight for x in self.contents]
        total_weight = sum(weights)
        new_weight = total_weight / float(self.size)

        probs = np.array(weights) / total_weight
        sampled_indices = np.random.choice(current_size, self.size, p=probs)

        self.contents = [StructEntry(self.contents[i].index, new_weight) for i in sampled_indices]
