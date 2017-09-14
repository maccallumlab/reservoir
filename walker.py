import random
import math
import numpy as np


class Walker1D(object):
    def __init__(self, start, energy_landscape, steps):
        self.energy_landscape = np.array(energy_landscape)
        self.max_x = self.energy_landscape.shape[0]
        self.steps = steps

        assert start >= 0 and start < self.max_x
        self.x = start

    def update(self):
        for _ in range(self.steps):
            delta = random.choice([-1, 1])
            trial_x = self.x + delta
            if trial_x < 0 or trial_x >= self.max_x:
                continue

            delta_e = self.energy_landscape[trial_x] - self.energy_landscape[self.x]
            if delta_e < 0:
                accept = True
            else:
                metrop = math.exp(-delta_e)
                if random.random() < metrop:
                    accept = True
                else:
                    accept = False
            if accept:
                self.x = trial_x

    def get_energy(self, x):
        return self.energy_landscape[x]
