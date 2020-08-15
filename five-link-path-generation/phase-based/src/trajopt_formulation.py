import numpy as np
import casadi as ca
from five_link_model import Biped

class NLP():
    def __init__(self, knot_points_per_phase, steps, total_duration, model):
        super().__init__()
        self.knot_points_per_phase = knot_points_per_phase
        self.num_phases            = steps
        self.total_duration        = total_duration

        if model == 'biped':
            self.model = Biped()

problem = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')            