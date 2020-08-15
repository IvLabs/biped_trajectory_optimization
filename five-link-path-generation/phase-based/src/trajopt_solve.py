import numpy as np
import casadi as ca

from trajopt_formulation import NLP

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        self.formulation.opti.solver("ipopt",p_opts,s_opts)
    
problem = TrajOptSolve()
problem.formulation.opti.solve_limited()
