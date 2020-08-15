import numpy as np
import casadi as ca
from matplotlib import pyplot as plt
from trajopt_formulation import NLP

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        self.formulation.opti.solver("ipopt",p_opts,s_opts)

    def solve(self):
        sol = self.formulation.opti.solve_limited()
        self.sol_q = sol.value(self.formulation.q['2_2'])
        self.sol_qdot = sol.value(self.formulation.qdot['2_2'])

    def plot(self):
        print(self.sol_q)
        print(self.sol_qdot)

problem = TrajOptSolve()
problem.solve()
problem.plot()
