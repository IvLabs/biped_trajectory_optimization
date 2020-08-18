import numpy as np
import casadi as ca
from matplotlib import pyplot as plt
from trajopt_formulation import NLP1
from trajopt_formulation import NLP2

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP2(knot_points_per_phase=40, steps=2, total_duration=1, model='biped')
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        self.formulation.opti.solver("ipopt",p_opts,s_opts)

    def solve(self):
        sol = self.formulation.opti.solve_limited()
        for keys in self.formulation.lq:
            print(sol.value(self.formulation.lq[keys]))
            print(sol.value(self.formulation.lqdot[keys]))

    # def plot(self):
    #     print(self.sol_q)
    #     print(self.sol_qdot)

######################################
# Remodel the dynamics according to h#
######################################

problem = TrajOptSolve()
problem.solve()
# problem.plot()
