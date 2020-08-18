import numpy as np
import casadi as ca
from matplotlib import pyplot as plt
from trajopt_formulation import NLP1
from trajopt_formulation import NLP2

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP1(knot_points_per_phase=20, steps=2, total_duration=1, model='biped')
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        self.formulation.opti.solver("ipopt",p_opts,s_opts)

    def solve(self):
        sol = self.formulation.opti.solve_limited()
        self.sol_lq1  = [] 
        self.sol_dlq1 = []
        
        self.sol_lq2  = [] 
        self.sol_dlq2 = []
        
        self.sol_rq1  = [] 
        self.sol_drq1 = []

        self.sol_rq2  = [] 
        self.sol_drq2 = []

        self.sol_tq  = [] 
        self.sol_dtq = []

        self.sol_lp  = [] 
        self.sol_rp  = []

        self.sol_lf  = [] 
        self.sol_rf  = []


        for keys in self.formulation.lq:
            self.sol_lq1.append(sol.value(self.formulation.lq[keys][0]))
            self.sol_rq1.append(sol.value(self.formulation.rq[keys][0]))

            self.sol_lq2.append(sol.value(self.formulation.lq[keys][1]))
            self.sol_rq2.append(sol.value(self.formulation.rq[keys][1]))

            self.sol_tq.append(sol.value(self.formulation.tq[keys][0]))

            self.sol_dlq1.append(sol.value(self.formulation.lqdot[keys][0]))
            self.sol_drq1.append(sol.value(self.formulation.rqdot[keys][0]))

            self.sol_dlq2.append(sol.value(self.formulation.lqdot[keys][1]))
            self.sol_drq2.append(sol.value(self.formulation.rqdot[keys][1]))

            self.sol_dtq.append(sol.value(self.formulation.tqdot[keys][0]))

            self.sol_lf.append(sol.value(self.formulation.lforce[keys]))
            self.sol_rf.append(sol.value(self.formulation.rforce[keys]))

            self.sol_lp.append(sol.value(self.formulation.lpos[keys]))
            self.sol_rp.append(sol.value(self.formulation.rpos[keys]))


        self.time = np.linspace(0.0, self.formulation.total_duration, len(self.sol_lq1))
        
    def plot(self):

        plt.subplot(211)
        plt.plot(self.time, self.sol_lq1, label='lq1')
        plt.plot(self.time, self.sol_lq2, label='lq2')
        plt.plot(self.time, self.sol_rq1, label='rq1')
        plt.plot(self.time, self.sol_rq2, label='rq2')
        plt.plot(self.time, self.sol_tq , label='tq ')
        plt.legend()

        plt.subplot(212)
        plt.plot(self.time, self.sol_dlq1, label='dlq1')
        plt.plot(self.time, self.sol_dlq2, label='dlq2')
        plt.plot(self.time, self.sol_drq1, label='drq1')
        plt.plot(self.time, self.sol_drq2, label='drq2')
        plt.plot(self.time, self.sol_dtq , label='dtq ')
        plt.legend()

        plt.show()
        # print(self.sol_q)
        # print(self.sol_qdot)

######################################
# Remodel the dynamics according to h#
######################################

problem = TrajOptSolve()
problem.solve()
problem.plot()
