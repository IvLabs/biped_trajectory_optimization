import numpy as np
import casadi as ca
from matplotlib import pyplot as plt
from matplotlib import animation

from trajopt_formulation import NLP1
from trajopt_formulation import NLP2

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP1(knot_points_per_phase=40, steps=2, total_duration=1, model='biped')
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

        self.sol_lpx  = [] 
        self.sol_rpx  = []

        self.sol_lpy  = [] 
        self.sol_rpy  = []

        self.sol_lf  = [] 
        self.sol_rf  = []

        self.sol_u1 = []
        self.sol_u2 = []
        self.sol_u3 = []
        self.sol_u4 = []


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

            self.sol_lf.append(np.linalg.norm(sol.value(self.formulation.lforce[keys])))
            self.sol_rf.append(np.linalg.norm(sol.value(self.formulation.rforce[keys])))

            self.sol_lpx.append(sol.value(self.formulation.lpos[keys][0]))
            self.sol_rpx.append(sol.value(self.formulation.rpos[keys][0]))
       
            self.sol_lpy.append(sol.value(self.formulation.lpos[keys][1]))
            self.sol_rpy.append(sol.value(self.formulation.rpos[keys][1]))

            self.sol_u1.append(sol.value(self.formulation.u[keys][0]))
            self.sol_u2.append(sol.value(self.formulation.u[keys][1]))
            self.sol_u3.append(sol.value(self.formulation.u[keys][2]))
            self.sol_u4.append(sol.value(self.formulation.u[keys][3]))

        self.time = np.linspace(0.0, self.formulation.total_duration, len(self.sol_lq1))
        
    def plot(self):
        fig = plt.figure()

        ax1 = fig.add_subplot(521, autoscale_on=False)
        ax1.plot(self.time, self.sol_lq1, label='lq1')
        ax1.plot(self.time, self.sol_lq2, label='lq2')
        ax1.plot(self.time, self.sol_rq1, label='rq1')
        ax1.plot(self.time, self.sol_rq2, label='rq2')
        ax1.plot(self.time, self.sol_tq , label='tq ')
        ax1.legend()

        ax2 = fig.add_subplot(523, autoscale_on=False)
        ax2.plot(self.time, self.sol_dlq1, label='dlq1')
        ax2.plot(self.time, self.sol_dlq2, label='dlq2')
        ax2.plot(self.time, self.sol_drq1, label='drq1')
        ax2.plot(self.time, self.sol_drq2, label='drq2')
        ax2.plot(self.time, self.sol_dtq , label='dtq ')
        ax2.legend()

        ax3 = fig.add_subplot(525, autoscale_on=False)
        ax3.plot(self.time, self.sol_u1, label='u1')
        ax3.plot(self.time, self.sol_u2, label='u2')
        ax3.plot(self.time, self.sol_u3, label='u3')
        ax3.plot(self.time, self.sol_u4, label='u4')
        ax3.legend()

        ax4 = fig.add_subplot(522, autoscale_on=False)
        ax4.plot(self.time, self.sol_lf, label='lf')
        ax4.plot(self.time, self.sol_rf, label='rf')
        ax4.legend()

        ax5 = fig.add_subplot(524, autoscale_on=False)
        ax5.plot(self.time, self.sol_lpx, label='lpx')
        ax5.plot(self.time, self.sol_lpy, label='lpy')
        ax5.plot(self.time, self.sol_rpx, label='rpx')
        ax5.plot(self.time, self.sol_rpy, label='rpy')
        ax5.legend()

        plt.show()

    def visualize(self):
        l   = self.formulation.model.length
        m   = self.formulation.model.mass
        
        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=False)
        ax.grid()

        lp1 = ax.plot([],[],'r')
        lp2 = ax.plot([],[],'g')
        tp  = ax.plot([],[],'b')
        rp1 = ax.plot([],[],'c')
        rp2 = ax.plot([],[],'y')

        def init():
            lp1.set_data([], [])
            lp2.set_data([], [])
            tp.set_data( [], [])
            rp1.set_data([], [])
            rp2.set_data([], [])
            return lp1, lp2, tp, rp1, rp2

        def animate(i):
            lp0 = np.array([self.sol_lpx[i], self.sol_lpy[i]])
            rp0 = np.array([self.sol_rpx[i], self.sol_rpy[i]])

            l_link1_x = [lp0[0], l[0]*np.sin(self.sol_lq1[i]) + lp0[0]]
            l_link1_y = [lp0[1], l[0]*np.cos(self.sol_lq1[i]) + lp0[1]]

            l_link2_x = [l_link1_x[1], l[1]*np.sin(self.sol_lq2[i])]
            l_link2_y = [l_link1_y[1], l[1]*np.cos(self.sol_lq2[i])]

            r_link1_x = [rp0[0], l[4]*np.sin(self.sol_rq1[i]) + rp0[0]]
            r_link1_y = [rp0[1], l[4]*np.cos(self.sol_rq1[i]) + rp0[1]]

            r_link2_x = [r_link1_x[1], l[3]*np.sin(self.sol_rq2[i])]
            r_link2_y = [r_link1_y[1], l[3]*np.cos(self.sol_rq2[i])]

            t_link_x = [l_link2_x[1], l[2]*np.sin(self.sol_tq[i])]
            t_link_y = [l_link2_y[1], l[2]*np.cos(self.sol_tq[i])]

            lp1.set_data(l_link1_x, l_link1_y)
            lp2.set_data(l_link2_x, l_link2_y)
            tp.set_data(  t_link_x,  t_link_y)
            rp1.set_data(r_link1_x, r_link1_y)
            rp2.set_data(r_link2_x, r_link2_y)

            return lp1, lp2, tp, rp1, rp2

        ani = animation.FuncAnimation(fig, animate, np.arange(0, len(self.time)), init_func=init,
                               interval=20, blit=True)

        plt.show()

problem = TrajOptSolve()
problem.solve()
problem.plot()
problem.visualize()
