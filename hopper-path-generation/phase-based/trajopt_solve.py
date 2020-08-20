import numpy as np
import casadi as ca
from matplotlib import pyplot as plt
from matplotlib import animation

from trajopt_formulation import NLP1
from trajopt_formulation import NLP2

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP1(knot_points_per_phase=20, steps=3, total_duration=1, model='hopper')
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        self.formulation.opti.solver("ipopt",p_opts,s_opts)

    def solve(self):
        sol = self.formulation.opti.solve_limited()
        self.sol_q1  = [] 
        self.sol_dq1 = []
        
        self.sol_q2  = [] 
        self.sol_dq2 = []
        
        self.sol_q3  = [] 
        self.sol_dq3 = []

        self.sol_p0x  = [] 
        self.sol_p0y  = [] 

        self.sol_fx  = [] 
        self.sol_fy  = []

        self.sol_u1 = []
        self.sol_u2 = []

        for step in (self.formulation.num_phases):

            for knot_point in range(self.formulation.knot_points_per_phase):
                self.sol_q1.append(sol.value(self.formulation.q[str(step)][knot_point][0]))
                self.sol_q2.append(sol.value(self.formulation.q[str(step)][knot_point][1]))
                self.sol_q3.append(sol.value(self.formulation.q[str(step)][knot_point][2]))

                self.sol_dq1.append(sol.value(self.formulation.dq[str(step)][knot_point][0]))
                self.sol_dq2.append(sol.value(self.formulation.dq[str(step)][knot_point][1]))
                self.sol_dq3.append(sol.value(self.formulation.dq[str(step)][knot_point][2]))

                self.sol_fx.append(sol.value(self.formulation.f10[str(step)][knot_point][0]))
                self.sol_fy.append(sol.value(self.formulation.f10[str(step)][knot_point][1]))

                self.sol_p0x.append(sol.value(self.formulation.p10[str(step)][knot_point][0]))
                self.sol_p0x.append(sol.value(self.formulation.p10[str(step)][knot_point][1]))
        
                self.sol_u1.append(sol.value(self.formulation.u[str(step)][knot_point][0]))
                self.sol_u2.append(sol.value(self.formulation.u[str(step)][knot_point][1]))

        self.time = np.linspace(0.0, self.formulation.total_duration, len(self.sol_q1))

    def plot(self):
        fig = plt.figure()
        fig.tight_layout()

        ax1 = fig.add_subplot(521)
        ax1.plot(self.time, self.sol_lq1, label='lq1')
        ax1.plot(self.time, self.sol_lq2, label='lq2')
        ax1.plot(self.time, self.sol_rq1, label='rq1')
        ax1.plot(self.time, self.sol_rq2, label='rq2')
        ax1.plot(self.time, self.sol_tq , label='tq ')
        ax1.legend()

        ax2 = fig.add_subplot(523)
        ax2.plot(self.time, self.sol_dlq1, label='dlq1')
        ax2.plot(self.time, self.sol_dlq2, label='dlq2')
        ax2.plot(self.time, self.sol_drq1, label='drq1')
        ax2.plot(self.time, self.sol_drq2, label='drq2')
        ax2.plot(self.time, self.sol_dtq , label='dtq ')
        ax2.legend()

        ax3 = fig.add_subplot(525)
        ax3.plot(self.time, self.sol_u1, label='u1')
        ax3.plot(self.time, self.sol_u2, label='u2')
        ax3.plot(self.time, self.sol_u3, label='u3')
        ax3.plot(self.time, self.sol_u4, label='u4')
        ax3.legend()

        ax4 = fig.add_subplot(522)
        ax4.plot(self.time, self.sol_lf, label='lf')
        ax4.plot(self.time, self.sol_rf, label='rf')
        ax4.legend()

        ax5 = fig.add_subplot(524)
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
        ax = fig.add_subplot(111)
        ax.grid()

        # lp0 = [self.sol_lpx[0], self.sol_lpy[0]]
        # rp0 = [self.sol_rpx[0], self.sol_rpy[0]]

        # l_link1_x = [lp0[0], l[0]*np.sin(self.sol_lq1[0]) + lp0[0]]
        # l_link1_y = [lp0[1], l[0]*np.cos(self.sol_lq1[0]) + lp0[1]]

        # print(l_link1_x, l_link1_y)

        lp1, = ax.plot([],[],'r')
        lp2, = ax.plot([],[],'g')
        tp,  = ax.plot([],[],'blue')
        rp1, = ax.plot([],[],'c')
        rp2, = ax.plot([],[],'y')

        terrain, = ax.plot([],[],'black')

        def init():
            lp1.set_data([], [])
            lp2.set_data([], [])
            tp.set_data( [], [])
            rp1.set_data([], [])
            rp2.set_data([], [])
            terrain.set_data([],[])
            return lp1, lp2, tp, rp1, rp2, terrain

        def animate(i):

            lp0 = [self.sol_lpx[i], self.sol_lpy[i]]
            rp0 = [self.sol_rpx[i], self.sol_rpy[i]]

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

            # ax.set_xlim([-2+l_link2_x[1], 2+l_link2_x[1]])
            # ax.set_ylim([-2+l_link2_y[1], 2+l_link2_y[1]])
            
            # terrain.set_data([-2+l_link2_x[1], 2+l_link2_x[1]], [0, 0])

            ax.set_xlim([-2, 2])
            ax.set_ylim([-2, 2])
        
            terrain.set_data([-2, 2], [0, 0])


            lp1.set_data(l_link1_x, l_link1_y)
            lp2.set_data(l_link2_x, l_link2_y)
            tp.set_data(  t_link_x,  t_link_y)
            rp1.set_data(r_link1_x, r_link1_y)
            rp2.set_data(r_link2_x, r_link2_y)

            return lp1, lp2, tp, rp1, rp2, terrain

        ani = animation.FuncAnimation(fig, animate, np.arange(0, len(self.time)), init_func=init,
                               interval=60, blit=True)

        plt.show()

problem = TrajOptSolve()
problem.solve()
problem.plot()
# problem.visualize()
