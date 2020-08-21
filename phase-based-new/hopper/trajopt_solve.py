import numpy as np
import casadi as ca
import math
from matplotlib import pyplot as plt
from matplotlib import animation

from trajopt_formulation import NLP

class TrajOptSolve():
    def __init__(self):
        super().__init__()
        self.formulation = NLP(knot_points_per_phase=40, steps=3, total_duration=1, model='hopper')
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

        for step in range(self.formulation.num_phases):

            for knot_point in range(self.formulation.knot_points_per_phase):
                self.sol_q1.append(sol.value(self.formulation.q[str(step)][knot_point][0]))
                self.sol_q2.append(sol.value(self.formulation.q[str(step)][knot_point][1]))
                self.sol_q3.append(sol.value(self.formulation.q[str(step)][knot_point][2]))

                self.sol_dq1.append(sol.value(self.formulation.qdot[str(step)][knot_point][0]))
                self.sol_dq2.append(sol.value(self.formulation.qdot[str(step)][knot_point][1]))
                self.sol_dq3.append(sol.value(self.formulation.qdot[str(step)][knot_point][2]))

                self.sol_fx.append(sol.value(self.formulation.f10[str(step)][knot_point][0]))
                self.sol_fy.append(sol.value(self.formulation.f10[str(step)][knot_point][1]))

                self.sol_p0x.append(sol.value(self.formulation.p0[str(step)][knot_point][0]))
                self.sol_p0y.append(sol.value(self.formulation.p0[str(step)][knot_point][1]))
        
                self.sol_u1.append(sol.value(self.formulation.u[str(step)][knot_point][0]))
                self.sol_u2.append(sol.value(self.formulation.u[str(step)][knot_point][1]))

        self.time = np.linspace(0.0, self.formulation.total_duration, len(self.sol_q1))

    def plot(self):
        fig = plt.figure()
        fig.tight_layout()

        ax1 = fig.add_subplot(521)
        ax1.plot(self.time, self.sol_q1, label='q1')
        ax1.plot(self.time, self.sol_q2, label='q2')
        ax1.plot(self.time, self.sol_q3, label='q3')
        ax1.legend()

        ax2 = fig.add_subplot(523)
        ax2.plot(self.time, self.sol_dq1, label='dq1')
        ax2.plot(self.time, self.sol_dq2, label='dq2')
        ax2.plot(self.time, self.sol_dq3, label='dq3')
        ax2.legend()

        ax3 = fig.add_subplot(525)
        ax3.plot(self.time, self.sol_u1, label='u1')
        ax3.plot(self.time, self.sol_u2, label='u2')
        ax3.legend()

        ax4 = fig.add_subplot(522)
        ax4.plot(self.time, self.sol_fx, label='fx')
        ax4.plot(self.time, self.sol_fy, label='fy')
        ax4.legend()

        ax5 = fig.add_subplot(524)
        ax5.plot(self.time, self.sol_p0x, label='p0x')
        ax5.plot(self.time, self.sol_p0y, label='p0y')
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

        p1, = ax.plot([],[],'r')
        p2, = ax.plot([],[],'g')
        p3,  = ax.plot([],[],'blue')
        feet, = ax.plot([],[],'yo')
        terrain, = ax.plot([],[],'black')
        force, = ax.plot([],[],'c') 
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        
        def init():
            p1.set_data([], [])
            p2.set_data([], [])
            p3.set_data( [], [])
            feet.set_data([],[])
            force.set_data([],[])
            terrain.set_data([],[])
            return p1, p2, p3, feet, force, terrain

        def animate(i):

            p0 = [self.sol_p0x[i], self.sol_p0y[i]]
            
            link1_x = [p0[0], l[0]*np.sin(self.sol_q1[i]) + p0[0]]
            link1_y = [p0[1], l[0]*np.cos(self.sol_q1[i]) + p0[1]]

            link2_x = [link1_x[1], link1_x[1] + l[1]*np.sin(self.sol_q2[i])]
            link2_y = [link1_y[1], link1_y[1] + l[1]*np.cos(self.sol_q2[i])]

            link3_x = [link2_x[1],link2_x[1] + l[2]*np.sin(self.sol_q3[i])]
            link3_y = [link2_y[1],link2_y[1] + l[2]*np.cos(self.sol_q3[i])]
            force_x = [p0[0], self.sol_fx[i]/math.sqrt(self.sol_fx[i]**2 + self.sol_fy[i]**2) + p0[0]]
            force_y = [p0[1], self.sol_fy[i]/math.sqrt(self.sol_fx[i]**2 + self.sol_fy[i]**2) + p0[1]]
            
            # ax.set_xlim([-2+l_link2_x[1], 2+l_link2_x[1]])
            # ax.set_ylim([-2+l_link2_y[1], 2+l_link2_y[1]])
            
            # terrain.set_data([-2+l_link2_x[1], 2+l_link2_x[1]], [0, 0])

            # ax.set_xlim([-2, 2])
            # ax.set_ylim([-2, 2])
        
            terrain.set_data([-2, 2], [0, 0])


            p1.set_data(link1_x, link1_y)
            p2.set_data(link2_x, link2_y)
            p3.set_data(link3_x, link3_y)
            feet.set_data(p0[0], p0[1])
            force.set_data(force_x,force_y)
            return p1, p2, p3, feet, force, terrain

        ani = animation.FuncAnimation(fig, animate, np.arange(0, len(self.time)), init_func=init,
                               interval=60, blit=True)

        plt.show()

problem = TrajOptSolve()
problem.solve()
problem.plot()
problem.visualize()
