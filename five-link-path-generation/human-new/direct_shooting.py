import casadi as ca
import numpy as np

class walker():
    def __init__(self):
        
        #################################
        # Optimization hyper-parameters #
        #################################
        
        self.T = 1. # Time Horizon
        self.N = 50 # Number of Control Intervals
        self.tau_max = 1. # Max effort
        self.h = self.T/self.N
        gl = np.array([-0.3,0.7,0.0,-0.5,-0.6]).reshape(5, 1)
        self.goal = np.array([gl, gl[::-1]]).reshape(5, 2) # Goal list
        self.p0 = self.p0 = np.array([0., 0.]).reshape(2,1) # Tip of stance leg

        ####################
        # Model parameters #
        ####################
        
        self.m = np.array([0.5,0.5,0.5,0.5,0.5]).reshape(5, 1)
        self.l = np.array([0.5,0.5,0.5,0.5,0.5]).reshape(5, 1)
        self.i = self.m*(self.l**2)/12
        self.gravity = -9.81
        
        ###################
        # Model variables #
        ###################

        q1 = ca.SX.sym('q1')
        q2 = ca.SX.sym('q2')
        q3 = ca.SX.sym('q3')
        q4 = ca.SX.sym('q4')
        q5 = ca.SX.sym('q5')
        self.q = ca.vertcat(q1, q2, q3, q4, q5)
        
        dq1 = ca.SX.sym('dq1')
        dq2 = ca.SX.sym('dq2')
        dq3 = ca.SX.sym('dq3')
        dq4 = ca.SX.sym('dq4')
        dq5 = ca.SX.sym('dq5')
        self.dq = ca.vertcat(dq1, dq2, dq3, dq4, dq5)
        
        self.x = ca.vertcat(self.q, self.dq)

        u1 = ca.SX.sym('u1')
        u2 = ca.SX.sym('u2')
        u3 = ca.SX.sym('u3')
        u4 = ca.SX.sym('u4')
        self.u = ca.vertcat(u1, u2, u3, u4)

        ###################
        # Model equations #
        ###################

        self.p, self.g, self.dg, self.ddq = self.getModel(self.q, self.dq, self.u)
        self.xdot = ca.vertcat(self.dq, self.ddq)
        
        ######################
        # Objective Function #
        ######################

        self.J = self.u**2 

        ####################################
        # Formulate Discrete Time Dynamics #
        ####################################

        self.formulateDTD()
        
        ############################
        # Evaluate at a test point #
        ############################

        print('##############')
        print('# Test Point #')
        print('##############')
        Fk = self.F(x0=[-0.3,0.7,0.0,-0.5,-0.6,0.1,0.1,0.1,0.1,0.1],u=[0.,0.,0.,0.])
        print(Fk['xf'])
        print(Fk['j'])

    def getModel(self, q, dq, u):
        q = ca.reshape(q, 5, 1)
        dq = ca.reshape(dq, 5, 1)
        p0 = ca.SX(self.p0[0, 0])

        p1 = self.l[0, 0]*ca.sin(q[0, 0]) + p0
        p2 = self.l[1, 0]*ca.sin(q[1, 0]) + p1
        p3 = self.l[2, 0]*ca.sin(q[2, 0]) + p2
        p4 = self.l[3, 0]*ca.sin(q[3, 0]) + p2
        p5 = self.l[4, 0]*ca.sin(q[4, 0]) + p4

        c1 = self.l[0, 0]*ca.sin(q[0, 0])/2 + p0
        c2 = self.l[1, 0]*ca.sin(q[1, 0])/2 + p1
        c3 = self.l[2, 0]*ca.sin(q[2, 0])/2 + p2
        c4 = self.l[3, 0]*ca.sin(q[3, 0])/2 + p2
        c5 = self.l[4, 0]*ca.sin(q[4, 0])/2 + p4


        dc1 = self.l[0, 0]*ca.cos(q[0, 0])*dq[0, 0]/2
        dc2 = self.l[1, 0]*ca.cos(q[1, 0])*dq[1, 0]/2 + self.l[0, 0]*ca.cos(q[0, 0])*dq[0, 0]
        dc3 = self.l[2, 0]*ca.cos(q[2, 0])*dq[2, 0]/2 + self.l[1, 0]*ca.cos(q[1, 0])*dq[1, 0] + self.l[0, 0]*ca.cos(q[0, 0])*dq[0, 0]
        dc4 = self.l[3, 0]*ca.cos(q[3, 0])*dq[3, 0]/2 + self.l[1, 0]*ca.cos(q[1, 0])*dq[1, 0] + self.l[0, 0]*ca.cos(q[0, 0])*dq[0, 0]
        dc5 = (self.l[4, 0]*ca.cos(q[4, 0])*dq[4, 0]/2 + self.l[3, 0]*ca.cos(q[3, 0])*dq[3, 0] 
                + self.l[1, 0]*ca.cos(q[1, 0])*dq[1, 0] + self.l[0, 0]*ca.cos(q[0, 0])*dq[0, 0])

        dC = ca.reshape(ca.vertcat(dc1, dc2, dc3, dc4, dc5), 5, 1)

        ddc1 = (self.l[0, 0]*ca.sin(q[0, 0])*(-dq[0, 0]**2)/2) 
        ddc2 = (self.l[1, 0]*ca.sin(q[1, 0])*(-dq[1, 0]**2)/2) + (self.l[0, 0]*ca.sin(q[0, 0])*(-dq[0, 0]**2))
        ddc3 = ((self.l[2, 0]*ca.sin(q[2, 0])*(-dq[2, 0]**2)/2) + (self.l[1, 0]*ca.sin(q[1, 0])*(-dq[1, 0]**2)) 
                + (self.l[0, 0]*ca.sin(q[0, 0])*(-dq[0, 0]**2)))
        ddc4 = ((self.l[3, 0]*ca.sin(q[3, 0])*(-dq[3, 0]**2)/2) + (self.l[1, 0]*ca.sin(q[1, 0])*(-dq[1, 0]**2)) 
                + (self.l[0, 0]*ca.sin(q[0, 0])*(-dq[0, 0]**2)))
        ddc5 = ((self.l[4, 0]*ca.sin(q[0, 0])*(-dq[4, 0]**2)/2) + (self.l[3, 0]*ca.sin(q[3, 0])*(-dq[3, 0]**2)) 
                + (self.l[1, 0]*ca.sin(q[1, 0])*(-dq[1, 0]**2)) + (self.l[0, 0]*ca.sin(q[0, 0])*(-dq[0, 0]**2)))

        P = ca.SX(5, 5)
        P[0, :] = p0
        P[1, 1:4:1] = p1
        P[2, 2:4:1] = p2
        P[3, 3:4:1] = p3
        P[4, 4:4:1] = p4
        # print(P.shape)

        G = ca.SX(5, 5)
        G[0, 0] = c1
        G[0:1:1, 1] = c2
        G[0:2:1, 2] = c3
        G[0:3:1, 3] = c4
        G[:, 4] = c5
        # print(G.shape)

        ddC = ca.reshape(ca.vertcat(ddc1, ddc2, ddc3, ddc4, ddc5), 5, 1)
        # print(ddC.shape)

        U = ca.reshape(ca.vertcat(0., self.u), 5, 1)
        # print(U.shape)

        M = ca.SX(self.m.reshape(5, 1))
        # print(M)

        I = ca.DM([
            [self.i[0], self.i[1], self.i[2], self.i[3], self.i[4]],
            [       0., self.i[1], self.i[2], self.i[3], self.i[4]],
            [       0.,        0., self.i[2], self.i[3], self.i[4]],
            [       0.,        0.,        0., self.i[3], self.i[4]],
            [       0.,        0.,        0.,        0., self.i[4]]
            ])
        # print(I)
        iI = ca.inv(I)
        # print(iI) 
        Q = U + ca.mtimes((G - P), M*(self.gravity - ddC))
        ddq = ca.mtimes(ca.SX(iI), Q)
        # print(ddq)
        return [p0, p1, p2, p3, p4, p5], [c1, c2, c3, c4, c5], dC, ddq   

    def formulateDTD(self):
        # Fixed step Runge-Kutta 4 integrator
        M = 4 # RK4 steps per interval
        dt = self.T/self.N/M
        f = ca.Function('f', [self.x, self.u], [self.xdot, self.J])
        X0 = ca.SX.sym('X0', 10, 1)
        U = ca.SX.sym('U', 4, 1)
        X = X0
        J = 0
        for j in range(M):
            k1, k1_q = f(X, U)
            k2, k2_q = f(X + dt/2 * k1, U)
            k3, k3_q = f(X + dt/2 * k2, U)
            k4, k4_q = f(X + dt * k3, U)
            X = X + dt/6*(k1 +2*k2 +2*k3 +k4)
            J = J + dt/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
        self.F = ca.Function('F', [X0, U], [X, J],['x0','u'],['xf','j'])


biped = walker()
# print(biped.ddq)        