import casadi as ca
import numpy as np

class walker():
    def __init__(self):
        # super().__init__()

        #################################
        # Optimization hyper-parameters #
        #################################
        
        self.T = 0.2 # Time Horizon
        self.N = 50 # Number of Control Intervals
        self.tau_max = 1.5 # Max effort
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

        q1 = ca.MX.sym('q1')
        q2 = ca.MX.sym('q2')
        q3 = ca.MX.sym('q3')
        q4 = ca.MX.sym('q4')
        q5 = ca.MX.sym('q5')
        self.q = ca.vertcat(q1, q2, q3, q4, q5)

        dq1 = ca.MX.sym('dq1')
        dq2 = ca.MX.sym('dq2')
        dq3 = ca.MX.sym('dq3')
        dq4 = ca.MX.sym('dq4')
        dq5 = ca.MX.sym('dq5')
        self.dq = ca.vertcat(dq1, dq2, dq3, dq4, dq5)
        
        self.x = ca.vertcat(self.q, self.dq)

        u1 = ca.MX.sym('u1')
        u2 = ca.MX.sym('u2')
        u3 = ca.MX.sym('u3')
        u4 = ca.MX.sym('u4')
        self.u = ca.vertcat(u1, u2, u3, u4)

        ###################
        # Model equations #
        ###################

        self.p, self.g, self.dg, self.ddq = self.getModel(self.q, self.dq, self.u)
        self.xdot = ca.vertcat(self.dq, self.ddq)
        
        ######################
        # Objective Function #
        ######################

        self.J = ca.sumsqr(self.u) 

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
        Fk = self.F(x0=[-0.3,0.7,0.0,-0.5,-0.6,0.1,0.1,0.1,0.1,0.1], u=[0., 0., 0., 0.])
        print(Fk['xf'])
        print(Fk['j'])

    def getModel(self, q, dq, u):
        q = ca.reshape(q, 5, 1)
        dq = ca.reshape(dq, 5, 1)
        p0 = ca.MX(self.p0[0, 0])

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

        P = ca.MX(5, 5)
        P[0, :] = p0
        P[1, 1:4:1] = p1
        P[2, 2:4:1] = p2
        P[3, 3:4:1] = p3
        P[4, 4:4:1] = p4
        # print(P.shape)

        G = ca.MX(5, 5)
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

        M = ca.MX(self.m.reshape(5, 1))
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
        ddq = ca.mtimes(ca.MX(iI), Q)
        # print(ddq)
        return [p0, p1, p2, p3, p4, p5], [c1, c2, c3, c4, c5], dC, ddq   

    def formulateDTD(self):
        # Fixed step Runge-Kutta 4 integrator
        M = 4 # RK4 steps per interval
        dt = self.T/self.N/M
        self.f = ca.Function('f', [self.x, self.u], [self.xdot, self.J])
        X0 = ca.MX.sym('X0', 10, 1)
        U = ca.MX.sym('U', 4, 1)
        X = X0
        Q = 0
        for j in range(M):
            k1, k1_q = self.f(X, U)
            k2, k2_q = self.f(X + dt/2 * k1, U)
            k3, k3_q = self.f(X + dt/2 * k2, U)
            k4, k4_q = self.f(X + dt * k3, U)
            X = X + dt/6*(k1 +2*k2 +2*k3 +k4)
            Q = Q + dt/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
        self.F = ca.Function('F', [X0, U], [X, Q],['x0','u'],['xf','j'])

class nlp(walker):
    def __init__(self, walker):
        # super().__init__()

        # Start with an empty NLP
        self.w=[]
        self.w0 = []
        self.lbw = []
        self.ubw = []
        self.J = 0.
        self.g=[]
        self.lbg = []
        self.ubg = []

        self.opti = ca.Opti()

        #####################
        # Formulate the NLP #
        #####################

        ## choose a any one

        # self.formulateNLP(walker)
        self.formulateOptiStack(walker)

        ########################
        # Create an NLP solver #
        ########################

        # print(self.w.is_symbolic())
        # self.prob = {'f': self.J, 'x': ca.vertcat(*self.w), 'g': ca.vertcat(*self.g)}
        # options = {'ignore_check_vec': True}
        # self.solver = ca.nlpsol('solver', 'ipopt', self.prob, options)

    def formulateNLP(self, walker):
        Xk = ca.MX([-0.3, 0.7, 0.0, -0.5, -0.6, 0., 0., 0., 0., 0.])
        Xk = ca.reshape(Xk, 10, 1)
        self.g += [Xk]
        self.lbg += [-np.pi/2]*10
        self.ubg += [np.pi/2]*10

        for k in range(walker.N):
            if k == walker.N - 1:
                Xk = ca.MX([-0.6, -0.5, 0.0, 0.7, -0.3, 0., 0., 0., 0., 0.])
                Xk = ca.reshape(Xk, 10, 1)
                self.g += [Xk]
                self.lbg += [-np.pi/2]*10
                self.ubg += [np.pi/2]*10

            # New NLP variable for the control
            Uk = ca.MX.sym('U_' + str(k), 4, 1)
            # Uk = ca.MX([0.])
            self.w += [Uk]
            self.lbw += [-walker.tau_max]*4
            self.ubw += [walker.tau_max]*4
            self.w0 += [0.]*4

            # Integrate till the end of the interval
            Fk = walker.F(x0=Xk, u=Uk)
            Xk = Fk['xf']
            self.J = self.J + Fk['j']
            # print(Xk)
            # print(Fk['j'])
            # Add inequality constraint
            # self.g += [Xk]
            # self.lbg += [-np.pi/2]*10
            # self.ubg += [np.pi/2]*10

    def formulateOptiStack(self, walker):
        Xk = ca.MX([-0.3, 0.7, 0.0, -0.5, -0.6, 0., 0., 0., 0., 0.])

        tau_bound = ca.MX([walker.tau_max, walker.tau_max, walker.tau_max, walker.tau_max])

        J = 0

        for k in range(walker.N):
            Uk = self.opti.variable(4, 1)
            self.opti.subject_to(self.opti.bounded(-tau_bound, Uk, tau_bound))
            Fk = walker.F(x0=Xk, u=Uk)
            Xk = Fk['xf']
            J += Fk['j']

            # self.opti.subject_to() 

        self.opti.minimize(J)
        self.opti.solver('ipopt')

# Create biped and NLP 
biped = walker()
problem = nlp(biped)
sol = problem.opti.solve()
# # Solve the NLP
# sol = problem.solver(x0=problem.w0, lbx=problem.lbw, ubx=problem.ubw, lbg=problem.lbg, ubg=problem.ubg)
# w_opt = sol['x']

# # Plot the solution
# u_opt = w_opt

x_opt = [[-0.3, 0.7, 0.0, -0.5, -0.6, 0, 0, 0, 0, 0]]
for k in range(biped.N):
    Fk = biped.F(x0=x_opt[-1], u=u_opt[k])
    x_opt += [Fk['xf'].full()]
x1_opt = [r[0] for r in x_opt]
x2_opt = [r[1] for r in x_opt]
x3_opt = [r[2] for r in x_opt]
x4_opt = [r[3] for r in x_opt]
x5_opt = [r[4] for r in x_opt]

tgrid = [biped.T/(biped.N)*k for k in range(biped.N+1)]
import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.subplot(221)
plt.plot(tgrid, x1_opt, '-')
plt.plot(tgrid, x2_opt, '-')
plt.plot(tgrid, x3_opt, '-')
plt.plot(tgrid, x4_opt, '-')
plt.plot(tgrid, x5_opt, '-')
# plt.legend(['x1','x2','x3','x4','x5'])
plt.xlabel('t')
plt.ylabel('q')

plt.subplot(222)
plt.plot(tgrid, ca.vertcat(ca.DM.nan(1), u_opt[0:u_opt.numel():4]), '-'
        ,tgrid, ca.vertcat(ca.DM.nan(1), u_opt[1:u_opt.numel():4]), '-'
        ,tgrid, ca.vertcat(ca.DM.nan(1), u_opt[2:u_opt.numel():4]), '-'
        ,tgrid, ca.vertcat(ca.DM.nan(1), u_opt[3:u_opt.numel():4]), '-')
plt.xlabel('t')
plt.ylabel('u')
# plt.legend(['u1','u2','u3','u4',])
plt.show()  