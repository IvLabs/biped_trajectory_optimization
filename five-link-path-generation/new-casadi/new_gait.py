import numpy as np
import casadi as ca

class walker():
    def __init__(self, start_angles, start_angular_vel, start_pos):
        # set our parameters of optimization
        self.opti = ca.Opti()
        self.terrain_factor = 0.5
        self.terrain = ['sin','wedge'][1]
        self.N = 80; self.T = .3
        self.step_max = 0.5; self.tauMax = 50
        self.pi = np.pi; 
        self.length = ca.MX([0.5,0.5,0.5,0.5,0.5])
        self.mass = ca.MX([0.25,0.25,0.25,0.25,0.25])
        self.inertia = self.mass * (self.length**2) /12
        self.gravity = 10
        self.h = self.T/self.N
        self.goal = [start_angles, start_angular_vel]
        self.ini_goal = self.goal[0].to_DM()
        self.fin_goal = ca.DM([self.ini_goal[4],self.ini_goal[3],self.ini_goal[2],self.ini_goal[1],self.ini_goal[0]])
        self.p0 = ca.MX(start_pos)     
        self.comh = self.length[0]*0.5
        #set our optimization variables
        self.x = []
        self.xdot = []
        self.u = []
        for i in range(self.N): 
            self.x.append(self.opti.variable(5))
            self.xdot.append(self.opti.variable(5))
            self.u.append(self.opti.variable(4))

        self.pos = [];self.com = [];self.ddq = []
        self.dpos = []

        for n in range(self.N):
            p,dp,ddp,c,dc,ddc = self.getKinematics(self.x[n], self.xdot[n]) 
            self.pos.append(p);self.dpos.append(dp);self.com.append(c)
            ddq = self.getDynamics(self.x[n], self.xdot[n], self.u[n], p, ddp, c, ddc)
            self.ddq.append(ddq)

            if n == self.N - 1:
               self.x_impact, self.xdot_impact = self.impactMap(self.x[n], self.xdot[n], p, dp, c, dc)

    def getDynamics(self, q, dq, u, p, ddp, c, ddc):
        p0 = ca.reshape(self.p0, 1, 2)

        ddq = ca.MX.zeros(5)

        ddq[4] = (u[3] + self.mass[4]*(c[4, 0] - p[3, 0])*self.gravity
                       - self.mass[4]*self.crossProduct2D(c[4, :] - p[3, :], ddc[4, :]) - ddq[4])
        
        ddq[3] = (u[2] + self.mass[4]*(c[4, 0] - p[2, 0])*self.gravity 
                       + self.mass[3]*(c[3, 0] - p[2, 0])*self.gravity
                       - self.mass[4]*self.crossProduct2D(c[4, :] - p[2, :], ddc[4, :])
                       - self.mass[3]*self.crossProduct2D(c[3, :] - p[2, :], ddc[3, :]) - ddq[4])
        
        ddq[2] = (u[1] + self.mass[4]*(c[4, 0] - p[1, 0])*self.gravity 
                       + self.mass[3]*(c[3, 0] - p[1, 0])*self.gravity
                       + self.mass[2]*(c[2, 0] - p[1, 0])*self.gravity 
                       - self.mass[4]*self.crossProduct2D(c[4, :] - p[1, :], ddc[4, :])
                       - self.mass[3]*self.crossProduct2D(c[3, :] - p[1, :], ddc[3, :]) 
                       - self.mass[2]*self.crossProduct2D(c[2, :] - p[1, :], ddc[2, :]) - ddq[3])

        ddq[1] = (u[0] + self.mass[4]*(c[4, 0] - p[0, 0])*self.gravity 
                       + self.mass[3]*(c[3, 0] - p[0, 0])*self.gravity
                       + self.mass[2]*(c[2, 0] - p[0, 0])*self.gravity
                       + self.mass[1]*(c[1, 0] - p[0, 0])*self.gravity 
                       - self.mass[4]*self.crossProduct2D(c[4, :] - p[0, :], ddc[4, :])
                       - self.mass[3]*self.crossProduct2D(c[3, :] - p[0, :], ddc[3, :]) 
                       - self.mass[2]*self.crossProduct2D(c[2, :] - p[0, :], ddc[2, :])
                       - self.mass[1]*self.crossProduct2D(c[1, :] - p[0, :], ddc[1, :]) - ddq[2])

        ddq[0] = (       self.mass[4]*(c[4, 0] - p0[0, 0])*self.gravity 
                       + self.mass[3]*(c[3, 0] - p0[0, 0])*self.gravity
                       + self.mass[2]*(c[2, 0] - p0[0, 0])*self.gravity
                       + self.mass[1]*(c[1, 0] - p0[0, 0])*self.gravity 
                       + self.mass[0]*(c[0, 0] - p0[0, 0])*self.gravity 
                       - self.mass[4]*self.crossProduct2D(c[4, :] - p0[0, :], ddc[4, :])
                       - self.mass[3]*self.crossProduct2D(c[3, :] - p0[0, :], ddc[3, :]) 
                       - self.mass[2]*self.crossProduct2D(c[2, :] - p0[0, :], ddc[2, :])
                       - self.mass[1]*self.crossProduct2D(c[1, :] - p0[0, :], ddc[1, :]) 
                       - self.mass[0]*self.crossProduct2D(c[0, :] - p0[0, :], ddc[2, :]) - ddq[2])

        return ddq/self.inertia

    def getKinematics(self, q, dq):
        p0 = [self.p0[0],self.p0[1]]
        
        p = ca.MX.zeros(5, 2)
        c = ca.MX.zeros(5, 2)
        p0 = [self.p0[0],self.p0[1]]
        
        p[0,0],p[0,1] = self.length[0]*ca.sin(q[0]) + p0[0] , self.length[0]*ca.cos(q[0]) + p0[1]
        p[1,0],p[1,1] = self.length[1]*ca.sin(q[1]) + p[0,0], self.length[1]*ca.cos(q[1]) + p[0,1]
        p[2,0],p[2,1] = self.length[2]*ca.sin(q[2]) + p[1,0], self.length[2]*ca.cos(q[2]) + p[1,1]
        p[3,0],p[3,1] = self.length[3]*ca.sin(q[3]) + p[1,0],-self.length[3]*ca.cos(q[3]) + p[1,1]
        p[4,0],p[4,1] = self.length[4]*ca.sin(q[4]) + p[3,0],-self.length[4]*ca.cos(q[4]) + p[3,1]
        
        c[0,0],c[0,1] = self.length[0]*ca.sin(q[0])/2 + p0[0] , self.length[0]*ca.cos(q[0])/2 + p0[1]
        c[1,0],c[1,1] = self.length[1]*ca.sin(q[1])/2 + p[0,0], self.length[1]*ca.cos(q[1])/2 + p[0,1]
        c[2,0],c[2,1] = self.length[2]*ca.sin(q[2])/2 + p[1,0], self.length[2]*ca.cos(q[2])/2 + p[1,1]
        c[3,0],c[3,1] = self.length[3]*ca.sin(q[3])/2 + p[1,0],-self.length[3]*ca.cos(q[3])/2 + p[1,1]
        c[4,0],c[4,1] = self.length[4]*ca.sin(q[4])/2 + p[3,0],-self.length[4]*ca.cos(q[4])/2 + p[3,1]
        
        dp = ca.jtimes(p,q,dq)
        dc = ca.jtimes(c,q,dq)
        ddp = ca.jtimes(dp,q,dq)
        ddc = ca.jtimes(dc,q,dq)

        return p,dp,ddp,c,dc,ddc        

    def impactMap(self, q, dq, p, dp, g, dg):
        q_plus = q[::-1]

        p_plus = ca.MX.zeros(5, 2)

        p_plus[0:4, 0], p_plus[0:4, 1] = p[0:4, 0][::-1], p[0:4, 1][::-1]
        p_plus[4, 0], p_plus[4, 1] = self.p0[0], self.p0[1]
        p_plus = ca.MX.zeros(5, 2)


        g_plus = ca.MX.zeros(5, 2)
        g_plus[0:5, 0], g_plus[0:5, 1] = g[0:5, 0][::-1], g[0:5, 1][::-1]

        dg_plus = ca.MX.zeros(5, 2)
        dg_plus[0:5, 0], dg_plus[0:5, 1] = dg[0:5, 0][::-1], dg[0:5, 1][::-1]


        dq_plus = ca.MX.zeros(5, 1)

        dq_plus[0] = ( self.inertia[0]*dq[0] 
                        + self.mass[0]*self.crossProduct2D(g[0, :]-p[0, :], dg[0, :]) 
                        - self.mass[4]*self.crossProduct2D(g_plus[0, :]-p_plus[0, :], dg_plus[0, :])
                        )/self.inertia[4]
        
        dq_plus[1] = ( self.inertia[0]*dq[0]
                        + self.inertia[1]*dq[1] 
                        + self.mass[0]*self.crossProduct2D(g[0, :]-p[0, :], dg[0, :])
                        + self.mass[1]*self.crossProduct2D(g[1, :]-p[1, :], dg[1, :]) 
                        - self.mass[4]*self.crossProduct2D(g_plus[0, :]-p_plus[0, :], dg_plus[0, :])
                        - self.mass[3]*self.crossProduct2D(g_plus[1, :]-p_plus[1, :], dg_plus[1, :])
                        - self.inertia[4]*dq_plus[0] 
                        )/self.inertia[3]
        
        dq_plus[2] = ( self.inertia[0]*dq[0]
                        + self.inertia[1]*dq[1] 
                        + self.inertia[2]*dq[2] 
                        + self.mass[0]*self.crossProduct2D(g[0, :]-p[0, :], dg[0, :])
                        + self.mass[1]*self.crossProduct2D(g[1, :]-p[1, :], dg[1, :]) 
                        + self.mass[2]*self.crossProduct2D(g[2, :]-p[2, :], dg[2, :]) 
                        - self.mass[4]*self.crossProduct2D(g_plus[0, :]-p_plus[0, :], dg_plus[0, :])
                        - self.mass[3]*self.crossProduct2D(g_plus[1, :]-p_plus[1, :], dg_plus[1, :])
                        - self.mass[2]*self.crossProduct2D(g_plus[2, :]-p_plus[2, :], dg_plus[2, :])
                        - self.inertia[4]*dq_plus[0] 
                        - self.inertia[3]*dq_plus[1] 
                        )/self.inertia[2]

        dq_plus[3] = ( self.inertia[0]*dq[0]
                        + self.inertia[1]*dq[1] 
                        + self.inertia[2]*dq[2] 
                        + self.inertia[3]*dq[3] 
                        + self.mass[0]*self.crossProduct2D(g[0, :]-p[0, :], dg[0, :])
                        + self.mass[1]*self.crossProduct2D(g[1, :]-p[1, :], dg[1, :]) 
                        + self.mass[2]*self.crossProduct2D(g[2, :]-p[2, :], dg[2, :]) 
                        + self.mass[3]*self.crossProduct2D(g[3, :]-p[3, :], dg[3, :]) 
                        - self.mass[4]*self.crossProduct2D(g_plus[0, :]-p_plus[0, :], dg_plus[0, :])
                        - self.mass[3]*self.crossProduct2D(g_plus[1, :]-p_plus[1, :], dg_plus[1, :])
                        - self.mass[2]*self.crossProduct2D(g_plus[2, :]-p_plus[2, :], dg_plus[2, :])
                        - self.mass[1]*self.crossProduct2D(g_plus[3, :]-p_plus[3, :], dg_plus[3, :])
                        - self.inertia[4]*dq_plus[0] 
                        - self.inertia[3]*dq_plus[1] 
                        - self.inertia[2]*dq_plus[2] 
                        )/self.inertia[1]        

        dq_plus[4] = ( self.inertia[0]*dq[0]
                        + self.inertia[1]*dq[1] 
                        + self.inertia[2]*dq[2] 
                        + self.inertia[3]*dq[3] 
                        + self.inertia[4]*dq[4] 
                        + self.mass[0]*self.crossProduct2D(g[0, :]-p[0, :], dg[0, :])
                        + self.mass[1]*self.crossProduct2D(g[1, :]-p[1, :], dg[1, :]) 
                        + self.mass[2]*self.crossProduct2D(g[2, :]-p[2, :], dg[2, :]) 
                        + self.mass[3]*self.crossProduct2D(g[3, :]-p[3, :], dg[3, :]) 
                        + self.mass[4]*self.crossProduct2D(g[4, :]-p[4, :], dg[4, :]) 
                        - self.mass[4]*self.crossProduct2D(g_plus[0, :]-p_plus[0, :], dg_plus[0, :])
                        - self.mass[3]*self.crossProduct2D(g_plus[1, :]-p_plus[1, :], dg_plus[1, :])
                        - self.mass[2]*self.crossProduct2D(g_plus[2, :]-p_plus[2, :], dg_plus[2, :])
                        - self.mass[1]*self.crossProduct2D(g_plus[3, :]-p_plus[3, :], dg_plus[3, :])
                        - self.mass[0]*self.crossProduct2D(g_plus[4, :]-p_plus[4, :], dg_plus[4, :])
                        - self.inertia[4]*dq_plus[0] 
                        - self.inertia[3]*dq_plus[1] 
                        - self.inertia[2]*dq_plus[2] 
                        - self.inertia[1]*dq_plus[3] 
                        )/self.inertia[0]
        return q_plus, dq_plus

    def crossProduct2D(self, a, b):
        return (a[0]*b[1]) - (b[1]*a[1])

    def heightMap(self, x):
        if self.terrain == 'sin':
            return self.terrain_factor*ca.sin(x)
        if self.terrain == 'wedge':
            return self.terrain_factor*x

    def heightMapNormalVector(self, x):
        if self.terrain == 'sin':
            tangent_vector = ca.MX.ones(2, 1)
            tangent_vector[1, 0] = self.terrain_factor*ca.cos(x)
            normal_vector = ca.MX.ones(2, 1)
            normal_vector[0, 0] = -tangent_vector[1, 0]
            normal_vector = normal_vector/ca.norm_2(normal_vector)
            return normal_vector
        if self.terrain == 'wedge':    
            tangent_vector = ca.MX.ones(2, 1)
            tangent_vector[1, 0] = self.terrain_factor
            normal_vector = ca.MX.ones(2, 1)
            normal_vector[0, 0] = -tangent_vector[1, 0]
            normal_vector = normal_vector/ca.norm_2(normal_vector)
            return normal_vector
        # m_tangent = self.terrain_factor*ca.cos(x)
        # m_normal = -1/m_tangent
        # theta_normal = ca.atan(m_normal)
        # return theta_normal


class nlp(walker):
    def __init__(self, walker):
        self.cost = self.getCost(walker.u,walker.N,walker.h)
        walker.opti.minimize(self.cost)
        self.ceq = self.getConstraints(walker)
        walker.opti.subject_to(self.ceq)
        self.bounds = self.getBounds(walker)
        walker.opti.subject_to(self.bounds)
        p_opts = {"expand":True}
        s_opts = {"max_iter": 3000}
        walker.opti.solver("ipopt",p_opts,s_opts)
        self.initial = self.initalGuess(walker)

    def initalGuess(self,walker):
        iniq = ca.DM.zeros((5,walker.N))
        inidq = ca.DM.zeros((5,walker.N))
        iniu = ca.DM.zeros((4,walker.N))
        for j in range(5):
            for i in range(walker.N):
                walker.opti.set_initial(walker.x[i][j],
                            (walker.ini_goal[j] + (i/(walker.N - 1))*(walker.fin_goal[j] - walker.ini_goal[j])))
                iniq[j,i] = (walker.ini_goal[j] + (i/(walker.N - 1))*(walker.fin_goal[j] - walker.ini_goal[j]))
                
                walker.opti.set_initial(walker.xdot[i][j],
                             (walker.fin_goal[j] - walker.ini_goal[j])/(walker.N-1))
                inidq[j,i] = (walker.fin_goal[j] - walker.ini_goal[j])/(walker.N-1)                
                
                if j < 4:
                    walker.opti.set_initial(walker.u[i][j], 0)
                
        return [iniq,inidq,iniu]

    def getCost(self,u,N,h):
        result = 0
        for i in range(N-1): 
            result += (h/2)*(ca.sumsqr(u[i]) + ca.sumsqr(u[i+1]))
        return result

    def getConstraints(self,walker):
        ceq = []
        for i in range(walker.N-1):
            q1 = walker.x[i]
            q2 = walker.x[i+1]
            dq1 = walker.xdot[i]
            dq2 = walker.xdot[i+1]
            ddq1 = walker.ddq[i]
            ddq2 = walker.ddq[i+1]
            ceq.extend(self.getCollocation(q1,q2,
                                        dq1,dq2,ddq1,ddq2,
                                        walker.h))
        
        q0 = walker.x[0]
        dq0 = walker.xdot[0]
        qf = walker.x[-1]
        dqf = walker.x[-1]
        ceq.extend(self.getBoundaryConstrainsts(q0, dq0, qf, dqf, walker.goal, walker.x_impact, walker.xdot_impact))
        # ceq.extend([walker.pos[0][4, 0] <= walker.p0[0]])
        # ceq.extend([walker.pos[0][4, 1] == walker.heightMap(walker.pos[0][4, 0])])

        ceq.extend([    
                    # (ca.dot(walker.dpos[0][4,:].T, walker.heightMapNormalVector(walker.pos[0][4, 0])) > 0.),
                    # (ca.dot(walker.dpos[-1][4,:].T, walker.heightMapNormalVector(walker.pos[-1][4, 0])) < 0.)
                    ])

        # ceq.extend([
        #             (walker.dpos[0][4, 1]*ca.sin(walker.heightMapNormal(walker.pos[0][4, 0], walker.pos[0][4, 1])) > 0),
        #             (walker.dpos[0][4, 0]*ca.cos(walker.heightMapNormal(walker.pos[0][4, 0], walker.pos[0][4, 1])) > 0),
        #             (walker.dpos[-1][4, 0]*ca.cos(walker.heightMapNormal(walker.pos[-1][4, 0], walker.pos[-1][4, 1])) < 0),
        #             (walker.dpos[-1][4, 1]*ca.sin(walker.heightMapNormal(walker.pos[-1][4, 0], walker.pos[-1][4, 1])) < 0),
        #             ])

        for i in range(1,walker.N-1):
            ceq.extend([
                        # (walker.pos[i][4, 0] <=  walker.step_max + walker.p0[0, 0]),
                        # (walker.pos[i][4, 0] >= -walker.step_max - walker.p0[0, 0]),

                        (walker.pos[i][0, 1] > walker.heightMap(walker.pos[i][0, 0])),
                        (walker.pos[i][1, 1] > walker.heightMap(walker.pos[i][1, 0])),
                        (walker.pos[i][2, 1] > walker.heightMap(walker.pos[i][2, 0])),
                        (walker.pos[i][3, 1] > walker.heightMap(walker.pos[i][3, 0])),
                        (walker.pos[i][4, 1] > walker.heightMap(walker.pos[i][4, 0])),
                        # (walker.pos[i][4, 1] < walker.heightMap(walker.pos[i][4, 0]) + walker.comh),
                        
                        # (walker.pos[i][3, 1] > walker.heightMap(walker.pos[i][3, 0])),
                        # (walker.pos[i][2, 1] > walker.heightMap(walker.pos[i][2, 0])),
                        # (walker.pos[i][1, 1] > walker.heightMap(walker.pos[i][1, 0])),
                        # (walker.pos[i][0, 1] > walker.heightMap(walker.pos[i][0, 0])),

                        (walker.pos[i][4, 1] < walker.pos[i][3, 1]),
                        # (walker.pos[i][1, 1] < walker.pos[i][2, 1]),
                        (walker.pos[i][4, 1] < walker.pos[i][0, 1]),
                        
                        # (walker.pos[i][4, 1] - walker.heightMap(walker.pos[i][4, 0] + walker.p0[0, 0]) > 0),
                        # (walker.pos[i][3, 1] - walker.heightMap(walker.pos[i][3, 0]) >= walker.comh),
                        # (walker.pos[i][2, 1] - walker.heightMap(walker.pos[i][2, 0]) >= walker.comh),
                        # (walker.pos[i][1, 1] - walker.heightMap(walker.pos[i][1, 0]) >= walker.comh),
                        # (walker.pos[i][0, 1] - walker.heightMap(walker.pos[i][0, 0]) >= walker.comh)
                        ])
            # ceq.extend([((walker.x[i][0, 0] + walker.x[i][1, 0]) < 0)])
            # ceq.extend([((-walker.x[i][0, 0] - walker.x[i][1, 0]) + ca.pi > ca.pi*3/2)])
            # ceq.extend([
            #     (ca.acos(ca.norm_2(walker.pos[i][1,:]-walker.pos[i][0,:])*ca.norm_2(walker.pos[i][0,:])/(ca.dot(walker.pos[i][1,:]-walker.pos[i][0,:],walker.pos[i][0,:]))) > 0)
            # ])

            ceq.extend([
                        # (ca.fabs(walker.x[i][0, 0]) < ca.fabs(walker.x[i][1, 0])),
                        # (ca.fabs(walker.x[i][4, 0]) < ca.fabs(walker.x[i][3, 0])),
                        # (ca.fabs(walker.x[i][0, 0]-walker.x[i][1, 0]) > 0),
                        # (ca.fabs(walker.x[i][4, 0]-walker.x[i][3, 0]) > 0),
                        
                        # (walker.x[i][0, 0] > walker.x[i][1, 0]), # human
                        # (walker.x[i][4, 0] < walker.x[i][3, 0]), # human
                        
                        (walker.x[i][0, 0] < walker.x[i][1, 0]), # ostrich
                        (walker.x[i][4, 0] > walker.x[i][3, 0]), # ostrich
                        ])

            ceq.extend([(walker.x[i][2, 0] <= walker.pi/5)])
            ceq.extend([(walker.x[i][2, 0] >= -walker.pi/5)])
            # ceq.extend([(walker.x[i][1, 0] >= 0 )])
            # ceq.extend([(walker.x[i][1, 0] >= 0 )])


            # ceq.extend([((walker.state[i][0][4, 0] - walker.state[i][0][3]) <= walker.pi/2)])

            # ceq.extend([((walker.pos[i][1] - walker.p0[1]) >= walker.comh)])
            # ceq.extend([((walker.pos[i][1][1] - walker.p0[1]) >= walker.comh)])
            # ceq.extend([((walker.pos[i][2][1] - walker.p0[1]) >= walker.comh)])
            # if i == 0:
            #     ceq.extend([((walker.pos[i][4][1]) == walker.p0[1])])
            # comy = 0
            # for j in range(4):
            #     comy += walker.com[i][j][1]
            # ceq.extend([comy - walker.p0[1] >= walker.comh])

        ceq.extend([walker.pos[-1][4, 0] >= 0.5*walker.step_max + walker.p0[0]])
        ceq.extend([walker.pos[-1][4, 0] <= 1.5*walker.step_max + walker.p0[0]])
        # ceq.extend([walker.pos[int(len(walker.pos)/2)][4, 0] >= walker.p0[0]])
        ceq.extend([walker.pos[-1][4, 1] == walker.heightMap(walker.pos[-1][4, 0])])
        
        return ceq

    def getCollocation(self,q1,q2,dq1,dq2,ddq1,ddq2,h):
        cc = []
        cc.extend([(((h/2)*(dq2 + dq1)) - (q2 - q1)==0)])
        cc.extend([(((h/2)*(ddq2 + ddq1)) - (dq2 - dq1)==0)])

        return cc

    def getBoundaryConstrainsts(self, q0, dq0, qf, dqf, goal, q_plus, dq_plus):
        c = []
        c.extend([  
                    # (q0 == q_plus),
                    (dq0 == dq_plus),
                    (q0 - goal[0] == 0),
                    # (dq0 - goal[1] == 0),
                    # (qf - q_plus == 0),
                    # (dqf - dq_plus == 0)
        ])
        return c
    
    def getBounds(self,walker):
        c = []
        # f = 3
        for i in range(walker.N):
            q = walker.x[i]
            dq = walker.xdot[i]
            u = walker.u[i]
            c.extend([  walker.opti.bounded(ca.MX([-walker.pi/2]*5),q,ca.MX([walker.pi/2]*5)),
                        # walker.opti.bounded(ca.MX([-f*walker.pi]*5),dq,ca.MX([f*walker.pi]*5)),
                        walker.opti.bounded(ca.MX([-walker.tauMax]*4),u,ca.MX([walker.tauMax]*4)),
            ])
                    
        return c


