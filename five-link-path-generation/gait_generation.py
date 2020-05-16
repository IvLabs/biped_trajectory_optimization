import casadi as ca
import numpy as np

class walker():
    def __init__(self,start,start_pos):
        # set our parameters of optimization
        self.opti = ca.Opti()
        self.N = 100; self.T = 0.1
        self.step_max = 1.5; self.tauMax = 1.0
        self.pi = np.pi; 
        self.l1 = 0.5; self.l2 = 0.5; self.l3 = 0.6
        self.m = [0.5,0.5,0.5,0.5,0.5]#; self.m2 = 0.5; self.m3 = 0.5
        self.i1 = self.m[0]*(self.l1**2)/12; self.i2 = self.m[1]*(self.l2**2)/12; self.i3 = self.m[2]*(self.l3**2)/12
        self.i = [self.i1,self.i2,self.i3,self.i2,self.i1]
        self.g = -9.81
        self.h = self.T/self.N
        self.comh = 0.1
        goali = start; goalf = goali[::-1]
        self.goal = [goali,goalf]
        self.p0 = start_pos
        #set our optimization variables
        self.state = []
        self.u = []
        for i in range(self.N): 
            rowu = []; rowq = []
            rowq.append(self.opti.variable(5))    
            rowq.append(self.opti.variable(5))
            rowu.append(self.opti.variable(4))
            self.state.append(rowq)
            self.u.append(rowu)

        self.pos = [];self.com = [];self.ddq = []
        for i in range(self.N):
            p,dp,g,dg,ddq = self.getModel(self.state[i],self.u[i])
            self.pos.append(p); self.com.append(g);self.ddq.append(ddq)
            if i == 0:
                self.dp0 = dp
            if i == self.N - 1:
                self.dpN = dp
                self.impactmap = self.heelStrike(self.state[i][0],self.state[i][1],p,dp,g,dg)
        
    def getModel(self,state,u):
        q = state[0]
        dq = state[1]
        u = u[0]
        p0 = [self.p0[0],self.p0[1]]
        p,dp,g,dc = self.getKinematics(q,dq)
    
        ddc10,ddc11 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc20,ddc21 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc30,ddc31 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc40,ddc41 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc50,ddc51 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc = [[ddc10,ddc11],[ddc20,ddc21],[ddc30,ddc31],[ddc40,ddc41],[ddc50,ddc51]]
        s = [0,0,0,0,0]
        for i in range(5):
            s[0] += (((g[i][0]-p0[0])*self.m[i]*self.g)) + (((g[i][0]-p0[0])*self.m[i]*ddc[i][1])-((g[i][1])*self.m[i]*ddc[i][0]))
            if i > 0:
                s[1] += (((g[i][0]-p[0][0])*self.m[i]*self.g)) + (((g[i][0]-p[0][0])*self.m[i]*ddc[i][1])-((g[i][1]-p[0][1])*self.m[i]*ddc[i][0]))
                if i > 1:    
                    s[2] += (((g[i][0]-p[1][0])*self.m[i]*self.g)) + (((g[i][0]-p[1][0])*self.m[i]*ddc[i][1])-((g[i][1]-p[1][1])*self.m[i]*ddc[i][0]))
                    if i > 2:
                        s[3] += (((g[i][0]-p[1][0])*self.m[i]*self.g)) + (((g[i][0]-p[1][0])*self.m[i]*ddc[i][1])-((g[i][1]-p[1][1])*self.m[i]*ddc[i][0]))
                        if i > 3:
                            s[4] += (((g[i][0]-p[3][0])*self.m[i]*self.g)) + (((g[i][0]-p[3][0])*self.m[i]*ddc[i][1])-((g[i][1]-p[3][1])*self.m[i]*ddc[i][0]))
        
        ddq5 = (u[3]- s[4]) / self.i1
        ddq4 = ((u[2]-s[3]) - (self.i1*ddq5))/self.i2  
        ddq3 = ((u[1]-s[2]) - (self.i1*ddq5) - (self.i2*ddq4))/self.i3  
        ddq2 = ((u[0]-s[1]) - (self.i1*ddq5) - (self.i2*ddq4) - (self.i3*ddq3))/self.i2  
        ddq1 = ((-s[0]) - (self.i1*ddq5) - (self.i2*ddq4) - (self.i3*ddq3) - (self.i2*ddq4))/self.i1  
        ddq = [ddq1,ddq2,ddq3,ddq4,ddq5]
                
        return p,dp,g,dc,ddq
        
    def getKinematics(self,q,dq):    
        p10 = ca.MX.sym('p10',1); c10 = ca.MX.sym('g10',1)
        p11 = ca.MX.sym('p11',1); c11 = ca.MX.sym('g11',1)
        p20 = ca.MX.sym('p20',1); c20 = ca.MX.sym('g20',1)
        p21 = ca.MX.sym('p21',1); c21 = ca.MX.sym('g21',1)
        p30 = ca.MX.sym('p30',1); c30 = ca.MX.sym('g30',1)
        p31 = ca.MX.sym('p31',1); c31 = ca.MX.sym('g31',1)
        p40 = ca.MX.sym('p40',1); c40 = ca.MX.sym('g40',1)
        p41 = ca.MX.sym('p41',1); c41 = ca.MX.sym('g41',1)
        p50 = ca.MX.sym('p50',1); c50 = ca.MX.sym('g50',1)
        p51 = ca.MX.sym('p51',1); c51 = ca.MX.sym('g51',1)
        
        p0 = [self.p0[0],self.p0[1]]
        
        p10,p11 = -self.l1*ca.sin(q[0]) + p0[0],self.l1*ca.cos(q[0]) + p0[1]
        p20,p21 = -self.l2*ca.sin(q[1]) + p10,self.l2*ca.cos(q[1]) + p11
        p30,p31 = self.l3*ca.sin(q[2]) + p20,self.l3*ca.cos(q[2]) + p21
        p40,p41 = (self.l2*ca.sin(q[3])) + p20,(-self.l2*ca.cos(q[3])) + p21
        p50,p51 = (self.l1*ca.sin(q[4])) + p40,(-self.l1*ca.cos(q[4])) + p41

        dp10,dp11 = ca.jtimes(p10,ca.MX(q),dq) , ca.jtimes(p11,ca.MX(q),dq)
        dp20,dp21 = ca.jtimes(p20,ca.MX(q),dq) , ca.jtimes(p21,ca.MX(q),dq)
        dp30,dp31 = ca.jtimes(p30,ca.MX(q),dq) , ca.jtimes(p31,ca.MX(q),dq)
        dp40,dp41 = ca.jtimes(p40,ca.MX(q),dq) , ca.jtimes(p41,ca.MX(q),dq)
        dp50,dp51 = ca.jtimes(p50,ca.MX(q),dq) , ca.jtimes(p51,ca.MX(q),dq) 
        p = [[p10,p11],[p20,p21],[p30,p31],[p40,p41],[p50,p51]]
        dp = [[dp10,dp11],[dp20,dp21],[dp30,dp31],[dp40,dp41],[dp50,dp51]]

        ##########################################
        c10,c11 = -self.l1*ca.sin(q[0])/2 + p0[0],self.l1*ca.cos(q[0])/2 + p0[0]
        c20,c21 = -self.l2*ca.sin(q[1])/2 + p10,self.l2*ca.cos(q[1])/2 + p11
        c30,c31 = self.l3*ca.sin(q[2])/2 + p20,self.l3*ca.cos(q[2])/2 + p21
        c40,c41 = (self.l2*ca.sin(q[3])/2) + p20,(-self.l2*ca.cos(q[3])/2) + p21
        c50,c51 = (self.l1*ca.sin(q[4])/2) + p40,(-self.l1*ca.cos(q[4])/2) + p41
        # print(c1)
        # print(c2)
        dc10,dc11 = ca.jtimes(c10,ca.MX(q),dq) , ca.jtimes(c11,ca.MX(q),dq)
        dc20,dc21 = ca.jtimes(c20,ca.MX(q),dq) , ca.jtimes(c21,ca.MX(q),dq)
        dc30,dc31 = ca.jtimes(c30,ca.MX(q),dq) , ca.jtimes(c31,ca.MX(q),dq)
        dc40,dc41 = ca.jtimes(c40,ca.MX(q),dq) , ca.jtimes(c41,ca.MX(q),dq)
        dc50,dc51 = ca.jtimes(c50,ca.MX(q),dq) , ca.jtimes(c51,ca.MX(q),dq)
        g = [[c10,c11],[c20,c21],[c30,c31],[c40,c41],[c50,c51]]
        dc = [[dc10,dc11],[dc20,dc21],[dc30,dc31],[dc40,dc41],[dc50,dc51]]
        # print(dc1)
        return p,dp,g,dc

    def heelStrike(self,q,dq,p,dp,g,dg):
        qi = [q[4],q[3],q[2],q[1],q[0]]
        pi = [p[4],p[3],p[2],p[1],p[0]]
        dpi = [dp[4],dp[3],dp[2],dp[1],dp[0]]
        gi = [g[4],g[3],g[2],g[1],g[0]]
        dgi = [dg[4],dg[3],dg[2],dg[1],dg[0]]
        p0 = [self.p0[0],self.p0[1]]
        s = [0,0,0,0,0]
        for i in range(5):
            s[0] += ((self.m[i]*(((g[i][0] - p[4][0])*dg[i][1])-((g[i][1] - p[4][1])*dg[i][0]))) + (self.i[i]*dq[i])
                    - (self.m[i]*(((gi[i][0] - p0[0])*dgi[i][1])-((gi[i][1] - p0[1])*dgi[i][0]))))
            if i < 4:
                s[1] += ((self.m[i]*(((g[i][0] - p[3][0])*dg[i][1])-((g[i][1] - p[3][1])*dg[i][0]))) + (self.i[i]*dq[i])
                        - (self.m[i]*(((gi[i][0] - pi[1][0])*dgi[i][1])-((gi[i][1] - pi[1][1])*dgi[i][0]))))
                if i < 3:
                    s[2] += ((self.m[i]*(((g[i][0] - p[1][0])*dg[i][1])-((g[i][1] - p[1][1])*dg[i][0]))) + (self.i[i]*dq[i])
                            - (self.m[i]*(((gi[i][0] - pi[1][0])*dgi[i][1])-((gi[i][1] - pi[1][1])*dgi[i][0]))))
                    if i < 2:
                        s[3] += ((self.m[i]*(((g[i][0] - p[1][0])*dg[i][1])-((g[i][1] - p[1][1])*dg[i][0]))) + (self.i[i]*dq[i])
                                - (self.m[i]*(((gi[i][0] - pi[1][0])*dgi[i][1])-((gi[i][1] - pi[1][1])*dgi[i][0]))))
                        if i < 1:
                            s[4] += ((self.m[i]*(((g[i][0] - p[0][0])*dg[i][1])-((g[i][1] - p[0][1])*dg[i][0]))) + (self.i[i]*dq[i])
                                    - (self.m[i]*(((gi[i][0] - pi[3][0])*dgi[i][1])-((gi[i][1] - pi[3][1])*dgi[i][0]))))
        dqi = [0,0,0,0,0]
        dqi[4] = s[4] / self.i[4]
        dqi[3] = s[3] - s[4] / self.i[3]
        dqi[2] = s[2] - s[3] / self.i[2]
        dqi[1] = s[1] - s[2]/ self.i[1]
        dqi[0] = s[0] - s[1] / self.i[0]

        return [qi,dqi]

class nlp(walker):
    def __init__(self, walker):
        self.cost = self.getCost(walker.u,walker.N,walker.h)
        walker.opti.minimize(self.cost)
        self.ceq = self.getConstraints(walker)
        walker.opti.subject_to(self.ceq)
        self.bounds = self.getBounds(walker)
        walker.opti.subject_to(self.bounds)
        p_opts = {"expand":True}
        s_opts = {"max_iter": 1000}
        walker.opti.solver("ipopt",p_opts,s_opts)
        self.initial = self.initalGuess(walker)

    def initalGuess(self,walker):
        iniq = np.zeros((5,walker.N))
        inidq = np.zeros((5,walker.N))
        iniu = np.zeros((4,walker.N))
        for j in range(5):
            for i in range(walker.N):
                walker.opti.set_initial(walker.state[i][0][j],
                            (walker.goal[0][j] + (i/(walker.N - 1))*(walker.goal[1][j] - walker.goal[0][j])))
                iniq[j,i] = (walker.goal[0][j] + (i/(walker.N - 1))*(walker.goal[1][j] - walker.goal[0][j]))
                
                walker.opti.set_initial(walker.state[i][1][j],
                            (walker.goal[1][j] - walker.goal[0][j])/(walker.N-1))
                inidq[j,i] = (walker.goal[1][j] - walker.goal[0][j])/(walker.N-1)                
                
                if j < 4:
                    walker.opti.set_initial(walker.u[i][0][j],0)
                
        return [iniq,inidq,iniu]

    def getCost(self,u,N,h):
        result = 0
        for i in range(N-1): 
            for j in range(4):
                result += (h/2)*(u[i][0][j]**2 + u[i+1][0][j]**2)
        return result

    def getConstraints(self,walker):
        ceq = []
        for i in range(walker.N-1):
            q1 = (walker.state[i][0])
            q2 = (walker.state[i+1][0])
            dq1 = (walker.state[i][1])
            dq2 = (walker.state[i+1][1])
            ddq1 = walker.ddq[i]
            ddq2 = walker.ddq[i+1]
            ceq.extend(self.getCollocation(q1,q2,
                                        dq1,dq2,ddq1,ddq2,
                                        walker.h))
        
        q0 = (walker.state[0][0])
        dq0 = (walker.state[0][1])
        qf = (walker.state[-1][0])
        dqf = (walker.state[-1][1])
        ceq.extend(self.getBoundaryConstrainsts(q0,dq0,qf,walker.goal,walker.impactmap))
        ceq.extend([(walker.dp0[4][1] > 0),(walker.dpN[4][1] < 0)])
        for i in range(walker.N):
            ceq.extend([((walker.pos[i][4][0]) <= (walker.step_max + walker.p0[0]))])
            ceq.extend([((walker.pos[i][4][0]) >= -(walker.step_max + walker.p0[0]))])
            ceq.extend([((walker.pos[i][4][1]) >= walker.p0[1])])
            ceq.extend([((walker.pos[i][3][1] - walker.p0[1]) >= walker.comh)])
            # ceq.extend([((walker.pos[i][0][1] - walker.p0[1]) >= walker.comh)])
            # ceq.extend([((walker.pos[i][1][1] - walker.p0[1]) >= walker.comh)])
            # ceq.extend([((walker.pos[i][2][1] - walker.p0[1]) >= walker.comh)])
            if i == 0:
                ceq.extend([((walker.pos[i][4][1]) == walker.p0[1])])
        # ceq.extend([walker.pos[-1][4][1] == walker.p0[1]])
            comy = 0
            for j in range(4):
                comy += walker.com[i][j][1]
            ceq.extend([comy - walker.p0[1] >= walker.comh])

        return ceq

    def getCollocation(self,q1,q2,dq1,dq2,ddq1,ddq2,h):
        cc = []
        for i in range(4):
            cc.extend([(((h/2)*(ddq2[i] + ddq1[i])) - (dq2[i] - dq1[i])==0)])
        cc.extend([(((h/2)*(dq2 + dq1)) - (q2 - q1)==0)])
        return cc

    def getBoundaryConstrainsts(self,state1,dstate1,state2,goal,impact):
        c = []
        for i in range(4): 
            c.extend([(state1[i] - impact[0][i] == 0),(dstate1[i] - impact[1][i] == 0)
                                    ,((state2[i] - goal[1][i]) ==0)]) 

        return c
    
    def getBounds(self,walker):
        c = []
        f = 15
        for i in range(walker.N):
            q = (walker.state[i][0])
            dq = (walker.state[i][1])
            u = (walker.u[i][0])
            c.extend([walker.opti.bounded(-walker.pi,q[0],walker.pi),
                    walker.opti.bounded(-walker.pi,q[1],walker.pi),
                    walker.opti.bounded(-walker.pi,q[2],walker.pi),
                    walker.opti.bounded(-walker.pi,q[3],walker.pi),
                    walker.opti.bounded(-walker.pi,q[4],walker.pi),
                    walker.opti.bounded(-f*walker.pi,dq[0],f*walker.pi),
                    walker.opti.bounded(-f*walker.pi,dq[1],f*walker.pi),
                    walker.opti.bounded(-f*walker.pi,dq[2],f*walker.pi),
                    walker.opti.bounded(-f*walker.pi,dq[3],f*walker.pi),
                    walker.opti.bounded(-f*walker.pi,dq[4],f*walker.pi),
                    # walker.opti.bounded(0,u[0],0),
                    walker.opti.bounded(-walker.tauMax,u[0],walker.tauMax),
                    walker.opti.bounded(-walker.tauMax,u[1],walker.tauMax),
                    walker.opti.bounded(-walker.tauMax,u[2],walker.tauMax),
                    walker.opti.bounded(-walker.tauMax,u[3],walker.tauMax)])
        return c
