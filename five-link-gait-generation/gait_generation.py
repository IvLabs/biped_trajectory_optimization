import casadi as ca
import numpy as np
from matplotlib import pyplot as plt
import cv2
from PIL import Image

class walker():
    def __init__(self):
        # set our parameters of optimization
        self.opti = ca.Opti()
        self.N = 100; self.T = 0.1
        self.step_max = 1; self.tauMax = 1.5
        self.pi = np.pi; 
        self.l1 = 0.5; self.l2 = 0.5; self.l3 = 0.6
        self.m = [0.5,0.5,0.5,0.5,0.5]#; self.m2 = 0.5; self.m3 = 0.5
        self.i1 = self.m[0]*(self.l1**2)/12; self.i2 = self.m[1]*(self.l2**2)/12; self.i3 = self.m[2]*(self.l3**2)/12
        self.i = [self.i1,self.i2,self.i3,self.i2,self.i1]
        self.g = -9.81
        self.h = self.T/self.N
        goali = [-0.3,0.7,0.0,-0.5,-0.3]; goalf = goali[::-1]
        self.goal = [goali,goalf]
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
        p,dp,g,dc = self.getKinematics(q,dq)
    
        ddc10,ddc11 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc20,ddc21 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc30,ddc31 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc40,ddc41 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc50,ddc51 = ca.jtimes(dc[0][0],dq,dq) , ca.jtimes(dc[0][1],dq,dq)
        ddc = [[ddc10,ddc11],[ddc20,ddc21],[ddc30,ddc31],[ddc40,ddc41],[ddc50,ddc51]]
        s = [0,0,0,0,0]
        for i in range(5):
            s[0] += (((g[i][0])*self.m[i]*self.g)) + (((g[i][0])*self.m[i]*ddc[i][1])-((g[i][1])*self.m[i]*ddc[i][0]))
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
        
        p0 = [0,0]
        
        p10,p11 = -self.l1*ca.sin(q[0]) + p0[0],self.l1*ca.cos(q[0]) + p0[0]
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
        s = [0,0,0,0,0]
        for i in range(5):
            s[0] += ((self.m[i]*(((g[i][0] - p[4][0])*dg[i][1])-((g[i][1] - p[4][1])*dg[i][0]))) + (self.i[i]*dq[i])
                    - (self.m[i]*(((gi[i][0] - pi[0][0])*dgi[i][1])-((gi[i][1] - pi[0][1])*dgi[i][0]))))
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
        ceq.extend(self.getBoundaryConstrainsts(q0,dq0,qf,walker.goal,walker.impactmap))
        ceq.extend([(walker.dp0[4][1] > 0),(walker.dpN[4][1] < 0)])
        for i in range(walker.N):
            ceq.extend([((walker.pos[i][4][0]) <= walker.step_max)])
            ceq.extend([((walker.pos[i][4][0]) >= -walker.step_max)])
            ceq.extend([((walker.pos[i][4][1]) > 0)])
            ceq.extend([((walker.pos[i][4][1]) <= 0.1)])
        return ceq

    def getCollocation(self,q1,q2,dq1,dq2,ddq1,ddq2,h):
        cc = []
        for i in range(4):
            cc.extend([(((h/2)*(ddq2[i] + ddq1[i])) - (dq2[i] - dq1[i])==0)])
        cc.extend([(((h/2)*(dq2 + dq1)) - (q2 - q1)==0)])
        return cc

    def getBoundaryConstrainsts(self,state1,dstate1,state2,goal,impact):
        c = []
        for i in range(4): c.extend([(state1[i] - impact[0][i] == 0),(dstate1[i] - impact[1][i] == 0)
                                    ,((state2[i] - goal[1][i]) ==0)]) 
        return c
    
    def getBounds(self,walker):
        c = []
        f = 20
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

me = walker()
n1 = nlp(me)

sol1 = me.opti.solve()

q = []; dq = []; u = []; pos = []
for j in range(5):
    tempq = [];tempdq = [];tempu = [];temp = []
    for i in range(me.N):
        tempq.append(sol1.value(me.state[i][0][j]))    
        tempdq.append(sol1.value(me.state[i][1][j]))
        temp.append([sol1.value(me.pos[i][j][0]),sol1.value(me.pos[i][j][1])]) 
        if j < 4:
            tempu.append(sol1.value(me.u[i][0][j]))
    q.append(tempq);pos.append(temp)        
    dq.append(tempdq)        
    u.append(tempu)        
time = np.arange(0.0, me.T, me.h)

from matplotlib import animation
from celluloid import Camera

fig = plt.figure()
camera = Camera(fig)
for i in range(me.N):
    p1 = [pos[0][i][1],pos[0][i][0]]
    p2 = [pos[1][i][1],pos[1][i][0]]
    p3 = [pos[2][i][1],pos[2][i][0]]
    p4 = [pos[3][i][1],pos[3][i][0]]
    p5 = [pos[4][i][1],pos[4][i][0]]
    # plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    # plt.plot([0,-p1[1]], [0,p1[0]],'r',[-p1[1],-p2[1]], [p1[0],p2[0]],'b',
    #     [-p2[1],-p3[1]], [p2[0],p3[0]],'c',
    #     [-p2[1],p4[1] - 2*p2[1]], [p2[0],2*p2[0]-p4[0]],'b',
    #     [p4[1] - 2*p2[1],p5[1]], [2*p2[0]-p4[0],(p5[0] - 2*p2[0])],'r')
    plt.plot([0,p1[1]], [0,p1[0]],'r',[p1[1],p2[1]], [p1[0],p2[0]],'b',
            [p2[1],p3[1]], [p2[0],p3[0]],'c', [p2[1],p4[1]], [p2[0],p4[0]],'b',
            [p4[1],p5[1]], [p4[0],p5[0]],'r')
    
    plt.plot([-2,2],[0,0],'g')    
    # if cv2.waitKey(0) & 0xFF == ord("q"):
    # #     break
    camera.snap()
animation = camera.animate(interval=60)
# animation.save('animation2.mp4')
plt.show()
plt.close()

name = ['q','dq','u']

plt.subplot(322)
plt.title('Optimised Solution')
plt.plot(time,q[:][0],'r',time,q[:][1],'g',time,q[:][2],'b',
        time,q[:][3],'y',time,q[:][4],'c')

plt.subplot(321)
plt.title('Initial Guess')
iniq = n1.initial[0]
plt.plot(time,iniq[:][0],'r',time,iniq[:][1],'g',time,iniq[:][2],'b',
        time,iniq[:][3],'y',time,iniq[:][4],'c')
plt.ylabel(name[0])

plt.subplot(324)
plt.plot(time,dq[:][0],'r',time,dq[:][1],'g',time,dq[:][2],'b',
        time,dq[:][3],'y',time,dq[:][4],'c')

plt.subplot(323)
inidq = n1.initial[1]
plt.plot(time,inidq[:][0],'r',time,inidq[:][1],'g',time,inidq[:][2],'b',
        time,inidq[:][3],'y',time,inidq[:][4],'c')
plt.ylabel(name[1])

plt.subplot(326)
plt.plot(time,u[:][0],'g',time,u[:][1],'b',time,u[:][2],'y',
        time,u[:][3],'c')

plt.subplot(325)
iniu = n1.initial[2]
plt.plot(time,iniu[:][0],'r',time,iniu[:][1],'g',time,iniu[:][2],'b',
        time,iniu[:][3],'y')
plt.ylabel(name[2])

plt.suptitle('Five-Link')
plt.show()






























































# print(pos)
# femur = 3
# tibia = 1
# torso = 2
# d = {1: (255,175,0), #BGR format
#      2: (0,255,0),
#      3: (0,0,255)}
# SIZE = 300
# # fourcc = cv2.VideoWriter_fourcc('M','J','P','G')
# # out = cv2.VideoWriter('cartpole.mp4', fourcc, 60.0, (300, 300))
# for i in range(me.N):
#     env = np.ones((SIZE,SIZE,3),  dtype=np.uint8)

#     p1 = [pos[0][i][0]*50,pos[0][i][1]*50];p1 = list(map(np.int, p1)) 
#     p2 = [pos[1][i][0]*50,pos[1][i][1]*50];p2 = list(map(np.int, p2))
#     p3 = [pos[2][i][0]*50,pos[2][i][1]*50];p3 = list(map(np.int, p3))
#     p4 = [pos[3][i][0]*50,pos[3][i][1]*50];p4 = list(map(np.int, p4))
#     p5 = [pos[4][i][0]*50,pos[4][i][1]*50];p5 = list(map(np.int, p5))
#     # print([p5[0]-p4[0],p5[1]-p4[1]])
#     # ax= mycart.l*np.sin(ang[i])*100
#     # ay = mycart.l*np.cos(ang[i])*100
#     # ax = -ax.astype(np.int)
#     # ay = ay.astype(np.int)
#     # # print(ang[i])
#     stance1 = cv2.line(env, (SIZE//2,SIZE//2), (-p1[1] + SIZE//2,-p1[0] + SIZE//2), d[femur], 1)   
#     stance2 = cv2.line(stance1, (-p1[1]+SIZE//2,-p1[0]+SIZE//2), (-p2[1]+SIZE//2,-p2[0]+SIZE//2), d[tibia], 1)   
#     tor = cv2.line(stance2, (-p2[1]+SIZE//2,-p2[0]+SIZE//2), (-p3[1]+SIZE//2,-p3[0]+SIZE//2), d[torso], 1)   
#     swing1 = cv2.line(tor, (-p2[1]+SIZE//2,-p2[0]+SIZE//2), (-p4[1]+SIZE//2,-p4[0]+SIZE//2), d[tibia], 1)   
#     swing2 = cv2.line(swing1, (-p4[1]+SIZE//2,-p4[0]+SIZE//2), (-p5[1]+SIZE//2,-p5[0]+SIZE//2), d[femur], 1)      
    
#     img = Image.fromarray(swing2, "RGB")
#     # img = img.resize((300,300))
#     cv2.imshow("", np.array(img))
#     # out.write(image2)
    
#     # status = cv2.imwrite(f'/home/shitwalker/PycharmProjects/Biped/animate/image-{i}.png', image2)
#     # print(status)
#     if cv2.waitKey(50) & 0xFF == ord("q"):
#         break





    # def getCollocation(self,walker):
    #     i=0# print(walker.state[0,0])
    #     # for i in range(walker.N-1):
    #     #     P,dP = self.getModel(walker.state[i,0:walker.SIZE//2:1],
    #     #                             walker.state[i,walker.SIZE//2:walker.SIZE:1],
    #     #                             walker.l1,walker.l2,walker.l3)
    #     P,dP = self.getModel(walker.state[i,0:walker.SIZE//2:1],
    #                                 walker.state[i,walker.SIZE//2:walker.SIZE:1],
    #                                 walker.l1,walker.l2,walker.l3)
    #     return [P,dP]

    # def getModel(self,q,dq,l1,l2,l3):
        # p = ca.MX.sym('p',ca.Sparsity.dense(5,2))
        # P[0,0] = 0
        # P[0,1] = 0
        # P[1,0] = l1*ca.cos(q[0]) + P[0,0]
        # P[1,1] = l1*ca.sin(q[0]) + P[0,1]
        # P[2,0] = l2*ca.cos(q[1]) + P[1,0]
        # P[2,1] = l2*ca.sin(q[1]) + P[1,1]
        # P[3,0] = l3*ca.cos(q[2]) + P[2,0]
        # P[3,1] = l3*ca.sin(q[2]) + P[2,1]
        # P[4,0] = l2*ca.cos(q[3]) + P[2,0]
        # P[4,1] = l2*ca.sin(q[3]) + P[2,1]
        # P[5,1] = l1*ca.sin(q[3]) + P[4,1]
        # P[5,1] = l1*ca.sin(q[3]) + P[4,1]

        # G = ca.MX.sym('g',6,2)
        # G[0,0] = ca.MX(0)
        # G[0,1] = ca.MX(0) 
        # G[1,0] = (l1*ca.cos(q[0])/2) + G[0,0]
        # G[1,1] = (l1*ca.sin(q[0])/2) + G[0,1]
        # G[2,0] = (l2*ca.cos(q[1])/2) + G[1,0]
        # G[2,1] = (l2*ca.sin(q[1])/2) + G[1,1]
        # G[3,0] = (l3*ca.cos(q[2])/2) + G[2,0]
        # G[3,1] = (l3*ca.sin(q[2])/2) + G[2,1]
        # G[4,0] = (l2*ca.cos(q[3])/2) + G[2,0]
        # G[4,1] = (l2*ca.sin(q[3])/2) + G[2,1]
        # G[5,1] = (l1*ca.sin(q[3])/2) + G[4,1]
        # G[5,1] = (l1*ca.sin(q[3])/2) + G[4,1]
        # p = ca.Sparsity.dense(5,2)
        # p[0,1] = p1[1]
        # p[0,0] = p1[0]
        # p[0,1] = p1[1]
        # p[1,0] = p2[0]
        # p[1,1] = p2[1]
        # p[2,0] = p3[0]
        # p[2,1] = p3[1]
        # p[3,0] = p4[0]
        # p[3,1] = p4[1]
        # p[4,0] = p5[0]
        # p[4,1] = p5[1]
        # print(p1)
        # qd = ca.MX.sym('q',1)
        # dqd = ca.MX.sym('dq',1)
        # x = ca.MX.sym('x',1,5)
        # # v = ca.MX.sym('v',5)
        # dp1 = ca.jacobian(q@q.T,q)
        # p = p1
        # dp = dp1
        # # print(p)
        # return p,dp

# print(n1.cost)
# print(n1.ceq)
        # p0 = [0,0]
        # p1 = [self.l1*ca.cos(self.state[:,0]) + p0[0],self.l1*ca.sin(self.state[:,0]) + p0[0]]
        # p2 = [self.l2*ca.cos(self.state[:,1]) + p1[0],self.l2*ca.sin(self.state[:,1]) + p1[1]]
        # p3 = [self.l3*ca.cos(self.state[:,2]) + p2[0],self.l3*ca.sin(self.state[:,2]) + p2[1]]
        # p4 = [self.l2*ca.cos(self.state[:,3]) + p2[0],self.l2*ca.sin(self.state[:,3]) + p2[1]]
        # p5 = [self.l1*ca.sin(self.state[:,4]) + p4[0],self.l1*ca.sin(self.state[:,4]) + p4[1]]

        # g1 = [self.l1*ca.cos(self.state[:,0])/2 + p0[0],self.l1*ca.sin(self.state[:,0])/2 + p0[0]]
        # g2 = [self.l2*ca.cos(self.state[:,1])/2 + p1[0],self.l2*ca.sin(self.state[:,1])/2 + p1[1]]
        # g3 = [self.l3*ca.cos(self.state[:,2])/2 + p2[0],self.l3*ca.sin(self.state[:,2])/2 + p2[1]]
        # g4 = [self.l2*ca.cos(self.state[:,3])/2 + p2[0],self.l2*ca.sin(self.state[:,3])/2 + p2[1]]
        # g5 = [self.l1*ca.sin(self.state[:,4])/2 + p4[0],self.l1*ca.sin(self.state[:,4])/2 + p4[1]]

        # p1 = ca.MX.sym('p1',self.N,2)
        # p1[:,0],p1[:,1] = self.l1*ca.cos(self.state[:,0]) + p0[0],self.l1*ca.sin(self.state[:,0]) + p0[0]
        # dp1 = ca.MX([ca.jacobian(p1,self.q1),ca.jacobian(p1,self.q2),ca.jacobian(p1,self.q3),ca.jacobian(p1,self.q4),ca.jacobian(p1,self.q5)])
        # # dp1 = ca.jtimes(p1,self.state[:,0:self.SIZE//2:1],self.state[:,self.SIZE//2:self.SIZE:1])
        # print(dp1)
        # return [p1,p2,p3,p4,p5],[g1,g2,g3,g4,g5]
        # A = ca.sparsify(ca.MX([[self.i1,self.i2,self.i3,self.i2,self.i1],
        #             [0,self.i2,self.i3,self.i2,self.i1],
        #             [0,0,self.i3,self.i2,self.i1],
        #             [0,0,0,self.i2,self.i1],
        #             [0,0,0,0,self.i1]]))
