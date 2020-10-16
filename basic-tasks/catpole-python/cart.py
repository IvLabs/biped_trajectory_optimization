
import casadi as ca
import numpy as np
from matplotlib import pyplot as plt
import cv2 
from PIL import Image

class cart():

    def __init__(self):
        self.opti = ca.Opti()
        self.N = 200
        self.d = 2.0
        self.dmax = 20.0
        self.umax = 20.0
        self.pi = -np.pi
        self.l = 0.5
        self.m1 = 0.5
        self.m2 = 0.1
        self.g = -9.81
        self.T = 2.0
        self.h = self.T/self.N

        goal1 = self.opti.parameter()
        goal2 = self.opti.parameter()
        goal3 = self.opti.parameter()
        goal4 = self.opti.parameter()
        self.opti.set_value(goal1,self.d)
        self.opti.set_value(goal2,self.pi)
        self.opti.set_value(goal3,0)
        self.opti.set_value(goal4,0)
        self.goal = [goal1,goal2,goal3,goal4]

        self.q1 = self.opti.variable(self.N)
        self.q2 = self.opti.variable(self.N)    
        self.q1_dot = self.opti.variable(self.N)
        self.q2_dot = self.opti.variable(self.N)
        
        self.u = self.opti.variable(self.N)

        self.q1_ddot = (((self.l*self.m2*ca.sin(self.q2)*(self.q2**2))
                    + (self.u) + (self.m2*self.g*ca.cos(self.q2)*ca.sin(self.q2)))
                     /(self.m1 + self.m2*(ca.sin(self.q2)**2)))
        self.q2_ddot = (((self.l*self.m2*ca.cos(self.q2)*ca.sin(self.q2)*(self.q2_dot**2)) 
                    + (self.u*ca.cos(self.q2)) + ((self.m1+self.m2)*self.g*ca.sin(self.q2)))
                    /(self.m1*self.l + self.l*self.m2*(ca.sin(self.q2)**2)))

        self.state = ca.horzcat(self.q1,self.q2,self.q1_dot,self.q2_dot)
        self.model = ca.horzcat(self.q1_dot,self.q2_dot,self.q1_ddot,self.q2_ddot)

class nlp_cart(cart):
    def __init__(self,cart):
        # self.x = ca.horzcat(cart.state,cart.u)#ca.vertcat(cart.state,cart.u)
        cart.opti.minimize(self.getCost(cart.u,cart.N,cart.h))
        cart.opti.subject_to(self.getConstraints(cart))
        cart.opti.subject_to(self.getBounds(cart))
        p_opts = {"expand":True}
        s_opts = {"max_iter": 1000}
        cart.opti.solver("ipopt",p_opts,s_opts)
        self.initial = self.initalGuess(cart)

    def initalGuess(self,cart):
        ini = np.zeros((cart.N,5))
        for i in range(cart.N):
            cart.opti.set_initial(cart.state[i,0],(i/(cart.N - 1))*cart.d)
            ini[i,0] = (i/(cart.N - 1))*cart.d
            cart.opti.set_initial(cart.state[i,1],(i/(cart.N - 1))*cart.pi)
            ini[i,1]=(i/(cart.N - 1))*cart.pi
            cart.opti.set_initial(cart.state[i,2],cart.d/(cart.N - 1))
            ini[i,2] = cart.d/(cart.N - 1)
            cart.opti.set_initial(cart.state[i,3],cart.pi/(cart.N - 1))
            ini[i,3] = cart.pi/(cart.N - 1)
            cart.opti.set_initial(cart.u[i],0)
        return ini

    def getCost(self,u,N,h):
        result = 0
        for i in range(N-1): result += (h/2)*(u[i]**2 + u[i+1]**2)  
        return result

    def getConstraints(self,cart):
        ceq = []
        for i in range(cart.N - 1): 
            ceq.append(self.getCollocationConstraints(cart.state[i,:],cart.state[i+1,:]
                                              ,cart.model[i,:],cart.model[i+1,:]
                                              ,cart.h))
        ceq.extend(self.getBoundaryConstrainsts(cart.goal,cart.state[0,:],cart.state[-1,:]))
        return ceq

    def getCollocationConstraints(self,state1,state2,model1,model2,h):
        return (((h/2)*(model2 + model1)) - (state2 - state1)==0)

    def getBoundaryConstrainsts(self,goal,state1,state2):
        ceq = []
        for i in range(4): ceq.extend([(state1[i] == 0),((state2[i] - goal[i]) ==0)]) 
        return ceq

    def getBounds(self,cart):
        c = [] 
        c.extend([cart.opti.bounded(-cart.dmax,cart.state[:,0],cart.dmax),
                cart.opti.bounded(cart.pi,cart.state[:,1],-cart.pi),
                # cart.opti.bounded(-cart.d,cart.state[:,0],cart.d),
                # cart.opti.bounded(-cart.d,cart.state[:,0],cart.d),
                cart.opti.bounded(-cart.umax,cart.u[:],cart.umax)])
        return c

mycart = cart()
ini = nlp_cart(mycart)

sol1 = mycart.opti.solve()
# # print(sol1.stats()["iter_count"])

pos = sol1.value(mycart.state[:,0])    
ang = sol1.value(mycart.state[:,1])    
u = sol1.value(mycart.u[:])
time = np.arange(0.0, mycart.T, mycart.h)

ball = 3
box = 1
rod = 2
d = {1: (255,175,0), #BGR format
     2: (0,255,0),
     3: (0,0,255)}
SIZE = 300

for i in range(len(pos)):
    env = np.ones((SIZE,SIZE,3),  dtype=np.uint8)

    p1= pos[i]*50
    ax= mycart.l*np.sin(ang[i])*100
    ay = mycart.l*np.cos(ang[i])*100
    ax = -ax.astype(np.int)
    ay = ay.astype(np.int)
    p = p1.astype(np.int)
    # print(ang[i])
    image = cv2.rectangle(env, (p+50 - 10,SIZE//2 - 3), (p+50+10 ,SIZE//2 + 3), d[box], -1)   
    image1 = cv2.line(image, (p+50,SIZE//2), (p+50 + ax, SIZE//2 + ay), d[rod], 1)
    image2 = cv2.circle(image1, (p+50+ax,SIZE//2 +ay), 5, d[ball], -1)
    img = Image.fromarray(image2, "RGB")
    # img = img.resize((300,300))
    cv2.imshow("", np.array(img))
    status = cv2.imwrite(f'/home/shitwalker/PycharmProjects/Biped/animate/image-{i}.png', image2)
    print(status)
    if cv2.waitKey(10) & 0xFF == ord("q"):
        break


fourcc = cv2.VideoWriter_fourcc('M','J','P','G')
out = cv2.VideoWriter('cartpole.avi', fourcc, 60.0, (300, 300))

for i in range(len(pos)):
    img_path = f"/home/shitwalker/PycharmProjects/Biped/animate/image-{i}.png"
    print(img_path)
    frame = cv2.imread(img_path)
    out.write(frame)

out.release()

name = ['position : x', 'angle : theta', 'controlForce : u']

plt.subplot(322)
plt.title('Optimised Solution')
plt.plot(time,pos)


plt.subplot(321)
plt.title('Initial Guess')
plt.plot(time,ini.initial[:,0],'r')
plt.ylabel(name[0])

plt.subplot(324)
plt.plot(time,ang)

plt.subplot(323)
plt.plot(time,ini.initial[:,1],'r')
plt.ylabel(name[1])

plt.subplot(326)
plt.plot(time,u)
plt.subplot(325)
plt.plot(time,ini.initial[:,4],'r')
plt.ylabel(name[2])

plt.suptitle('CartPole')
plt.show()




