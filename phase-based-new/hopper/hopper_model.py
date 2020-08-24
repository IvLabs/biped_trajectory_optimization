import numpy as np
import casadi as ca
import helper_functions as hf

class Hopper():
    def __init__(self):
        super().__init__()
        self.num_ee    = 1
        self.length    = np.array([0.5,0.5,0.5])
        self.mass      = np.array([0.25,0.25,0.25])
        self.initial_q = np.zeros((3,1))
        self.final_q   = np.zeros((3,1))
        self.inertia   = self.mass * (self.length**2) /12
        self.gravity   = -10

        self.gravity_vector = ca.MX.zeros(2)
        self.gravity_vector[1] = self.gravity
        
        self.b = ca.mmax(self.length*2)
        self.a = ca.mmax(self.length)

    def setFullState(self, q, dq, c3, dc3, u, f10):
        self.q   = ca.reshape( q, 3, 1)
        self.dq  = ca.reshape(dq, 3, 1)

        self.c3   = ca.reshape(  c3, 2, 1)
        self.dc3  = ca.reshape( dc3, 2, 1)
        self.ddc3 = ca.reshape(ddc3, 2, 1)

        self.u   = ca.reshape(  u, 2, 1)
        self.f10 = ca.reshape(f10, 2, 1)

        self.p, self.dp, self.ddp, self.c, self.dc, self.ddc = self.getKinematics()
        self.dynamics = self.getDynamics()

    def getKinematics(self):        
        # Assume p,c as 2x1

        ##################
        #####--Leg--######
        ##################
        p = ca.MX.zeros(2, 4)
        c = ca.MX.zeros(2, 3)

        c3 = self.c3

        p[0,3],p[1,3] =  self.length[2]*ca.sin(self.q[2])/2 + c3[0,0],  self.length[2]*ca.cos(self.q[2])/2 + c3[1,0]
        p[0,2],p[1,2] = -self.length[2]*ca.sin(self.q[2])/2 + c3[0,0], -self.length[2]*ca.cos(self.q[2])/2 + c3[1,0]
        p[0,1],p[1,1] = -self.length[1]*ca.sin(self.q[1]) + p[0,2]   , -self.length[1]*ca.cos(self.q[1]) + p[1,2]
        p[0,0],p[1,0] = -self.length[0]*ca.sin(self.q[0]) + p[0,1]   , -self.length[1]*ca.cos(self.q[0]) + p[1,1]

        c[0,2],c[1,2] = c3[0,0], c3[1,0]
        c[0,1],c[1,1] = -self.length[1]*ca.sin(self.q[1])/2 + p[0,2], -self.length[1]*ca.cos(self.q[1])/2 + p[1,2]
        c[0,0],c[1,0] = -self.length[0]*ca.sin(self.q[0])/2 + p[0,1], -self.length[0]*ca.cos(self.q[0])/2 + p[1,1]
 
        ###--Derivatives--###
        dp  = ca.jtimes( p, self.q, self.dq) + self.dc3
        dc  = ca.jtimes( c, self.q, self.dq) + self.dc3

        ddp = ca.jtimes(dp, self.q, self.dq**2) + self.ddc3
        ddc = ca.jtimes(dc, self.q, self.dq**2) + self.ddc3


        ###--Form a dictionary--###
        p   = {'Leg' :   p, 'Constraint' : (self.q[0]<=self.q[1])}
        dp  = {'Leg' :  dp}
        ddp = {'Leg' : ddp}

        c   = {'Leg' :   c}
        dc  = {'Leg' :  dc}
        ddc = {'Leg' : ddc}

        return p,dp,ddp,c,dc,ddc      

    def getDynamics(self):

        ###--Leg--###
        ddq  = ca.MX.zeros(3)
        c3  = self.c3
        p   = self.p  ['Leg']
        c   = self.c  ['Leg']
        ddc = self.ddc['Leg']
        f10 = self.f10

        f32 = self.mass[2]*(ddc[:,2] - self.gravity_vector)
        ddq[2] = self.u[1] + hf.crossProduct2D(p[:,2]-c[:,2], f32)

        f12 = self.mass[0]*(ddc[:,0] - self.gravity_vector) - f10
        ddq[0] = -self.u[0] + hf.crossProduct2D(p[:,1]-c[:,0], f12) + hf.crossProduct2D(p[:,0]-c[:,0], f10)

        f21 = -f12
        f23 = -f32
        ddq[1] = self.u[0]-self.u[1] + hf.crossProduct2D(p[:,2]-c[:,1], f23) + hf.crossProduct2D(p[:,1]-c[:,1], f21)

        ddq /= self.inertia

        dynamics = {'Leg': ddq}

        return dynamics

# test check for sanity

test_hopper = Hopper()

q    = ca.MX.sym(   'q', 3, 1)
dq   = ca.MX.sym(  'dq', 3, 1)

p0   = ca.MX.sym(  'p0', 2, 1)

dp0  = ca.MX.sym( 'dp0', 2, 1)

u    = ca.MX.sym(   'u', 2, 1)
f10  = ca.MX.sym( 'f10', 2, 1)

test_hopper.setFullState(q, dq, p0, dp0, u, f10)

print(test_hopper.dynamics)