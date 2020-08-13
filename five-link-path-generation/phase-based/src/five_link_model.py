import numpy as np
import casadi as ca
import helper_functions as hf

class Biped():
    def __init__(self):
        super().__init__()
        self.length = ca.MX([0.5,0.5,0.5,0.5,0.5])
        self.mass = ca.MX([0.25,0.25,0.25,0.25,0.25])
        self.inertia = self.mass * (self.length**2) /12
        self.gravity = 10
        self.b = 1.
        self.a = 0.5

    def setFullState(self, q, dq, lp0, dlp0, rp0, drp0, u, lf10, rf10):
        self.q   = ca.reshape(  q, 5, 1)
        self.dq  = ca.reshape( dq, 5, 1)

        self.lp0  = ca.reshape( lp0, 2, 1) 
        self.rp5  = ca.reshape( rp0, 2, 1)
        self.dlp0 = ca.reshape(dlp0, 2, 1)
        self.drp0 = ca.reshape(drp0, 2, 1)
        
        self.u    = ca.reshape(   u, 4, 1)
        self.lf10 = ca.reshape(lf10, 2, 1)
        self.rf10 = ca.reshape(rf10, 2, 1)

        self.p, self.dp, self.ddp, self.c, self.dc, self.ddc = self.getKinematics()

    def getKinematics(self):        
        # Assume p,c as 2x1

        ###--Left Leg--###
        lp = ca.MX.zeros(2, 2)
        lc = ca.MX.zeros(2, 2)

        lp0 = self.lp0

        lp[0,0],lp[1,0] = self.length[0]*ca.sin(self.q[0]) + lp0[0,0], self.length[0]*ca.cos(self.q[0]) + lp0[1,0]
        lp[0,1],lp[1,1] = self.length[1]*ca.sin(self.q[1]) +  lp[0,0], self.length[1]*ca.cos(self.q[1]) +  lp[1,0]

        lc[0,0],lc[1,0] = self.length[0]*ca.sin(self.q[0])/2 + lp0[0,0], self.length[0]*ca.cos(self.q[0])/2 + lp0[1,0]
        lc[0,1],lc[1,1] = self.length[1]*ca.sin(self.q[1])/2 +  lp[0,0], self.length[1]*ca.cos(self.q[1])/2 +  lp[1,0]
                
        ###--Right Leg--###
        rp = ca.MX.zeros(2, 2)
        rc = ca.MX.zeros(2, 2)

        rp0 = self.rp0

        rp[0,0],rp[1,0] = self.length[4]*ca.sin(self.q[4]) + rp0[0,0], self.length[4]*ca.cos(self.q[4]) + rp0[1,0]
        rp[0,1],rp[1,1] = self.length[3]*ca.sin(self.q[3]) +  rp[0,0], self.length[3]*ca.cos(self.q[3]) +  rp[1,0]

        rc[0,0],rc[1,0] = self.length[4]*ca.sin(self.q[4])/2 + rp0[0,0], self.length[4]*ca.cos(self.q[4])/2 + rp0[1,0]
        rc[1,1],rc[1,1] = self.length[3]*ca.sin(self.q[3])/2 +  rp[0,0], self.length[3]*ca.cos(self.q[3])/2 +  rp[1,0]

        ###--Torso--###
        p[0,2],p[1,2] = self.length[2]*ca.sin(self.q[2]) +   p[0,1], self.length[2]*ca.cos(self.q[2]) +   p[0,1]
        
        c[0,2],c[1,2] = self.length[2]*ca.sin(self.q[2])/2 + p[0,2], self.length[2]*ca.cos(self.q[2])/2 + p[1,2]
        
        dp  = ca.jtimes( p, self.q, self.dq)
        dc  = ca.jtimes( c, self.q, self.dq)
        ddp = ca.jtimes(dp, self.q, self.dq)
        ddc = ca.jtimes(dc, self.q, self.dq)



        return p,dp,ddp,c,dc,ddc      

    def getDynamics(self):
        ddq = ca.MX.zeros(5)

        f10 = self.f10
        f50 = self.f50

        f12 = self.mass[0]*ddc[:, 0] - f10
        ddq[0] = -u[0] + hf.crossProduct2D(f12, p[:, 1]-c[0, :]) + hf.crossProduct2D(f10, p[0, :]-c[0, :])

        f21 = -f12
        f23_24 = self.mass[1]*ddc[0] - f21
        ddq[1] = u[0]-u[1] + self.crossProduct2D(f21, p[2, :]-c[1, :]) + self.crossProduct2D(f23_24, p[1, :]-c[1, :])


        f54 = self.mass[4]*ddc[4] - f50
        ddq[4] = -u[3] + self.crossProduct2D(f54, p[4, :]-c[4, :]) + self.crossProduct2D(f50, p[5, :]-c[4, :])

        f45 = -f54
        f43_42 = self.mass[3]*ddc[3] - f45
        ddq[3] = u[3]-u[2] + self.crossProduct2D(f43_42, p[3, :]-c[3, :]) + self.crossProduct2D(f45, p[4, :]-c[3, :])


        f32_34 = self.mass[2]*ddc[2, :]
        ddq[2] = u[2]+u[1] + self.crossProduct2D(f32_34, p[2, :]-c[2, :])

        return ddq/self.inertia

# test_biped = Biped()

# q   = ca.MX.sym(  'q', 5, 1)
# dq  = ca.MX.sym( 'dq', 5, 1)
 
# p0  = ca.MX.sym( 'p0', 2, 1)
# p5  = ca.MX.sym( 'p5', 2, 1)
# dp0 = ca.MX.sym('dp0', 2, 1)
# dp5 = ca.MX.sym('dp5', 2, 1)

# u   = ca.MX.sym(  'u', 4, 1)
# f10 = ca.MX.sym('f10', 2, 1)
# f50 = ca.MX.sym('f50', 2, 1)

# test_biped.setFullState(q, dq, p0, dp0, p5, dp5, u, f10, f50)