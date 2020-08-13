import numpy as np
import casadi as ca

class Biped():
    def __init__(self):
        super().__init__()
        self.length = ca.MX([0.5,0.5,0.5,0.5,0.5])
        self.mass = ca.MX([0.25,0.25,0.25,0.25,0.25])
        self.inertia = self.mass * (self.length**2) /12
        self.gravity = 10

    def setFullState(self, q, dq, p0, dp0, p5, dp5, u, f10, f50):
        self.q   = ca.reshape(  q, 5, 1)
        self.dq  = ca.reshape( dq, 5, 1)

        self.p0  = ca.reshape( p0, 2, 1) 
        self.p5  = ca.reshape( p5, 2, 1)
        self.dp0 = ca.reshape(dp0, 2, 1)
        self.dp5 = ca.reshape(dp5, 2, 1)
        
        self.u   = ca.reshape(  u, 4, 1)
        self.f10 = ca.reshape(f10, 2, 1)
        self.f50 = ca.reshape(f50, 2, 1)

        self.p, self.dp, self.ddp, self.c, self.dc, self.ddc = self.getKinematics()

    def getKinematics(self):        
        p = ca.MX.zeros(2, 6)
        c = ca.MX.zeros(2, 5)
        

        p[0,0],p[1,0] = self.p0[0,0] , self.p0[1,0]
        p[0,1],p[1,1] = self.length[0]*ca.sin(self.q[0]) + p[0,0], self.length[0]*ca.cos(self.q[0]) + p[1,0]
        p[0,2],p[1,2] = self.length[1]*ca.sin(self.q[1]) + p[0,1], self.length[1]*ca.cos(self.q[1]) + p[1,1]
        p[0,3],p[1,3] = self.length[2]*ca.sin(self.q[2]) + p[0,2], self.length[2]*ca.cos(self.q[2]) + p[1,2]
        p[0,4],p[1,4] = self.length[3]*ca.sin(self.q[3]) + p[0,2],-self.length[3]*ca.cos(self.q[3]) + p[1,2]
        p[0,5],p[1,5] = self.length[4]*ca.sin(self.q[4]) + p[0,4],-self.length[4]*ca.cos(self.q[4]) + p[1,4]
        
        c[0,0],c[1,0] = self.length[0]*ca.sin(self.q[0])/2 + p[0,0], self.length[0]*ca.cos(self.q[0])/2 + p[1,0]
        c[0,1],c[1,1] = self.length[1]*ca.sin(self.q[1])/2 + p[0,1], self.length[1]*ca.cos(self.q[1])/2 + p[1,1]
        c[0,2],c[1,2] = self.length[2]*ca.sin(self.q[2])/2 + p[0,2], self.length[2]*ca.cos(self.q[2])/2 + p[1,2]
        c[0,3],c[1,3] = self.length[3]*ca.sin(self.q[3])/2 + p[0,2],-self.length[3]*ca.cos(self.q[3])/2 + p[1,2]
        c[0,4],c[1,4] = self.length[4]*ca.sin(self.q[4])/2 + p[0,4],-self.length[4]*ca.cos(self.q[4])/2 + p[1,4]
        
        dp  = ca.jtimes( p, self.q, self.dq)
        dc  = ca.jtimes( c, self.q, self.dq)
        ddp = ca.jtimes(dp, self.q, self.dq)
        ddc = ca.jtimes(dc, self.q, self.dq)

        return p,dp,ddp,c,dc,ddc      


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