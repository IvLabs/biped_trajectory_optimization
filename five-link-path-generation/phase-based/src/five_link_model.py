import numpy as np
import casadi as ca
import helper_functions as hf

class Biped():
    def __init__(self):
        super().__init__()
        self.length = ca.MX([0.5,0.5,0.5,0.5,0.5])
        self.mass = ca.MX([0.25,0.25,0.25,0.25,0.25])
        self.inertia = self.mass * (self.length**2) /12
        self.gravity = -10

        self.gravity_vector = ca.MX.zeros(2)
        self.gravity_vector[1] = self.gravity

        self.b = ca.mmax(self.length*2)
        self.a = ca.mmax(self.length)

    def setFullState(self, q, dq, lp0, dlp0, ddlp0, rp0, drp0, ddrp0, u, lf10, rf10):
        self.q   = ca.reshape(  q, 5, 1)
        self.dq  = ca.reshape( dq, 5, 1)

        self.lp0   = ca.reshape(  lp0, 2, 1) 
        self.rp0   = ca.reshape(  rp0, 2, 1)
        self.dlp0  = ca.reshape( dlp0, 2, 1)
        self.drp0  = ca.reshape( drp0, 2, 1)
        self.ddlp0 = ca.reshape(ddlp0, 2, 1)
        self.ddrp0 = ca.reshape(ddrp0, 2, 1)

        self.u    = ca.reshape(   u, 4, 1)
        self.lf10 = ca.reshape(lf10, 2, 1)
        self.rf10 = ca.reshape(rf10, 2, 1)

        self.p, self.dp, self.ddp, self.c, self.dc, self.ddc = self.getKinematics()
        self.dynamics = self.getDynamics()

    def getKinematics(self):        
        # Assume p,c as 2x1

        ##################
        ###--Left Leg--###
        ##################
        lp = ca.MX.zeros(2, 2)
        lc = ca.MX.zeros(2, 2)

        lp0 = self.lp0

        lp[0,0],lp[1,0] = self.length[0]*ca.sin(self.q[0]) + lp0[0,0], self.length[0]*ca.cos(self.q[0]) + lp0[1,0]
        lp[0,1],lp[1,1] = self.length[1]*ca.sin(self.q[1]) +  lp[0,0], self.length[1]*ca.cos(self.q[1]) +  lp[1,0]

        lc[0,0],lc[1,0] = self.length[0]*ca.sin(self.q[0])/2 + lp0[0,0], self.length[0]*ca.cos(self.q[0])/2 + lp0[1,0]
        lc[0,1],lc[1,1] = self.length[1]*ca.sin(self.q[1])/2 +  lp[0,0], self.length[1]*ca.cos(self.q[1])/2 +  lp[1,0]

        ###--Derivatives--###
        dlp  = ca.jtimes( lp, self.q, self.dq) + self.dlp0
        dlc  = ca.jtimes( lc, self.q, self.dq) + self.dlp0

        ddlp = ca.jtimes(dlp, self.q, self.dq) + self.ddlp0
        ddlc = ca.jtimes(dlc, self.q, self.dq) + self.ddlp0

        ###################
        ###--Right Leg--###
        ###################
        rp = ca.MX.zeros(2, 2)
        rc = ca.MX.zeros(2, 2)

        rp0 = self.rp0

        rp[0,0],rp[1,0] = self.length[4]*ca.sin(self.q[4]) + rp0[0,0], self.length[4]*ca.cos(self.q[4]) + rp0[1,0]
        rp[0,1],rp[1,1] = self.length[3]*ca.sin(self.q[3]) +  rp[0,0], self.length[3]*ca.cos(self.q[3]) +  rp[1,0]

        rc[0,0],rc[1,0] = self.length[4]*ca.sin(self.q[4])/2 + rp0[0,0], self.length[4]*ca.cos(self.q[4])/2 + rp0[1,0]
        rc[0,1],rc[1,1] = self.length[3]*ca.sin(self.q[3])/2 +  rp[0,0], self.length[3]*ca.cos(self.q[3])/2 +  rp[1,0]

        ###--Derivatives--###
        drp  = ca.jtimes( rp, self.q, self.dq) + self.drp0
        drc  = ca.jtimes( rc, self.q, self.dq) + self.drp0

        ddrp = ca.jtimes(drp, self.q, self.dq) + self.ddrp0
        ddrc = ca.jtimes(drc, self.q, self.dq) + self.ddrp0

        ###############
        ###--Torso--###
        ###############
        tp = ca.MX.zeros(2)
        tc = ca.MX.zeros(2)

        tp[0],tp[1] = self.length[2]*ca.sin(self.q[2]) + lp[0,1], self.length[2]*ca.cos(self.q[2]) + lp[1,1]
        
        tc[0],tc[1] = self.length[2]*ca.sin(self.q[2])/2 + lp[0,1], self.length[2]*ca.cos(self.q[2])/2 + lp[1,1]
        
        ###--Derivatives--###
        dtp  = ca.jtimes( tp, self.q, self.dq) + self.dlp0
        dtc  = ca.jtimes( tc, self.q, self.dq) + self.dlp0

        ddtp = ca.jtimes(dtp, self.q, self.dq) + self.ddlp0
        ddtc = ca.jtimes(dtc, self.q, self.dq) + self.ddlp0

        ###--Form a dictionary--###
        p   = {'Left Leg' :   lp, 'Right Leg':   rp, 'Torso' :   tp, 'Constraint' : (  lp[:,1] ==   rp[:,1])}
        dp  = {'Left Leg' :  dlp, 'Right Leg':  drp, 'Torso' :  dtp, 'Constraint' : ( dlp[:,1] ==  drp[:,1])}
        ddp = {'Left Leg' : ddlp, 'Right Leg': ddrp, 'Torso' : ddtp, 'Constraint' : (ddlp[:,1] == ddrp[:,1])}

        c   = {'Left Leg' :   lc, 'Right Leg':   rc, 'Torso' :   tc}
        dc  = {'Left Leg' :  dlc, 'Right Leg':  drc, 'Torso' :  dtc}
        ddc = {'Left Leg' : ddlc, 'Right Leg': ddrc, 'Torso' : ddtc}

        return p,dp,ddp,c,dc,ddc      

    def getDynamics(self):
        ddq  = ca.MX.zeros(5)

        ###--Left Leg--###

        lp0  = self.lp0
        lp   = self.p['Left Leg']
        lc   = self.c['Left Leg']
        ddlc = self.ddc['Left Leg']
        f10  = self.f10

        f12 = self.mass[0]*(ddlc[:,0] - self.gravity_vector) - f10
        ddq[0] = -u[0] + hf.crossProduct2D(lp[:,0]-lc[:,0], f12) + hf.crossProduct2D(lp0-lc[:,0], f10)

        f21 = -f12
        f23_24 = self.mass[1]*(ddlc[:,1] - self.gravity_vector) - f21
        ddq[1] = u[0]-u[1] + hf.crossProduct2D(lp[:,1]-lc[:,1], f23_24) + hf.crossProduct2D(lp[:,0]-lc[:,1], f21)

        ###--Right Leg--###

        rp0  = self.rp0
        rp   = self.p['Right Leg']
        rc   = self.c['Right Leg']
        ddrc = self.ddc['Right Leg']
        f50  = self.f50
        
        f54 = self.mass[4]*(ddrc[:,0] - self.gravity_vector) - f50
        ddq[4] = -u[3] + hf.crossProduct2D(rp[:,0]-rc[:,0],f54) + hf.crossProduct2D(rp[:,0]-rc[:,0],f50)

        f45 = -f54
        f43_42 = self.mass[3]*(ddrc[:,1] - self.gravity_vector) - f45
        ddq[3] = u[3]-u[2] + hf.crossProduct2D(rp[:,1]-rc[:,1], f43_42) + hf.crossProduct2D(rp[:,0]-rc[:,1], f45)

        ###--Torso--###
        tp   = self.p['Torso']
        tc   = self.c['Torso']
        ddtc = self.ddc['Torso']

        f32_34 = self.mass[2]*(ddtc - self.gravity_vector)
        ddq[2] = u[2]+u[1] + self.crossProduct2D(lp[:,1]-tc, f32_34)

        ddq /= self.inertia
        
        dynamics = {'ddq':ddq, 'Contraints': (f32_34 + f23_24 + f43_42 == 0)}

        return dynamics

test_biped = Biped()

q    = ca.MX.sym(   'q', 5, 1)
dq   = ca.MX.sym(  'dq', 5, 1)
 
p0   = ca.MX.sym(  'p0', 2, 1)
p5   = ca.MX.sym(  'p5', 2, 1)
dp0  = ca.MX.sym( 'dp0', 2, 1)
dp5  = ca.MX.sym( 'dp5', 2, 1)
ddp0 = ca.MX.sym('ddp0', 2, 1)
ddp5 = ca.MX.sym('ddp5', 2, 1)
 
u    = ca.MX.sym(   'u', 4, 1)
f10  = ca.MX.sym( 'f10', 2, 1)
f50  = ca.MX.sym( 'f50', 2, 1)

test_biped.setFullState(q, dq, p0, dp0, ddp0, p5, dp5, ddp5, u, f10, f50)