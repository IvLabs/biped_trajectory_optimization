import numpy as np
import casadi as ca
import helper_functions as hf

class Biped():
    def __init__(self):
        super().__init__()
        self.num_ee = 2
        self.length = ca.MX([0.5,0.5,0.5,0.5,0.5])
        self.mass = ca.MX([0.25,0.25,0.25,0.25,0.25])
        self.inertia = self.mass * (self.length**2) /12
        self.gravity = -10

        self.gravity_vector = ca.MX.zeros(2)
        self.gravity_vector[1] = self.gravity
        
        self.b = ca.mmax(self.length*2)
        self.a = ca.mmax(self.length)

    def setFullState(self, lq, dlq, rq, drq, tq, dtq, lp0, dlp0, ddlp0, rp0, drp0, ddrp0, u, lf10, rf10):
        self.lq   = ca.reshape(  lq, 2, 1)
        self.dlq  = ca.reshape( dlq, 2, 1)

        self.rq   = ca.reshape(  rq, 2, 1)
        self.drq  = ca.reshape( drq, 2, 1)
        
        self.tq   = ca.reshape(  tq, 1, 1)
        self.dtq  = ca.reshape( dtq, 1, 1)

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

        lp[0,0],lp[1,0] = self.length[0]*ca.sin(self.lq[0]) + lp0[0,0]  , self.length[0]*ca.cos(self.lq[0]) + lp0[1,0]
        lp[0,1],lp[1,1] = self.length[1]*ca.sin(self.lq[1]) +  lp[0,0]  , self.length[1]*ca.cos(self.lq[1]) +  lp[1,0]
        
        lc[0,0],lc[1,0] = self.length[0]*ca.sin(self.lq[0])/2 + lp0[0,0], self.length[0]*ca.cos(self.lq[0])/2 + lp0[1,0]
        lc[0,1],lc[1,1] = self.length[1]*ca.sin(self.lq[1])/2 +  lp[0,0], self.length[1]*ca.cos(self.lq[1])/2 +  lp[1,0]

        ###--Derivatives--###
        dlp  = ca.jtimes( lp, self.lq, self.dlq) + self.dlp0
        dlc  = ca.jtimes( lc, self.lq, self.dlq) + self.dlp0

        ddlp = ca.jtimes(dlp, self.lq, self.dlq) + self.ddlp0
        ddlc = ca.jtimes(dlc, self.lq, self.dlq) + self.ddlp0

        ###################
        ###--Right Leg--###
        ###################
        rp = ca.MX.zeros(2, 2)
        rc = ca.MX.zeros(2, 2)

        rp0 = self.rp0

        rp[0,0],rp[1,0] = self.length[4]*ca.sin(self.rq[0]) + rp0[0,0]  , self.length[4]*ca.cos(self.rq[0]) + rp0[1,0]
        rp[0,1],rp[1,1] = self.length[3]*ca.sin(self.rq[1]) +  rp[0,0]  , self.length[3]*ca.cos(self.rq[1]) +  rp[1,0]

        rc[0,0],rc[1,0] = self.length[4]*ca.sin(self.rq[0])/2 + rp0[0,0], self.length[4]*ca.cos(self.rq[0])/2 + rp0[1,0]
        rc[0,1],rc[1,1] = self.length[3]*ca.sin(self.rq[1])/2 +  rp[0,0], self.length[3]*ca.cos(self.rq[1])/2 +  rp[1,0]

        ###--Derivatives--###
        drp  = ca.jtimes( rp, self.rq, self.drq) + self.drp0
        drc  = ca.jtimes( rc, self.rq, self.drq) + self.drp0

        ddrp = ca.jtimes(drp, self.rq, self.drq) + self.ddrp0
        ddrc = ca.jtimes(drc, self.rq, self.drq) + self.ddrp0

        ###############
        ###--Torso--###
        ###############
        tp = ca.MX.zeros(2)
        tc = ca.MX.zeros(2)

        tp[0],tp[1] = self.length[2]*ca.sin(self.tq[0]) + lp[0,1], self.length[2]*ca.cos(self.tq[0]) + lp[1,1]
        
        tc[0],tc[1] = self.length[2]*ca.sin(self.tq[0])/2 + lp[0,1], self.length[2]*ca.cos(self.tq[0])/2 + lp[1,1]
        
        ###--Derivatives--###
        dtp  = ca.jtimes( tp, self.tq, self.dtq) + self.dlp0
        dtc  = ca.jtimes( tc, self.tq, self.dtq) + self.dlp0

        ddtp = ca.jtimes(dtp, self.tq, self.dtq) + self.ddlp0
        ddtc = ca.jtimes(dtc, self.tq, self.dtq) + self.ddlp0

        ###--Form a dictionary--###
        p   = {'Left Leg' :   lp, 'Right Leg':   rp, 'Torso' :   tp, 'Constraint' : (  lp[:,1] ==   rp[:,1])}
        dp  = {'Left Leg' :  dlp, 'Right Leg':  drp, 'Torso' :  dtp, 'Constraint' : ( dlp[:,1] ==  drp[:,1])}
        ddp = {'Left Leg' : ddlp, 'Right Leg': ddrp, 'Torso' : ddtp, 'Constraint' : (ddlp[:,1] == ddrp[:,1])}

        c   = {'Left Leg' :   lc, 'Right Leg':   rc, 'Torso' :   tc}
        dc  = {'Left Leg' :  dlc, 'Right Leg':  drc, 'Torso' :  dtc}
        ddc = {'Left Leg' : ddlc, 'Right Leg': ddrc, 'Torso' : ddtc}

        return p,dp,ddp,c,dc,ddc      

    def getDynamics(self):

        ###--Left Leg--###
        ddlq  = ca.MX.zeros(2)
        lp0  = self.lp0
        lp   = self.p['Left Leg']
        lc   = self.c['Left Leg']
        ddlc = self.ddc['Left Leg']
        f10  = self.lf10

        f12 = self.mass[0]*(ddlc[:,0] - self.gravity_vector) - f10
        ddlq[0] = -self.u[0] + hf.crossProduct2D(lp[:,0]-lc[:,0], f12) + hf.crossProduct2D(lp0-lc[:,0], f10)

        f21 = -f12
        f23_24 = self.mass[1]*(ddlc[:,1] - self.gravity_vector) - f21
        ddlq[1] = self.u[0]-self.u[1] + hf.crossProduct2D(lp[:,1]-lc[:,1], f23_24) + hf.crossProduct2D(lp[:,0]-lc[:,1], f21)

        ###--Right Leg--###
        ddrq  = ca.MX.zeros(2)
        rp0  = self.rp0
        rp   = self.p['Right Leg']
        rc   = self.c['Right Leg']
        ddrc = self.ddc['Right Leg']
        f50  = self.rf10
        
        f54 = self.mass[4]*(ddrc[:,0] - self.gravity_vector) - f50
        ddrq[0] = -self.u[3] + hf.crossProduct2D(rp[:,0]-rc[:,0],f54) + hf.crossProduct2D(rp[:,0]-rc[:,0],f50)

        f45 = -f54
        f43_42 = self.mass[3]*(ddrc[:,1] - self.gravity_vector) - f45
        ddrq[1] = self.u[3]-self.u[2] + hf.crossProduct2D(rp[:,1]-rc[:,1], f43_42) + hf.crossProduct2D(rp[:,0]-rc[:,1], f45)

        ###--Torso--###
        ddtq  = ca.MX.zeros(1)
        tp   = self.p['Torso']
        tc   = self.c['Torso']
        ddtc = self.ddc['Torso']

        f32_34 = self.mass[2]*(ddtc - self.gravity_vector)
        ddtq[0] = self.u[2]+self.u[1] + hf.crossProduct2D(lp[:,1]-tc, f32_34)

        ddlq /= self.inertia[0:2]
        ddrq /= self.inertia[3:5]
        ddtq /= self.inertia[2]

        
        dynamics = {'Left Leg': ddlq, 'Right Leg': ddrq, 'Torso': ddtq, 'Constraint': (f32_34 + f23_24 + f43_42 == 0)}

        return dynamics


# test check for sanity

# test_biped = Biped()

# lq    = ca.MX.sym(   'q', 2, 1)
# ldq   = ca.MX.sym(  'dq', 2, 1)

# rq    = ca.MX.sym(   'q', 2, 1)
# rdq   = ca.MX.sym(  'dq', 2, 1)

# tq    = ca.MX.sym(   'q', 1, 1)
# tdq   = ca.MX.sym(  'dq', 1, 1)


# p0   = ca.MX.sym(  'p0', 2, 1)
# p5   = ca.MX.sym(  'p5', 2, 1)
# dp0  = ca.MX.sym( 'dp0', 2, 1)
# dp5  = ca.MX.sym( 'dp5', 2, 1)
# ddp0 = ca.MX.sym('ddp0', 2, 1)
# ddp5 = ca.MX.sym('ddp5', 2, 1)
 
# u    = ca.MX.sym(   'u', 4, 1)
# f10  = ca.MX.sym( 'f10', 2, 1)
# f50  = ca.MX.sym( 'f50', 2, 1)

# test_biped.setFullState(lq, ldq, rq, rdq, tq, tdq, p0, dp0, ddp0, p5, dp5, ddp5, u, f10, f50)
# print(test_biped.dynamics)