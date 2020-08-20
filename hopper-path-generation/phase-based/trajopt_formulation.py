import numpy as np
import casadi as ca

# from five_link_model import Biped1
# from five_link_model import Biped2

from hopper_model import Hopper

from terrain import Terrain
# import helper_functions as hf 

##################################
#### Remove ddpos ###############
#################################

class NLP():
    def __init__(self, knot_points_per_phase, steps, total_duration, model='hopper', terrain='flat'):
        super().__init__()
        self.knot_points_per_phase = knot_points_per_phase
        self.num_phases            = steps
        self.total_duration        = total_duration

        if model == 'hopper':
            self.model = Hopper()
            self.time_phases = {}
        
        self.terrain = Terrain(type=terrain)

        self.opti   = ca.Opti()

        self.q    = {}
        self.qdot = {}
        self.u    = {}
        self.p0   = {} # i = 0
        self.dp0  = {} # i = 0
        self.f10  = {} # i = 0

        self.setVariables()
        self.setConstraints()

    def setVariables(self):    
        for step in range(self.num_phases):
            
            self.time_phases.update({str(step) : self.opti.variable(1)})

            self.q.update(   {str(step) : []})
            self.qdot.update({str(step) : []})

            self.p0.update(  {str(step) : []})
            self.dp0.update( {str(step) : []})
            
            self.u.update(   {str(step) : []})
            self.f10.update( {str(step) : []})

            for knot_point in range(self.knot_points_per_phase):
                self.q   [str(step)].append(self.opti.variable(3))
                self.qdot[str(step)].append(self.opti.variable(3))

                self.p0  [str(step)].append(self.opti.variable(2))
                self.dp0 [str(step)].append(self.opti.variable(2))
                
                self.u   [str(step)].append(self.opti.variable(2))
                self.f10 [str(step)].append(self.opti.variable(2))
                
                q    = self.q   [str(step)][-1]  
                qdot = self.qdot[str(step)][-1]
                u    = self.u   [str(step)][-1]
                p0   = self.p0  [str(step)][-1]
                dp0  = self.dp0 [str(step)][-1]
                f10  = self.f10 [str(step)][-1]
            
                self.model.setFullState(q=q, dq=qdot, 
                                        p0=p0, dp0=dp0, 
                                        u=u, f10=f10)

                ###--Add Model Constraints--###
                self.opti.subject_to(self.model.p       ['Constraint'])                                                 
                self.opti.subject_to(self.model.dynamics['Constraint'])

                # if j==0 and n==0:
                #     self.setInitialBoundary()
                #     com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                #     self.opti.subject_to(com_x==0)
                
                # if j==self.num_phases-1 and n==self.num_phases-1:
                #     com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                #     self.opti.subject_to(com_x==2)
        
    def setConstraints(self):
        self.setBoundary()
        self.setContactConstraints()
        self.setCollocationContraints()

    def setBoundary(self):
        self.opti.subject_to(self.q [str(0)][0] == self.model.initial_q)
        self.opti.subject_to(self.p0[str(0)][0] == np.zeros((2,1)))
        self.opti.subject_to(self.q[str(self.num_phases-1)][-1] == self.model.final_q)

    def setContactConstraints(self):

        self.opti.subject_to(sum(self.time_phases) == self.total_duration)

        for step in self.num_phases:
            
            # if step is even ==> Contact
            # if step is odd  ==> No Contact
            
            if step%2 != 0:
                for knot_point in range(self.knot_points_per_phase):
                    f10 = self.f10[str(step)][knot_point]
                    self.opti.subject_to(f10==0)
            else: 
                for knot_point in range(self.knot_points_per_phase):
                    p0  = self.p0 [knot_point]
                    dp0 = self.dp0[knot_point]
                    f10 = self.f10[knot_point]
                    
                    self.opti.subject_to((self.terrain.mu*f10[0])**2 - f10[1]**2 >= 0)
                    self.opti.subject_to(ca.dot(f10,p0) >= 0)
                    self.opti.subject_to(p0[1]==self.terrain.heightMap(p0[0]))
                    self.opti.subject_to(dp0==0)

    def setCollocationContraints(self):
        for step in range(self.num_phases):
            
            h = self.time_phases[str(step)]/self.knot_points_per_phase

            for n in range(self.knot_points_per_phase-1):
                q1   , q2     = self.q   [str(step)][n], self.q   [str(step)][n+1]
                qdot1, qdot2  = self.qdot[str(step)][n], self.qdot[str(step)][n+1]
                
                lpos1  , lpos2   = self.lpos  [str(j)+ '_' +str(n)], self.lpos  [str(j)+ '_' +str(n+1)]
                dlpos1 , dlpos2  = self.dlpos [str(j)+ '_' +str(n)], self.dlpos [str(j)+ '_' +str(n+1)]
                
                u1     , u2      = self.u     [str(j)+ '_' +str(n)], self.u     [str(j)+ '_' +str(n+1)]
                rforce1, rforce2 = self.rforce[str(j)+ '_' +str(n)], self.rforce[str(j)+ '_' +str(n+1)]

                m1 = Biped1()
                m1.setFullState(        lq=lq1,   dlq=lqdot1,        rq=rq1, 
                                    drq=rqdot1,       tq=tq1,    dtq=tqdot1, 
                                     lp0=lpos1,  dlp0=dlpos1, ddlp0=ddlpos1, 
                                     rp0=rpos1,  drp0=drpos1, ddrp0=ddrpos1, 
                                          u=u1, lf10=lforce1,  rf10=rforce1)
                
                lqddot1 = m1.dynamics['Left Leg']
                rqddot1 = m1.dynamics['Right Leg']
                tqddot1 = m1.dynamics['Torso']
                
                m2 = Biped1()
                m2.setFullState(        lq=lq2,   dlq=lqdot2,        rq=rq2, 
                                    drq=rqdot2,       tq=tq2,    dtq=tqdot2, 
                                     lp0=lpos2,  dlp0=dlpos2, ddlp0=ddlpos2, 
                                     rp0=rpos2,  drp0=drpos2, ddrp0=ddrpos2, 
                                          u=u2, lf10=lforce2,  rf10=rforce2)
                
                lqddot2 = m2.dynamics['Left Leg']
                rqddot2 = m2.dynamics['Right Leg']
                tqddot2 = m2.dynamics['Torso']
                
                self.opti.subject_to( (l_h/2) * (lqdot2  + lqdot1) == (lq2 - lq1) )
                self.opti.subject_to( (l_h/2) * (lqddot2 + lqddot1) == (lqdot2 - lqdot1) )

                self.opti.subject_to( (l_h/2) * (dlpos2  + dlpos1) == (lpos2 - lpos1) )
                # self.opti.subject_to( (l_h/2) * (ddlpos2 + ddlpos1) == (dlpos2 - dlpos1) )

                self.opti.subject_to( (r_h/2) * (rqdot2 + rqdot1) == (rq2 - rq1) )
                self.opti.subject_to( (r_h/2) * (rqddot2 + rqddot1) == (rqdot2 - rqdot1) )

                self.opti.subject_to( (r_h/2) * (drpos2  + drpos1) == (rpos2 - rpos1) )
                # self.opti.subject_to( (r_h/2) * (ddrpos2 + ddrpos1) == (drpos2 - drpos1) )

                self.opti.subject_to( (l_h/2) * (tqdot2 + tqdot1) == (tq2 - tq1) )
                self.opti.subject_to( (l_h/2) * (tqddot2 + tqddot1) == (tqdot2 - tqdot1) )



                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )
                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )

                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )
                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )

class NLP2():
    def __init__(self, knot_points_per_phase, steps, total_duration, model='biped', terrain='flat'):
        super().__init__()
        self.knot_points_per_phase = knot_points_per_phase
        self.num_phases            = steps
        self.total_duration        = total_duration

        if model == 'biped':
            self.model = Biped2()
            self.time_phases = {'Left Leg' : {}, 'Right Leg': {}}
        
        self.terrain = Terrain(type=terrain)

        self.opti   = ca.Opti()

        self.lq     = {}
        self.lqdot  = {}

        self.rq     = {}
        self.rqdot  = {}

        self.tq     = {}
        self.tqdot  = {}

        self.u      = {}

        self.lpos   = {} # i = 0
        self.dlpos  = {} # i = 0
        self.lforce = {} # i = 0
        
        self.rpos   = {} # i = 1
        self.drpos  = {} # i = 1
        self.rforce = {} # i = 1
        
        self.setVariable()
        self.setConstraints()

    def setVariable(self):
        for j in range(self.num_phases):
            
            for leg in self.time_phases:
                if j%2 == 0:
                    self.time_phases[leg].update({'Contact @ ' + str(j) : self.opti.variable(1)})
                else:
                    self.time_phases[leg].update({'No Contact @ ' + str(j): self.opti.variable(1)})

            for n in range(self.knot_points_per_phase):
                self.lq.update(    {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.lqdot.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})
                
                self.rq.update(    {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.rqdot.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})

                self.tq.update(    {str(j)+ '_' +str(n) : self.opti.variable(1)})
                self.tqdot.update( {str(j)+ '_' +str(n) : self.opti.variable(1)})
                
                self.u.update(     {str(j)+ '_' +str(n) : self.opti.variable(4)})

                self.lpos.update(  {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.dlpos.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.lforce.update({str(j)+ '_' +str(n) : self.opti.variable(2)})

                self.rpos.update(  {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.drpos.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.rforce.update({str(j)+ '_' +str(n) : self.opti.variable(2)})
                
                lq     = self.lq    [str(j)+ '_' +str(n)]  
                lqdot  = self.lqdot [str(j)+ '_' +str(n)]

                rq     = self.rq    [str(j)+ '_' +str(n)]
                rqdot  = self.rqdot [str(j)+ '_' +str(n)]

                tq     = self.tq    [str(j)+ '_' +str(n)]
                tqdot  = self.tqdot [str(j)+ '_' +str(n)]

                u      = self.u     [str(j)+ '_' +str(n)]

                lpos   = self.lpos  [str(j)+ '_' +str(n)]
                dlpos  = self.dlpos [str(j)+ '_' +str(n)]
                lforce = self.lforce[str(j)+ '_' +str(n)]

                rpos   = self.rpos  [str(j)+ '_' +str(n)]
                drpos  = self.drpos [str(j)+ '_' +str(n)]
                rforce = self.rforce[str(j)+ '_' +str(n)]

                self.model.setFullState(lq=lq, dlq=lqdot, rq=rq, drq=rqdot, 
                                        tq=tq, dtq=tqdot, 
                                        lp0=lpos, dlp0=dlpos, 
                                        rp0=rpos, drp0=drpos,
                                        u=u, lf10=lforce, rf10=rforce)

                ###--Add Model Constraints--###
                self.opti.subject_to(self.model.p  ['Constraint'])
                self.opti.subject_to(self.model.dp ['Constraint'])                                                 
                self.opti.subject_to(self.model.ddp['Constraint'])                                                 
                self.opti.subject_to(self.model.dynamics['Constraint'])

                # if j==0 and n==0:
                #     com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                #     self.opti.subject_to(com_x==0)
                
                # if j==self.num_phases-1 and n==self.num_phases-1:
                #     com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                #     self.opti.subject_to(com_x==2)
  
    def setConstraints(self):
        # pass
        self.setContactConstraints()
        # self.setCollocationContraints()

    def setContactConstraints(self):

        total_left  = 0
        total_right = 0
        
        ##--Left Leg--###
        j = 0
        for is_contact, del_T in self.time_phases['Left Leg'].items():
            
            if is_contact == 'No Contact @ ' + str(j):
                for n in range(self.knot_points_per_phase):
                    lforce = self.lforce[str(j)+ '_' +str(n)]
                    self.opti.subject_to(lforce==0)
            else: # Contact
                for knot_point in self.dlpos:
                    lpos   = self.lpos  [knot_point]
                    dlpos  = self.dlpos [knot_point]
                    lforce = self.lforce[knot_point]
                    
                    self.opti.subject_to((self.terrain.mu*lforce[0])**2 - lforce[1]**2 >= 0)
                    self.opti.subject_to(ca.dot(lforce,lpos) >= 0)
                    self.opti.subject_to(lpos[1]==self.terrain.heightMap(lpos[0]))
                    self.opti.subject_to(dlpos==0)

            total_left += del_T        
            j += 1

        self.opti.subject_to(total_left==self.total_duration)

        ##--Right Leg--###
        j = 0
        for is_contact, del_T in self.time_phases['Right Leg'].items():
            
            if is_contact == 'No Contact @ ' + str(j):
                for n in range(self.knot_points_per_phase):
                    rforce = self.rforce[str(j)+ '_' +str(n)]
                    self.opti.subject_to(rforce==0)
            else:
                for knot_point in self.drpos:
                    rpos   = self.rpos  [knot_point]
                    drpos  = self.drpos [knot_point]
                    rforce = self.rforce[knot_point]

                    self.opti.subject_to((self.terrain.mu*rforce[0])**2 - rforce[1]**2 >= 0)                    
                    self.opti.subject_to(ca.dot(rforce,rpos) >= 0)
                    self.opti.subject_to(rpos[1]==self.terrain.heightMap(rpos[0]))
                    self.opti.subject_to(drpos==0)                
            total_right += del_T        
            j += 1

        self.opti.subject_to(total_right==self.total_duration)

    def setCollocationContraints(self):
        for j in range(self.num_phases):
            
            if j%2 == 0:
                l_h = self.time_phases['Left Leg']['Contact @ ' + str(j)]
                r_h = self.time_phases['Right Leg']['Contact @ ' + str(j)]
            else:
                l_h = self.time_phases['Left Leg']['No Contact @ ' + str(j)]
                r_h = self.time_phases['Right Leg']['No Contact @ ' + str(j)]

            for n in range(self.knot_points_per_phase-1):
                lq1    , lq2     = self.lq    [str(j)+ '_' +str(n)], self.lq    [str(j)+ '_' +str(n+1)]
                lqdot1 , lqdot2  = self.lqdot [str(j)+ '_' +str(n)], self.lqdot [str(j)+ '_' +str(n+1)]
                lqddot1, lqddot2 = self.lqddot[str(j)+ '_' +str(n)], self.lqddot[str(j)+ '_' +str(n+1)]
                
                rq1    , rq2     = self.rq    [str(j)+ '_' +str(n)], self.rq    [str(j)+ '_' +str(n+1)]
                rqdot1 , rqdot2  = self.rqdot [str(j)+ '_' +str(n)], self.rqdot [str(j)+ '_' +str(n+1)]
                rqddot1, rqddot2 = self.rqddot[str(j)+ '_' +str(n)], self.rqddot[str(j)+ '_' +str(n+1)]
                
                tq1    , tq2     = self.tq    [str(j)+ '_' +str(n)], self.tq    [str(j)+ '_' +str(n+1)]
                tqdot1 , tqdot2  = self.tqdot [str(j)+ '_' +str(n)], self.tqdot [str(j)+ '_' +str(n+1)]
                tqddot1, tqddot2 = self.tqddot[str(j)+ '_' +str(n)], self.tqddot[str(j)+ '_' +str(n+1)]
                
                # u1     , u2      = self.u     [str(j)+ '_' +str(n)], self.u     [str(j)+ '_' +str(n+1)]
                lpos1  , lpos2   = self.lpos  [str(j)+ '_' +str(n)], self.lpos  [str(j)+ '_' +str(n+1)]
                dlpos1 , dlpos2  = self.dlpos [str(j)+ '_' +str(n)], self.dlpos [str(j)+ '_' +str(n+1)]
                ddlpos1, ddlpos2 = self.ddlpos[str(j)+ '_' +str(n)], self.ddlpos[str(j)+ '_' +str(n+1)]
                # lforce1, lforce2 = self.lforce[str(j)+ '_' +str(n)], self.lforce[str(j)+ '_' +str(n+1)]

                rpos1  , rpos2   = self.rpos  [str(j)+ '_' +str(n)], self.rpos  [str(j)+ '_' +str(n+1)]
                drpos1 , drpos2  = self.drpos [str(j)+ '_' +str(n)], self.drpos [str(j)+ '_' +str(n+1)]
                ddrpos1, ddrpos2 = self.ddrpos[str(j)+ '_' +str(n)], self.ddrpos[str(j)+ '_' +str(n+1)]
                # rforce1, rforce2 = self.rforce[str(j)+ '_' +str(n)], self.rforce[str(j)+ '_' +str(n+1)]

                
                self.opti.subject_to( (l_h/2) * (lqdot2  + lqdot1) == (lq2 - lq1) )
                self.opti.subject_to( (l_h/2) * (lqddot2 + lqddot1) == (lqdot2 - lqdot1) )

                # self.opti.subject_to( (l_h/2) * (dlpos2  + dlpos1) == (lpos2 - lpos1) )
                # self.opti.subject_to( (l_h/2) * (ddlpos2 + ddlpos1) == (dlpos2 - dlpos1) )

                # self.opti.subject_to( (r_h/2) * (rqdot2 + rqdot1) == (rq2 - rq1) )
                # self.opti.subject_to( (r_h/2) * (rqddot2 + rqddot1) == (rqdot2 - rqdot1) )

                # self.opti.subject_to( (r_h/2) * (drpos2  + drpos1) == (rpos2 - rpos1) )
                # self.opti.subject_to( (r_h/2) * (ddrpos2 + ddrpos1) == (drpos2 - drpos1) )

                # self.opti.subject_to( (l_h/2) * (tqdot2 + tqdot1) == (tq2 - tq1) )
                # self.opti.subject_to( (l_h/2) * (tqddot2 + tqddot1) == (tqdot2 - tqdot1) )



                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )
                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )

                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )
                # self.opti.subject_to( (l_h/2) * (lqdot2 + lqdot1) == (lq2 - lq1) )


# test check for sanity
# test_problem = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')            
# print((test_problem.time_phases))

#    def getConstraints(self):
        # for j in range(self.num_phases):
        #     for knot_point in self.q:
        #         q      = self.q     [knot_point]  
        #         qdot   = self.qdot  [knot_point]
        #         u      = self.u     [knot_point]
        #         lpos   = self.lpos  [knot_point]
        #         dlpos  = self.dlpos [knot_point]
        #         ddlpos = self.ddlpos[knot_point]
        #         lforce = self.lforce[knot_point]
        #         rpos   = self.rpos  [knot_point]
        #         drpos  = self.drpos [knot_point]
        #         ddrpos = self.ddrpos[knot_point]
        #         rforce = self.rforce[knot_point]

        #         self.model.setFullState(q=q, dq=qdot, lp0=lpos, dlp0=dlpos, ddlp0=ddlpos, 
        #                                               rp0=rpos, drp0=drpos, ddrp0=ddrpos, 
        #                                               u=u, lf10=lforce, rf10=rforce)

        #         self.opti.subject_to(self.model.p  ['Constraint'])
        #         self.opti.subject_to(self.model.dp ['Constraint'])                                                 
        #         self.opti.subject_to(self.model.ddp['Constraint'])                                                 

        #         self.opti.subject_to(self.model.dynamics['Constraint'])
        #         print(knot_point)
