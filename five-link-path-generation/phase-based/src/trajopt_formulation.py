import numpy as np
import casadi as ca

from five_link_model import Biped
from terrain import Terrain
# import helper_functions as hf 

class NLP():
    def __init__(self, knot_points_per_phase, steps, total_duration, model='biped', terrain='flat'):
        super().__init__()
        self.knot_points_per_phase = knot_points_per_phase
        self.num_phases            = steps
        self.total_duration        = total_duration

        if model == 'biped':
            self.model = Biped()
            self.time_phases = {'Left Leg' : {}, 'Right Leg': {}}
        
        self.terrain = Terrain(type=terrain)

        self.opti   = ca.Opti()

        self.q      = {}
        self.qdot   = {}
        self.u      = {}
        self.lpos   = {} # i = 0
        self.dlpos  = {} # i = 0
        self.ddlpos = {} # i = 0
        self.lforce = {} # i = 0
        self.rpos   = {} # i = 1
        self.drpos  = {} # i = 1
        self.ddrpos = {} # i = 1
        self.rforce = {} # i = 1
        
        for j in range(self.num_phases):
            
            for leg in self.time_phases:
                if j%2 == 0:
                    self.time_phases[leg].update({'Contact @ ' + str(j) : self.opti.variable(1)})
                else:
                    self.time_phases[leg].update({'No Contact @ ' + str(j): self.opti.variable(1)})

            for n in range(self.knot_points_per_phase):
                self.q.update(     {str(j)+ '_' +str(n) : self.opti.variable(5)})
                self.qdot.update(  {str(j)+ '_' +str(n) : self.opti.variable(5)})
                self.u.update(     {str(j)+ '_' +str(n) : self.opti.variable(4)})
                self.lpos.update(  {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.dlpos.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.ddlpos.update({str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.lforce.update({str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.rpos.update(  {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.drpos.update( {str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.ddrpos.update({str(j)+ '_' +str(n) : self.opti.variable(2)})
                self.rforce.update({str(j)+ '_' +str(n) : self.opti.variable(2)})
                
                q      = self.q     [str(j)+ '_' +str(n)]  
                qdot   = self.qdot  [str(j)+ '_' +str(n)]
                u      = self.u     [str(j)+ '_' +str(n)]
                lpos   = self.lpos  [str(j)+ '_' +str(n)]
                dlpos  = self.dlpos [str(j)+ '_' +str(n)]
                ddlpos = self.ddlpos[str(j)+ '_' +str(n)]
                lforce = self.lforce[str(j)+ '_' +str(n)]
                rpos   = self.rpos  [str(j)+ '_' +str(n)]
                drpos  = self.drpos [str(j)+ '_' +str(n)]
                ddrpos = self.ddrpos[str(j)+ '_' +str(n)]
                rforce = self.rforce[str(j)+ '_' +str(n)]

                self.model.setFullState(q=q, dq=qdot, lp0=lpos, dlp0=dlpos, ddlp0=ddlpos, 
                                                      rp0=rpos, drp0=drpos, ddrp0=ddrpos, 
                                                      u=u, lf10=lforce, rf10=rforce)

                ###--Add Model Constraints--###
                self.opti.subject_to(self.model.p  ['Constraint'])
                self.opti.subject_to(self.model.dp ['Constraint'])                                                 
                self.opti.subject_to(self.model.ddp['Constraint'])                                                 

                self.opti.subject_to(self.model.dynamics['Constraint'])

                if j==0 and n==0:
                    com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                    self.opti.subject_to(com_x==0)
                
                if j==self.num_phases-1 and n==self.num_phases-1:
                    com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                    self.opti.subject_to(com_x==2)

        self.getConstraints()
        
    def getConstraints(self):
        self.setContactConstraints()
        self.setCollocationContraints()

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
                    lpos  = self.lpos   [knot_point]
                    dlpos = self.dlpos  [knot_point]
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
                    rforce = self.lforce[str(j)+ '_' +str(n)]
                    self.opti.subject_to(rforce==0)
            else:
                for knot_point in self.drpos:
                    rpos  = self.rpos   [knot_point]
                    drpos = self.drpos  [knot_point]
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
            for n in range(self.knot_points_per_phase-1):
                q1     , q2      = self.q     [str(j)+ '_' +str(n)], self.q     [str(j)+ '_' +str(n+1)]
                qdot1  , qdot2   = self.qdot  [str(j)+ '_' +str(n)], self.qdot  [str(j)+ '_' +str(n+1)]
                u1     , u2      = self.u     [str(j)+ '_' +str(n)], self.u     [str(j)+ '_' +str(n+1)]
                lpos1  , lpos2   = self.lpos  [str(j)+ '_' +str(n)], self.lpos  [str(j)+ '_' +str(n+1)]
                dlpos1 , dlpos2  = self.dlpos [str(j)+ '_' +str(n)], self.dlpos [str(j)+ '_' +str(n+1)]
                ddlpos1, ddlpos2 = self.ddlpos[str(j)+ '_' +str(n)], self.ddlpos[str(j)+ '_' +str(n+1)]
                lforce1, lforce2 = self.lforce[str(j)+ '_' +str(n)], self.lforce[str(j)+ '_' +str(n+1)]
                rpos1  , rpos2   = self.rpos  [str(j)+ '_' +str(n)], self.rpos  [str(j)+ '_' +str(n+1)]
                drpos1 , drpos2  = self.drpos [str(j)+ '_' +str(n)], self.drpos [str(j)+ '_' +str(n+1)]
                ddrpos1, ddrpos2 = self.ddrpos[str(j)+ '_' +str(n)], self.ddrpos[str(j)+ '_' +str(n+1)]
                rforce1, rforce2 = self.rforce[str(j)+ '_' +str(n)], self.rforce[str(j)+ '_' +str(n+1)]

                self.opti.subject_to()
# test check for sanity
test_problem = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')            
# print(len(test_problem.q))

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
