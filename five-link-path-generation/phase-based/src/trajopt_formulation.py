import numpy as np
import casadi as ca
from five_link_model import Biped

class NLP():
    def __init__(self, knot_points_per_phase, steps, total_duration, model='biped'):
        super().__init__()
        self.knot_points_per_phase = knot_points_per_phase
        self.num_phases            = steps
        self.total_duration        = total_duration

        if model == 'biped':
            self.model = Biped()
            self.time_phases = {'Left Leg' : {}, 'Right Leg': {}}
        
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

                ###--Add Model Constraint--###
                self.opti.subject_to(self.model.p  ['Constraint'])
                self.opti.subject_to(self.model.dp ['Constraint'])                                                 
                self.opti.subject_to(self.model.ddp['Constraint'])                                                 

                self.opti.subject_to(self.model.dynamics['Constraint'])

        # self.getConstraints()

    def getConstraints(self):
        for j in range(self.num_phases):
            for knot_point in self.q:
                q      = self.q     [knot_point]  
                qdot   = self.qdot  [knot_point]
                u      = self.u     [knot_point]
                lpos   = self.lpos  [knot_point]
                dlpos  = self.dlpos [knot_point]
                ddlpos = self.ddlpos[knot_point]
                lforce = self.lforce[knot_point]
                rpos   = self.rpos  [knot_point]
                drpos  = self.drpos [knot_point]
                ddrpos = self.ddrpos[knot_point]
                rforce = self.rforce[knot_point]

                self.model.setFullState(q=q, dq=qdot, lp0=lpos, dlp0=dlpos, ddlp0=ddlpos, 
                                                      rp0=rpos, drp0=drpos, ddrp0=ddrpos, 
                                                      u=u, lf10=lforce, rf10=rforce)

                self.opti.subject_to(self.model.p  ['Constraint'])
                self.opti.subject_to(self.model.dp ['Constraint'])                                                 
                self.opti.subject_to(self.model.ddp['Constraint'])                                                 

                self.opti.subject_to(self.model.dynamics['Constraint'])
                print(knot_point)

# test check for sanity

test_problem = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')            
# print(len(test_problem.q))