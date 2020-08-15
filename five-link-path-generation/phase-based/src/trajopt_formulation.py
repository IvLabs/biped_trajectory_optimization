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
        
        for leg in self.time_phases:
            for j in range(self.num_phases):
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
                    

# test check for sanity

test_problem = NLP(knot_points_per_phase=40, steps=4, total_duration=3, model='biped')            
