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
        self.opti.subject_to(self.p0[str(self.num_phases-1)][-1][0] == 1)

    def setContactConstraints(self):

        self.opti.subject_to(sum(self.time_phases.values()) == self.total_duration)

        for step in range(self.num_phases):
            
            # if step is even ==> Contact
            # if step is odd  ==> No Contact
            
            if step%2 != 0:
                for knot_point in range(self.knot_points_per_phase):
                    f10 = self.f10[str(step)][knot_point]
                    self.opti.subject_to(f10==0)
            else: 
                for knot_point in range(self.knot_points_per_phase):
                    p0  = self.p0 [str(step)][knot_point]
                    dp0 = self.dp0[str(step)][knot_point]
                    f10 = self.f10[str(step)][knot_point]
                    
                    self.opti.subject_to((self.terrain.mu*f10[0])**2 - f10[1]**2 >= 0)
                    self.opti.subject_to(ca.dot(f10,p0) >= 0)
                    self.opti.subject_to(p0[1]==self.terrain.heightMap(p0[0]))
                    self.opti.subject_to(dp0==0)

    def setCollocationContraints(self):
        for step in range(self.num_phases):
            
            h = self.time_phases[str(step)]/self.knot_points_per_phase

            for n in range(self.knot_points_per_phase-1):
                q_1   , q_2    = self.q   [str(step)][n], self.q   [str(step)][n+1]
                qdot_1, qdot_2 = self.qdot[str(step)][n], self.qdot[str(step)][n+1]
                
                p0_1 , p0_2  = self.p0 [str(step)][n], self.p0 [str(step)][n+1]
                dp0_1, dp0_2 = self.dp0[str(step)][n], self.dp0[str(step)][n+1]
                
                u_1  , u_2   = self.u  [str(step)][n], self.u  [str(step)][n+1]
                f10_1, f10_2 = self.f10[str(step)][n], self.f10[str(step)][n+1]

                m_1 = Hopper()
                m_1.setFullState(q=q_1, dq=qdot_1, 
                                 p0=p0_1, dp0=dp0_1, 
                                 u=u_1, f10=f10_1)
                
                qddot_1 = m_1.dynamics['Leg']

                m_2 = Hopper()
                m_2.setFullState(q=q_2, dq=qdot_2, 
                                 p0=p0_2, dp0=dp0_2, 
                                 u=u_2, f10=f10_2)
                
                qddot_2 = m_2.dynamics['Leg']

                
                self.opti.subject_to( (h/2) * (qdot_2  + qdot_1)  == (q_2 - q_1) )
                self.opti.subject_to( (h/2) * (qddot_2 + qddot_1) == (qdot_2 - qdot_1) )

                self.opti.subject_to( (h/2) * (dp0_2  + dp0_1) == (p0_2 - dp0_1) )


# test check for sanity
test_problem = NLP(knot_points_per_phase=40, steps=3, total_duration=2, model='hopper')            
print((test_problem.time_phases))
