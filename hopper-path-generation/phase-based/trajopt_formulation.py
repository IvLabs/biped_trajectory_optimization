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

        self.ceq = []
        self.ciq = []

        self.setVariables()
        self.setConstraints()
        self.setBounds()

    def setVariables(self):    
        for step in range(self.num_phases):
            
            self.time_phases.update({str(step) : self.opti.variable(1)})
            self.ceq.append(0<=self.time_phases[str(step)])

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
                # self.opti.subject_to(self.model.p       ['Constraint'])                                                 
                # self.opti.subject_to(self.model.dynamics['Constraint'])
                
                self.ceq.append(self.model.p['Constraint'])
                # self.ceq.append(self.model.dynamics['Constraint'])

                # self.opti.minimize(p0[1]**2)
                # self.opti.minimize(ca.sumsqr(self.model.p['Leg'][1,:]))

                for i in range(3):
                    self.opti.subject_to((self.model.p['Leg'][1,i]) >= self.terrain.heightMap(self.model.p['Leg'][0,:]))                                                 

                # if j==0 and n==0:
                #     self.setInitialBoundary()
                #     com_x = (self.model.c['Left Leg'][0,:]) + (self.model.c['Right Leg'][0,:]) + (self.model.c['Torso'][0,:]) 
                #     self.opti.subject_to(com_x==0)
                
                # if j==self.num_phases-1 and n==self.num_phases-1:
                # com_y = ca.sum1(self.model.c['Leg'][1,:]) 
                # self.opti.subject_to(com_y<=self.model.b)
        
    def setConstraints(self):
        # self.setBoundary()
        self.setContactConstraints()
        self.setCollocationContraints()

        self.opti.subject_to(self.ceq)

    # def setBoundary(self):
        # self.opti.subject_to(self.q [str(0)][0] == self.model.initial_q)
        # self.opti.subject_to(self.p0[str(0)][0] == np.zeros((2,1)))
        # self.opti.subject_to(self.q[str(self.num_phases-1)][-1] == self.model.final_q)
        # self.opti.subject_to(self.p0[str(self.num_phases-1)][-1][0] == 1)
        
        # self.ceq.append(self.q [str(0)][0] == self.model.initial_q)
        # self.ceq.append(self.p0[str(0)][0] == np.zeros((2,1)))

    def setContactConstraints(self):

        # self.opti.subject_to(sum(self.time_phases.values()) == self.total_duration)
        self.ceq.append(sum(self.time_phases.values()) == self.total_duration)

        for step in range(self.num_phases):
            
            # if step is even ==> Contact
            # if step is odd  ==> No Contact
            for knot_point in range(self.knot_points_per_phase):
                
                # if step > 0:
                #     # self.opti.subject_to(self.p0 [str(step-1)][-1] - self. p0[str(step)][0] == 0)
                #     # self.opti.subject_to(self.f10[str(step-1)][-1] - self.f10[str(step)][0] == 0)
                    
                #     self.ceq.append(self.p0 [str(step-1)][-1] - self. p0[str(step)][0] == 0)
                #     self.ceq.append(self.f10[str(step-1)][-1] - self.f10[str(step)][0] == 0)

                if step%2 != 0:
                    f10 = self.f10[str(step)][knot_point]
                    p0  = self.p0 [str(step)][knot_point]

                    # self.opti.bounded(0,self.terrain.heightMap(p0[0])-p0[1],ca.inf)
                    # self.opti.subject_to(f10==0)
                    
                    self.ceq.append(self.terrain.heightMap(p0[0])<=p0[1])   
                    self.ceq.append(f10==0)
                else: 
                    p0  = self.p0 [str(step)][knot_point]
                    dp0 = self.dp0[str(step)][knot_point]
                    f10 = self.f10[str(step)][knot_point]
                    
                    # self.opti.subject_to((self.terrain.mu*f10[0])**2 - f10[1]**2 >= 0)
                    # self.opti.subject_to(ca.dot(f10,p0) >= 0)
                    # self.opti.subject_to(p0[1]==self.terrain.heightMap(p0[0]))
                    # self.opti.subject_to(dp0==0)

                    self.ceq.append((self.terrain.mu*f10[0])**2 - f10[1]**2 >= 0)
                    self.ceq.append(ca.dot(f10,p0) >= 0)
                    self.ceq.append(p0[1]==self.terrain.heightMap(p0[0]))
                    self.ceq.append(dp0==0)

    def setCollocationContraints(self):
        for step in range(self.num_phases):
            
            h = self.time_phases[str(step)]/self.knot_points_per_phase    
            # h = 0.01
            for n in range(self.knot_points_per_phase-1):

                if step>0 and n==0:
                    q_1   , q_2    = self.q   [str(step-1)][-1], self.q   [str(step)][n]
                    qdot_1, qdot_2 = self.qdot[str(step-1)][-1], self.qdot[str(step)][n]
                    
                    p0_1 , p0_2  = self.p0 [str(step-1)][-1], self.p0 [str(step)][n]
                    dp0_1, dp0_2 = self.dp0[str(step-1)][-1], self.dp0[str(step)][n]
                    
                    u_1  , u_2   = self.u  [str(step-1)][-1], self.u  [str(step)][n]
                    f10_1, f10_2 = self.f10[str(step-1)][-1], self.f10[str(step)][n]

                    # A = ca.MX.zeros(6)
                    # A[0,0] = q_2 - q_1;A[1,0] = qdot_2 - qdot_1;A[2,0] = p0_2-p0_1
                    # A[3,0] = dp0_2-dp0_1;A[4,0] = u_2-u_1;A[5,0] = f10_2-f10_1
    
                    self.ceq.append(q_2 - q_1== 0);self.ceq.append(qdot_2 - qdot_1== 0)
                    self.ceq.append(p0_2-p0_1== 0);self.ceq.append(dp0_2-dp0_1== 0)
                    self.ceq.append(u_2-u_1== 0);self.ceq.append(f10_2-f10_1== 0)


                    # m_1 = Hopper()
                    # m_1.setFullState(q=q_1, dq=qdot_1, 
                    #                 p0=p0_1, dp0=dp0_1, 
                    #                 u=u_1, f10=f10_1)
                    
                    # qddot_1 = m_1.dynamics['Leg']

                    # m_2 = Hopper()
                    # m_2.setFullState(q=q_2, dq=qdot_2, 
                    #                 p0=p0_2, dp0=dp0_2, 
                    #                 u=u_2, f10=f10_2)
                    
                    # qddot_2 = m_2.dynamics['Leg']

                    # self.ceq.append((h/2) * (qdot_2  + qdot_1)  == (q_2 - q_1))
                    # self.ceq.append((h/2) * (qddot_2 + qddot_1) == (qdot_2 - qdot_1))
                    # self.ceq.append((h/2) * (dp0_2  + dp0_1) == (p0_2 - dp0_1))
                    # self.ceq.append((h)*(self.model.gravity) == (dp0_2[1] - dp0_1[1]))
                # else:
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

                self.ceq.append((h/2) * (qdot_2  + qdot_1)  == (q_2 - q_1))
                self.ceq.append((h/2) * (qddot_2 + qddot_1) == (qdot_2 - qdot_1))
                self.ceq.append((h/2) * (dp0_2  + dp0_1) == (p0_2 - dp0_1))
                self.ceq.append((h)*(self.model.gravity) == (dp0_2[1] - dp0_1[1]))

    def setBounds(self):
        # traj=0
        for step in range(self.num_phases):
            for knot_point in range(self.knot_points_per_phase):
                q    = self.q   [str(step)][knot_point]
                qdot = self.qdot[str(step)][knot_point]
                
                p0  = self.p0 [str(step)][knot_point]
                dp0 = self.dp0[str(step)][knot_point]
                
                u   = self.u  [str(step)][knot_point]
                f10 = self.f10[str(step)][knot_point]

                self.opti.bounded([-np.pi/2]*3, q,[np.pi/2]*3)
                self.opti.bounded([-1.5]*2, u,[1.5]*2)
                
                # traj+=(p0[1]**2)
                

        # self.opti.minimize(traj)
# test check for sanity
# test_problem = NLP(knot_points_per_phase=40, steps=3, total_duration=2, model='hopper')            
# print((test_problem.time_phases))
