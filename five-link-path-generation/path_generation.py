from matplotlib import pyplot as plt
from gait_generation import *

start = [-0.3,0.7,0.0,-0.5,-0.3]

me = walker(start)
problem = nlp(me)
sol1 = me.opti.solve()

q = []; dq = []; u = []; pos = []
for j in range(5):
    tempq = [];tempdq = [];tempu = [];temp = []
    for i in range(me.N):
        tempq.append(sol1.value(me.state[i][0][j]))    
        tempdq.append(sol1.value(me.state[i][1][j]))
        temp.append([sol1.value(me.pos[i][j][0]),sol1.value(me.pos[i][j][1])]) 
        if j < 4:
            tempu.append(sol1.value(me.u[i][0][j]))
    q.append(tempq);pos.append(temp)        
    dq.append(tempdq)        
    u.append(tempu)        

time = np.arange(0.0, me.T, me.h)

print(q[0][-1],q[1][-1],q[2][-1],q[3][-1],q[4][-1])




