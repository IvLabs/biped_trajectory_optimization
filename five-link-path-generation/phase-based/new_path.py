from new_gait import *
from matplotlib import pyplot as plt
from time import sleep
# start_angles = ca.MX([-0.3,0.3,0.2,-0.3,0.3])
# start_angles = ca.MX([-0.2,0.1,0.2,-0.5,1.5]) # ostrich 2
# start_angles = ca.MX([-0.6,0.2,0.0,-0.7,0.9]) # ostrich 3
start_angles = ca.MX.zeros(5)
# start_angles = ca.MX([-0.6,0.7,0.0,-0.5,0.9]) # ostrich
# start_angles = ca.MX([0.3,-0.3,0.1,-0.3,-0.45]) # human

start_pos_left = [[0,0]]
start_pos_right = [[0,0]]

start_angular_vel = ca.MX.zeros(5)
q = []; dq = []; u = []; pos = []; time = []
f = 1
for k in range(f):
    model = walker(start_angles, start_angular_vel, start_pos_left[-1], start_pos_right[-1])  
    # exit()
    problem = nlp(model)
    sol = model.opti.solve_limited()

    for j in range(5):
        tempq = [];tempdq = [];tempu = []
        for i in range(model.N):
            tempq.append(sol.value(model.x[i][j]))    
            tempdq.append(sol.value(model.xdot[i][j]))
            if j < 4:
                tempu.append(sol.value(model.u[i][j]))
            if j==0:
                print(sol.value(model.l_force[i]))
                print(sol.value(model.r_force[i]))
        if k > 0:
            q[j].extend(tempq)
            dq[j].extend(tempdq)
            if j < 4:    
                u[j].extend(tempu)
        else : 
            q.append(tempq)
            dq.append(tempdq)
            u.append(tempu)

    for j in range(6):
        temp = []
        for i in range(model.N):
            temp.append([sol.value(model.pos[i][j, 0]),sol.value(model.pos[i][j, 1])]) 
            if i==0:
                print('p'+str(j)+'_'+str(i)+' == '+str(temp[-1]))
        pos.append(temp)   
    
    start_angles = ca.MX([q[0][-1],q[1][-1],q[2][-1],q[3][-1],q[4][-1]])
    start_angular_vel = ca.MX([dq[0][-1],dq[1][-1],dq[2][-1],dq[3][-1],dq[4][-1]])

    print('##################')
    print('step number = ', k)
    print('step time = ', model.T)
    print('##################')

    start_pos_left.append([pos[0][-1][0], pos[0][-1][1]])
    start_pos_right.append([pos[5][-1][0], pos[5][-1][1]])

timodel = np.linspace(0.0, f*model.T, len(q[:][0]))

from matplotlib import animation
from celluloid import Camera

fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
camodelra = Camera(fig)
terrain = np.linspace(-2,2,1000*f)
if model.terrain == 'sin':
    terrain_y = model.terrain_factor*ca.sin(terrain)
elif model.terrain == 'wedge':
    terrain_y = model.terrain_factor*terrain
elif model.terrain == 'smooth_stair':
    k = -50
    terrain_y = terrain*k - ca.sin(terrain*k) - ca.sin(terrain*k - ca.sin(terrain*k)) - ca.sin(terrain*k - ca.sin(terrain*k) - ca.sin(terrain*k - ca.sin(terrain*k))) - ca.sin(terrain*k - ca.sin(terrain*k) - ca.sin(terrain*k - ca.sin(terrain*k)) - ca.sin(terrain*k - ca.sin(terrain*k) - ca.sin(terrain*k - ca.sin(terrain*k))))
    terrain_y /= abs(k)
k = 0
ax.set_xlim([-1., 7]) # sin
# ax.set_xlim([-1., 5]) # wedge
# ax.set_ylim([-1, 5]) # wedge
ax.set_ylim([-1, 3]) # sin


# ax.set_ylim([-1, 2]) # ss
# ax.set_xlim([-1., 2]) # ss

for i in range(f*model.N):
    # print(i)    
    p0 = [pos[0][i][0],pos[0][i][1]]
    p1 = [pos[1][i][0],pos[1][i][1]]
    p2 = [pos[2][i][0],pos[2][i][1]]
    p3 = [pos[3][i][0],pos[3][i][1]]
    p4 = [pos[4][i][0],pos[4][i][1]]
    p5 = [pos[5][i][0],pos[5][i][1]]
    # plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    # plt.axes(xlim=(-2., 2.), ylim=(-0.1, 3))
    # plt.plot([0,-p1[1]], [0,p1[0]],'r',[-p1[1],-p2[1]], [p1[0],p2[0]],'b',
    #     [-p2[1],-p3[1]], [p2[0],p3[0]],'c',
    #     [-p2[1],p4[1] - 2*p2[1]], [p2[0],2*p2[0]-p4[0]],'b',
    #     [p4[1] - 2*p2[1],p5[1]], [2*p2[0]-p4[0],(p5[0] - 2*p2[0])],'r')
    

    # plt.plot([p0[0],p1[1]], [p0[1],p1[0]],'r',[p1[1],p2[1]], [p1[0],p2[0]],'g',
    #         [p2[1],p3[1]], [p2[0],p3[0]],'b', [p2[1],p4[1]], [p2[0],p4[0]],'y',
    #         [p4[1],p5[1]], [p4[0],p5[0]],'c')
    # ax.set_xlim([-2.+p2[0], 2.+p2[0]])
    # ax.set_ylim([-2+p2[1], 2+p2[1]])
    ax.plot([p0[0],p1[0]], [p0[1],p1[1]],'r',[p1[0],p2[0]], [p1[1],p2[1]],'g',
        [p2[0],p3[0]], [p2[1],p3[1]],'b', [p2[0],p4[0]], [p2[1],p4[1]],'y',
        [p4[0],p5[0]], [p4[1],p5[1]],'c')

    # plt.plot([-2,6],[0,0],'g')  

    if model.terrain == 'sin':
        ax.plot(terrain, terrain_y,'black')   # sin
    if model.terrain == 'wedge':
        ax.plot(terrain, terrain_y,'black')   # wedge 
    if model.terrain == 'smooth_stair':
        ax.plot(terrain, terrain_y,'black')   # smooth stair 


    camodelra.snap()
    ax.grid()

    # plt.draw() 
    # plt.pause(1e-5)
    # ax.cla()
animation = camodelra.animate(interval=60)
# animation.save('path_sstairs_down_N_40_human.mp4')
plt.show()


# print(model.heightMapNormalVector(model.pos[0][0, 0], model.pos[0][1, 0]))

namodel = ['q','dq','u']

# plt.plot(test['q1'])
# plt.plot(test['q2'])
# plt.plot(test['q3'])
# plt.plot(test['q4'])

# plt.show()
plt.subplot(311)
plt.title('Optimised Solution')
plt.plot(timodel, q[:][0],'r', timodel, q[:][1],'g', timodel,q[:][2],'b',
        timodel,q[:][3],'y', timodel, q[:][4],'c')
plt.grid()

# plt.subplot(321)
# plt.title('Initial Guess')
# iniq = n1.initial[0]
# plt.plot(timodel,np.append(iniq[:][0],iniq[:][0]),'r',timodel,np.append(iniq[:][1],iniq[:][1]),'g',timodel,np.append(iniq[:][2],iniq[:][2]),'b',
#         timodel,np.append(iniq[:][3],iniq[:][3]),'y',timodel,np.append(iniq[:][4],iniq[:][4]),'c')
plt.ylabel(namodel[0])

plt.subplot(312)
plt.plot(timodel,dq[:][0],'r',timodel,dq[:][1],'g',timodel,dq[:][2],'b',
        timodel,dq[:][3],'y',timodel,dq[:][4],'c')
plt.grid()

# plt.subplot(323)
# inidq = n1.initial[1]
# plt.plot(timodel,np.append(inidq[:][0],inidq[:][0]),'r',timodel,np.append(inidq[:][1],inidq[:][1]),'g',timodel,np.append(inidq[:][2],inidq[:][2]),'b',
#         timodel,np.append(inidq[:][3],inidq[:][3]),'y',timodel,np.append(inidq[:][4],inidq[:][4]),'c')
plt.ylabel(namodel[1])

plt.subplot(313)
plt.plot(timodel,u[:][0],'g',timodel,u[:][1],'b',timodel,u[:][2],'y',
        timodel,u[:][3],'c')
plt.grid()

# plt.subplot(325)
# iniu = n1.initial[2]
# plt.plot(timodel,np.append(iniu[:][0],iniu[:][0]),'r',timodel,np.append(iniu[:][1],iniu[:][1]),'g',timodel,np.append(iniu[:][2],iniu[:][2]),'b',
#         timodel,np.append(iniu[:][3],iniu[:][3]),'y',timodel,np.append(iniu[:][4],iniu[:][4]),'c')
plt.ylabel(namodel[2])

plt.suptitle('Five-Link')
plt.show()