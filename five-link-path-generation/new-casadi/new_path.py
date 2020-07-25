from new_gait import *
from matplotlib import pyplot as plt
from time import sleep
# start_angles = ca.MX([-0.3,0.3,0.2,-0.3,0.3])
# start_angles = ca.MX([-0.2,0.1,0.2,-0.5,1.5]) # ostrich 2
start_angles = ca.MX([-0.6,0.2,0.0,-0.7,0.9]) # ostrich 3
# start_angles = ca.MX.zeros(5)
# start_angles = ca.MX([-0.6,0.7,0.0,-0.5,0.9]) # ostrich
# start_angles = ca.MX([0.3,-0.3,0.1,-0.3,-0.45]) # human

start_pos = [[0,0]]
start_angular_vel = ca.MX.zeros(5)
q = []; dq = []; u = []; pos = []; time = []

f = 6
for k in range(f):
    # try:
    model = walker(start_angles, start_angular_vel, start_pos[-1])        
    problem = nlp(model)
    sol = model.opti.solve()
    # except:
    #     pass
    # p = [0, 0, 0, 0, 0]
    for j in range(5):
        # tq = []; tdq = []; tu = []; tpos = []
        tempq = [];tempdq = [];tempu = [];temp = []
        for i in range(model.N):
            tempq.append(sol.value(model.x[i][j]))    
            tempdq.append(sol.value(model.x[i][j]))
            temp.append([sol.value(model.pos[i][j, 0]),sol.value(model.pos[i][j, 1])]) 

            # if k != 0:
            #     p0x, p0y = start_pos[-1][0].to_DM(), start_pos[-1][1].to_DM()
            # else:
            #     p0x, p0y = start_pos[-1][0], start_pos[-1][1]
            # p[0] = np.array([tempq[-1]]),  np.cos(tempq[-1]) + np.full((2,), [p0x, p0y])
            # p[1] = np.array([tempq[-1]]),  np.cos(tempq[-1]) + p[0]
            # p[2] = np.array([tempq[-1]]),  np.cos(tempq[-1]) + p[1]
            # p[3] = np.array([tempq[-1]]), -np.cos(tempq[-1]) + p[1]
            # p[4] = np.array([tempq[-1]]), -np.cos(tempq[-1]) + p[3]

            # temp.append([p[j][0],p[j][1]]) 
            if j < 4:
                tempu.append(sol.value(model.u[i][j]))

        # if k%2 != 0:
        #     tempq = tempq[::-1] 
        #     tempdq = tempdq[::-1]
        #     tempu = tempu[::-1]
        #     # temp = temp[::-1]
        # test['q' + str(k)] = tempq
        if k > 0:
            q[j].extend(tempq)
            dq[j].extend(tempdq)
            pos[j].extend(temp)
            if j < 4:    
                u[j].extend(tempu)
            # print(j)    
        else : 
            q.append(tempq)
            dq.append(tempdq)
            pos.append(temp)   
            u.append(tempu)

    # start = [q[4][-1],q[3][-1],q[2][-1],q[1][-1],q[0][-1]]
    # print(-sol.value(model.x_impact))
    start_angles, start_angular_vel = ca.MX(-sol.value(model.x_impact)), ca.MX(-sol.value(model.xdot_impact))
    # start_angles = ca.MX([q[4][-1],q[3][-1],q[2][-1],q[1][-1],q[0][-1]])
    # start_angular_vel = ca.MX([dq[4][-1],dq[3][-1],dq[2][-1],dq[1][-1],dq[0][-1]])
    # start_angles = ca.MX([q[0][-1],q[1][-1],q[2][-1],q[3][-1],q[4][-1]])
    # start_angular_vel = ca.MX([dq[0][-1],dq[1][-1],dq[2][-1],dq[3][-1],dq[4][-1]])
    # start_angles = ca.MX(sol.value(model.x[-1])[::-1])
    # start_angular_vel = ca.MX(sol.value(model.xdot[-1])[::-1])

    start_pos.append([pos[4][-1][0], pos[4][-1][1]])

#     if k > 1:
        # print(start_pos[k][1])

timodel = np.arange(0.0, f*model.T, model.h)

from matplotlib import animation
from celluloid import Camera

# fig = plt.figure()
# ax = fig.add_subplot(111, autoscale_on=False)

# ax.grid()

# time_template = 'time = %.1fs'
# time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
# p1, = ax.plot([], [], '-', color='r')
# p2, = ax.plot([], [], '-', color='g')
# p3, = ax.plot([], [], '-', color='b')
# p4, = ax.plot([], [], '-', color='y')
# p5, = ax.plot([], [], '-', color='c')

# def init():
#     p1.set_data([], [])
#     p2.set_data([], [])
#     p3.set_data([], [])
#     p4.set_data([], [])
#     p5.set_data([], [])
#     time_text.set_text('')
#     return p1, p2, p3, p4, p5, time_text

# k = 0
# def animate(i):
#     global k
#     if i%model.N == 0:
#         k += 1
#         # print(k)
#     if i == len(q[:])+1:
#         k = 0    
#     ax.set_xlim([-2.+pos[4][i][0], 2.+pos[4][i][0]])
#     ax.set_ylim([-1+pos[4][i][1], 3+pos[4][i][1]])
#     p0 = start_pos[k]
#     p1x = [       p0[0],pos[0][i][0]]
#     p2x = [pos[0][i][0],pos[1][i][0]]
#     p3x = [pos[1][i][0],pos[2][i][0]]
#     p4x = [pos[1][i][0],pos[3][i][0]]
#     p5x = [pos[3][i][0],pos[4][i][0]]

#     p1y = [       p0[1],pos[0][i][1]]
#     p2y = [pos[0][i][1],pos[1][i][1]]
#     p3y = [pos[1][i][1],pos[2][i][1]]
#     p4y = [pos[1][i][1],pos[3][i][1]]
#     p5y = [pos[3][i][1],pos[4][i][1]]

    
#     p1.set_data(p1x, p1y)
#     p2.set_data(p2x, p2y)
#     p3.set_data(p3x, p3y)
#     p4.set_data(p4x, p4y)
#     p5.set_data(p5x, p5y)
#     time_text.set_text(time_template % (i*20/len(timodel)))

#     return p1, p2, p3, p4, p5, time_text

# ani = animation.FuncAnimation(fig, animate, np.arange(0, len(q[:][0])), init_func=init,
#                                interval=20, blit=True)

# # ani.save('spoiler.mp4')
# plt.show()

# print(len(pos[0]))

# from matplotlib.animation import FuncAnimation
# plt.style.use('seaborn-pastel')


# fig = plt.figure()
# ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
# line, = ax.plot([], [], lw=3)

# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
#     x = np.linspace(0, 4, 1000)
#     y = np.sin(2 * np.pi * (x - 0.01 * i))
#     line.set_data(x, y)
#     return line,

# anim = FuncAnimation(fig, animate, init_func=init,
#                                frames=200, interval=20, blit=True)


# anim.save('sine_wave.gif', writer='imagemagick')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
camodelra = Camera(fig)
k = 0
terrain = np.linspace(-2,f*6,int(f*100/2))
ax.set_xlim([-2., (f/2)*2])
ax.set_ylim([-2, 2])

for i in range(f*model.N):
    # print(i)    
    p1 = [pos[0][i][0],pos[0][i][1]]
    p2 = [pos[1][i][0],pos[1][i][1]]
    p3 = [pos[2][i][0],pos[2][i][1]]
    p4 = [pos[3][i][0],pos[3][i][1]]
    p5 = [pos[4][i][0],pos[4][i][1]]
    # plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    # plt.axes(xlim=(-2., 2.), ylim=(-0.1, 3))
    # plt.plot([0,-p1[1]], [0,p1[0]],'r',[-p1[1],-p2[1]], [p1[0],p2[0]],'b',
    #     [-p2[1],-p3[1]], [p2[0],p3[0]],'c',
    #     [-p2[1],p4[1] - 2*p2[1]], [p2[0],2*p2[0]-p4[0]],'b',
    #     [p4[1] - 2*p2[1],p5[1]], [2*p2[0]-p4[0],(p5[0] - 2*p2[0])],'r')
    if i%model.N == 0:
            p0 = start_pos[k]
            print(p0)
            k += 1 

    # plt.plot([p0[0],p1[1]], [p0[1],p1[0]],'r',[p1[1],p2[1]], [p1[0],p2[0]],'g',
    #         [p2[1],p3[1]], [p2[0],p3[0]],'b', [p2[1],p4[1]], [p2[0],p4[0]],'y',
    #         [p4[1],p5[1]], [p4[0],p5[0]],'c')
    # ax.set_xlim([-2.+p2[0], 2.+p2[0]])
    # ax.set_ylim([-4+p2[1], 4+p2[1]])
    ax.plot([p0[0],p1[0]], [p0[1],p1[1]],'r',[p1[0],p2[0]], [p1[1],p2[1]],'g',
        [p2[0],p3[0]], [p2[1],p3[1]],'b', [p2[0],p4[0]], [p2[1],p4[1]],'y',
        [p4[0],p5[0]], [p4[1],p5[1]],'c')
    # plt.plot([-2,6],[0,0],'g')  
    ax.plot(terrain, model.terrain_factor*np.sin(terrain),'g')    
    camodelra.snap()
    # ax.grid()

    # plt.draw() 
    # plt.pause(1e-1) 
    # if cv2.waitKey(0) & 0xFF == ord("q"):
    # #     break
    # ax.cla()
animation = camodelra.animate(interval=30)
# animation.save('path_test_human_2.mp4')
plt.show()
plt.close()

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

# plt.subplot(321)
# plt.title('Initial Guess')
# iniq = n1.initial[0]
# plt.plot(timodel,np.append(iniq[:][0],iniq[:][0]),'r',timodel,np.append(iniq[:][1],iniq[:][1]),'g',timodel,np.append(iniq[:][2],iniq[:][2]),'b',
#         timodel,np.append(iniq[:][3],iniq[:][3]),'y',timodel,np.append(iniq[:][4],iniq[:][4]),'c')
plt.ylabel(namodel[0])

plt.subplot(312)
plt.plot(timodel,dq[:][0],'r',timodel,dq[:][1],'g',timodel,dq[:][2],'b',
        timodel,dq[:][3],'y',timodel,dq[:][4],'c')

# plt.subplot(323)
# inidq = n1.initial[1]
# plt.plot(timodel,np.append(inidq[:][0],inidq[:][0]),'r',timodel,np.append(inidq[:][1],inidq[:][1]),'g',timodel,np.append(inidq[:][2],inidq[:][2]),'b',
#         timodel,np.append(inidq[:][3],inidq[:][3]),'y',timodel,np.append(inidq[:][4],inidq[:][4]),'c')
plt.ylabel(namodel[1])

plt.subplot(313)
plt.plot(timodel,u[:][0],'g',timodel,u[:][1],'b',timodel,u[:][2],'y',
        timodel,u[:][3],'c')

# plt.subplot(325)
# iniu = n1.initial[2]
# plt.plot(timodel,np.append(iniu[:][0],iniu[:][0]),'r',timodel,np.append(iniu[:][1],iniu[:][1]),'g',timodel,np.append(iniu[:][2],iniu[:][2]),'b',
#         timodel,np.append(iniu[:][3],iniu[:][3]),'y',timodel,np.append(iniu[:][4],iniu[:][4]),'c')
plt.ylabel(namodel[2])

plt.suptitle('Five-Link')
plt.show()