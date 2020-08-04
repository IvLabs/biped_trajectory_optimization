from gait import * 

q = []; dq = []; u = []; pos = []; time = []
f = 1
test = {'q1': [], 'q2': [], 'q3': [], 'q4': []}
from matplotlib import pyplot as plt
start = [-0.6,0.7,0.0,-0.5,-0.3]
# start = [0.3,-0.3,0.1,-0.3,-0.3]
# start = [0,-0,0,-0,-0]

start_vel = [0., 0., 0., 0., 0.]
start_pos = [[0,0]]
for k in range(f):
    if k == 0:
        me = walker(start,start_pos[k])
        n1 = nlp(me)
        debug = me.opti.debug.value
        print(debug)
        sol1 = me.opti.solve()
        
    else :
        me = walker(start,start_pos[k], start_vel)
        n1 = nlp(me)
        debug = me.opti.debug.value
        print(debug)
        sol1 = me.opti.solve()

    for j in range(5):
        # tq = []; tdq = []; tu = []; tpos = []
        tempq = [];tempdq = [];tempu = [];temp = []
        for i in range(me.N):
            tempq.append(sol1.value(me.state[i][0][j]))    
            tempdq.append(sol1.value(me.state[i][1][j]))
            temp.append([sol1.value(me.pos[i][j][0]),sol1.value(me.pos[i][j][1])]) 
            if j < 4:
                tempu.append(sol1.value(me.u[i][0][j]))

        # if k%2 != 0:
        #     tempq = tempq[::-1] 
        #     tempdq = tempdq[::-1]
        #     tempu = tempu[::-1]
        #     # temp = temp[::-1]
        test['q' + str(k)] = tempq
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
    start = [q[0][-1],q[1][-1],q[2][-1],q[3][-1],q[4][-1]]
    start_vel = [dq[0][-1],dq[1][-1],dq[2][-1],dq[3][-1],dq[4][-1]]
    start_pos.append([pos[4][-1][0],pos[4][-1][1]])

#     if k > 1:
        # print(start_pos[k][1])

time = np.arange(0.0, f*me.T, me.h)
terrain = np.linspace(-2,6,100)
from matplotlib import animation
from celluloid import Camera
# print(len(pos[0]))
fig = plt.figure()
camera = Camera(fig)
k = 0
for i in range(f*me.N):
    # print(i)    
    p1 = [pos[0][i][1],pos[0][i][0]]
    p2 = [pos[1][i][1],pos[1][i][0]]
    p3 = [pos[2][i][1],pos[2][i][0]]
    p4 = [pos[3][i][1],pos[3][i][0]]
    p5 = [pos[4][i][1],pos[4][i][0]]
    # plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    plt.axes(xlim=(-1, 6), ylim=(-2, 2))
    # plt.plot([0,-p1[1]], [0,p1[0]],'r',[-p1[1],-p2[1]], [p1[0],p2[0]],'b',
    #     [-p2[1],-p3[1]], [p2[0],p3[0]],'c',
    #     [-p2[1],p4[1] - 2*p2[1]], [p2[0],2*p2[0]-p4[0]],'b',
    #     [p4[1] - 2*p2[1],p5[1]], [2*p2[0]-p4[0],(p5[0] - 2*p2[0])],'r')
    if i%me.N == 0:
            p0 = start_pos[k]
            print(p0)
            k += 1 

    plt.plot([p0[0],p1[1]], [p0[1],p1[0]],'r',[p1[1],p2[1]], [p1[0],p2[0]],'g',
            [p2[1],p3[1]], [p2[0],p3[0]],'b', [p2[1],p4[1]], [p2[0],p4[0]],'y',
            [p4[1],p5[1]], [p4[0],p5[0]],'c')
    
    plt.plot(terrain,me.terrain_factor*np.sin(terrain),'g')    
    # if cv2.waitKey(0) & 0xFF == ord("q"):
    # #     break
    camera.snap()
animation = camera.animate(interval=60)
# animation.save('path_.mp4')
plt.show()
plt.close()

name = ['q','dq','u']

# plt.plot(test['q1'])
# plt.plot(test['q2'])
# plt.plot(test['q3'])
# plt.plot(test['q4'])

# plt.show()
plt.subplot(311)
plt.title('Optimised Solution')
plt.plot(time, q[:][0],'r', time, q[:][1],'g', time,q[:][2],'b',
        time,q[:][3],'y', time, q[:][4],'c')

# plt.subplot(321)
# plt.title('Initial Guess')
# iniq = n1.initial[0]
# plt.plot(time,np.append(iniq[:][0],iniq[:][0]),'r',time,np.append(iniq[:][1],iniq[:][1]),'g',time,np.append(iniq[:][2],iniq[:][2]),'b',
#         time,np.append(iniq[:][3],iniq[:][3]),'y',time,np.append(iniq[:][4],iniq[:][4]),'c')
plt.ylabel(name[0])

plt.subplot(312)
plt.plot(time,dq[:][0],'r',time,dq[:][1],'g',time,dq[:][2],'b',
        time,dq[:][3],'y',time,dq[:][4],'c')

# plt.subplot(323)
# inidq = n1.initial[1]
# plt.plot(time,np.append(inidq[:][0],inidq[:][0]),'r',time,np.append(inidq[:][1],inidq[:][1]),'g',time,np.append(inidq[:][2],inidq[:][2]),'b',
#         time,np.append(inidq[:][3],inidq[:][3]),'y',time,np.append(inidq[:][4],inidq[:][4]),'c')
plt.ylabel(name[1])

plt.subplot(313)
plt.plot(time,u[:][0],'g',time,u[:][1],'b',time,u[:][2],'y',
        time,u[:][3],'c')

# plt.subplot(325)
# iniu = n1.initial[2]
# plt.plot(time,np.append(iniu[:][0],iniu[:][0]),'r',time,np.append(iniu[:][1],iniu[:][1]),'g',time,np.append(iniu[:][2],iniu[:][2]),'b',
#         time,np.append(iniu[:][3],iniu[:][3]),'y',time,np.append(iniu[:][4],iniu[:][4]),'c')
plt.ylabel(name[2])

plt.suptitle('Five-Link')
plt.show()