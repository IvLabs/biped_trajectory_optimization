import numpy as np
from matplotlib import pyplot as plt
# start = [0,-0,0,-0,-0]
start = [-0.6,0.2,0.0,-0.7,0.9]
start = np.array([0.2422276,  0.26576328, 0.04551796, 0.21421479, 0.28891085])
start2 = -start[::-1]
# start = [0.3,-0.3,0.1,-0.3,-0.5]
# start = start[::-1]
l = np.array([0.5,0.5,0.5,0.5,0.5])

origin = np.array([0,0])
origin2 = np.array([0.5000001147989219, 6.969980148596733e-13])
p1 = l[0]*np.array([np.sin(start[0]), l[0]*np.cos(start[0])]) + origin
p2 = l[0]*np.array([np.sin(start[1]), l[0]*np.cos(start[1])]) + p1
p3 = l[0]*np.array([np.sin(start[2]), l[0]*np.cos(start[2])]) + p2
p4 = l[0]*np.array([np.sin(start[3]), l[0]*-np.cos(start[3])]) + p2
p5 = l[0]*np.array([np.sin(start[4]), l[0]*-np.cos(start[4])]) + p4

p12 = l[0]*np.array([np.sin(start2[0]), l[0]*np.cos(start2[0])]) + origin2
p22 = l[0]*np.array([np.sin(start2[1]), l[0]*np.cos(start2[1])]) + p12
p32 = l[0]*np.array([np.sin(start2[2]), l[0]*np.cos(start2[2])]) + p22
p42 = l[0]*np.array([np.sin(start2[3]), l[0]*-np.cos(start2[3])]) + p22
p52 = l[0]*np.array([np.sin(start2[4]), l[0]*-np.cos(start2[4])]) + p42

plt.subplot(211)
plt.grid()
plt.title('0th step')
plt.plot([origin[0],p1[0]], [origin[1],p1[1]],'r',[p1[0],p2[0]], [p1[1],p2[1]],'g',
        [p2[0],p3[0]], [p2[1],p3[1]],'b', [p2[0],p4[0]], [p2[1],p4[1]],'y',
        [p4[0],p5[0]], [p4[1],p5[1]],'c')
plt.plot([-4, 4],[0, 0], 'black')        
# plt.axes(xlim=(-4, 4), ylim=(-4, 4))

plt.subplot(212)
plt.grid()
plt.title('1st step')
plt.plot([origin2[0],p12[0]], [origin2[1],p12[1]],'r',[p12[0],p22[0]], [p12[1],p22[1]],'g',
        [p22[0],p32[0]], [p22[1],p32[1]],'b', [p22[0],p42[0]], [p22[1],p42[1]],'y',
        [p42[0],p52[0]], [p42[1],p52[1]],'c')
plt.plot([-4, 4],[0, 0], 'black')        
# plt.axes(xlim=(-4, 4), ylim=(-4, 4))

plt.show()

