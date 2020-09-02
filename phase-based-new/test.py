import casadi as ca
import numpy as np
import matplotlib.pyplot as plt

# p = ca.MX.sym('t', 2)
opti = ca.Opti()

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# xgrid = np.linspace(0,6,4)
# V = [-1,-1,-2,-3, 100]
# lut = ca.interpolant('LUT','bspline',[xgrid],V)
# print(lut(2.5))

# x = np.linspace(0,6,50)

# plt.plot(x, lut(x))
# plt.plot(xgrid, V, 'ro')
# plt.show()

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# Splines


delta_T = ca.MX.sym('delta_T',1)
t = ca.MX.sym('t', 1)

x0  = ca.MX.sym( 'x0_1', 1)
x1  = ca.MX.sym( 'x1_1', 1)
dx0 = ca.MX.sym('dx0_1', 1)
dx1 = ca.MX.sym('dx1_1', 1)

a0 = x0
a1 = dx0
a2 = -(delta_T**(-2))*(3*(x0 - x1) + delta_T*(2*dx0 + dx1))
a3 = (delta_T**(-3))*(2*(x0 - x1) + delta_T*(dx0 + dx1))
a4 = -(delta_T**(-4))*((x0 - x1) + delta_T*(dx1))

x = a0 + a1*t + a2*(t**2) + a3*(t**3) 
# x = a0 + a1*t + a2*(t**2) + a3*(t**3) + a4*(t**4)

f = ca.Function('f', [delta_T, t, x0, x1, dx0, dx1], [x], ['delta_T', 't', 'x0', 'x1', 'dx0', 'dx1'], ['x'])

# y = []
# y_super = []
# T = ca.MX.sym('T',1)
# delta_T = T/3
# x0_1  = x0  
# x1_1  = x1  
# dx0_1 = dx0 
# dx1_1 = dx1 

# x0_2  = x1_1
# x1_2  = ca.MX.sym( 'x1_2', 1)
# dx0_2 = dx1_1
# dx1_2 = ca.MX.sym('dx1_2', 1)

# x0_3  = x1_2
# x1_3  = ca.MX.sym( 'x1_3', 1)
# dx0_3 = dx1_2
# dx1_3 = ca.MX.sym('dx1_3', 1)

# # for j in range(2):
# for i in range(3):
#     if i == 0:
#         x0  =  x0_1
#         x1  =  x1_1
#         dx0 = dx0_1
#         dx1 = dx1_1
#         y.append(f(delta_T, t, x0, x1, dx0, dx1))
#     elif i == 1:
#         x0  =  x0_2 
#         x1  =  x1_2 
#         dx0 = dx0_2
#         dx1 = dx1_2
#         y.append(f(delta_T, t, x0, x1, dx0, dx1))
#     else:
#         x0  =  x0_3 
#         x1  =  x1_3 
#         dx0 = dx0_3
#         dx1 = dx1_3
#         y.append(f(delta_T, t, x0, x1, dx0, dx1))

#     # T = ca.MX.sym('T' + str(j+1), 1)
#     # delta_T = T/3
#     # x0_1  =  x0_3 
#     # x1_1  = dx1_3 
#     # dx0_1 =  x0_3 
#     # dx1_1 = dx1_3 

#     # x0_2  = x1_1
#     # x1_2  = ca.MX.sym( 'x1_2' + str(j+1), 1)
#     # dx0_2 = dx1_1
#     # dx1_2 = ca.MX.sym('dx1_2' + str(j+1), 1)

#     # x0_3  = x1_2
#     # x1_3  = ca.MX.sym( 'x1_30' + str(j+1), 1)
#     # dx0_3 = dx1_2
#     # dx1_3 = ca.MX.sym('dx1_30' + str(j+1), 1)

# Y = ca.vcat(y)
# F = ca.Function("F", [x0_1, x1_1, dx0_1, dx1_1,
#                     x1_2, dx1_2, x1_3, dx1_3, T, t], [Y])
# dphase_spline = F.jacobian()
# print(dphase_spline()['jac'][:,-1])

# T = 0.5
# del_t = T/3
# x0_1, x1_1, dx0_1, dx1_1, x1_2, dx1_2, x1_3, dx1_3 = 0, 3, 0, 0, 2, 5, 0, 0  
# N = 10
# x = ca.linspace(0, T, N)
# # x = x.reshape(N, 1) 
# y_plot = []

# print(dphase_spline(x0_1, x1_1, dx0_1, dx1_1,
#                     x1_2, dx1_2, x1_3, dx1_3, T, x[2], True))
# a0 = x0_1
# a1 = dx0_1
# a2 = -(del_t**(-2))*(3*(x0_1 - x1_1) + del_t*(2*dx0_1 + dx1_1))
# a3 = (del_t**(-3))*(2*(x0_1 - x1_1) + del_t*(dx0_1 + dx1_1))

# print(a1 + 2*a2*x[2] + 3*a3*(x[2]**2))
# # for j in range(2):
# for i in range(x.numel()):
#     if i <= (N-1)/3:
#         y_plot.append(F(x0_1, x1_1, dx0_1, dx1_1,
#                     x1_2, dx1_2, x1_3, dx1_3, T, x[i])[0])
#     elif (N-1)/3 < i <= 2*(N-1)/3:
#         y_plot.append(F(x0_1, x1_1, dx0_1, dx1_1,
#                     x1_2, dx1_2, x1_3, dx1_3, T, x[i-int((N-1)/3)])[1])
#     else:
#         y_plot.append(F(x0_1, x1_1, dx0_1, dx1_1,
#                     x1_2, dx1_2, x1_3, dx1_3, T, x[i-int(2*(N-1)/3)])[2])
    

#     # x0_1  =  x0_3 
#     # x1_1  = dx1_3 
#     # dx0_1 =  x0_3 
#     # dx1_1 = dx1_3 

#     # x0_2  = x1_1
#     # x1_2  = 3
#     # dx0_2 = dx1_1
#     # dx1_2 = 0

#     # x0_3  = x1_2
#     # x1_3  = 0
#     # dx0_3 = dx1_2
#     # dx1_3 = 0

# plt.plot(x, ca.vcat(y_plot), 'ro')
# plt.show()


######################
# delta_T = opti.variable(1)
# t = 0.2
# x0 = opti.variable(1)
# x1 = opti.variable(1)
# dx0 = opti.variable(1)
# dx1 = opti.variable(1)  

# print(f(delta_T=delta_T, t=t, x0=x0, x1=x1, dx0=dx0, dx1=dx1)['x'])

######################


t = 0.0
x0 = 0
x1 = 3
dx0 = 0
dx1 = 0  
# x0 = 2
# x1 = 0
# dx0 = 0
# dx1 = 0  
y = []
N = 100
T = 0.5
delta_T = T/3
count = 0
while len(y) < N:
    y.append(f(delta_T=delta_T, t=t, x0=x0, x1=x1, dx0=dx0, dx1=dx1)['x'])
    t += T/N
    if t >= T/3 and count == 0:
        # delta_T = T - delta_T
        x0 = x1
        x1 = 2
        dx0 = dx0
        dx1 = 5
        t = 0
        count += 1
    if t >= T/3 and count == 1:
        # delta_T = T - delta_T
        x0 = x1
        x1 = 0
        dx0 = dx0
        dx1 = 0
        t = 0
        count += 1    
        # print(f(delta_T=delta_T, t=t, x0=x0, x1=x1, dx0=dx0, dx1=dx1)['x'])
    # print("x0 = " + str(x0) + ", x1 = " + str(x1) + ", dx0 = " + str(dx0) + ", dx1 = " + str(dx1) + ", y = "+ str(y[-1]))
time = ca.linspace(0, T, N)
plt.plot(time, y, label='cubic')
plt.legend()


delta_T = ca.MX.sym('delta_T',1)
t = ca.MX.sym('t', 1)

x0  = ca.MX.sym( 'x0_1', 1)
x1  = ca.MX.sym( 'x1_1', 1)
dx0 = ca.MX.sym('dx0_1', 1)
dx1 = ca.MX.sym('dx1_1', 1)

a0 = x0
a1 = dx0
a2 = -(delta_T**(-2))*(3*(x0 - x1) + delta_T*(2*dx0 + dx1))
a3 = (delta_T**(-3))*(2*(x0 - x1) + delta_T*(dx0 + dx1))
a4 = -(delta_T**(-2))*((3*x0 - x1) + delta_T*(2*dx0 + dx1))

x = a0 + a1*t + a2*(t**2) + a3*(t**3) + a4*(t**4)

f = ca.Function('f', [delta_T, t, x0, x1, dx0, dx1], [x], ['delta_T', 't', 'x0', 'x1', 'dx0', 'dx1'], ['x'])

t = 0.0
x0 = 0
x1 = 3
dx0 = 0
dx1 = 0  
# x0 = 2
# x1 = 0
# dx0 = 0
# dx1 = 0  
y = []
N = 100
T = 0.5
delta_T = T/3
count = 0
while len(y) < N:
    y.append(f(delta_T=delta_T, t=t, x0=x0, x1=x1, dx0=dx0, dx1=dx1)['x'])
    t += T/N
    if t >= T/3 and count == 0:
        # delta_T = T - delta_T
        x0 = x1
        x1 = 2
        dx0 = dx0
        dx1 = 5
        t = 0
        count += 1
    if t >= T/3 and count == 1:
        # delta_T = T - delta_T
        x0 = x1
        x1 = 0
        dx0 = dx0
        dx1 = 0
        t = 0
        count += 1    
        # print(f(delta_T=delta_T, t=t, x0=x0, x1=x1, dx0=dx0, dx1=dx1)['x'])
    # print("x0 = " + str(x0) + ", x1 = " + str(x1) + ", dx0 = " + str(dx0) + ", dx1 = " + str(dx1) + ", y = "+ str(y[-1]))
time = ca.linspace(0, T, N)
plt.plot(time, y, label='bi_quad')
plt.legend()


plt.show()

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# # Inverse Kinnematics

# q = opti.variable(3)
# c_real = [0.5, 0.8]
# length = ca.DM([0.5,0.5,0.5])

# p0 = ca.DM([0, 0])
# p = ca.MX.zeros(2,3)
# p[0,0],p[1,0] = length[0]*ca.sin(q[0]) + p0[0,0]  , length[0]*ca.cos(q[0]) + p0[1,0]
# p[0,1],p[1,1] = length[1]*ca.sin(q[1]) +  p[0,0]  , length[1]*ca.cos(q[1]) +   p[1,0]
# p[0,2],p[1,2] = length[2]*ca.sin(q[2]) +  p[0,1]  , length[2]*ca.cos(q[2]) +   p[1,1]

# c = ca.MX.zeros(2,3)
# c[0,0],c[1,0] =  length[0]*ca.sin(q[0])/2 + p0[0,0], length[0]*ca.cos(q[0])/2 + p0[1,0]
# c[0,1],c[1,1] =  length[1]*ca.sin(q[1])/2 +  p[0,0], length[1]*ca.cos(q[1])/2 +  p[1,0]
# c[0,2],c[1,2] =  length[2]*ca.sin(q[2])/2 +  p[0,1], length[2]*ca.cos(q[2])/2 +  p[1,1]

# opti.minimize(ca.sumsqr(c_real - c[:,2]) + ca.sumsqr(c[0,2] - p0[0,0]))

# opti.subject_to(q[1]>=q[0])
# opti.subject_to(ca.vec(p[:,1])>=0)

# opti.subject_to(ca.fabs(q[2])<=np.pi/3)

# opti.bounded([-np.pi/2]*3,q,[np.pi/2]*3)

# p_opts = {"expand":True}
# s_opts = {"max_iter": 100}
# opti.solver("ipopt",p_opts,
#                     s_opts)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.grid()

# def callplot(q0):
#     plt.cla()
#     ax.grid()
#     ax.set_xlim([-2, 2])
#     ax.set_ylim([-2, 2])

#     temp = ca.DM.zeros(2,3)

#     temp[0,0],temp[1,0] = length[0]*ca.sin(q0[0]) + p0[0,0]  , length[0]*ca.cos(q0[0]) + p0[1,0]
#     temp[0,1],temp[1,1] = length[1]*ca.sin(q0[1]) + temp[0,0], length[1]*ca.cos(q0[1]) + temp[1,0]
#     temp[0,2],temp[1,2] = length[2]*ca.sin(q0[2]) + temp[0,1], length[2]*ca.cos(q0[2]) + temp[1,1]

#     ax.plot([p0[0,0], temp[0,0]], [p0[1,0], temp[1,0]], 'r', lw=4)
#     ax.plot([temp[0,0], temp[0,1]], [temp[1,0], temp[1,1]], 'g', lw=4)
#     ax.plot([temp[0,1], temp[0,2]], [temp[1,1], temp[1,2]], 'b', lw=10)
#     plt.draw()
#     plt.pause(1e-5)

# opti.callback(lambda x : callplot(opti.debug.value(q)))

# sol = opti.solve()

# print(sol.value(q))

# plt.show()

############################################################################################
############################################################################################
############################################################################################
############################################################################################


































# random things
# y = a0*t

# f = ca.Function('f',[a0,t],[y],['t','a0'],['y'])

# d = ca.MX.sym('d',5)
# atest = ca.MX.sym('atest',1)
# e = 1.

# print(f(t=e,a0=atest))
# # atest = 2.
# d = np.linspace(1,2,20)
# print(f(t=d,a0=atest)['y'].shape)
# from new_gait import *
# from matplotlib import pyplot as plt
# from time import sleep
# # start_angles = ca.MX([-0.3,0.3,0.2,-0.3,0.3])
# # start_angles = ca.MX([-0.2,0.1,0.2,-0.5,1.5]) # ostrich 2
# # start_angles = ca.MX([-0.6,0.2,0.0,-0.7,0.9]) # ostrich 3
# start_angles = ca.MX.zeros(5)
# # start_angles = ca.MX([-0.6,0.7,0.0,-0.5,0.9]) # ostrich
# # start_angles = ca.MX([0.3,-0.3,0.1,-0.3,-0.45]) # human

# start_pos = [[0,0]]
# start_angular_vel = ca.MX.zeros(5)
# q = []; dq = []; u = []; pos = []; time = []
# f = 1

# # def callplot(x, y):
# #     plt.plot(x[0], y[0], 'o')


# for k in range(f):
#     # try:
#     model = walker(start_angles, start_angular_vel, start_pos[-1])        
#     problem = nlp(model)
#     # model.opti.callback(lambda x : callplot(opti.debug.value(model.x[-1]), opti.debug.value(model.xdot[-1])))

#     sol = model.opti.solve_limited()
    
# # import casadi as ca
# # from matplotlib import pyplot as plt

# # opti = ca.Opti()

# # x = opti.variable()
# # y = opti.variable()
# # opti.minimize(  (y-x**2)**2   )
# # opti.subject_to( x**2+y**2==1 )
# # opti.subject_to(       x+y>=1 )

# # p_opts = {"expand":True}
# # s_opts = {"max_iter": 100, "print_level":0}
# # opti.solver("ipopt",p_opts,
# #                     s_opts)
                    
# # def callplot(x, y):
# #     plt.plot(x, y, 'ro')


# # opti.callback(lambda x : callplot(opti.debug.value(x),opti.debug.value(y)))

# # sol = opti.solve()

# # print(sol.value(x))
# # print(sol.value(y))

# # plt.show()