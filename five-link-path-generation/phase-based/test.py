from new_gait import *
from matplotlib import pyplot as plt
from time import sleep
# start_angles = ca.MX([-0.3,0.3,0.2,-0.3,0.3])
# start_angles = ca.MX([-0.2,0.1,0.2,-0.5,1.5]) # ostrich 2
# start_angles = ca.MX([-0.6,0.2,0.0,-0.7,0.9]) # ostrich 3
start_angles = ca.MX.zeros(5)
# start_angles = ca.MX([-0.6,0.7,0.0,-0.5,0.9]) # ostrich
# start_angles = ca.MX([0.3,-0.3,0.1,-0.3,-0.45]) # human

start_pos = [[0,0]]
start_angular_vel = ca.MX.zeros(5)
q = []; dq = []; u = []; pos = []; time = []
f = 1

# def callplot(x, y):
#     plt.plot(x[0], y[0], 'o')


for k in range(f):
    # try:
    model = walker(start_angles, start_angular_vel, start_pos[-1])        
    problem = nlp(model)
    # model.opti.callback(lambda x : callplot(opti.debug.value(model.x[-1]), opti.debug.value(model.xdot[-1])))

    sol = model.opti.solve_limited()
    
# import casadi as ca
# from matplotlib import pyplot as plt

# opti = ca.Opti()

# x = opti.variable()
# y = opti.variable()
# opti.minimize(  (y-x**2)**2   )
# opti.subject_to( x**2+y**2==1 )
# opti.subject_to(       x+y>=1 )

# p_opts = {"expand":True}
# s_opts = {"max_iter": 100, "print_level":0}
# opti.solver("ipopt",p_opts,
#                     s_opts)
                    
# def callplot(x, y):
#     plt.plot(x, y, 'ro')


# opti.callback(lambda x : callplot(opti.debug.value(x),opti.debug.value(y)))

# sol = opti.solve()

# print(sol.value(x))
# print(sol.value(y))

# plt.show()