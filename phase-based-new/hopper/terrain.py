import numpy as np
import casadi as ca

class Terrain():
    def __init__(self, name='flat'):
        super().__init__()

        x_pos = ca.MX.sym('x_pos', 1)
        self.mu = 1.
        if name == 'flat':
            self.terrain_factor = 0.
            y_pos = self.terrain_factor*x_pos
            self.f = ca.Function('flat',[x_pos],[y_pos],['x'],['y'])
            self.df = self.f.jacobian()
        elif name == 'sin':
            self.terrain_factor = 1.
            y_pos = self.terrain_factor*ca.sin(x_pos)
            self.f = ca.Function('sin',[x_pos],[y_pos],['x'],['y'])
            self.df = self.f.jacobian()
        elif name == 'wedge':
            self.terrain_factor = 1.
            self.f = ca.Function('wedge',[x_pos],[y_pos],['x'],['y'])
            self.df = self.f.jacobian()
        elif name == 'stairs': # smooth
            self.terrain_factor = 50
            self.order = 4
            f = ca.Function('smooth_stair', [x_pos],[(self.terrain_factor*x_pos - ca.sin(self.terrain_factor*x_pos))/self.terrain_factor], ['x'],['y'])
            # y_pos = x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor)) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor))) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor)) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor) - ca.sin(x_pos*self.terrain_factor - ca.sin(x_pos*self.terrain_factor))))
            # y_pos /= abs(self.terrain_factor)
            # self.f = ca.Function('smooth_stair',[x_pos],[y_pos],['x'],['y'])
            self.f = f.fold(self.order)
            self.df = self.f.jacobian()

        self.name = name

    def heightMap(self, x):
        return self.f(x=x)['y']    

    def heightMapNumericalSlope(self, x):
        return self.df(x=x)['jac']

    def heightMapNormalVector(self, x):
        tangent_vector = ca.MX.ones(2, 1)
        tangent_vector[1, 0] = self.df(x=x)['jac']
        normal_vector = ca.MX.ones(2, 1)
        normal_vector[0, 0] = -tangent_vector[1, 0]
        normal_vector /= ca.norm_2(normal_vector)
        return normal_vector

    def heightMapTangentVector(self, x):
        tangent_vector = ca.MX.ones(2, 1)
        tangent_vector[1, 0] = self.df(x=x)['jac']
        tangent_vector /= ca.norm_2(tangent_vector)
        return tangent_vector

# test check for sanity

# test_terrain = Terrain(name='stairs')

# print(test_terrain.heightMap(2), test_terrain.heightMapNormalVector(2))
