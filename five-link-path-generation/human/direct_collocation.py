import casadi as ca
import numpy as np
import matplotlib.pyplot as plt

class walker():
    def __init__(self):
        self.makePolynomial()

    def makePolynomial(self):        
        # Degree of interpolating polynomial
        self.d = 3

        # Get collocation points
        self.tau_root = np.append(0, ca.collocation_points(self.d, 'legendre'))

        # Coefficients of the collocation equation
        self.C = np.zeros((self.d+1, self.d+1))

        # Coefficients of the continuity equation
        self.D = np.zeros(self.d+1)

        # Coefficients of the quadrature function
        self.B = np.zeros(self.d+1)
    
biped = walker()
print(biped.tau_root)    