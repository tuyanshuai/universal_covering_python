import numpy as np
def hyperbolic_deck_transorm(s1,t1,s2,t2):
    theta1 = np.angle((t1 - s1) / (1 - np.conj(s1) * t1))
    theta2 = np.angle((t2 - s2) / (1 - np.conj(s2) * t2))

    def f1(z):
        return (z-s1) / (1-np.conj(s1)*z)

    def f2(z):
        return np.exp(-theta1*1j)*f1(z)

    def f3(z):
        return np.exp(theta2*1j)*f2(z)

    def f4(z):

        return (f3(z)+s2) / (1+np.conj(s2)*f3(z))

    return f4