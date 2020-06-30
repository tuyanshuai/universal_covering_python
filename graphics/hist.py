import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def hist(x):
    n, bins, patches = plt.hist(x, 50, facecolor='green', alpha=0.75)
    plt.xlabel('x')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.show()