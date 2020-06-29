import numpy as np
def euclidean_deck_transform(s1,t1,s2,t2):
    def f(z):
        return z + s2 - s1

    return f
