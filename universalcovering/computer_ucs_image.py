'''
comptuer the feature image based on the universal covering space
'''
import scipy.misc
import numpy as np
def computer_ucs_image(ucs, feature, shape = (640,640)):

    image_array = np.zeros(shape)


    scipy.misc.toimage(image_array, cmin=0.0, cmax=...)
    return image