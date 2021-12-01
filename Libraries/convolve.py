# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:55:31 2021
@author: Hooman Ayat

This function create object masks from a 2D precipitation map
for more information please refer to https://doi.org/10.21203/rs.3.rs-783979/v1
"""
import numpy as np
from scipy import ndimage
def convolve(two_d_arrr,R,Th):
    two_d_arr=np.copy(two_d_arrr)
    two_d_arr[two_d_arr<0]=0
    two_d_arr[np.isnan(two_d_arr)]=0
    convolved=ndimage.uniform_filter(two_d_arr, size=2*R+1, mode='constant')
    convolved[convolved<Th]=0
    convolved[convolved>0]=1
    return convolved