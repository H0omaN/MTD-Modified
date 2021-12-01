"""

Created on Wed Dec  1 12:55:31 2021
@author: Hooman Ayat

This function create object masks from a 2D precipitation map
Note: The term "object" is a group of connected pixels higher than a specified threshold in a convolved precipitation map
file_address: is the address of NetCDF input file 
variable: is the varaible of the input file based on which the objects are defined. For example fo the samplefile the variable is "rainrate"
label0: object numbers start from this number when are being tracked
th: is the selected threshold to indetify the objects
r: is width/height of the smoothing window
for more information please refer to https://doi.org/10.21203/rs.3.rs-783979/v1

"""

import os
from Libraries.Track import Track

file_names=os.listdir(os.getcwd()+"/Data/")
for f in range(len(file_names)):
    print('Tracking ... '+str(file_names[f]))
    file_address,variable,label_0,th,r=os.getcwd()+"/Data/"+file_names[f],'rainrate',0,3,3    
    xr_tracked,connected_objs=Track(f,file_address,variable,label_0,th,r)
    



