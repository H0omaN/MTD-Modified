# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:59:48 2021
@author: Hooman Ayat: Hooman.ayat@gmail.com

This function detect and track the objects in the input data
for more information please refer to https://doi.org/10.21203/rs.3.rs-783979/v1
"""
from scipy import ndimage
import numpy as np
from skimage.measure import label
from skimage.measure import regionprops
from collections import OrderedDict
import os
import xarray as xr
from Libraries.convolve import convolve
def Track(fileno,file_address,variable,label_0,th,r):
    
    
    ### opening the input netcdf file
    xr_data=xr.open_dataset(file_address)
    ### creating object masks 
    MTD_mask_data=list()
    for i in range(len(xr_data[variable])):
        MTD_mask_data.append(convolve(xr_data[variable][i],3,3)[:,:])   
    ### adding masks to the original data    
    xr_data["MTD_masks"]=([ 'time','y', 'x'],  MTD_mask_data)    
    MTD_Cube=np.copy(xr_data["MTD_masks"])

    ### Tracked_maps: time series of maps of the tracked objects 
    Tracked_maps=list()
    
    ### IDNo: Current id for the detected object              
    IDNo=label_0
    
    ### appendedpoints: this is the lists of merged or sperated objects
    appendedpoints=list()
    
    for step in range(len(MTD_Cube)): 
        if step==0:
            ### for this step we just name the objects randomly starting from IDNo (label_image is the map of detected objects at the current time step)
            label_image = np.asarray(label(MTD_Cube[step]))+IDNo
            
            ### take the list of randomly named detected objects in the current precipitation map
            obj_list_label_image=np.unique(label_image)
            obj_list_label_image=obj_list_label_image[obj_list_label_image>0]
            Tracked_maps.append(label_image)
            IDNo=IDNo+int(len(obj_list_label_image))     
        if step>0:
            
            ### In this section we find the objects (at the present time step) that are connected to the objects at previous time step in order to name after them
            ### If the merge or split is occuring we name them with a different id number (IDno)            
            label_image = label(MTD_Cube[step])  
            
            ### object properties at current time steps
            label_image_prop=regionprops(label_image)
            
            ### object properties at previous time step
            Mod_Objects_prop=regionprops(Tracked_maps[-1])
    
            ### we use temp3dobj to connect object maps at previous and current time step
            Temp3Dobj=list()
            Temp3Dobj.append(np.copy(Tracked_maps[-1]))
            Temp3Dobj.append(np.copy(label_image))
            Temp3Dobj[0][Temp3Dobj[0]>0]=1
            Temp3Dobj[1][Temp3Dobj[1]>0]=1
     
            ### identifying the connected objects in time
            label_image3D=label(np.asarray(Temp3Dobj)) 
            label_image3D_prop=regionprops(label_image3D)
            for label_image3D_member in label_image3D_prop:
                
                ### label_image3D_member.bbox shows the 3d box fits the object element 3 and zero is for time dimension for examle :(0, 166, 79, 1, 173, 82)
                if label_image3D_member.bbox[3]==2 and label_image3D_member.bbox[0]==0:   ### this condition is for the objects connected in time
                    
                    ### omitting all other objects except the selected one ### slice0 and slice1 are for selected connected object at two different time steps
                    imagevaset0=np.copy(label_image3D_member._label_image[0])
                    imagevaset0[imagevaset0 != label_image3D_member.label]=0                  
                    label_image3D_slice0=regionprops(np.asarray(imagevaset0))
                    imagevaset1=np.copy(label_image3D_member._label_image[1])
                    imagevaset1[imagevaset1 != label_image3D_member.label]=0
                    label_image3D_slice1=regionprops(np.asarray(imagevaset1))
      
                    ### Handling split/merge events: in this section we are using the centroid of the connected objects in time at each time step to detect split/merge events
                    MergFinder=1
                    for Mod_Objects_prop_member in Mod_Objects_prop:
                        ### ckecking if the centorid of slice 0 from a label 3d output is equivalent with the object's cetroid from previuos step or not
                        if label_image3D_slice0[0].centroid == Mod_Objects_prop_member.centroid:
                             newlabel=Mod_Objects_prop_member.label
                             MergFinder=0 
                             
                    SeperationFinder=1          
                    for label_image_prop_member in label_image_prop:
                         
                         if label_image3D_slice1[0].centroid == label_image_prop_member.centroid:
                             oldlabel=label_image_prop_member.label                           
                             SeperationFinder=0 
                           
                    ### Handling merging events                
                    if MergFinder==1 and SeperationFinder!=1:
                        IDNo=IDNo+1
                        newlabel=IDNo
                        
                        ### find the connected objects (with different ids) via merge events ... labels that are connected after each merging event    
                        Devided_Merged_Last_Step_prop=regionprops(label(np.asarray(imagevaset0)))
                        for obj0 in Devided_Merged_Last_Step_prop:
                            for obj0ref in Mod_Objects_prop:         
                                if obj0.centroid == obj0ref.centroid:
                                    if (obj0ref.label,newlabel) not in appendedpoints and (newlabel,obj0ref.label) not in appendedpoints:
                                        appendedpoints.append((obj0ref.label,newlabel))
    
                    if SeperationFinder==0:         
                        label_image[label_image==oldlabel]=-1*newlabel  
                                        
                    ### handling the spliting events  : Note that if split and merge occurring at the same time, it is considered as a simple split.                  
                    if SeperationFinder==1:
                        oldlabellist=list()
                        newlabellist=list()
                        centroids1=list()
                        Area1=list()
                        plist=list()
                        ID1=list()    
                        Devided_Merged_Current_Step_prop=regionprops(label(np.asarray(imagevaset1)))
                        for objj in Devided_Merged_Current_Step_prop:
                            for label_image_prop_member in label_image_prop:
                                if objj.centroid == label_image_prop_member.centroid:
                                    plist.append(label_image_prop_member.centroid)
                                    oldlabellist.append(label_image_prop_member.label )
                                    IDNo=IDNo+1
                                    newlabellist.append(IDNo)  
                                    centroids1.append(label_image_prop_member.centroid)
                                    Area1.append(label_image_prop_member.area)
                                    ID1.append(IDNo)                             
                        iiii=0
                        for labels in oldlabellist: 
                            label_image[label_image==labels]=-1*newlabellist[iiii]                       
                            iiii=iiii+1       
                        Area0=list()
                        centroids0=list()
                        ID0=list()
                        Devided_Merged_Last_Step_prop=regionprops(label(np.asarray(imagevaset0)))
                        for obj0 in Devided_Merged_Last_Step_prop:
                            for obj0ref in Mod_Objects_prop:         
                                if obj0.centroid == obj0ref.centroid:
                                    Area0.append(obj0ref.area)
                                    centroids0.append(obj0ref.centroid)
                                    ID0.append(obj0ref.label)  
                        for r0 in range(len(Area0)):
                            for r1 in range(len(Area1)):
                                p0=ID0[r0]
                                p1=ID1[r1]
                                if (p0,p1) not in appendedpoints and (p1,p0) not in appendedpoints:
                                    appendedpoints.append((p0,p1))
    
                ### Handling newly appeared objects:
                if label_image3D_member.bbox[3]==2 and label_image3D_member.bbox[0]==1: ### this is for new objects appear in the current time step
                    IDNo=IDNo+1
                    for label_image_prop_member in label_image_prop:
                        x=label_image3D_member.centroid[1]
                        y=label_image3D_member.centroid[2]
                        if (x,y) == label_image_prop_member.centroid:
                            oldlabel=label_image_prop_member.label 
                            newlabel=IDNo
                    label_image[label_image==oldlabel]=-1*newlabel
            label_image=np.absolute(label_image)  
            Tracked_maps.append(label_image)
            
    ### list of connected objects in time        
    Connections= list(OrderedDict.fromkeys(appendedpoints))    
    
    ### Adding the tracked objects into the original input file
    xr_data["MTD_tracked"]=([ 'time','y', 'x'],  Tracked_maps)

    Connected_objects = xr.DataArray(data=Connections,dims=["x", "y"],coords=dict(x=(["x"], np.arange(len(Connections))),y=(["y"],np.arange(2) )),attrs=dict(description="Connected objects via split/merge events", units="object no.", ),)
    
    ### saving the outputs
    if not os.path.exists(os.getcwd()+"/Output/"):
        os.makedirs(os.getcwd()+"/Output/")
    xr_data.to_netcdf(path=os.getcwd()+"/Output/"+str(fileno)+"_Tracked.nc")
    Connected_objects.to_netcdf(path=os.getcwd()+"/Output/"+str(fileno)+"_Connected_objects.nc")
    
    ### returning the outputs
    return xr_data,Connected_objects