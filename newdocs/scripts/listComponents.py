# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:37:28 2024

@author: ame
"""

import oimodeler as oim
import csv

list_components_fourier=[]
list_components_image=[]
list_components_radial=[]
list_interpolators=[]
list_fitters=[]
list_filters=[]
list_class=[]

for obj in oim.__dict__:
    try:
        if issubclass(oim.__dict__[obj], oim.oimComponentFourier):
            list_components_fourier.append(obj)  
        if issubclass(oim.__dict__[obj], oim.oimComponentImage):
            list_components_image.append(obj)      
        if issubclass(oim.__dict__[obj], oim.oimComponentRadialProfile):
            list_components_radial.append(obj)  
        if issubclass(oim.__dict__[obj], oim.oimParamInterpolator):
            list_interpolators.append(obj) 
        if issubclass(oim.__dict__[obj], oim.oimFitter):
            list_fitters.append(obj)    
        if issubclass(oim.__dict__[obj], oim.oimDataFilterComponent):
                list_filters.append(obj) 
        if issubclass(oim.__dict__[obj], object):
                list_class.append(obj)                 
    except:
        
        pass
#%%  
print("*"*10,"components Fourier","*"*10)
table_components_fourier =[]
for cname in list_components_fourier:
    try:
        print(cname)
        ti=[cname]
        c = oim.__dict__[cname]()
        p = c.params
        name = c.name
        
        ti.append(name)
        txt=""
        for pname in p:
            txt+=":abbr:`"
            txt+=pname
            txt+="("
            txt+=p[pname].description
            
            txt+=")`, "
        txt=txt[:-2]
        ti.append(txt)
        table_components_fourier.append(ti)
    except:
        pass
f=open("table_components_fourier.csv","w")
w=csv.writer(f,delimiter="|")
w.writerows(table_components_fourier)
f.close()
#%%
print("*"*10,"interpolators","*"*10)
param0=oim.oimParam()
table_param_interpolators =[]
for cname in list_interpolators:
    print(cname)
    try:
        ti=[cname]
        #c = oim.__dict__[cname]()
        
        
        table_param_interpolators.append(ti)
    except:
        pass
f=open("table_param_interpolators.csv","w")
w=csv.writer(f,delimiter="|")
w.writerows(table_param_interpolators)
f.close()
#%%
print("*"*10,"filters","*"*10)

table_filters =[["Filter Name","Short description","Class Keywords"]]
for cname in list_filters:
    print(cname)
    try:
        
        filt = oim.__dict__[cname]()
        boo=[cname]
        boo.append(filt.description)

        txt=""
        for pname in filt.params:
            txt+=pname+", "
        boo.append(txt[:-2])
        table_filters.append(boo)
    except:
        pass
f=open("table_dataFilter.csv","w",newline='')
w=csv.writer(f,delimiter="|")
w.writerows(table_filters)
f.close()