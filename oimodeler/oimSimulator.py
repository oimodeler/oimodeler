# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:26:42 2021

@author: Ame
"""

import numpy as np
from astropy.io import fits
import os
import oimodeler as oim


def corrFlux2Vis2(vcompl):
    nB=vcompl.shape[0]
    norm=np.outer(np.ones(nB-1),vcompl[0,:])
    return np.abs(vcompl[1:,:]/norm)**2

def corrFlux2VisAmpAbs(vcompl):
    nB=vcompl.shape[0]
    norm=np.outer(np.ones(nB-1),vcompl[0,:])
    return np.abs(vcompl[1:,:]/norm)

# TODO : not real formula for diff Vis 
def corrFlux2VisAmpDif(vcompl):
    nlam=vcompl.shape[1]
    norm=np.outer(np.mean(vcompl[1:,:],axis=1),np.ones(nlam))
    return np.abs(vcompl[1:,:]/norm)
        
def corrFlux2VisAmpCor(vcompl):
    return np.abs(vcompl[1:,:])

def corrFlux2VisPhiAbs(vcompl):
    return np.rad2deg(np.angle(vcompl[1:,:]))

# TODO : not real formula for diff phase + should do it in complex space
def corrFlux2VisPhiDif(vcompl):
    nlam=vcompl.shape[1]
    norm=np.outer(np.mean(vcompl[1:,:],axis=1),np.ones(nlam))
    
    phi= np.rad2deg(np.angle(vcompl[1:,:]*np.conjugate(norm)))
    #norm=np.outer(np.mean(phi,axis=1),np.ones(nlam))
    return phi#-norm

#TODO special function doing T3Amp and T3Phi simultaneously
def corrFlux2T3Amp(vcompl):
    nB=vcompl.shape[0]
    nCP=(nB-1)//3
    norm=np.outer(np.ones(nCP),vcompl[0,:])
    BS=vcompl[1:nCP+1,:]*vcompl[nCP+1:2*nCP+1,:]*np.conjugate(vcompl[2*nCP+1:,:])/norm**3
    return np.abs(BS)

def corrFlux2T3Phi(vcompl):
    nB=vcompl.shape[0]
    nCP=(nB-1)//3
    norm=np.outer(np.ones(nCP),vcompl[0,:])   
    BS=vcompl[1:nCP+1,:]*vcompl[nCP+1:2*nCP+1,:]*np.conjugate(vcompl[2*nCP+1:,:])/norm**3
    return np.rad2deg(np.angle(BS))

def corrFlux2Flux(vcompl):
    return np.abs(vcompl)



class OImSimulator(object):
    """
    contains 
    """
    def __init__(self,data=None,model=None,fitter=None):
        self.data=oim.OImData()
        self.simulatedData=oim.OImData()
        self.model=None
        
        if data!=None:
            self.addData(data)
            
        if model!=None:        
            self.setModel(model)
    
    def setModel(self,model):
        self.model=model
    
    def addData(self,data):
            self.data.addData(data)
            self.simulatedData.addData(data)
         
    def prepareData(self):
        self.data.prepareData()
    

    def compute(self,computeChi2=True,computeSimulatedData=True):
        self.vcompl=self.model.getComplexCoherentFlux(self.data.vect_u,self.data.vect_v,self.data.vect_wl)
       
        nelChi2=0
        chi2=0
        chi2List=[]
       
          
        if (computeChi2==True)|(computeSimulatedData==True):
            
            idx=0
            nfiles=len(self.data.struct_u)
            for ifile in range(nfiles):
                #print("Data {}".format(ifile))
                narr=len(self.data.struct_arrType[ifile])
                for iarr in range(narr):
                    
                    arrNum=self.data.struct_arrNum[ifile][iarr]
                    arrType=self.data.struct_arrType[ifile][iarr]
                    dataType=self.data.struct_dataType[ifile][iarr]
                    nB=self.data.struct_nB[ifile][iarr]
                    nwl=self.data.struct_nwl[ifile][iarr]
                    vcompli=self.vcompl[idx:idx+nB*nwl]              
                    vcompli=np.reshape(vcompli,[nB,nwl])
                    
                    
                    dataVal=self.data.struct_val[ifile][iarr]
                    dataErr=self.data.struct_err[ifile][iarr]
                        
                    idx+=nB*nwl
                    quantities=[]
                    val=[]
                    
                    #print("Arr {} idx={}".format(arrType,arrNum))
                    if arrType=="OI_VIS2":  
                        val.append(corrFlux2Vis2(vcompli))
                        quantities.append("VIS2DATA")
                          
                    #Computing observables from complex Coherent Flux
                    elif arrType=="OI_VIS":
                        if dataType&oim.OImDataType.VISAMP_ABS:
                            val.append(corrFlux2VisAmpAbs(vcompli))
                            quantities.append("VISAMP")
                        elif dataType&oim.OImDataType.VISAMP_DIF:
                            val.append(corrFlux2VisAmpDif(vcompli))
                            quantities.append("VISAMP")
                        elif dataType&oim.OImDataType.VISAMP_COR:
                            val.append(corrFlux2VisAmpCor(vcompli))
                            quantities.append("VISAMP")
                        
                        if dataType&oim.OImDataType.VISPHI_ABS:
                            val.append(corrFlux2VisPhiAbs(vcompli))
                            quantities.append("VISPHI")                        
                        elif dataType&oim.OImDataType.VISPHI_DIF:
                            val.append(corrFlux2VisPhiDif(vcompli))
                            quantities.append("VISPHI") 
                    elif arrType=="OI_T3":
                        if dataType&oim.OImDataType.T3AMP:
                            val.append(corrFlux2T3Amp(vcompli))
                            quantities.append("T3AMP") 
                        if  dataType&oim.OImDataType.T3PHI:
                            val.append(corrFlux2T3Phi(vcompli))
                            quantities.append("T3PHI")                         
                    elif arrType=="OI_FLUX":
                        val.append(corrFlux2Flux(vcompli))
                        quantities.append("FLUXDATA")
                    
                    
                    #Filling the simulatedData astropy array with the computed values
                    if computeSimulatedData==True:
                        for ival in range(len(val)):
                            self.simulatedData.data[ifile][arrNum].data[quantities[ival]]=val[ival]
                            
                    #Computing the chi2
                    if computeChi2==True:
                        for ival in range(len(val)):
                            chi2+=np.sum(((dataVal[ival]-val[ival])/dataErr[ival])**2)
                            nelChi2+=np.size(dataVal[ival])
                            chi2List.append(((dataVal[ival]-val[ival])/dataErr[ival])**2)
    
       
        if computeChi2==True: 
             self.chi2=chi2
             self.chi2r=chi2/nelChi2
             self.chi2List=chi2List
                
                    

