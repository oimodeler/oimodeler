# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 11:21:49 2022

@author: Ame
"""
import matplotlib.pyplot as plt
import matplotlib.projections as proj
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection


import numpy as np
import os
from astropy.io import fits
import oimodeler as oim

#TODO replace functions imported from external library oifitstools
import oifitstools


###############################################################################
def _errorplot(axe,X,Y,dY,smooth=1,line_kw=None, **kwargs ):
    Ys=Y
    if not("alpha" in kwargs):
        kwargs["alpha"]=0.4
    if smooth!=1:      
        ker=np.ones(smooth)/smooth
        Ys=np.convolve(Y,ker,mode="same")#[smooth//2:-smooth//2-1]
    XX=np.concatenate([X,np.flip(X)])
    YY=np.concatenate([Ys-dY,np.flip(Ys+dY)])
    axe.fill(XX,YY,**kwargs)
    if line_kw:
        axe.plot(X,Y, **line_kw)
        
        
      
###############################################################################
        
def _colorPlot(axe,x,y,z,**kwargs):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    if not("cmap" in kwargs):
        kwargs['cmap']='plasma'

    maxi=[np.max(z)]
    mini=[np.min(z)]
    for ci in axe.collections:
        maxi.append(np.max(ci.get_array()))
        mini.append(np.min(ci.get_array()))
     
    maxi=np.max(maxi)
    mini=np.min(mini)
    
    if not("norm" in kwargs):
        norm=plt.Normalize(mini,maxi)
    else:
        norm=kwargs["norm"]
        
    
    lc = LineCollection(segments, **kwargs)
    lc.set_array(z) 
    line = axe.add_collection(lc)
    
    """
    if 'label' in kwargs:
        label=kwargs.pop('label')
        print(label)
        legend=axe.get_legend_handles_labels()
        legend[0].append(lc[0])
        legend[1].append=label
    """
    
    for ci in axe.collections:
        ci.set_norm(norm)
    
    xmax=[]
    ymax=[]
    xmin=[]
    ymin=[]    
    for ci in axe.collections:
        xy=np.array(ci.get_segments())
        
        if len(np.shape(xy))==3:  
            xx=xy[:,:,0]  
            yy=xy[:,:,1]
            
            xmax.append(np.max(xx))
            ymax.append(np.max(yy))
            
            xmin.append(np.min(xx))
            ymin.append(np.min(yy))
            
    xmax=np.max(xmax)
    xmin=np.min(xmin)    
    ymax=np.max(ymax)
    ymin=np.min(ymin)           
    
    axe.set_xlim(xmin,xmax)
    axe.set_ylim(ymin,ymax)
    
    
    
    
    return line
        
        
        
###############################################################################

oimPlotParamName=np.array(["B","PA","UCOORD","VCOORD","SPAFREQ","EFF_WAVE",
                    "VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI","FLUXDATA"])
oimPlotParamError=np.array(["","","","","","EFF_BAND","VIS2ERR","VISAMPERR",
                         "VISPHIERR","T3AMPERR","T3PHIERR","FLUXERR"])
oimPlotParamArr=np.array(["","","","","","OI_WAVELENGTH","OI_VIS2","OI_VIS",
                       "OI_VIS","OI_T3","OI_T3","OI_FLUX"])
oimPlotParamLabel=np.array(["Baseline Length","Baseline Orientation","U","V",
                         "Spatial Frequency","Wavelength","Square Visibility",
                          "Differential Visibility","Differential Phase",
                          "CLosure Amplitude","Closure Phase","Flux"])
oimPlotParamLabelShort=np.array(["B","PA","U","V","B/$\lambda$","$\lambda$",
                "V$^2$","V$_{diff}$","$\phi_{diff}$","Clos. Amp.","CP","Flux"])
oimPlotParamUnit0=np.array(["m","$^o$","m","m","cycles/rad","$\mu$m","","","$^o$",
                         "","$^o$",""])
oimPlotParamIsUVcoord=np.array([1,1,1,1,1,0,0,0,0,0,0,0])

oimPlotParamColorCycle=plt.rcParams['axes.prop_cycle'].by_key()['color']

###############################################################################

def oimPlot(oifitsList,xname,yname,axe=None,xunitname=None,xunitmultiplier=1,
            yunitname=None,yunitmultiplier=1,xlim=None,ylim=None,xscale=None,
            yscale=None,shortLabel=True,color="byFile",colorTab=None,errorbar=False,
            showFlagged=False,kwargs_error={},**kwargs):

    
    res= None
    
    if type(oifitsList)!=type([]):
        oifitsList=[oifitsList]
        
    #TODO colors with lam, baselines, ...
    icol=0
    if colorTab==None:
        colorTab=oimPlotParamColorCycle
        
    if not(color in ["byBaseline","byFile","byWavelength"]):
        colorTab=[color]
   
    ncol=len(colorTab)    
        
    ndata=len(oifitsList)
 
    idxX=np.where(oimPlotParamName == xname)[0][0]
    idxY=np.where(oimPlotParamName == yname)[0][0]
    
    xerrname=oimPlotParamError[idxX]
    yerrname=oimPlotParamError[idxY]
    
    xarr=oimPlotParamArr[idxX]
    yarr=oimPlotParamArr[idxY]

    if not(shortLabel):
        xlabel=oimPlotParamLabel[idxX]
        ylabel=oimPlotParamLabel[idxY]
    else:
        xlabel=oimPlotParamLabelShort[idxX]
        ylabel=oimPlotParamLabelShort[idxY]    
        
        
    #TODO implement units with astropy
    if xunitname:
        xunit0=xunitname
    else:
        xunit0=oimPlotParamUnit0[idxX]
           
    if yunitname:
         yunit0=yunitname  
    else:
         yunit0=oimPlotParamUnit0[idxY]
    
    if xunit0!="": xlabel+=" ("+xunit0+")"
    if yunit0!="": ylabel+=" ("+yunit0+")" 
        

    xIsUVcoord=oimPlotParamIsUVcoord[idxX]
    yIsUVcoord=oimPlotParamIsUVcoord[idxY]    

    if not(xIsUVcoord):
        xdata=[d[xarr].data[xname] for d in oifitsList]
    elif xname=="SPAFREQ":
        xdata=[oim.getSpaFreq(d,arr=yarr,unit=xunitname) for d in oifitsList]
    elif xname=="UCOORD" and yname!="VCOORD":
        pass
        #TODO
    elif xname=="B":
        xdata=[np.transpose(np.tile(oim.getBaselineLengthAndPA(d,arr=yarr,unit=xunitname)[0],
                (np.shape(d[yarr].data[yname])[1],1))) for d in oifitsList]
        
    elif xname=="PA":
        xdata=[np.transpose(np.tile(oim.getBaselineLengthAndPA(d,arr=yarr,unit=xunitname)[1],
                (np.shape(d[yarr].data[yname])[1],1))) for d in oifitsList]
         

    if not(yIsUVcoord):
        ydata=[d[yarr].data[yname] for d in oifitsList]
        ydataerr=[d[yarr].data[yerrname] for d in oifitsList]
    elif yname=="SPAFREQ":
        ydata=[oim.getSpaFreq(d,arr=yarr,unit=yunitname) for d in oifitsList]
    elif yname=="VCOORD" and yname!="UCOORD":
        pass
        #TODO
    elif yname=="B":
        ydata=[np.transpose(np.tile(oim.getBaselineLengthAndPA(d,arr=yarr,unit=yunitname)[0],
                (np.shape(d[yarr].data[yname])[1],1))) for d in oifitsList]
        
    elif yname=="PA":
        ydata=[np.transpose(np.tile(oim.getBaselineLengthAndPA(d,arr=yarr,unit=yunitname)[1],
                (np.shape(d[yarr].data[yname])[1],1))) for d in oifitsList]


    if 'label' in kwargs:
        label=kwargs.pop('label')
    else:
        label=""
    if not(axe):axe=plt.axes()
       
    
    #looping through oifits files
    for idata in range(ndata):
        
        
        if color=="byWavelength":
            wl=oifitsList[idata]["OI_WAVELENGTH"].data["EFF_WAVE"]

        
        shapex=np.shape(xdata[idata])
        shapey=np.shape(ydata[idata])
        
        #Dealing with the data dimensions
        if (np.size(shapex)==np.size(shapey)):
            if np.size(shapex)==1: # if 1 baseline only just change array dim
                nlam=np.size(xdata)
                xdata=np.reshape((1,nlam))
                ydata=np.reshape((1,nlam))    
        elif (np.size(shapex))==1:# if x=1D and y=2D
            xdata[idata]=np.outer(np.ones(shapey[0]),xdata[idata])
            
        elif (np.size(shapey))==1:# if x=2D and y=1D
            ydata[idata]=np.outer(np.ones(shapex[0]),ydata[idata])

        shapex=np.shape(xdata[idata])
        shapey=np.shape(ydata[idata])
                
        # separate multiples baselines
        nB=shapex[0]
        for iB in range(nB):
            if showFlagged==False: 
                flags=oifitsList[idata][yarr].data["FLAG"][iB,:]
                nflags=len(flags)
                flag0=True
                ilam0=0
                for ilam,flagi in enumerate(flags):
                    doPlot=False              
                    if flag0!=flagi:
                        if flagi==False:
                            ilam0=ilam
                        else:
                            doPlot=True
                        flag0=flagi 
                    elif ilam==(nflags-1) and flagi==False:
                            doPlot=True
                    
                    if doPlot==True:  
                        if color!="byWavelength":
                            axe.plot(xdata[idata][iB,ilam0:ilam]*
                                 xunitmultiplier,ydata[idata][iB,ilam0:ilam],
                                 color=colorTab[icol%ncol],label=label,**kwargs)
                            if errorbar==True:
                                _errorplot(axe,xdata[idata][iB,ilam0:ilam]*xunitmultiplier,
                                            ydata[idata][iB,ilam0:ilam],
                                            ydataerr[idata][iB,ilam0:ilam],color=colorTab[icol%ncol],
                                            **kwargs_error)

                            label=None
                        else:
                            
                            res=_colorPlot(axe, xdata[idata][iB,ilam0:ilam]*
                                 xunitmultiplier, ydata[idata][iB,ilam0:ilam], wl[ilam0:ilam],
                                 label=label,**kwargs)
                            if errorbar==True:
                                _errorplot(axe,xdata[idata][iB,ilam0:ilam]*xunitmultiplier,
                                            ydata[idata][iB,ilam0:ilam],
                                            ydataerr[idata][iB,ilam0:ilam],color="gray",alpha=0.2,
                                            **kwargs_error)
                            label=None
    
                        
            else:
                axe.plot(xdata[idata][iB,:]*xunitmultiplier,
                         ydata[idata][iB,:],color=colorTab[icol%ncol])
                if errorbar==True:
                    _errorplot(axe,xdata[idata][iB,:]*xunitmultiplier,ydata[idata][iB,:],
                               ydataerr[idata][iB,:],color=colorTab[icol%ncol],**kwargs_error)
                 
           
            if color=="byBaseline":icol+=1
        if color=="byFile":icol+=1

       
    if yscale!=None:
        axe.set_yscale(yscale)

    if xscale!=None:        
        axe.set_xscale(xscale)
    
    axe.set_xlim(xlim)
    axe.set_ylim(ylim)
    axe.set_xlabel(xlabel)
    axe.set_ylabel(ylabel)
    
    return res

class HandlerColorLineCollection(HandlerLineCollection):
    def create_artists(self, legend, artist ,xdescent, ydescent,
                        width, height, fontsize,trans):
        x = np.linspace(0,width,self.get_numpoints(legend)+1)
        y = np.zeros(self.get_numpoints(legend)+1)+height/2.-ydescent
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=artist.cmap,
                     transform=trans)
        lc.set_array(x)
        lc.set_linewidth(artist.get_linewidth())
        return [lc]



###############################################################################
class oimAxes(plt.Axes):
    name = 'oimAxes'
    
    ytype=None
    xtype=None
            
    def uvplot(self,oifits, **kwargs ):
        oifitstools.uvplot(oifits,axe=self)
        
    def oiplot(self,oifitsList,xname,yname, xunit=None,**kwargs ):
          res=oimPlot(oifitsList,xname,yname,axe=self,xunitname=xunit,**kwargs) 
          self.ytype=yname
          self.xtype=xname   
          return res

    def autolim(self):
        if self.ytype in ["VIS2DATA","VISAMP","T3AMP"]:
            if self.get_yscale()=='linear':
                self.set_ylim(0,1)
            else:
                lines=self.get_lines()
                mini=np.Infinity
                for li in lines: 
                    mi= np.min(li.get_ydata())
                    if mi<mini:
                        mini=mi  
                self.set_ylim(mini,1)
                
        elif self.ytype in ["VISPHI","T3PHI"]:
            self.set_ylim(-180,180)
            
    def legend(self,**kwargs):
        h,l=self.get_legend_handles_labels()
        hmap={}
        for hi in h:
            if isinstance(hi,LineCollection):
                hmap[hi]=HandlerColorLineCollection(numpoints=100)
        super().legend(h,l,handler_map=hmap,**kwargs)            

###############################################################################

