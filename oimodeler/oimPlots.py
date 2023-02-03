# -*- coding: utf-8 -*-
"""
various plotting function and classes
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.legend_handler import HandlerLineCollection
from astropy.io import fits

from .oimData import oimData
from .oimUtils import getBaselineName, getBaselineLengthAndPA, getConfigName, getSpaFreq


def _errorplot(axe,X,Y,dY,smooth=1, **kwargs ):
    Ys=Y
    if not ("alpha" in kwargs):
        kwargs["alpha"]=0.4
    if not ("color" in kwargs):
        kwargs["color"]="grey"      
    if smooth!=1:      
        ker=np.ones(smooth)/smooth
        Ys=np.convolve(Y,ker,mode="same")#[smooth//2:-smooth//2-1]
    XX=np.concatenate([X,np.flip(X)])
    YY=np.concatenate([Ys-dY,np.flip(Ys+dY)])
    axe.fill(XX,YY,**kwargs)

        
        
      
        
def _colorPlot(axe,x,y,z,setlim=False,**kwargs):
    
    if not ("cmap" in kwargs):
        kwargs['cmap']='plasma'

    maxi=[np.max(z)]
    mini=[np.min(z)]
    for ci in axe.collections:
        maxii=np.max(ci.get_array())
        minii=np.min(ci.get_array())
        if maxii is not None and minii is not None:
            maxi.append(maxii)
            mini.append(minii)

    maxi=np.max(maxi)
    mini=np.min(mini)

    
    if not ("norm" in kwargs):
        norm=plt.Normalize(mini,maxi)
    else:
        norm=kwargs["norm"]
    
    if x.size>1:
    
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        
        

        if "marker" in kwargs:
            kwargs.pop("marker")
        
        lc = LineCollection(segments, **kwargs)
        lc.set_array(z) 
        res = axe.add_collection(lc)
        
        
    if x.size==1:
        res =axe.scatter(x,y,c=z,**kwargs)

    for ci in axe.collections:
        ci.set_norm(norm)
        

    if setlim:
        axe.autoscale_view()
    
    return res
        
###############################################################################


def uvPlot(oifits,extension="OI_VIS2",marker="o", facecolors='red',
           edgecolors='k',size=10,axe=None,maxi=None,xytitle=[True,True],
           title=None,gridcolor="k",grid=True,fontsize=None,**kwargs):

    
    if isinstance(oifits,oimData):
        oifits=oifits.data
    
    kwargs2={}
    for key in kwargs:
        if key!="label":
            kwargs2[key]=kwargs[key]
        
    data=[]
    u=np.array([])
    v=np.array([])
    if isinstance(oifits, str):
        data.append(fits.open(oifits))
    elif isinstance(oifits, list):
        for item in oifits:
            if isinstance(item, str):
                data.append(fits.open(item))
            else:
                data.append(item)
    else:
        data.append(oifits)

    for datai in data:

        extnames=np.array([datai[i].name for i in range(len(datai))])
        idx=np.where(extnames==extension)[0]
        for j in idx:
            u=np.append(u,datai[j].data['UCOORD'])
            v=np.append(v,datai[j].data['VCOORD'])


    if not axe:
        fig,axe=plt.subplots(nrows=1,ncols=1)
    axe.scatter(u,v,marker=marker,facecolors=facecolors, edgecolors=edgecolors,
                s=size,zorder=10,lw=1,**kwargs)

    axe.scatter(-u,-v,marker=marker,facecolors=facecolors, edgecolors=edgecolors,
                s=size,zorder=10,lw=1,**kwargs2)
    if not maxi:
        maxi=1.1*np.max(np.abs(np.array([u,v])))
    if grid:
        axe.plot([-maxi,maxi],[0,0],linewidth=1,color=gridcolor,zorder=5)
        axe.plot([0,0],[-maxi,maxi],linewidth=1,color=gridcolor,zorder=5)
    axe.set_aspect('equal','box')
    axe.set_xlim([maxi,-maxi])
    axe.set_ylim([-maxi,maxi])
    if xytitle[0]:
        axe.set_xlabel('u (m)',fontsize=fontsize)
    if xytitle[1]:
        axe.set_ylabel('v (m)',fontsize=fontsize)
    if title:
        axe.set_title(title)
    return [axe]

###############################################################################

oimPlotParamName=np.array(["LENGTH","PA","UCOORD","VCOORD","SPAFREQ","EFF_WAVE",
                    "VIS2DATA","VISAMP","VISPHI","T3AMP","T3PHI","FLUXDATA","MJD"])
oimPlotParamError=np.array(["","","","","","EFF_BAND","VIS2ERR","VISAMPERR",
                         "VISPHIERR","T3AMPERR","T3PHIERR","FLUXERR",""])
oimPlotParamArr=np.array(["","","","","","OI_WAVELENGTH","OI_VIS2","OI_VIS",
                       "OI_VIS","OI_T3","OI_T3","OI_FLUX",""])
oimPlotParamLabel=np.array(["Baseline Length","Baseline Orientation","U","V",
                         "Spatial Frequency","Wavelength","Square Visibility",
                          "Differential Visibility","Differential Phase",
                          "CLosure Amplitude","Closure Phase","Flux","MJD"])
oimPlotParamLabelShort=np.array(["B","PA","U","V","B/$\lambda$","$\lambda$",
                "V$^2$","V$_{diff}$","$\phi_{diff}$","Clos. Amp.","CP","Flux","MJD"])
oimPlotParamUnit0=np.array(["m","$^o$","m","m","cycles/rad","$\mu$m","","","$^o$",
                         "","$^o$","","day"])
oimPlotParamIsUVcoord=np.array([1,1,1,1,1,0,0,0,0,0,0,0,0])

oimPlotParamColorCycle=plt.rcParams['axes.prop_cycle'].by_key()['color']


###############################################################################
def getColorIndices(oifitsList,color,yarr,yname):
    idx=[]
    names=[]
    for idata,datai in enumerate(oifitsList):
        extnames=np.array([dj.name for dj in datai])
        idx_yext=np.where(extnames==yarr)[0]
        idxi=[]
        for j,jdata in enumerate(idx_yext):
            val=datai[jdata].data[yname]
            nB=val.shape[0]
            
            if color=="byFile":
                fname=datai.filename()
                if fname is None:
                    fname="File {}".format(idata)
                else:
                    fname=os.path.basename(fname)
                    
                if fname in names:
                    ifname = names.index(fname)
                else:
                    ifname=len(names)
                    names.append(fname)
                    
                idxi.append(np.zeros(nB,dtype=int)+ifname)

                
            elif color=="byArrname":
                array=datai[jdata].header['ARRNAME']
                if array in names:
                    iarr = names.index(array)
                else:
                    iarr=len(names)
                    names.append(array)
                idxi.append(np.zeros(nB,dtype=int)+iarr)    
                
                
            elif color=="byBaseline":
                    bnames=getBaselineName(datai,yarr,squeeze=False)[j]
                    idxj=[]
                    for bname in bnames:
                        if bname in names:
                            iB = names.index(bname)
                        else:
                            iB=len(names)
                            names.append(bname)   
                        idxj.append(iB)  
                    idxi.append(idxj)
                    
       
            elif color=="byConfiguration":
                    conf=getConfigName(datai,yarr,squeeze=False)[j]
                    if conf in names:
                        iconf = names.index(conf)
                    else:
                        iconf=len(names)
                        names.append(conf)
                    idxi.append(np.zeros(nB,dtype=int)+iconf)

            else:
                idxi.append(np.zeros(nB,dtype=int))
                names.append("")
        idx.append(idxi)
    return idx,names
###############################################################################

def oimPlot(oifitsList,xname,yname,axe=None,xunit=None,xunitmultiplier=1,
            yunit=None,yunitmultiplier=1,cname=None,cunit=None,cunitmultiplier=1,
            xlim=None,ylim=None,xscale=None,yscale=None,shortLabel=True,
            color=None,colorTab=None,errorbar=False,showFlagged=False,
            kwargs_error={},**kwargs):
    """
    

    Parameters
    ----------
    oifitsList : TYPE
        DESCRIPTION.
    xname : TYPE
        DESCRIPTION.
    yname : TYPE
        DESCRIPTION.
    axe : TYPE, optional
        DESCRIPTION. The default is None.
    xunit : TYPE, optional
        DESCRIPTION. The default is None.
    xunitmultiplier : TYPE, optional
        DESCRIPTION. The default is 1.
    yunit : TYPE, optional
        DESCRIPTION. The default is None.
    yunitmultiplier : TYPE, optional
        DESCRIPTION. The default is 1.
    cname : TYPE, optional
        DESCRIPTION. The default is None.
    cunit : TYPE, optional
        DESCRIPTION. The default is None.
    cunitmultiplier : TYPE, optional
        DESCRIPTION. The default is 1.
    xlim : TYPE, optional
        DESCRIPTION. The default is None.
    ylim : TYPE, optional
        DESCRIPTION. The default is None.
    xscale : TYPE, optional
        DESCRIPTION. The default is None.
    yscale : TYPE, optional
        DESCRIPTION. The default is None.
    shortLabel : TYPE, optional
        DESCRIPTION. The default is True.
    color : TYPE, optional
        DESCRIPTION. The default is None.
    colorTab : TYPE, optional
        DESCRIPTION. The default is None.
    errorbar : TYPE, optional
        DESCRIPTION. The default is False.
    showFlagged : TYPE, optional
        DESCRIPTION. The default is False.
    kwargs_error : TYPE, optional
        DESCRIPTION. The default is {}.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    
    res= None

    
    if isinstance(oifitsList,oimData):
        oifitsList=oifitsList.data
    
    if not isinstance(oifitsList, list):
        oifitsList=[oifitsList]
         
    ndata0=len(oifitsList)
 
    idxX=np.where(oimPlotParamName == xname)[0][0]
    idxY=np.where(oimPlotParamName == yname)[0][0]

    
    xerrname=oimPlotParamError[idxX]
    yerrname=oimPlotParamError[idxY]
    
    xarr=oimPlotParamArr[idxX]
    yarr=oimPlotParamArr[idxY]

    if not shortLabel:
        xlabel=oimPlotParamLabel[idxX]
        ylabel=oimPlotParamLabel[idxY]
    else:
        xlabel=oimPlotParamLabelShort[idxX]
        ylabel=oimPlotParamLabelShort[idxY]    
        
        
    # TODO implement units with astropy
    if xunit:
        xunit0=xunit
    else:
        xunit0=oimPlotParamUnit0[idxX]
           
    if yunit:
         yunit0=yunit  
    else:
         yunit0=oimPlotParamUnit0[idxY]
    
    if xunit0!="": xlabel+=" ("+xunit0+")"
    if yunit0!="": ylabel+=" ("+yunit0+")" 
        

    xIsUVcoord=oimPlotParamIsUVcoord[idxX]
    yIsUVcoord=oimPlotParamIsUVcoord[idxY]
    
    if not xIsUVcoord and xname!="EFF_WAVE":
        raise TypeError("X should be LENGTH, SPAFREQ, PA or EFF_WAVE")
  
    if yIsUVcoord:
        raise TypeError("Y shouldn't be UCOORD,VCOORD, SPAFREQ, or PA")
        
        
    
    if colorTab is None:
        colorTab=oimPlotParamColorCycle
    
    try :
        if not("by" in color):
            colorTab=[color]
    except Exception:
        pass
    
    ncol=len(colorTab)    
      
    colorIdx,ColorNames=getColorIndices(oifitsList,color,yarr,yname)


    if 'label' in kwargs:
        label=kwargs.pop('label')
    else:
        label=""
    if not axe:
        axe=plt.axes()
    
        
    for ifile,data in enumerate(oifitsList):
        # NOTE: yname can be anything but  UCOORD, VCOORD, LENGTH, SPAFREQ, PA or EFF_WAVE
        extnames=np.array([di.name for di in data])
        
        idx_yext=np.where(extnames==yarr)[0]
        yinsname=np.array([data[i].header['INSNAME'] for i in idx_yext])
        ydata=[data[i].data[yname] for i in idx_yext]
        ydataerr=[data[i].data[yerrname] for i in idx_yext]
        yflag=[data[i].data["FLAG"] for i in idx_yext] 
        # NOTE: xname can be LENGTH, SPAFREQ, PA or EFF_WAVE
        if xname=="EFF_WAVE":
            idx_xext=np.where(extnames==xarr)[0]
            xinsname=np.array([data[i].header['INSNAME'] for i in idx_xext])
            xdata=[]
            for idata in range(len(idx_yext)):
                iwlarr=idx_xext[np.where(xinsname==yinsname[idata])[0][0]]
                wl=data[iwlarr].data["EFF_WAVE"]
                nB=ydata[idata].shape[0]
                xdata.append(np.tile(wl[None,:], (nB,1)))
                
        elif xname=="SPAFREQ":
            xdata=getSpaFreq(data,arr=yarr,unit=xunit,squeeze=False) 
            
        elif xname=="LENGTH":
            B=getBaselineLengthAndPA(data,arr=yarr,squeeze=False)[0]
            xdata=[]
            for idata in range(len(idx_yext)):
                xdata.append(np.transpose(
                    np.tile(B[idata],(np.shape(data[idata].data[yname])[1],1))))

        elif xname=="PA":
            PA=getBaselineLengthAndPA(data,arr=yarr,squeeze=False)[1]
            xdata=[]
            for idata in range(len(idx_yext)):
                xdata.append(np.transpose(
                    np.tile(PA[idata],(np.shape(data[idata].data[yname])[1],1))))
         

        if cname is not None:
            
            idxC=np.where(oimPlotParamName == cname)[0][0]
            cIsUVcoord=oimPlotParamIsUVcoord[idxC] 
            
            if cunit:
                cunit0=cunit
            else:
                cunit0=oimPlotParamUnit0[idxC]
            
            if cname=="MJD":
                carr=yarr
                cdata=[data[i].data[cname] for i in idx_yext]
                
            elif cname=="EFF_WAVE":
                carr="OI_WAVELENGTH"
                idx_cext=np.where(extnames==carr)[0]
                cinsname=np.array([data[i].header['INSNAME'] for i in idx_cext])
                cdata=[]
                for idata in range(len(idx_yext)):
                    iwlarr=idx_cext[np.where(cinsname==yinsname[idata])[0][0]]
                    wl=data[iwlarr].data["EFF_WAVE"]
                    nB=ydata[idata].shape[0]
                    cdata.append(np.tile(wl[None,:], (nB,1)))

            elif not cIsUVcoord:
                try:
                    idxC=np.where(oimPlotParamName == cname)[0][0]
                    carr=oimPlotParamArr[idxC]
                    cdata=[data[i].data[cname] for i in idx_yext]
                except Exception:
                    idxC=np.where(oimPlotParamError == cname)[0][0]
                    carr=oimPlotParamArr[idxC]
                    cdata=[data[i].data[cname] for i in idx_yext]   
                
            elif cname=="SPAFREQ":
                cdata=getSpaFreq(data,arr=yarr,unit=cunit,squeez=False) 
                
            elif cname=="LENGTH":
                B=getBaselineLengthAndPA(data,arr=yarr,squeeze=False)[0]
                cdata=[]
                for idata in range(len(idx_yext)):
                    cdata.append(np.transpose(
                        np.tile(B[idata],(np.shape(data[idata].data[yname])[1],1))))
    
            elif cname=="PA":
                PA=getBaselineLengthAndPA(data,arr=yarr,squeeze=False)[1]
                cdata=[]
                for idata in range(len(idx_yext)):
                    cdata.append(np.transpose(
                        np.tile(PA[idata],(np.shape(data[idata].data[yname])[1],1))))
        
        # NOTE: looping through oifits files
        
        ndata=len(ydata)
        for idata in range(ndata):

            shapex=np.shape(xdata[idata])
            shapey=np.shape(ydata[idata])

            # NOTE: Dealing with the xy data dimensions
            if (np.size(shapex)==np.size(shapey)):
                if np.size(shapex)==1: # NOTE: if 1 baseline only just change array dim
                    nlam=np.size(xdata)
                    xdata=np.reshape((1,nlam))
                    ydata=np.reshape((1,nlam))    
            elif (np.size(shapex))==1:# NOTE: if x=1D and y=2D
                if shapex[0]==shapey[0]:
                    xdata[idata]=np.outer(xdata[idata],np.ones(shapey[1]))
                else:
                    xdata[idata]=np.outer(np.ones(shapey[0]),xdata[idata])
                
            elif (np.size(shapey))==1:# NOTE: if x=2D and y=1D
                if shapex[0]==shapey[0]:
                    ydata[idata]=np.outer(ydata[idata],np.ones(shapex[1]))
                    ydataerr[idata]=np.outer(ydataerr[idata],np.ones(shapex[1]))
                    
                else:
                    ydata[idata]=np.outer(np.ones(shapex[0]),ydata[idata])
                    ydataerr[idata]=np.outer(ydataerr[idata],np.ones(shapex[1]))
                    
    
            shapex=np.shape(xdata[idata])
            shapey=np.shape(ydata[idata])
                    
            if cname is not None:
                shapec=np.shape(cdata[idata])
                if (np.size(shapec)==1):
                    if shapec[0]==shapex[0]:
                        cdata[idata]=np.outer(cdata[idata],np.ones(shapex[1]))
                    else:
                        cdata[idata]=np.outer(np.ones(shapex[0]),cdata[idata])
                    shapec=np.shape(cdata[idata])
            # NOTE: separate multiples baselines
            nB=shapex[0]
            

            for iB in range(nB):
                if not showFlagged: 
                    flags=np.reshape(yflag[idata],shapey)[iB,:]
                    nflags=len(flags)
                    flag0=True
                    ilam0=0
                    for ilam,flagi in enumerate(flags):
                        doPlot=False
                        if np.isnan(ydata[idata][iB,ilam]):
                            flagi=True
                        if flag0!=flagi:
                            if not flagi:
                                ilam0=ilam
                            else:
                                doPlot=True
                            flag0=flagi 
                        if ilam==(nflags-1) and not flagi:
                                doPlot=True
     
                        if doPlot: 
                            labeli=label+ColorNames[colorIdx[ifile][idata][iB]]
                            
                            if cname is None:
                                if (xdata[idata][iB,ilam0:ilam+1]).size==1:
                                    axe.scatter(xdata[idata][iB,ilam0:ilam+1]*
                                         xunitmultiplier, ydata[idata][iB,ilam0:ilam+1],
                                         color=colorTab[colorIdx[ifile][idata][iB]%ncol],
                                         label=labeli,**kwargs)
                                else: 
                                
                                    axe.plot(xdata[idata][iB,ilam0:ilam+1]*
                                         xunitmultiplier,ydata[idata][iB,ilam0:ilam+1],
                                         color=colorTab[colorIdx[ifile][idata][iB]%ncol],
                                         label=labeli,**kwargs)
                                    if errorbar:
                                        
                                        if not('color' in kwargs_error):
                                            kwargs_errori=kwargs_error.copy()
                                            kwargs_errori['color']=colorTab[colorIdx[ifile][idata][iB]%ncol] 
                                            
                                           
                                        _errorplot(axe,xdata[idata][iB,ilam0:ilam+1]*xunitmultiplier,
                                                    ydata[idata][iB,ilam0:ilam+1],
                                                    ydataerr[idata][iB,ilam0:ilam+1],
                                                    **kwargs_errori)
                                    
                            else:
                               
                               
                                # NOTE: dummy plot with alpha=0 as _colorPLot works with collections
                                # NOTE: thus not updating the xlim and ylim automatically
                                axe.plot(xdata[idata][iB,ilam0:ilam+1]*
                                     xunitmultiplier,ydata[idata][iB,ilam0:ilam+1],
                                     color="k",alpha=0)
                                
                                res=_colorPlot(axe, xdata[idata][iB,ilam0:ilam+1]*
                                     xunitmultiplier, ydata[idata][iB,ilam0:ilam+1],
                                     cdata[idata][iB,ilam0:ilam+1]*cunitmultiplier,
                                     label=labeli,setlim=False,**kwargs)
                                
                                if errorbar==True:
                                    _errorplot(axe,xdata[idata][iB,ilam0:ilam+1]*xunitmultiplier,
                                                ydata[idata][iB,ilam0:ilam+1],
                                                ydataerr[idata][iB,ilam0:ilam+1],
                                                **kwargs_error)
                        
                else:
                    axe.plot(xdata[idata][iB,:]*xunitmultiplier,
                             ydata[idata][iB,:],color=colorTab[colorIdx[ifile][idata][iB]%ncol])
                    if errorbar:
                        _errorplot(axe,xdata[idata][iB,:]*xunitmultiplier,ydata[idata][iB,:],
                                   ydataerr[idata][iB,:],color=colorTab[colorIdx[ifile][idata][iB]%ncol],**kwargs_error)
                     

       
    if yscale is not None:
        axe.set_yscale(yscale)

    if xscale is not None:        
        axe.set_xscale(xscale)
        
    if not (xlim is None):
        axe.set_xlim(xlim)
    if not (ylim is None):
        axe.set_ylim(ylim)
        
    axe.set_xlabel(xlabel)
    axe.set_ylabel(ylabel)
    
    
    return res

class _HandlerColorLineCollection(HandlerLineCollection):
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



class oimAxes(plt.Axes):
    
    """
    Class derived from plt.Axes that allows easy plotting of oifits data
    """
    name = 'oimAxes'
    
    ytype=None
    xtype=None
            
    def uvplot(self,oifits, **kwargs ):
        uvPlot(oifits,axe=self,**kwargs)
        
    def oiplot(self,oifitsList,xname,yname,**kwargs ):
          res=oimPlot(oifitsList,xname,yname,axe=self,**kwargs) 
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
                    if mi<mini and mi>0:
                        mini=mi  
                
                self.set_ylim(mini,1)
                
        elif self.ytype in ["VISPHI","T3PHI"]:
            self.set_ylim(-180,180)
            
    def legend(self,**kwargs):
        h,l=self.get_legend_handles_labels()
        hmap={}
        for hi in h:
            if isinstance(hi,LineCollection):
                hmap[hi]=_HandlerColorLineCollection(numpoints=100)
        
        # NOTE: use to remove duplicate legend
        lh = dict(zip(l, h))
        super().legend(lh.values(),lh.keys(),handler_map=hmap,**kwargs)   


    def set_yscale(self,value,**kwargs):
        super().set_yscale(value,**kwargs)
        self.autoscale_view()

    def set_xscale(self,value,**kwargs):
        super().set_xscale(value,**kwargs)
        self.autoscale_view()


