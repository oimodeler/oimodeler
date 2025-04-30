# -*- coding: utf-8 -*-
"""
various plotting function and classes
"""
import os
from typing import Dict, List, Optional, Tuple, Union

import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.collections import Collection, LineCollection
from matplotlib.legend_handler import HandlerLineCollection

from .oimData import loadOifitsData, oimData
from .oimUtils import getBaselineName, getBaselineLengthAndPA, getConfigName, \
    getSpaFreq, getWlFromOifits, getDataArrname


# TODO: Move global variables somewhere else and make their name in capital letters
# (python standard)
oimPlotParamName = np.array(["LENGTH", "PA", "UCOORD", "VCOORD", "SPAFREQ", "EFF_WAVE",
                             "VIS2DATA", "VISAMP", "VISPHI", "T3AMP", "T3PHI", "FLUXDATA", "MJD","VISDATA"])
oimPlotParamError = np.array(["", "", "", "", "", "EFF_BAND", "VIS2ERR", "VISAMPERR",
                              "VISPHIERR", "T3AMPERR", "T3PHIERR", "FLUXERR", "",""])
oimPlotParamArr = np.array(["", "", "", "", "", "OI_WAVELENGTH", "OI_VIS2", "OI_VIS",
                            "OI_VIS", "OI_T3", "OI_T3", "OI_FLUX", "","OI_VIS2"])
oimPlotParamLabel = np.array(["Baseline Length", "Baseline Orientation", "U", "V",
                              "Spatial Frequency", "Wavelength", "Square Visibility",
                             "Visibility Amplitude", "Phase",
                              "Closure Amplitude", "Closure Phase", "Flux", "MJD","Visibility"])
oimPlotParamLabelShort = np.array(["B", "PA", "U", "V", r"B/$\lambda$", r"$\lambda$",
                                   "V$^2$", "Vis. Amp.", "Vis. Phi.", "Clos. Amp.", "CP", "Flux", "MJD","V"])
oimPlotParamUnit0 = np.array([u.m, u.deg, u.m, u.m, u.Unit("cycle/rad"), u.m, u.one, u.one, u.deg,
                              u.one, u.deg, u.one, u.day,u.one])
oimPlotParamIsUVcoord = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,0])

oimPlotParamColorCycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]


# NOTE: might be useful as unit cycles is not defined but cycle is
# u.add_enabled_units(u.def_unit("cycles",u.cycle))
def _errorplot(axe, x: np.ndarray, y: np.ndarray,
               dy: np.ndarray, smooth: Optional[int] = 1, **kwargs):
    """Creates a plot with error bars.

    Parameters
    ----------
    axe : matplotlib.axes.Axes
        The axes to plot on.
    x : numpy.ndarray
        The x data.
    y : numpy.ndarray
        The y data.
    dy : np.ndarray
        The error on the y data.
    smooth : int, optional
        The smoothing factor. The default is 1.
    """
    ys = y
    if "alpha" not in kwargs:
        kwargs["alpha"] = 0.4

    if "color" not in kwargs:
        kwargs["color"] = "grey"

    if smooth != 1:
        ker = np.ones(smooth)/smooth
        ys = np.convolve(y, ker, mode="same")  # [smooth//2:-smooth//2-1]


    axe.fill_between(x,ys-dy,ys+dy, **kwargs)

def _colorPlot(axe, x: np.ndarray,
               y: np.ndarray, z: np.ndarray, flag:np.ndarray=None,
               setlim: Optional[bool] = False, **kwargs) -> Collection:
    """Creates a plot of the x and y data with the z data as color.

    Parameters
    ----------
    axe : matplotlib.axes.Axes
        The axes to plot on.
    x : np.ndarray
        The x data.
    y : numpy.ndarray
        The y data.
    z : numpy.ndarray
        The colorbar data.
    flag : None or numpy.ndarray, optional
        Flagging of bad data (not to be plotted). The default is None. 
    setlim : bool, optional
        Automatically sets the plot's limit. The default is False.

    Returns
    -------
    res : matplotlib.collections.Collection
        The collection of the plot.
    """
    noline=False
    if "ls" in kwargs:
        if  kwargs["ls"]==None or kwargs["ls"]=="":
            noline=True
    if "linestyle" in kwargs:
        if kwargs["linestyle"]==None or kwargs["linestyle"]=="":
            noline=True
        
    if "cmap" not in kwargs:
        kwargs["cmap"] = "plasma"

    maxi = [np.max(z)]
    mini = [np.min(z)]
    for ci in axe.collections:
        maxii = np.max(ci.get_array())
        minii = np.min(ci.get_array())
        if maxii is not None and minii is not None:
            maxi.append(maxii)
            mini.append(minii)
            
    maxi = np.max(maxi)
    mini = np.min(mini)
    
    if type(flag) == type(None):
        flag=(0*x).astype(bool)

    yma=np.ma.masked_where(flag, y)


    if "norm" not in kwargs:
        norm = plt.Normalize(mini, maxi)
    else:
        norm = kwargs["norm"]

    res = None
    
    if x.size == 1 or  "marker" in kwargs:
        res = axe.scatter(x, yma, c=z, **kwargs)
    if x.size > 1 and noline==False:
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        flag_seg = flag[1:] | flag[0:-1]
        segments_flagged=segments[np.logical_not(flag_seg),:,:]
        z_seg=np.mean(np.array([z[:-1], z[1:]]),axis=0)
        z_flagged=z_seg[np.logical_not(flag_seg)]

        if "marker" in kwargs:
            kwargs.pop("marker")
        
        lc = LineCollection(segments_flagged, **kwargs)
        lc.set_array(z_flagged)
        res = axe.add_collection(lc)

    for ci in axe.collections:
        ci.set_norm(norm)

    if setlim:
        axe.autoscale_view()

    return res

def uvPlot(oifitsList: fits.HDUList, arrname: Optional[str] = "OI_VIS2",
           unit: Optional[u.Quantity] = u.m,
           stringunitformat: Optional[str] = "latex_inline", 
           color: Optional[str] = None, maxi: Optional[float] = None,
           grid: Optional[bool] = True, gridcolor: Optional[str] = "k",
           fontsize: Optional[int] = None, xytitle: Optional[List[bool]] = [True, True],
           showLegend: Optional[bool] = True, showColorbar: Optional[bool] = True,
           colorTab: Optional = None, axe: Optional = None,
           title: Optional[str] = None, cunit: Optional[u.Quantity] = u.m,
           legendkwargs: Optional[Dict] = {}, **kwargs) -> Axes:
    """Plot the uv coverage of the data.

    Parameters
    ----------
    oifitsList : astropy.io.fits.HDUList
        A list of oifits files.
    arrname : str, optional
        The name of the array. The default is "OI_VIS2".
    unit : astropy.units.Quantity, optional
        The unit of the data. The default is u.m.
    stringunitformat : str, optional
        The unit format. The default is "latex_inline".
    color : str, optional
        The color of the plot. The default is None.
    maxi : float, optional
        The maximum value of the plot. The default is None.
    grid : bool, optional
        Show the grid. The default is True.
    gridcolor : str, optional
        The color of the grid. The default is "k".
    fontsize : int, optional
        The fontsize of the labels. The default is None.
    xytitle : list of bool, optional
        Show the x and y labels. The default is [True, True].
    showLegend : bool, optional
        Show the legend. The default is True.
    showColorbar : bool, optional
        Show the colorbar. The default is True.

    Returns
    -------
    axe : matplotlib.axes.Axes
        The axes of the plot.
    """
    if axe is None:
        _, axe = plt.subplots(nrows=1, ncols=1)

    oifitsList = loadOifitsData(oifitsList)
    # nfiles = len(oifitsList)

    unit, cunit = u.Unit(unit), u.Unit(cunit)
    if colorTab is None:
        colorTab = oimPlotParamColorCycle

    try:
        if "by" not in color:
            colorTab = [color]
    except:
        pass

    try:
        label0 = kwargs.pop("label")
    except:
        label0 = ""

    wl, ucoord, vcoord = [], [], []

    try:
        if "by" not in color:
            colorTab = [color]
    except:
        pass

    if colorTab is None:
        colorTab = oimPlotParamColorCycle

    try:
        if "by" not in color:
            colorTab = [color]
    except:
        pass

    ncol = len(colorTab)
    colorIdx, ColorNames = getColorIndices(oifitsList, color, arrname, "UCOORD", flatten=True)

    for datai in oifitsList:
        _,_, ui, vi = getBaselineLengthAndPA(datai, arr=arrname, squeeze=False, returnUV=True)
        for iext in range(len(ui)):
            wlii = getWlFromOifits(datai, arr=arrname, extver=iext+1)
            nB=ui[iext].size
            ucoord.extend(ui[iext])
            vcoord.extend(vi[iext])
            for iB in range(nB):
                wl.append(wlii)

    ucoord = np.array(ucoord)
    vcoord = np.array(vcoord)

    if unit.is_equivalent(u.m):
        ucoord *= u.m.to(unit)
        vcoord *= u.m.to(unit)

        for icol in np.unique(colorIdx):
            idx = np.where(colorIdx == icol)
            axe.scatter(ucoord[idx], vcoord[idx], color=colorTab[icol%ncol], label=label0+ColorNames[icol], **kwargs)
            axe.scatter(-ucoord[idx], -vcoord[idx], color=colorTab[icol%ncol], **kwargs)

        if maxi is None:
            maxi = 1.1*np.max(np.abs(np.array([ucoord, vcoord])))

    elif unit.is_equivalent(u.Unit("cycle/rad")):
        for iB, ui, vi, wli in zip(range(len(ucoord)), ucoord, vcoord, wl):
            spfu = ui/wli*u.Unit("cycle/rad").to(unit)
            spfv = vi/wli*u.Unit("cycle/rad").to(unit)

            if color is None:
                res = _colorPlot(axe, spfu, spfv, wli*u.m.to(cunit), label=label0, setlim=True, **kwargs)
                _colorPlot(axe, -spfu, -spfv, wli*u.m.to(cunit), setlim=True, **kwargs)

            else:
                icol=colorIdx[iB]
                axe.plot(spfu, spfv,label=label0+ColorNames[icol], color=colorTab[icol%ncol], **kwargs)
                axe.plot(-spfu, -spfv,color=colorTab[icol%ncol], **kwargs)

    else:
        raise TypeError("invalid unit")
                   
    if maxi is None:
        xmax, ymax = -1e99, -1e99
        for line in axe.lines:
            x, y = line.get_data()
            x0, y0 = np.max(x), np.max(y)
            xmax, ymax = max(x0, xmax), max(y0, ymax)

        for collection in axe.collections:
            datalim = collection.get_datalim(axe.transData)
            x0 = np.max(np.asarray(datalim)[:,0])
            y0 = np.max(np.asarray(datalim)[:,1])
            xmax, ymax = max(x0, xmax), max(y0, ymax)

        maxi = 1.1*max(xmax,ymax)

    if grid:
        axe.plot([-maxi, maxi], [0, 0], linewidth=1, color=gridcolor, zorder=5)
        axe.plot([0, 0], [-maxi, maxi], linewidth=1, color=gridcolor, zorder=5)
    axe.set_aspect('equal', 'box')
    axe.set_xlim([maxi, -maxi])
    axe.set_ylim([-maxi, maxi])

    xyunittext = unit.to_string(stringunitformat)
    if xytitle[0]:
        axe.set_xlabel(f'u ({xyunittext})', fontsize=fontsize)
    if xytitle[1]:
        axe.set_ylabel(f'v ({xyunittext})', fontsize=fontsize)
    if title:
        axe.set_title(title)

    if unit.is_equivalent(u.Unit("cycle/rad")) and color is None and showColorbar:
        plt.colorbar(res, ax=axe,label=f"$\\lambda$ ({cunit:latex})")

    if showLegend:
        axe.legend(**legendkwargs)
    return axe


def getColorIndices(oifitsList: fits.HDUList, color: str,
                    yarr: np.ndarray, yname: str,
                    flatten: Optional[bool] = False) -> Tuple[np.ndarray, List[str]]:
    """Get the color indices of the data depending on the colorisation.

    Parameters
    ----------
    oifitsList : fits.HDUList
        A list of oifits files.
    color : str
        The color of the plot.
    yarr : numpy.ndarray
        The data.
    yname : str
        The data name.
    flatten : bool, optional
        If toggled, flattens the data. The default is False.

    Returns
    -------
    idx : numpy.ndarray
        The color indices.
    names : list of str
        The names of the colors.
    """
    idx, names = [], []
    for idata, datai in enumerate(oifitsList):
        extnames = np.array([dj.name for dj in datai])
        idx_yext = np.where(extnames == yarr)[0]
        idxi = []
        for j, jdata in enumerate(idx_yext):
            val = datai[jdata].data[yname]
            nB = val.shape[0]

            if color == "byFile":
                fname = datai.filename()
                if fname == None:
                    fname = "File {}".format(idata)
                else:
                    fname = os.path.basename(fname)

                if fname in names:
                    ifname = names.index(fname)
                else:
                    ifname = len(names)
                    names.append(fname)

                idxi.append(np.zeros(nB, dtype=int)+ifname)

            elif color == "byArrname":
                array = datai[jdata].header['ARRNAME']
                if array in names:
                    iarr = names.index(array)
                else:
                    iarr = len(names)
                    names.append(array)
                idxi.append(np.zeros(nB, dtype=int)+iarr)

            elif color == "byBaseline":
                bnames = getBaselineName(datai, yarr, squeeze=False)[j]
                idxj = []
                for bname in bnames:
                    if bname in names:
                        iB = names.index(bname)
                    else:
                        iB = len(names)
                        names.append(bname)
                    idxj.append(iB)
                idxi.append(idxj)

            elif color == "byConfiguration":
                conf = getConfigName(datai, yarr, squeeze=False)[j]
                if conf in names:
                    iconf = names.index(conf)
                else:
                    iconf = len(names)
                    names.append(conf)
                idxi.append(np.zeros(nB, dtype=int)+iconf)

            else:
                idxi.append(np.zeros(nB, dtype=int))
                names.append("")
        idx.append(idxi)

    #flatten version for the uvplot
    if flatten:
        idxf = []

        for ifile in range(len(idx)):
            for iext in range(len(idx[ifile])):
                    idxf.extend(idx[ifile][iext])
        idx = idxf

    return idx, names


def oimPlot(oifitsList: fits.HDUList,
            xname: str, yname: str, axe: Optional[Axes] = None,
            xunit: Optional[u.Quantity] = None,
            yunit: Optional[u.Quantity] = None,
            cname: Optional[str] = None, cunit: Optional[str] = None,
            xlim: Optional[float] = None, ylim: Optional[float] = None,
            xscale: Optional[str] = None, yscale: Optional[str] = None,
            shortLabel: Optional[bool] = True,
            color: Optional[str] = None,
            colorTab: Optional[str] = None,
            showColorbar: Optional[bool] = True,
            errorbar: Optional[bool] = False,
            showFlagged: Optional[bool] = False,
            colorbar: Optional[bool] = True,
            kwargs_error: Optional[Dict] = {},**kwargs):
    """Plot the data from the oifits files.

    Parameters
    ----------
    oifitsList : fits.HDUList
        A list of oifits files.
    xname : str
        The name of the x-axis parameter.
    yname : str
        The name of the y-axis parameter.
    axe : matplotlib.axes.Axes, optional
        The axes to plot on. The default is None.
    xunit : astropy.units.Quantity, optional
        The unit of the x-axis parameter. The default is None.
    yunit : astropy.units.Quantity, optional
        The unit of the y-axis parameter. The default is None.
    cname : str, optional
        The name of the colorbar. The default is None.
    cunit : astropy.units.Quantity, optional
        The unit of the colorbar. The default is None.
    xlim : float, optional
        The limit of the x-axis. The default is None.
    ylim : float, optional
        The limit of the y-axis. The default is None.
    xscale : str, optional
        The scale of the x-axis. The default is None.
    yscale : str, optional
        The scale of the y-axis. The default is None.
    shortLabel : str, optional
        The label format. The default is True.
    color : str, optional
        The color of the plot. The default is None.
    colorTab : str, optional
        The color table. The default is None.
    errorbar : bool, optional
        If toggled shows error bars. The default is False.
    showFlagged : bool, optional
        If toggled shows flagged data. The default is False.
    kwargs_error : dict, optional
        Additional keyword arguments for the error plot. The default is {}.
    **kwargs : dict
        Keyword arguments for the plot.

    Returns
    -------
    res : TYPE
        DESCRIPTION.
    """
    res = None
    oifitsList = loadOifitsData(oifitsList)
    # ndata0 = len(oifitsList)

    idxX = np.where(oimPlotParamName == xname)[0][0]
    idxY = np.where(oimPlotParamName == yname)[0][0]

    # xerrname = oimPlotParamError[idxX]
    yerrname = oimPlotParamError[idxY]

    xarr = oimPlotParamArr[idxX]
    yarr = oimPlotParamArr[idxY]

    if not shortLabel:
        xlabel = oimPlotParamLabel[idxX]
        ylabel = oimPlotParamLabel[idxY]
    else:
        xlabel = oimPlotParamLabelShort[idxX]
        ylabel = oimPlotParamLabelShort[idxY]

    # TODO: Implement units with astropy
    if xunit:
        xunit0 = u.Unit(xunit)
        xunitmultiplier = oimPlotParamUnit0[idxX].to(xunit0)
    else:
        xunit0 = oimPlotParamUnit0[idxX]
        xunitmultiplier = 1

    if yunit:
        yunit0 =u.Unit(yunit)
        yunitmultiplier = oimPlotParamUnit0[idxY].to(yunit0)
    else:
        yunit0 = oimPlotParamUnit0[idxY]
        yunitmultiplier = 1

    if xunit0 != u.one:
        xlabel += f" ({xunit0:latex_inline})"
    if yunit0 != u.one:
        ylabel += f" ({yunit0:latex_inline})"

    xIsUVcoord = oimPlotParamIsUVcoord[idxX]
    yIsUVcoord = oimPlotParamIsUVcoord[idxY]

    if not xIsUVcoord and xname != "EFF_WAVE":
        raise TypeError("X should be LENGTH, SPAFREQ, PA or EFF_WAVE")

    if yIsUVcoord:
        raise TypeError("Y shouldn't be UCOORD,VCOORD, SPAFREQ, or PA")

    if colorTab == None:
        colorTab = oimPlotParamColorCycle

    try:
        if "by" not in color:
            colorTab = [color]
    except:
        pass

    ncol = len(colorTab)
    
    yname0 = yname
    if yname == "VISDATA":
        yname0 = "VIS2DATA"
    colorIdx, ColorNames = getColorIndices(oifitsList, color, yarr, yname0)

    if "label" in kwargs:
        label = kwargs.pop("label")
    else:
        label = ""
    if not (axe):
        axe = plt.axes()

    for ifile, data in enumerate(oifitsList):
        # yname can be anything but  UCOORD, VCOORD, LENGTH, SPAFREQ, PA or EFF_WAVE
        extnames = np.array([di.name for di in data])
        idx_yext = np.where(extnames == yarr)[0]
        yinsname = np.array([data[i].header['INSNAME'] for i in idx_yext])
        
        #yname can be VISDATA =sqrt(VIS2DATA)
        if yname != "VISDATA":
            ydata = [data[i].data[yname] for i in idx_yext]
            ydataerr = [data[i].data[yerrname] for i in idx_yext]
            yflag = [data[i].data["FLAG"] for i in idx_yext]
        else:
            ydata = [np.sqrt(data[i].data["VIS2DATA"]) for i in idx_yext]
            ydataerr = [data[i].data["VIS2ERR"]/(2*np.sqrt(data[i].data["VIS2DATA"])) for i in idx_yext]
            yflag = [data[i].data["FLAG"] for i in idx_yext]
            
        # xname can be LENGTH, SPAFREQ, PA or EFF_WAVE
        if xname == "EFF_WAVE":
            idx_xext = np.where(extnames == xarr)[0]
            xinsname = np.array([data[i].header['INSNAME'] for i in idx_xext])
            xdata = []
            for idata in range(len(idx_yext)):
                iwlarr = idx_xext[np.where(xinsname == yinsname[idata])[0][0]]
                wl = data[iwlarr].data["EFF_WAVE"]
                nB = ydata[idata].shape[0]
                xdata.append(np.tile(wl[None, :], (nB, 1)))

        elif xname == "SPAFREQ":
            xdata = getSpaFreq(data, arr=yarr,  squeeze=False)

        elif xname == "LENGTH":
            B = getBaselineLengthAndPA(data, arr=yarr, squeeze=False,T3Max=True)[0]
            xdata = []
            for i,idata in enumerate(idx_yext):
                xdata.append(np.transpose(
                    np.tile(B[i], (np.shape(data[idata].data[yname])[1], 1))))

        elif xname == "PA":
            PA = getBaselineLengthAndPA(data, arr=yarr, squeeze=False,T3Max=True)[1]
            xdata = []
            for i,idata in enumerate(idx_yext):
                xdata.append(np.transpose(
                    np.tile(PA[i], (np.shape(data[idata].data[yname])[1], 1))))

        if cname is not None:
            idxC = np.where(oimPlotParamName == cname)[0][0]
            cIsUVcoord = oimPlotParamIsUVcoord[idxC]

            if not shortLabel:
                clabel = oimPlotParamLabel[idxC]
                cunitmultiplier = u.Unit(oimPlotParamUnit0[idxC]).to(u.Unit(cunit))
            else:
                clabel = oimPlotParamLabelShort[idxC]

            if cunit:
                cunit0 = u.Unit(cunit)
                cunitmultiplier = oimPlotParamUnit0[idxC].to(cunit0)
            else:
                cunit0 = oimPlotParamUnit0[idxC]
                cunitmultiplier = 1

            if cunit0 != u.one:
                clabel += " ("+f"{cunit0:latex}"+")"

            if cname == "MJD":
                carr = yarr
                cdata = [data[i].data[cname] for i in idx_yext]

            elif cname == "EFF_WAVE":
                carr = "OI_WAVELENGTH"
                idx_cext = np.where(extnames == carr)[0]
                cinsname = np.array([data[i].header['INSNAME']
                                    for i in idx_cext])
                cdata = []
                for idata in range(len(idx_yext)):
                    iwlarr = idx_cext[np.where(
                        cinsname == yinsname[idata])[0][0]]
                    wl = data[iwlarr].data["EFF_WAVE"]
                    nB = ydata[idata].shape[0]
                    cdata.append(np.tile(wl[None, :], (nB, 1)))

            elif not cIsUVcoord:
                try:
                    idxC = np.where(oimPlotParamName == cname)[0][0]
                    carr = oimPlotParamArr[idxC]
                    cdata = [data[i].data[cname] for i in idx_yext]
                except Exception:
                    idxC = np.where(oimPlotParamError == cname)[0][0]
                    carr = oimPlotParamArr[idxC]
                    cdata = [data[i].data[cname] for i in idx_yext]

            elif cname == "SPAFREQ":
                cdata = getSpaFreq(data, arr=yarr, squeeze=False)

            elif cname == "LENGTH":
                B = getBaselineLengthAndPA(data, arr=yarr, squeeze=False,T3Max=True)[0]
                cdata = []
                for i,idata in enumerate(idx_yext):
                    cdata.append(np.transpose(
                        np.tile(B[i], (np.shape(data[idata].data[yname])[1], 1))))

            elif cname == "PA":
                PA = getBaselineLengthAndPA(data, arr=yarr, squeeze=False,T3Max=True)[1]
                cdata = []
                for i,idata in enumerate(idx_yext):
                    cdata.append(np.transpose(
                        np.tile(PA[i], (np.shape(data[idata].data[yname])[1], 1))))

        for idata in range(len(xdata)):
            xdata[idata]=xdata[idata]*xunitmultiplier

        for idata in range(len(xdata)):
            ydata[idata]=ydata[idata]*yunitmultiplier

        if cunit:
            for idata in range(len(xdata)):
                cdata[idata]=cdata[idata]*cunitmultiplier

        # NOTE: Looping through oifits files
        for idata in range(len(ydata)):
            shapex = np.shape(xdata[idata])
            shapey = np.shape(ydata[idata])

            # NOTE: Dealing with the xy data dimensions
            if np.size(shapex) == np.size(shapey):
                if np.size(shapex) == 1:  # if 1 baseline only just change array dim
                    nlam = np.size(xdata)
                    xdata = np.reshape((1, nlam))
                    ydata = np.reshape((1, nlam))
            elif np.size(shapex) == 1:  # if x=1D and y=2D
                if shapex[0] == shapey[0]:
                    xdata[idata] = np.outer(xdata[idata], np.ones(shapey[1]))
                else:
                    xdata[idata] = np.outer(np.ones(shapey[0]), xdata[idata])

            elif np.size(shapey) == 1:  # if x=2D and y=1D
                if shapex[0] == shapey[0]:
                    ydata[idata] = np.outer(ydata[idata], np.ones(shapex[1]))
                    ydataerr[idata] = np.outer(
                        ydataerr[idata], np.ones(shapex[1]))
                else:
                    ydata[idata] = np.outer(np.ones(shapex[0]), ydata[idata])
                    ydataerr[idata] = np.outer(
                        ydataerr[idata], np.ones(shapex[1]))

            shapex = np.shape(xdata[idata])
            shapey = np.shape(ydata[idata])

            if cname is not None:
                shapec = np.shape(cdata[idata])
                if (np.size(shapec) == 1):
                    if shapec[0] == shapex[0]:
                        cdata[idata] = np.outer(
                            cdata[idata], np.ones(shapex[1]))
                    else:
                        cdata[idata] = np.outer(
                            np.ones(shapex[0]), cdata[idata])
                    shapec = np.shape(cdata[idata])

            # NOTE: Separate multiples baselines
            for iB in range(shapex[0]):
                if not showFlagged:
                    flags = np.reshape(yflag[idata], shapey)[iB, :]
                    nflags, flag0, ilam0 = len(flags), True, 0
                    for ilam, flagi in enumerate(flags):
                        doPlot = False
                        flagi = True if np.isnan(ydata[idata][iB, ilam]) else flagi
                        if flag0 != flagi:
                            if not flagi:
                                ilam0 = ilam
                            else:
                                doPlot = True

                            flag0 = flagi
                        if ilam == (nflags-1) and not flagi:
                            doPlot = True

                        if doPlot:
                            labeli = label + " " + \
                                ColorNames[colorIdx[ifile][idata][iB]]

                            if cname is None:
                                if xdata[idata][iB, ilam0:ilam+1].size == 1:
                                    axe.scatter(xdata[idata][iB, ilam0:ilam+1],
                                                ydata[idata][iB,ilam0:ilam+1],
                                                color=colorTab[colorIdx[ifile]
                                                               [idata][iB] % ncol],
                                                label=labeli, **kwargs)
                                else:

                                    axe.plot(xdata[idata][iB, ilam0:ilam+1],
                                             ydata[idata][iB, ilam0:ilam+1],
                                             color=colorTab[colorIdx[ifile]
                                                            [idata][iB] % ncol],
                                             label=labeli, **kwargs)
                                    if errorbar:
                                        if "color" not in kwargs_error:
                                            kwargs_errori = kwargs_error.copy()
                                            kwargs_errori["color"] = \
                                                colorTab[colorIdx[ifile][idata][iB] % ncol]

                                        _errorplot(axe, xdata[idata][iB, ilam0:ilam+1],
                                                   ydata[idata][iB, ilam0:ilam+1],
                                                   ydataerr[idata][iB, ilam0:ilam+1],
                                                   **kwargs_errori)

                            else:
                                # NOTE: dummy plot with alpha=0 as _colorPLot works
                                # with collections thus not updating the xlim and ylim
                                # automatically
                                axe.plot(xdata[idata][iB, ilam0:ilam+1],
                                         ydata[idata][iB, ilam0:ilam+1],
                                         color="k", alpha=0)

                                res = _colorPlot(axe, xdata[idata][iB, ilam0:ilam+1],
                                                 ydata[idata][iB, ilam0:ilam+1],
                                                 cdata[idata][iB, ilam0:ilam+1],
                                                 label=labeli, setlim=False, **kwargs)

                                if errorbar:
                                    _errorplot(axe, xdata[idata][iB, ilam0:ilam+1],
                                               ydata[idata][iB, ilam0:ilam+1],
                                               ydataerr[idata][iB, ilam0:ilam+1],
                                               **kwargs_error)

                else:
                    labeli = label + " " + ColorNames[colorIdx[ifile][idata][iB]]
                    axe.plot(xdata[idata][iB, :],ydata[idata][iB, :],
                             color=colorTab[colorIdx[ifile][idata][iB] % ncol],
                             label=labeli, **kwargs)
                    if errorbar:
                        _errorplot(axe, xdata[idata][iB, :], ydata[idata][iB, :],
                                   ydataerr[idata][iB, :], color=colorTab[
                                       colorIdx[ifile][idata][iB] % ncol],
                                   **kwargs_error)

    if yscale is not None:
        axe.set_yscale(yscale)

    if xscale is not None:
        axe.set_xscale(xscale)

    if xlim is not None:
        axe.set_xlim(xlim)

    if ylim is not None:
        axe.set_ylim(ylim)

    axe.set_xlabel(xlabel)
    axe.set_ylabel(ylabel)
    
    if cname and showColorbar:
        plt.colorbar(res, ax=axe, label=clabel) 
    return res


class _HandlerColorLineCollection(HandlerLineCollection):
    def create_artists(self, legend, artist, xdescent, ydescent,
                       width, height, fontsize, trans):
        """Create the artists for the legend."""
        x = np.linspace(0, width, self.get_numpoints(legend)+1)
        y = np.zeros(self.get_numpoints(legend)+1)+height/2.-ydescent
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=artist.cmap,
                            transform=trans)
        lc.set_array(x)
        lc.set_linewidth(artist.get_linewidth())
        return [lc]


class oimWlTemplatePlots(Figure):
    """A class to plot the wavelength of the data."""

    def __init__(self,*args, **kwargs):
        """The class initialisation."""
        super().__init__(*args, **kwargs)
        self.data = []
        self.xunit = u.m

    def autoShape(
            self, oifitsList: fits.HDUList,
            shape: Optional[List[List[str]]]=[["VIS2DATA", None], ["VISAMP", None], ["VISPHI", None],
                                              ["T3AMP", None], ["T3PHI", None], ["FLUXDATA", None]]
            ) -> None:
        """Automatically set the shape of the plot.

        Parameters
        ----------
        oifitsList : fits.HDUList
            A list of oifits files.
        shape : list of list of str, optional
            The shape of the plot.
        """
        if isinstance(oifitsList, oimData):
            oifitsList = oifitsList.data

        if type(oifitsList) != type([]):
            oifitsList = [oifitsList]

        self.nfiles = len(oifitsList)
        gs = gridspec.GridSpec(self.nfiles, 1, figure=self)
        self.datatype = []

        for filek in oifitsList:
            datatypek = []
            for shapei in shape:
                datatypei = []
                for shapeij in shapei:
                    if shapeij is not None:
                        idx = np.where(oimPlotParamName == shapeij)[0][0]
                        arrName = oimPlotParamArr[idx]
                        for hdui in filek:
                            if hdui.name == arrName:
                                extver = hdui.header["EXTVER"]
                                if np.any(hdui.data[shapeij]!=0):
                                    s = hdui.data[shapeij].shape[0]
                                    for iB in range(s):
                                        datatypei.append([extver,iB,shapeij])
 
                if len(datatypei) != 0:
                    datatypek.append(datatypei)

            self.datatype.append(datatypek)

        self.plotList = []
        self.columns_list, self.lines_list = [], []
        for k, datatypek in enumerate(self.datatype):
            ncs = np.array([len(c) for c in datatypek])
            nc, nl = np.max(ncs), np.size(ncs)

            self.columns_list.append(nc)
            self.lines_list.append(nl)

            plotListk = []
            gss = gridspec.GridSpecFromSubplotSpec(nl, nc, subplot_spec=gs[k])
            for il in range(nl):
                plotListkl = []
                datatypeil = self.datatype[k][il]
                nplots_col = len(datatypeil)
                showXlabel = (k == self.nfiles - 1) and (il == nl - 1 )

                datatype_previous = None
                for ic in range(nc):
                    self.add_subplot(gss[il,ic])
                    if ic < nplots_col:
                        extver, iB, dataType = datatypeil[ic]
                        if datatype_previous is None:
                            datatype_0 = dataType
                            showYlabel = 1
                        elif ic == nc-1 and datatype_0 != dataType:
                            showYlabel = 2
                        else:
                            showYlabel = 0

                        obj = [self.axes[-1], k, ic, il, extver, iB, dataType, showXlabel, showYlabel]
                    else:
                        obj = [self.axes[-1], k, ic, il, None, None, None, False, 0]
                        try:
                            plotListk[il-1][ic][7] = True
                        except:
                            pass

                    datatype_previous = dataType

                    plotListkl.append(obj)
                plotListk.append(plotListkl)
            self.plotList.append(plotListk)

        for ax in self.axes[1:]:
            ax.sharex(self.axes[0])

        yaxe_shares, yaxe_datatypes = [], []
        for k in range(len(self.plotList)):
            for l in range(len(self.plotList[k])):
                for ax, _, _, _, extver, iB, dataType, showXlabel, showYlabel in self.plotList[k][l]:
                    if dataType in yaxe_datatypes:
                        ax2share = yaxe_shares[yaxe_datatypes.index(dataType)]
                        ax.sharey(ax2share)
                    else:
                        yaxe_shares.append(ax)
                        yaxe_datatypes.append(dataType)

    def plot(self, oifitsList: fits.HDUList, add: Optional[bool] = True,
             plotFunction: Optional[Axes] = plt.Axes.plot,
             plotFunctionkwarg: Optional[Dict] = {}) -> None:
        """Plot the data from the oifits files."""
        if isinstance(oifitsList, oimData):
            oifitsList = oifitsList.data

        if type(oifitsList) != type([]):
            oifitsList = [oifitsList]

        if add:
            self.data.append(oifitsList)

        for k in range(len(self.plotList)):
            for l in range(len(self.plotList[k])):
                for ax, ifile, _, _, extver, iB, dataType, *_ in self.plotList[k][l]:
                    if dataType:
                        oimSimplePlotWavelength(
                                oifitsList, ifile, dataType, iB,
                                extver=extver, axe=ax, plotFunction=plotFunction,
                                xunit=self.xunit, plotFunctionkwarg=plotFunctionkwarg)
                    else:
                        ax.axis('off')

        self.set_xlabels("auto")
        self.set_ylabels("auto")

    def set_xunit(self, xunit: u.Quantity) -> None:
        """Set the unit of the x-axis."""
        self.xunit = u.Unit(xunit)

    def set_wllim(self,*arg) -> None:
        """Set the wavelength limits of the plot."""
        for k in range(len(self.plotList)):
            for l in range(len(self.plotList[k])):
                for ax, *_ in self.plotList[k][l]:
                    ax.set_xlim(*arg)

    def set_xlim(self,*arg) -> None:
        """Set the x-axis limits of the plot."""
        self.set_wllim(*arg)

    def set_ylim(self, dataTypes: Union[str, List[str]], *arg) -> None:
        """Set the y-axis limits of the plot."""
        if not isinstance(dataTypes,list):
            dataTypes = [dataTypes]

        for k in range(len(self.plotList)):
            for l in range(len(self.plotList[k])):
                for ax, _, _, _, _, _,dataType, *_ in self.plotList[k][l]:
                    if dataType in dataTypes:
                        ax.set_ylim(*arg)

    def set_xlabels(self, mode: str, *arg, **kwargs) -> None:
        """Set the x-axis label of the plot."""
        if mode == "auto":
            label = f"$\\lambda$ ({self.xunit:latex})"
            for k in range(len(self.plotList)):
                for l in range(len(self.plotList[k])):
                    for ax, _, _, _, _, _, dataType, showXlabel, _ in self.plotList[k][l]:
                        if dataType:
                            if showXlabel:
                                ax.set_xlabel(label, *arg,**kwargs)
                            else:
                                ax.get_xaxis().set_visible(False)

    def set_ylabels(self, mode: str, *arg) -> None:
        """Set the y-axis label of the plot."""
        if mode == "auto":
            for k in range(len(self.plotList)):
                for l in range(len(self.plotList[k])):
                    for ax, _, _, _, _, _, dataType, _, showYlabel in self.plotList[k][l]:
                        if dataType:
                            if showYlabel == 2:
                                ax.yaxis.set_label_position("right")
                                ax.yaxis.tick_right()
                            if showYlabel != 0:
                                idx = np.where(oimPlotParamName == dataType)[0][0]
                                label = oimPlotParamLabelShort[idx]
                                ax.set_ylabel(label,*arg)
                            else:
                                ax.get_yaxis().set_visible(False)


    def set_legends(self, *arg, **kwargs) -> None:
        """Set the legend of the plot."""
        if len(arg) == 2:
            text0, datatypes = arg
            x = y = 0.02
        elif len(arg) == 4:
             x, y, text0, datatypes = arg
        else :
            raise TypeError("Wrong number of argument : Should be 2 "\
                            "(text0, datatypes) or 4 (x, y, text0, datatypes)")

        if isinstance(datatypes,str):
            datatypes = [datatypes]

        for k in range(len(self.plotList)):
            for l in range(len(self.plotList[k])):
                for ax, _, _, _, _, iB, dataType, *_ in self.plotList[k][l]:
                    if  dataType in datatypes :
                        try :
                            LENGTH,PA = getBaselineLengthAndPA(
                                self.data[0][k], arr=getDataArrname(dataType))
                        except:
                            pass
                        BASELINE = getBaselineName(
                            self.data[0][k], hduname=getDataArrname(dataType))

                        text = text0
                        if ("$LENGTH$" in text0) :
                            text = text.replace("$LENGTH$",f"{LENGTH[iB]:.1f}")
                        if ("$PA$" in text0):
                            text = text.replace("$PA$",f"{PA[iB]:.1f}")
                        if ("$BASELINE$") in text0:
                            text =  text.replace("$BASELINE$",f"{BASELINE[iB]}")

                        ax.text(x, y, text, transform=ax.transAxes, **kwargs)


def oimSimplePlotWavelength(oifitsList: fits.HDUList,
                            ifile: int, dataType: str, iB: int,
                            extver = None, axe: Optional[Axes] = None,
                            plotFunction: Optional[Axes] = plt.Axes.plot,
                            xunit: Optional[u.Quantity] = u.m,
                            plotFunctionkwarg: Optional[Dict] = {}) -> None:
    """Set the x-axis to the wavelength of the data.

    Parameters
    ----------
    oifitsList : fits.HDUList
        A list of oifits files.
    ifile : int
    The index of the file.
    dataType : str
        The data type to plot.
    iB : int
    The baseline index.
    extver : int, optional
    The extension version. The default is None.
    """
    xmult = u.m.to(xunit)
    idx = np.where(oimPlotParamName == dataType)[0][0]
    arrName = oimPlotParamArr[idx]
    errorName = oimPlotParamError[idx]

    oifitsi=oifitsList[ifile]
    wl = getWlFromOifits(oifitsi, arrName, extver=extver)

    data = oifitsi[arrName, extver].data[dataType][iB, :]
    err = oifitsi[arrName, extver].data[errorName][iB, :]

    if axe is None:
        axe = plt.gca()

    if plotFunction.__name__ == "plot":
        plotFunction(axe, wl*xmult, data, **plotFunctionkwarg)
    else:
        plotFunction(axe, wl*xmult, data, err, **plotFunctionkwarg)


class oimAxes(plt.Axes):
    """Class derived from plt.Axes that allows easy plotting of oifits data"""
    name = "oimAxes"
    xtype, ytype = None, None
    colorbar = None

    def uvplot(self, oifits: fits.HDUList, **kwargs):
        """Plot the uv coverage of the oifits data."""
        uvPlot(oifits, axe=self, **kwargs)

    def oiplot(self, oifitsList: fits.HDUList, xname: str, yname: str, **kwargs):
        """Plot the data from the oifits files."""
        res = oimPlot(oifitsList, xname, yname, axe=self, **kwargs)
        self.xtype, self.ytype = xname, yname
        return res

    def autolim(self):
        """Set the limits of the axes to the data."""
        if self.ytype in ["VIS2DATA", "VISAMP", "T3AMP"]:
            if self.get_yscale() == 'linear':
                self.set_ylim(0, 1)
            else:
                lines = self.get_lines()
                mini = np.inf
                for li in lines:
                    mi = np.min(li.get_ydata())
                    if mi < mini and mi > 0:
                        mini = mi

                self.set_ylim(mini, 1)

        elif self.ytype in ["VISPHI", "T3PHI"]:
            self.set_ylim(-180, 180)

    def legend(self, **kwargs):
        """Add a legend to the plot."""
        handles, labels = self.get_legend_handles_labels()

        hmap = {}
        for hi in handles:
            if isinstance(hi, LineCollection):
                hmap[hi] = _HandlerColorLineCollection(numpoints=100)

        # NOTE: Use to remove duplicate legend. Is Set Maybe better for that?
        legend_hangle = dict(zip(labels, handles))
        super().legend(legend_hangle.values(), legend_hangle.keys(),
                       handler_map=hmap, **kwargs)

    def set_xscale(self, value, **kwargs):
        """Set the xscale of the plot."""
        super().set_xscale(value, **kwargs)
        self.autoscale_view()

    def set_yscale(self, value, **kwargs):
        """Set the yscale of the plot."""
        super().set_yscale(value, **kwargs)
        self.autoscale_view()
