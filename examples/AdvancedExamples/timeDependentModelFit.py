import oimodeler as oim
import matplotlib.pyplot as plt
import os
import numpy as np

path = os.path.dirname(oim.__file__)
pathData=os.path.join(path,os.pardir,"examples","testData","ASPRO_SPICA_GROWING_UD")

files=[os.path.abspath(os.path.join(pathData,fi)) for fi in os.listdir(pathData) if ".fits" in fi]
data=oim.oimData(files)

#%%
mjd_ref=[60076,60080]
ud=oim.oimUD(d=oim.oimInterp("time",mjd=mjd_ref,values=[0,0],extrapolate=True))
model=oim.oimModel(ud)
ud.params['d'].params[0].set(min=0,max=10)
ud.params['d'].params[1].set(min=0,max=10)
ud.params['f'].free=False

#%%
fit=oim.oimFitterEmcee(files,model,nwalkers=20)
fit.prepare(init="random")

f1=oim.oimRemoveArrayFilter(target="all",arr=["OI_VIS","OI_T3"])         
fit.data.setFilter(oim.oimDataFilter([f1]))

#%%
fit.run(nsteps=300,progress=True)

#%%
figWalkers,axeWalkers=fit.walkersPlot(cmap="plasma_r")
figCorner,axeCorner=fit.cornerPlot(discard=100)

#%%
median,err_l,err_u,err=fit.getResults(mode='median',discard=100)

#%%
fig, ax = plt.subplots(1,1, subplot_kw=dict(projection='oimAxes'),figsize=(8,8))
timecol=ax.oiplot(fit.simulator.data,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas",label="Data",cname="MJD",lw=4)
ax.oiplot(fit.simulator.simulatedData,"SPAFREQ","VIS2DATA" ,xunit="cycles/mas", color="k",ls=":",label="Model")
fig.colorbar(timecol, ax=ax,label="MJD (days)")
ax.legend()
#ax.set_yscale('log')
#ax.set_ylim(1e-5,1)

#%%
fig, ax = plt.subplots(1,1,figsize=(8,4))

diam_ref=ud.params['d'](0,mjd_ref)

mjd_Obs=np.unique(fit.data.vect_mjd)
diam_Obs=ud.params['d'](0,mjd_Obs)

mjd=np.linspace(60076,60080,num=100)
diam=ud.params['d'](0,mjd)
ax.plot(mjd,diam,color="k",lw=2,label="interpolated parameter")
ax.scatter(mjd_ref,diam_ref,marker="o",s=100,color="r",zorder=10,label="Ref values for model fitting")
ax.scatter(mjd_Obs,diam_Obs,marker="X",s=100,c=mjd_Obs,zorder=10,label="Observations")
ax.set_xlabel("Time (MJD)")
ax.set_ylabel("Diameter (mas)")
ax.legend()


