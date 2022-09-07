# oimodeler

A modular modelling tool for optical interferometry


>:warning: In early development!  
>The oimodel is partly working with gray and chromatic models defined in Fourier plane (as shown in the image below)  
>Partial implementation of oidata class: no filtering, no optimization of data)  
>Partial implementation of oiSimulator: simulated data and chi2 computation for all kind of data but no model-fitting  
 

Modules:
* **oimModel** : Create models with various components as bricks 
* **oimData** :  Handle interferometric, spectroscopic and photometric data
* **oimFitter** : Define and perform model-fitting   
* **oimSiumlulator** : Main class holding evertyhing together and producing final results :plots, tables...

Various example scripts are available in the examples directory.  
Here is the resulting plot from the createModelChromatic.py script.
![boo](./images/createModelChromatic.png)
