# oimodeler

A modular modelling tool for optical interferometry


>:warning: In early development!  
>The oimModel is partly working with gray and chromatic models defined in Fourier plane (as shown in the image below)  
>Partial implementation of oimData class: no filtering, no optimization of data)  
>Partial implementation of oimSiumlulator: simulated data and chi2 computation for all kind of data but no model-fitting  
 

Modules:
* **oimModel** : Create models with various components as bricks 
* **oimData** :  Handle interferometric, spectroscopic and photometric data
* **oimFilter** : Filtering and modifying data (wavlengths range cut, smoothing, removing flags...)  
* **oimSiumlulator** : Main class holding evertyhing together and producing final results :plots, tables...
* **oimFitter** : Define and perform model-fitting   
Various example scripts are available in the examples directory.  
Here is the resulting plot from the createModelChromatic.py script.
![boo](./images/createModelChromatic.png)
