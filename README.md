# oimodeler

A modular modelling tool for optical interferometry


> :warning: I early development!  
>  Only the oimodel is partly working with  gray and chromatic models (as shown in the image below)  

Modules:
####oimModel : Creation of Models with various components as bricks : from simple gray analytically defined ones (UD, Gaussian...) to complex ones with chromaticity or based on images computed by Radiative trasnfer codes
####oimData : Load oifits data and spectroscopic or photometric data and simulate data from oimodel.
####oimFitter : Fitting free parameters of an oimModel to data in oimData
####oimSiumlulator : Main class linking the oimModel, oimData and oimFitter together and produce final results (plots, tables...) from fit.


![boo](./images/createModelChromatic.png)
