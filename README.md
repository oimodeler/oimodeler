[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![Lifecycle:
Early Development](https://img.shields.io/badge/lifecycle-EarlyDevelopment-orange.svg)
![version](https://img.shields.io/badge/version-0.0.1-blue)
# oimodeler

A modular modelling tool for optical interferometry


>:warning: In early development!  
>Only the oimodel is partly working with  gray and chromatic models (as shown in the image below)  

Modules:
* **oimModel** : Create models with various components as bricks 
* **oimData** :  Handle interferometric, spectroscopic and photometric data
* **oimFitter** : Define and perform model-fitting   
* **oimSiumlulator** : Main class holding evertyhing together and producing final results :plots, tables...

Various example scripts are available in the examples directory.  
Here is the resulting plot from the createModelChromatic.py script.
![boo](./images/createModelChromatic.png)
