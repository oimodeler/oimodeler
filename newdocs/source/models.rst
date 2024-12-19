:tocdepth: 2

..  _models:

Building and using models
=========================


In the **oimodeler** framework, building models is based on three classes:

- :func:`oimModel <oimodeler.oimModel.oimModel>`: the model class that will be given to a :func:`oimSimulator <oimodeler.oimSimulator.oimSimulator>` or :func:`oimFitter <oimodeler.oimFitter.oimFitter>` 

- :func:`oimComponent <oimodeler.oimComponent.oimComponent>`: the semi-abstract class for each component of the model. :func:`oimModel <oimodeler.oimModel.oimModel>` is somehow a list of object derived from the <oimodeler.oimComponent.oimComponent>` class.

- :func:`oimParam <oimodeler.oimParam.oimParam>`: the parameter class. Each component has multiple parameters (for instance has the component position ``x`` and ``y``)


Models
------

Types of components
-------------------

List of basic Fourier components
--------------------------------

Image-plan components
---------------------

Importing fits images
---------------------

Using Radial-Profile
--------------------

Linking parameters
------------------

chromatic & time-dependent Interpolator
--------------------

