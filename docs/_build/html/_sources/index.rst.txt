.. CommunityFirnModel documentation master file, created by
   sphinx-quickstart on Wed Nov  6 14:46:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to documentation for the Community Firn Model!
======================================================

*This documentation is still a work in progress. Please inform me if you find glaring omissions!*

The Community Firn Model (CFM) is an open-source, modular firn-evolution model framework. Its most basic function is to predict the depth/density and depth/age profiles of firn (i.e. functioning as a firn-densification model), but it also includes modules to simulate heat transfer, meltwater percolation and refreezing, water isotope diffusion, firn-air diffusion and advection, and grain growth.

The modular nature of the CFM means that the user can easily choose to run the model using a firn-densification equation from any of a suite of published firn-densification equations. The user can also choose which of the optional physics modules (e.g. grain growth) to include with the model run. The model is designed so that different firn-densification equations or modules simulating different physics can be easily integrated into the model. 

The CFM is one dimensional and uses a Lagrangian (material-following) grid. The evolution of the density is calculated explicitly (e.g. :math:`\rho_{new}=\rho_{old}+d\rho/dt*dt)`. Heat and other (e.g. water isotope, firn air, enthalpy) diffusion is solved using a fully-implicit finite-volume method (Patankar, 1980). 

Requirements
------------

Python 3.6+, 
`numpy <http://www.numpy.org>`_, 
`scipy <http://www.scipy.org>`_, 
`h5py <https://www.h5py.org>`_,
`pandas <https://pandas.pydata.org>`_

The CFM should run on Windows, Linux, and OSX platforms. Additionally, plotting CFM results using Python requires the `matplotlib <http://matplotlib.org>`_ package. Installing Python 3 and the necessary packages is most easily done using a packaged Python distribution such as Anaconda. Explicit instructions on a scientific Python installation are beyond the scope of this document, but tutorials can be readily found online.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   running/index.rst
..   files/index.rst

.. toctree::
   :maxdepth: 1
   :caption: Files:

   files/index.rst

.. toctree::
   :maxdepth: 1
   :caption: Extras:
   
   extras/index.rst


Contributing
------------


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
