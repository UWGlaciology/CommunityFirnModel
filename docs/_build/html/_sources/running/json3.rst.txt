**************************************
The .json-formatted configuration file
**************************************

The CFM uses a .json-formatted file to configure individual model runs. JSON (JavaScript Object Notation) is a data-interchange file format. It consists of a number of names, each associated with a value. Values can be strings, Booleans, integers, floats, or arrays. Comments are not allowed, but can be added by considering the comment as a name/value pair. For the CFM, it provides a file format that is both easy to read and easy to alter in order to specify parameters for a particular model run. The configuration file is passed to the CFM, and the name/value pairs are read by the model and incorporated into the model run. The file format is editable in any text editor, and the name/value pairs are given by name: value, and different name/value pairs are separated by commas.

The specific names that are in the configuration .json file for the CFM are as follows. If any of the name/value pairs are missing, the model will generally return a message that that that name/value pair is missing and will use a default instead. For some name/value pairs the model run will fail. Note that in the .json file true/false are lowercase, but in the .py files they are True/False (first letter capitalized). The model automatically converts this. :math:`\rho_{s}` is the surface.

.. jsonschema:: 

    {
      "title": "InputFileFolder",
      "description": "Directory where the input csv files are located (usually a subdirectory of the\n\ndirectory that contains main.py, but user can specify an absolute paths as well.)\n\nUse '' if the input files are in the same directory as main.py.",
      "type": "string",
      "default": null,
      "examples": [inputdata]
    }

.. jsonschema:: 

    {
      "title": "InputFileNameXXXX",
      "description": "The names of the input files for temperature, accumulation/smb, water isotopes,\n\nsurface density, and melt. See 'Inputs for the CFM' section for more details.",
      "type": "string",
      "default": null,
      "examples": ["example_XXXX.csv"]
    }

.. jsonschema:: 

    {
      "title": "resultsFolder",
      "description": "Folder in which results will be stored.",
      "type": "string",
      "default": null,
      "examples": ["example_results"]
    }

.. jsonschema:: 

    {
      "title": "initfirnFile",
      "description": "File containing initial conditions if you are using firn measurements/data (e.g. temperature,\n\ndensity) to begin the model run. See 'Inputs for the CFM' section for more details.",
      "type": "string",
      "default": null,
      "examples": ["example_firndata.csv"]
    }

.. jsonschema:: 

    {
      "title": "initprofile",
      "description": "Whether or not the CFM should use the initfirnFile to generate an initial condition.",
      "type": "boolean",
      "default": false,
    }

.. jsonschema:: 

    {
      "title": "physRho",
      "description": "The firn-densification physics to use for the model run.",
      "type": "string",
      "default": null,
      "options": ["HLdynamic","HLSigfus","Li2004","Li2011","Helsen2008","Arthern2010S","Arthern2010T","Li2015","Goujon2003","Barnola1991","Morris2013","KuipersMunneke2015","Crocus","Ligtenberg2011"]
    }  

.. jsonschema:: 

    {
      "title": "MELT",
      "description": "Whether or not to include meltwater percolation physics in the model run.",
      "type": "boolean",
      "default": false,
    }

.. jsonschema:: 

    {
      "title": "ReehCorrectedT",
      "description": "If melt is enabled, (NEED TO FILL IN WHAT THIS DOES)",
      "type": "boolean",
      "default": false,
    }

.. jsonschema:: 

    {
      "title": "FirnAir",
      "description": "Whether or not to run the firn air module with the model run.",
      "type": "boolean",
      "default": false,
    }

.. jsonschema:: 

    {
      "title": "AirConfigName",
      "description": "Name of the .json configuration files that contains the parameters for the firn air module.",
      "type": "string",
      "default": "AirConfig.json",
    }

.. jsonschema:: 

    {
      "title": "TWriteInt",
      "description": "How often to write the results to file, relative to the time-step size, i.e. 1 will\n\nwrite at every time step, 10 will write at every 10th time step, etc.",
      "type": "int",
      "default": 1,
    }

.. jsonschema:: 

    {
      "title": "TWriteStart",
      "description": "The time at which to start saving model results. The time is model time, so must\n\ncorrespond to the time in the input forcing files.",
      "type": "float",
      "default": Null,
      "examples": ["If your model run is from 1958 to 2018, but you only want outputs from 2000\n\nonwards, 'TWriteStart' should be 2000."]
    }
    
.. jsonschema:: 

    {
      "title": "int_type",
      "description": "How to interpolate from the input file times to the model time. Use linear e.g. if\n\nyou have sparse ice core climate data. Use nearest e.g. if you have monthly climate data and want to\n\ntake monthly time steps (the model time may not be exactly the same as the input time).",
      "type": "string",
      "options": ["nearest","linear"]
    }

.. jsonschema:: 

    {
      "title": "SeasonalTcycle",
      "description": "Whether or not to add a seasonal temperature cycle (on top of the forcing data). Use this\n\nonly if you are using sub-annual time steps and your forcing data does not have a\n\nseasonal cycle already. Usually this would be if your forcing data is annual (or coarser resolution).",
      "type": "boolean",
      "default": false
    }

.. jsonschema:: 

    {
      "title": "SeasonalThemi",
      "description": "If 'SeasonalTCycle' is True, specify which hemisphere you are modeling to get the summer/winter\n\n timing correct.",
      "type": "string",
      "options": ['north','south']
    }

.. jsonschema:: 

    {
      "title": "coreless",
      "description": "If 'SeasonalTCycle' is True, add the coreless winter. ADD MORE INFO HERE",
      "type": "boolean",
    }

.. jsonschema:: 

    {
      "title": "TAmp",
      "description": "If 'SeasonalTCycle' is True, specify the amplitude of the cycle.",
      "type": "float",
      "default": 10,
      "units": ['K']
    }    

.. jsonschema:: 

    {
      "title": "other",
      "description": "If 'SeasonalTCycle' is True, specify the amplitude of the cycle. :math:`\rho_{s}` nn ",
      "type": "float",
      "default": 10,
      "units": [":math:`\rho_{s}`"]
    }  
























