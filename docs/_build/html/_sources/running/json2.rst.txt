**************************************
The .json-formatted configuration file
**************************************

The CFM uses a .json-formatted file to configure individual model runs. JSON (JavaScript Object Notation) is a data-interchange file format. It consists of a number of names, each associated with a value. Values can be strings, Booleans, integers, floats, or arrays. Comments are not allowed, but can be added by considering the comment as a name/value pair. For the CFM, it provides a file format that is both easy to read and easy to alter in order to specify parameters for a particular model run. The configuration file is passed to the CFM, and the name/value pairs are read by the model and incorporated into the model run. The file format is editable in any text editor, and the name/value pairs are given by name: value, and different name/value pairs are separated by commas.

The specific names that are in the configuration .json file for the CFM are as follows. If any of the name/value pairs are missing, the model will generally return a message that that that name/value pair is missing and will use a default instead. For some name/value pairs the model run will fail. Note that in the .json file true/false are lowercase, but in the .py files they are True/False (first letter capitalized). The model automatically converts this. 

.. jsonschema:: 

    {
      "$schema": "http://json-schema.org/draft-07/schema#",
      "title": "The CFM .json file",
      "id": "https://github.com/UWGlaciology/CommunityFirnModel/tree/master/CFM_main/example.json",
      "description": "Description of the fields in the CFM .json configuration file ",
      "type": "object",
      "properties": {
        
        "InputFileFolder": {
            "description": "The folder holding the input files.\n\nCan include a relative or absolute path",
            "type": "string",
            "default": null,
            "examples": [inputdata]},
        
        "InputFileNameTemp": {
            "description": "The name of the temperature forcing file",
            "type": "string",
            "default": null,
            "examples": ['example_tskin.csv']},

        "InputFileNamebdot": {
            "description": "The name of the accumulation or surface mass balance forcing file",
            "type": "string",
            "default": null,
            "examples": ['example_smb.csv']},

        "InputFileNameIso": {
            "description": "The name of the water isotope forcing file",
            "type": "string",
            "default": null,
            "examples": ['example_iso.csv']},

        "InputFileNamerho": {
            "description": "The name of the surface density forcing file",
            "type": "string",
            "default": null,
            "examples": ['example_rhos.csv']},

        "InputFileNamemelt": {
            "description": "The name of the surface melt forcing file",
            "type": "string",
            "default": null,
            "examples": ['example_melt.csv']},

        "resultsFolder": {
            "description": "folder in which results will be stored",
            "type": "string",
            "default": null,
            "examples": ['example_results']},                                            

        "yearSpin": {
            "description": "how many years to spinup for",
            "type": "float",
            "default": 400},

        "other_thing": {
            "description": "testng this",
            "type": boolean}
        }
    }
    
