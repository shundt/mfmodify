# mfmodify

## Description
mfmodify includes functions for reading and modifying existing MODFLOW 6 models, with the aim of performing analysis while making as much use of existing model files as possible. Current modules include:
- scenario_utils: 
    To quickly build scenarios for existing MODFLOW models with transient historic periods. Tools are available for extracting
     stress period data for specified years of the historic period and reassembling the years in any order. The resulting 
     scenario can be any length, can contain the same year multiple times, and can create synthetic years from weighted 
     averages of existing years.  The module also contains (simple) tools for weighted sampling from the 
     population of historic years, so as to make wet, dry, or average years (or whatever percentile range) more or less likely 
     than occured historically. These tools are envisioned to be used by a modeller who wants to develop baseline scenarios 
     and who wants to quickly approximate the uncertainty of unknown future (most likely weather-related) stresses.
- regrid_utils: 
    For converting a structured grid model to an unstructured grid model with local refinement. 

The functions of mfmodify make heavy use of the flopy library, but is not meant to extend the functionality of that library. To
the extent possible, mf_modify will avoid creating new classes. I want to claim that I made this choice for some philisophical 
reason, but it's probably just that I'm not a good enough programmer. Oh well, that's one less obscure data type for you to 
learn.


# mfmodify

## Description
`mfmodify` includes functions for reading and modifying existing MODFLOW 6 models, aimed at analysis while leveraging existing model files. 

### Modules
- **scenario_utils**: 
    - Quickly build scenarios for existing MODFLOW models with transient historic periods. 
    - Extract and reassemble stress period data for specified years.
    - Create synthetic years using weighted averages.
    - Tools for weighted sampling to adjust the likelihood of wet, dry, or average years.
  
- **regrid_utils**: 
    - Convert a structured grid model to an unstructured grid model with local refinement.

'mfmodify' makes use of the flopy library but does not extend it's functionality. The design will avoid creating new classes. 
I could pretend I have some strongly held opinion to justify this, but the truth is simple that I'm not a good enough
programmer.

## Examples
See the notebooks in the examples directory to see the intended use of mfmodify.