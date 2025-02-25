# mfmodify

## Description
`mfmodify` includes functions for reading and modifying existing MODFLOW 6 models, aimed at analysis while leveraging existing model files. 

### Modules
- **scenario_utils**: 
    - Quickly build scenarios for existing MODFLOW models with transient historic periods. 
    - Extract and reassemble stress period data for specified years into a new transient period.
    - Create synthetic years using weighted averages of years.
    - Tools for weighted sampling to adjust the likelihood of wet, dry, or average years.
  
- **regrid_utils (in progress)**: 
    - Convert a structured grid model to an unstructured grid model with local refinement.

'mfmodify' makes use of the flopy library but does not extend it's functionality. The design will avoid creating new classes. 
I could pretend I have some strongly held opinion to justify this, but the truth is simple that I'm not a good enough
programmer.

## Examples
See the notebooks in the examples directory to see the intended use of mfmodify.