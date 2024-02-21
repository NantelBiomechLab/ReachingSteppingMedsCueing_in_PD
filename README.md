# Reaching and stepping respond differently to medication and cueing in Parkinson's disease


This results of this research are fully reproducible using the source code, Jupyter
notebooks, and raw data contained within this repository. Please see the paper for a full
description of the data and methods.

Raw data included in this repository:

- raw motion capture trials, as `.c3d` files
- participant demographics

## Instructions

To run this analysis on your computer, both Julia and Jupyter (notebook or lab) must be available. A
version of Julia appropriate for your OS can be downloaded from the [Julia language
website](https://julialang.org/downloads/), and Jupyter can be installed from within Julia
(in the REPL) with

```julia
] add IJulia
```

Alternate instructions for installing Jupyter can be found on the [IJulia
github](https://github.com/JuliaLang/IJulia.jl) or the [Jupyter
homepage](https://jupyter.org/install) (advanced).

If you wish to use an external/pre-existing installation of Jupyter, before running `] add
IJulia`, set an environmental variable `"JUPYTER"` to the path of your Jupyter executable.

From within the directory of this repository, start Julia and then start Jupyter from the Julia
REPL, using this command:

```julia
using IJulia
notebook(;dir=pwd())
```

or if using a system Jupyter installation, start Jupyter from your favorite available shell
(e.g. Powershell on Windows, bash on any *nix variant, etc.).

## Description of notebooks:

- `Analysis (cleaned).ipynb`
    - Primary analysis for this study. Contains complete analysis workflow, from raw data
        (e.g. ".c3d" files) to final statistical results. Also includes ancillary
        demographic analyses and post-hoc tests.
- `Kinematic exploratory (cleaned).ipynb`
    - Initial exploration of motion capture data, and verification of event identification
        (reaches and steps). Detailed figures show important characteristics of signals
        during a trial, including "plausible freezing episodes" and supporting signals (e.g.
        the touch target signal for RRT trials, etc)
- `Randomization confirmation.ipynb`
    - Tests to confirm the block randomization of trial order.
- `Strip EMG and export.ipynb`
    - Original protocol included EMG data (see Cantú et al., 2018, 2019) that was not
        analyzed by the present study. This notebook records the modifications made to the
        original data to remove the EMG signals and then exported for publishing.

# Bibliography

- Cantú, H., Côté, J. N., & Nantel, J. (2018). A new method based on quiet stance baseline is more effective in identifying freezing in Parkinson’s disease. PLOS ONE, 13(11), e0207945. https://doi.org/10.1371/journal.pone.0207945

- Cantú, H., Nantel, J., Millán, M., Paquette, C., & Côté, J. N. (2019). Abnormal Muscle Activity and Variability Before, During, and After the Occurrence of Freezing in Parkinson’s Disease. Frontiers in Neurology, 10. https://doi.org/10.3389/fneur.2019.00951


