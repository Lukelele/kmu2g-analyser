# CERN NA62 Analysis Framework
This project is designed to select and study specific Kaon ($K^+$) decay events. It includes C++ modules for event selection based on kinematic and detector variables, and Jupyter notebooks for further analysis, machine learning (ML) model training, and data visualization.

The primary decay channel of interest is $K^+ \to \mu^+ \nu \gamma$ (Kmu2g).

## Project Structure

`/include`: Contains the C++ header files (.hh) defining the classes for different event selections.

`/src`: Contains the C++ source files (.cc) implementing the logic for the event selection classes.

`*.ipynb`: Jupyter notebooks used for data processing, analysis, model training, and visualization.

`*.log`: Log files, likely from data processing or batch jobs.


## C++ Selection Modules

`Kmu2gSelection`: Implements the selection for the $K^+ \to \mu^+ \nu \gamma$ decay channel.

`MyKmu2Selection`: Implements the selection for the $K^+ \to \mu^+ \nu$ decay channel.

`myK2piSelection`: Implements the selection for the $K^+ \to \pi^+ \pi^0$ decay channel.

`MyK3piSelection`: Implements the selection for the $K^+ \to \pi^+ \pi^+ \pi^-$ decay channel.


## Analysis Notebooks

`analysis.ipynb`: Main notebook for data analysis. Likely reads n-tuples, applies further cuts, and produces plots.

`visualisation.ipynb`: Notebook dedicated to creating plots and visualizing data.


Prerequisites

To use this project, you will likely need:

A C++ compiler (e.g., g++)

The ROOT data analysis framework
The NA62 experiment's official software framework (for compiling and running the C++ modules)

Python libraries:

`uproot` (for reading ROOT files in Python)
`pandas`
`numpy`
`matplotlib`
`scikit-learn` (for machine learning)
`xgboost`

