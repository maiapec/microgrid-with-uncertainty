# Incorporating uncertainty in Model Predictive Control for the optimal operation of a microgrid

This repository contains an implementation of two Model Predictive Control methods for energy management of an islanded microgrid, optimized to minimize unmet demand over the year:
- **Direct MPC**: optimization over a limited time horizon (one week) with error on weather forecasts.
- **Smooth MPC**: adaptive optimization similar to the Direct MPC, which penalizes deviations from previous actions to improve robustness against forecast errors.

Different types of error in the weather forecast are tested to analyze how they impact the optimization. The model uses the Gurobi solver.

**Acknowledgments:** The initial model was adapted from a passive model for the microgrid developed by Fred Fan (CEE, Stanford University)

### Operating the repository

A Julia environment is needed, containing the following packages: JuMP, CSV, DataFrames, PlotlyJS, Dates, XLSX, FileIO, Base, Random, Statistics, Gurobi, PyCall.

### Description of functions

- Main.jl contains the main implementation.
- MPCAlgorithms.jl contains the implementation of the two methods.
- ErrorFunctions.jl defines the different types of error added to weather data to simulate varying forecasts.
- MakeForecast.jl simulates the forecasts based on these error functions.
- MakeScenarios.jl defines the parameters and types of error for different scenarios tested.
- MPCValidation.jl runs different scenarios for the algorithms implemented, over chosen years.
- PassiveModel.jl contains the model describing how weather parameters translate to energy transfers for the microgrid structure.
- PlotResults.jl contains a script to plot results for the cumulative loss of load time series.
- Util.jl contains several basic functions used across the previous files. 
