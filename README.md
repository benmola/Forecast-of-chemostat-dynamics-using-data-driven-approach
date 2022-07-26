# Forecast-of-chemostat-dynamics-using-data-driven-approach

Using the Koopman Operator theory for the data-driven modeling of the Chemostat model.

In this repo I will share with you my work on the data-driven modeling of the Chemostat system using the Koopman operator theory. This work has been published in  the proccedings of the International Conference on Control, Automation and Diagnosis (ICCAD) 2021 in Gronoble France, you will find all the information in the following link : https://ieeexplore.ieee.org/document/9638749.


# Abstract
This paper deals with the forecast of chemostat
dynamics using a data-driven approach. We construct a datadriven model (predictor) based on the Koopman operator theory,
which can predict the future state of the nonlinear dynamical
system of the chemostat by only measuring the input and output
of the system. We are presenting a predictor with a linear
structure, that can be used for diagnostics, state estimation
and future state prediction and control of nonlinear chemostat.
Importantly, the method of generating such linear predictors is
entirely data-driven and extremely simple, leading to nonlinear
data transformation (embedding), and a linear least squares
problem in the embedded space which can be readily solved
for large data sets. We show in simulations that Koopman
approach best predicts the system trajectories compared to a
local linearization methods.

# Code

The code is an adaption/extension of the code associated with the paper "Linear predictors for nonlinear dynamical systems: Koopman operator meets model predictive control", Automatica 2018, by Milan Korda and Igor Mezic (available under: https://github.com/MilanKorda/KoopmanMPC).

The file named 'KoopmanChemostatModel.m' contains the main script for the work just run it and the results will show off. 
Others are just functions used for the main code.

