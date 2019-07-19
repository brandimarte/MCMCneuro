[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![mcmcneuro logo](/doc/images/brain.png "MCMCneuro code")

## About ##

The goal of __MCMCneuron__ is to find the most representative graph model for neuronal interactions vi Markov Chain Monte Carlo.

[https://brandimarte.github.io/coding/MCMC_Neuro.html](https://brandimarte.github.io/coding/MCMC_Neuro.html)

## Description ##

### graphPenalty ###

Computes a Markov Chain Monte Carlo on graphs representing neuronal connectivity to estimate the graph that best represents the observed data acquired from single neurons. Considers different values of the penalty constant within the 3 different following methods for computing the posterior probability:

 1. <Xi|Xj>=1 if Xi=1 and Xj=1
 2. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=0 (if Xi!=Xj)
 3. <Xi|Xj>=1 if Xi=Xj or <Xi|Xj>=-1 (if Xi!=Xj)

The output files are sent to 'out/outPenalty' folder.

### bestGraph ###

For a given penalty constant, chosen based on the previous executable, it computes a Markov Chain Monte Carlo on graphs representing neuronal connectivity to estimate the most representative graph for each set of observed neurons.

The output files are sent to 'out/outBestGraph' folder.
