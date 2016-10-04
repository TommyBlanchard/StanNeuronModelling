# StanNeuronModelling
A Bayesian model comparison of the tuning properties of a population of neurons, implemented in Stan

The purpose of this software is to take time-bins with spike rates that can be regressed against a task variable, and check for neurons that have no tuning to the input variables, mixed tuning (responsive to both variables), or pure tuning (only responsive to one variable).

Paper outlining this method and showing some results with it is currently under review.

Model code and general-purpose scripts are in the 'Model' folder. It also contains an example script that generates some fake data and tests the model on it.

Code used specifically for the Hayden Lab data sets is in HaydenDataScripts.

Expect updates to this repo to make it more user-friendly as the paper makes its way through the publication process
