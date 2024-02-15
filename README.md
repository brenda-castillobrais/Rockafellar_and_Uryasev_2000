# Rockafellar and Uryasev (2000)
An application of the paper "Rockafellar, R. and S. Uryasev (2000). “Optimization of conditional value-at risk.” Journal of Risk 3: 21-41".

I start by generating synthetic financial returns series assuming that they follow a multivariate Gaussian distribution. This is a simplifying assumption so that we focus in the optimization process. The most important thing here is to illustrate the optimization problem and how to solve it in R using either a linear or non-linear programming formulation.

Both formulations yield similar results; however, the non-linear programming is much faster than the linear programming.

The main output of the optimization problem is the vector of weigths that minimize portfolio's risk (measured by the expected shortfall). Then, those weights are used to compute one-day-ahead expected returns, which are required to plot the efficient frontier.

## Efficient frontier
I illustrate the efficient frontier obtained with both formulations.

## References
Rockafellar, R. and S. Uryasev. “Optimization of conditional value-at risk.” Journal of Risk 3 (2000): 21-41.
