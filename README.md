# Rockafellar and Uryasev (2000)
An application of the paper "Rockafellar, R. and S. Uryasev (2000). “Optimization of conditional value-at risk.” Journal of Risk 3: 21-41".

I start by generating synthetic financial returns assuming that they follow a multivariate Gaussian distribution. The assumption of this distribution is not realistic but it works to keep things easy, the most important thing here is the illustration of the optimization problem formulation and how to solve it using either a linnear or a non-linear programing formulation.

Both optimization formulations yield to similar results; however, the non-linear programming is much faster than the linear programming.

The main output of the optimization problem is the vector of weigths that minimize the portfolio risk (measured by the expected shortfall). Then, those weights are used to compute tomorrow's expected returns, which are required to plot the efficient frontier.

## Efficient frontier
I provide the efficient frontier using both formulations.

## References
Rockafellar, R. and S. Uryasev. “Optimization of conditional value-at risk.” Journal of Risk 3 (2000): 21-41.
