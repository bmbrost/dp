# mixtures

MCMC algorithms for implementing Dirichlet process mixture models

- **dp.utils.R**: Utility functions for Dirichlet proceses, e.g., functions for simulating data according to Dirichlet processes.
- **dp.mixture.sim.R**: Simulate Dirichlet process mixture data; fit Escobar, Neal, and blocked Gibbs sampler models; and examine model results. 
- **dp.mixture.escobar.1994.mcmc.R**: Fit Dirichlet mixture model based on Eq. 3, Escobar (1994). Note that this code needs to be vetted.
- **dp.mixture.neal.2000.algm.2.mcmc.R**: Fit Dirichlet mixture model using Algorithm 2, Neal (2000). Note that this code needs to be vetted.
- **dp.mixture.neal.2000.algm.5.mcmc.R**: Fit Dirichlet mixture model using Algorithm 5, Neal (2000). Note that this code needs to be vetted.
- **dp.mixture.neal.2000.algm.7.mcmc.R**: Fit Dirichlet mixture model using Algorithm 7, Neal (2000). Note that this code needs to be vetted.
- **dp.mixture.neal.2000.algm.8.mcmc.R**: Fit Dirichlet mixture model using Algorithm 8, Neal (2000). Note that this code needs to be vetted.
- **dp.mixture.blocked.mcmc.R**: Fit Dirichlet mixture model using blocked Gibbs sampler and the truncation approximation (Ishwaran and James 2001, Section 5; Gelman et al. 2014, BDA Section 23.3)
- **dp.mixture.blocked.2d.mcmc.R**: Fit Dirichlet mixture model using blocked Gibbs sampler to 2-dimensional data.
- **dp.mixture.blocked.2d.sim.R**: Simulate Dirichlet process mixture data in 2 dimensions, fit blocked Gibbs sampler models, and examine model results. 
