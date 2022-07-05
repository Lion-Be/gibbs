# Programming a Gibbs sampler
## :star: Purpose
This code consists of a programming exercise that I worked on in the first semester of my M.Sc. studies in [Methodology and Statistics](https://www.uu.nl/masters/en/methodology-and-statistics-behavioural-biomedical-and-social-sciences) at the [University of Utrecht](https://www.uu.nl/en). I manually compute a Gibbs sampler to construct posterior distributions in linear regression models. 

Specifically, I 
- compute the OLS solution using linear algebra to derive starting values for sampling
- run two separate MCMC chains sampling from the conditional posterior distribution of each parameter
- construct measures for autocorrelation and convergence

## :gem: How can you use it?
This exercise was primary a programming task as part of my studies. However, it did come in handy for me when I needed a quick Gibbs sampler for Bayesian inference that does not need any `Stan` compilation. 
