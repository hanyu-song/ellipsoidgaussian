##### ellipsoidgaussian
 Provides Bayesian estimation of the ellipsoid-Gaussian model
    introduced by Song and Dunson (2022). The model was inspired by the
    need for new models for characterizing dependence in multivariate
    data. The multivariate Gaussian distribution is routinely used, but
    cannot characterize nonlinear relationships in the data. Most
    non-linear extensions tend to be highly complex; for example,
    involving estimation of a non-linear regression model in latent
    variables. The relatively simple class of ellipsoid-Gaussian
    multivariate distributions, which are derived by using a Gaussian
    linear factor model involving latent variables having a von
    Mises-Fisher distribution on a unit hyper-sphere, can be used for
    modeling data distributed around a hyper-ellipsoid. In addition to
    fitting the model, the package also provides options for
    post-processing the posterior samples and generating samples from the
    posterior predictive distribution. Plotting functions are also
    included.
