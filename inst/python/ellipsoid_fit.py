# Ellispoid of best fit
import numpy as np
import numpy.matlib
import scipy
import scipy.stats
import sklearn

# Notation
from numpy.linalg import inv, norm
from numpy.random import multivariate_normal as multinormal
from numpy.random import normal
from scipy.linalg import null_space
from sklearn.covariance import empirical_covariance
from sklearn.decomposition import PCA

#---------- Loss function ----------#
def choose2(n):
    return int(n*(n-1)/2)

# Cayley transform
"""
Maps skew-symmetric k-by-k matrix S with entries s to a rotation matrix via
  : S --> (I - S)(I + S)^{-1}
"""
def cayley(s,k):
    S = np.zeros((k,k))
    for i in range(k):
        S[i,i+1:] = s[i*k - choose2(i+1):(i+1)*k - choose2(i+2)]
    S = S - S.T
    return (np.eye(k) - S) @ inv(np.eye(k) + S)

# Loss function for one point x
"""
Input: x in R^k
Output: ||diag(a) * R * (x - b)|| - 1 where
    * diag(a) is the diagonal matrix with entries a_1,...,a_k
    * b is the approximate center of the ellipsoid
    * R is the Cayley transform of the skew-symmetric matrix corresponding to s
"""
def loss_individual(theta,x):
    k = x.shape[0]
    a = theta[:k]
    b = theta[k:2*k]
    s = theta[2*k:]
    R = cayley(s,k)
    return norm(np.diag(a) @ R @ (x - b)) - 1

# n-dimensional array containing the loss function for each data point
def loss(theta,X):
    n = X.shape[0]
    output = np.empty(n)
    for i in range(n):
        output[i] = loss_individual(theta,X[i,:])
    return output


#---------- Parameters for trust region reflective (trr) algorithm ----------#
"""
Input: 
    * PCA-augmented data matrix X 
    * dimension k 
    * weight w
Output:
    * Initial value theta_0 = (a_0, b_0, s_0)
    * Parameter bounds (lb, ub)
"""
def trr_init(X,k,w):
    # Min and max along principal axis
    mins = np.empty(k)
    maxs = np.empty(k)
    for i in range(k):
        mins[i] = np.min(X[:,i])
        maxs[i] = np.max(X[:,i])
    M = maxs - mins

    a_lb = 1/(2 * M[0])
    a_ub = 2 * M[0]
    # a_lb = 1e-5
    # a_ub = np.inf
    b_lb = mins
    b_ub = maxs
    s_lb = -5
    s_ub = 5

    a_0 = M*np.ones(k)
    b_0 = np.average(X,axis=0)
    s_0 = np.zeros(choose2(k))

    # Set theta_0 = (a_0, b_0, s_0)
    theta_0 = np.concatenate((a_0, b_0, s_0))

    # Lower and upper bounds (lb and ub) for trust region reflective algorithm
    lb = np.empty(len(theta_0))
    for i in range(len(theta_0)):
        if i < k:
            lb[i] = a_lb
        elif i < 2*k:
            lb[i] = b_lb[i-k]
            # lb[i] = -np.inf
        else:
            lb[i] = s_lb
    ub = np.empty(len(theta_0))
    for i in range(len(theta_0)):
        if i < k:
            ub[i] = a_ub
        elif i < 2*k:
            ub[i] = b_ub[i-k]
            # ub[i] = np.inf
        else:
            ub[i] = s_ub

    return theta_0, w*lb, w*ub


#---------- Full ellipsoid fit algorithm ----------#
"""
Inputs: 
    * X: n-by-p data matrix
    * k: latent factor dimension
    * w (default 1.1): weight to apply to trust region bounds
    * ellipsoid_dimensions (default None): the k dimensions in which to fit ellipsoid
    * trr_params (default None): initial conditions for trust region algorithm. If 'None'
        applies trr_init function above. See that function for details.
Returns:
    * k: latent factor dimension
    * eigenvalues: eigenvalues of data covariance matrix
    * eigenvectors: matrix of eigenvectors of data covariance matrix
    * center: center of fitted ellipsoid
    * mu: estimated mu in vMF distribution
    * tau: estimated tau in vMF distribution
    * b: estimated center of k-dimensional ellipsoid
    * ax_lengths: axis lengths of fitted ellipsoid
    * s: s values corresponding to fitted ellipsoid
    * R: rotation of fitted ellipsoid (the Cayley transform of s)
    * Lambda: factor laoding matrix of fitted ellipsoid
    * Sigma: empirical covariance matrix for Gaussian noise
    * result: output of least squares algroithm, see 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
"""
from scipy.optimize import least_squares, root_scalar

def ellipsoid_fit(X, k, w=1.1, ellipsoid_dimensions=None, trr_params=None):
    n, p = X.shape
    avg = np.average(X, axis=0)
    X = X - avg

    #-- PCA --#
    pca = PCA(n_components = p)
    pca = pca.fit(X)
    eigenvals = pca.explained_variance_
    components = pca.components_
    V = components.T
    if ellipsoid_dimensions == None:
        Vk = V[:,:k]
        Wk = V[:,k:]
    else:
        Vk = V[:,ellipsoid_dimensions]
        perp_dims = []
        for i in range(p):
            if i in ellipsoid_dimensions:
                continue
            else:
                perp_dims.append(i)
        Wk = V[:,perp_dims]

    X_pca = X @ Vk
    perp_avg = np.average((X+avg) @ Wk @ Wk.T, axis=0)
    X_perp = X @ Wk @ Wk.T
    Q_perp = empirical_covariance(X_perp, assume_centered=True)

    #-- Trust region reflective algorithm --#
    if trr_params == None:
        theta_0, lb, ub = trr_init(X_pca,k,w)
    else:
        theta_0, lb, ub = trr_params
    result = least_squares(loss, theta_0, args=([X_pca]), method = 'trf', bounds=(lb, ub))
    a = result.x[:k]
    b = result.x[k:2*k]
    s = result.x[2*k:]
    R = cayley(s,k)

    residuals = result.fun
    cost = result.cost

    #-- MLE of von Mises-Fisher --#
    # Transform to sphere
    X_sphere = np.empty((n,k))
    for i in range(n):
        v = np.diag(a) @ R @ (X_pca[i,:] - b)
        X_sphere[i,:] = v/norm(v)

    # MLE of direction paprameter
    mu_hat = np.average(X_sphere, axis=0)
    r_hat = norm(mu_hat)
    mu_hat = mu_hat/r_hat

    # MLE of concentration parameter
    tau_hat = (r_hat * (k - r_hat**2)) / (1 - r_hat**2)

    # Error
    error = np.empty((n,p))
    M = np.diag(a) @ R
    N = Vk @ R.T @ np.diag(1/a)
    for i in range(n):
        y = Vk.T @ X[i,:]
        error[i,:] = Vk @ (y - b - (y-b)/norm(M @ (y-b)))
    Q = empirical_covariance(error)

    return {'k': k, 'eigenvalues': eigenvals, 'eigenvectors': V, 'center': avg + Vk @ b, 'mu': mu_hat, 'tau': tau_hat, 'b': b, 
                'ax_lengths': 1/a, 's': s, 'R': R, 'Lambda': Vk @ R.T @ np.diag(1/a), 'Sigma': Q + Q_perp, 'result': result}


