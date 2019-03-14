import numpy as np

from numpy       import pi, sqrt, log, exp
from scipy.stats import norm

#
# Normal cumulative distribution function, with extension
# for complex argument with small imaginary component
#

def norm_cdf(x):
    if not isinstance(x, np.ndarray):
        xr = x.real
        xi = x.imag
        if abs(xi) > 1.0e-10:
            raise ValueError('imag(x) too large in norm_cdf(x)')

        ncf = norm.cdf(xr)
        if abs(xi) > 0:
            ncf = ncf + 1.0j*xi*norm.pdf(xr)
    else:
        xr = np.real(x)
        xi = np.imag(x)
        if any(abs(xi) > 1.0e-10):
            raise ValueError('imag(x) too large in norm_cdf(x)')

        ncf = norm.cdf(xr)
        if any(abs(xi) > 0):
            ncf = ncf + 1.0j*xi*norm.pdf(xr)

    return ncf
# V = european_call(r,sigma,T,S,K,opt)
# r     - interest rate
# sigma - volatility
# T     - time interval
# S     - asset value(s)  (float or flattened numpy array)
# K     - strike price(s) (float or flattened numpy array)
# opt   - 'value'
# V     - option value(s) (float or flattened numpy array)


def european_call(r,D,sigma,T,S,K):

    S  = S + 1.0e-100     # avoids problems with S=0
    K  = K + 1.0e-100     # avoids problems with K=0

    d1 = ( log(S) - log(K) + (r-D +0.5*sigma**2)*T ) / (sigma*sqrt(T))
    d2 = ( log(S) - log(K) + (r-D - 0.5*sigma**2)*T ) / (sigma*sqrt(T))

    
    V = S*exp(-T*D)*norm_cdf(d1) - exp(-r*T)*K*norm_cdf(d2)

    return V


def european_put (r,D, sigma,T,S,K):

    S  = S + 1.0e-100    
    K  = K + 1.0e-100     

    d1 = ( log(S) - log(K) + (r- D+0.5*sigma**2)*T ) / (sigma*sqrt(T))
    d2 = ( log(S) - log(K) + (r- D-0.5*sigma**2)*T ) / (sigma*sqrt(T))

     
    V = - S*exp(-T*D)*norm_cdf(-d1) + exp(-r*T)*K*norm_cdf(-d2)

    return V


def binary_call (r,D,sigma,T,S,K):
    
    

    S  = S + 1.0e-100     
    K  = K + 1.0e-100     

    d1 = ( log(S) - log(K) + (r-D +0.5*sigma**2)*T ) / (sigma*sqrt(T))
    d2 = ( log(S) - log(K) + (r-D -0.5*sigma**2)*T ) / (sigma*sqrt(T))

     
    V = exp(-r*T)*norm_cdf(d2)

    return V


