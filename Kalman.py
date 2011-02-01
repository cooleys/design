
# Code from Chapter 15 of Machine Learning: An Algorithmic Perspective
# by Stephen Marsland (http://seat.massey.ac.nz/personal/s.r.marsland/MLBook.html)

# You are free to use, change, or redistribute the code in any way you wish for
# non-commercial purposes, but please maintain the name of the original author.
# This code comes with no warranty of any kind.

# Stephen Marsland, 2008

# The 1D Kalman filter

from pylab import *
from numpy import *
original_data = 0

def Kalman(obs=None,mu_init=array([0]),cov_init=0.1*ones((1)),nsteps=500,f=False):

    ndim = shape(mu_init)[0]

    if obs is None:
        mu_init = tile(mu_init,(1,nsteps))
        cov_init = tile(cov_init,(1,nsteps))
        obs = random.normal(mu_init,cov_init,(ndim,nsteps))
        for i in range(1,nsteps):
            obs[0][i] = obs[0][i] + (i*i*0.01)
        print obs
        global original_data
        original_data = obs[0]
        print original_data
    Sigma_x = eye(ndim)*5e-3
    A = eye(ndim)
    H = eye(ndim)
    mu_hat = 0
    cov = eye(ndim)
    R = eye(ndim)*0.01

    m = zeros((ndim,nsteps),dtype=float)
    ce = zeros((ndim,nsteps),dtype=float)

    for t in range(1,nsteps):
        # Make prediction
        mu_hat_est = dot(A,mu_hat)
        print "mu_hat ",mu_hat_est
        cov_est = dot(A,dot(cov,transpose(A))) + Sigma_x
        print "cov_est ",cov_est

        # Update estimate
        error_mu = obs[:,t] - dot(H,mu_hat_est)
        print "err_mu ",error_mu
        error_cov = dot(H,dot(cov,transpose(H))) + R
        print "err_cov ",error_cov
        K = dot(dot(cov_est,transpose(H)),linalg.inv(error_cov))
        print "K ",K
        mu_hat = mu_hat_est + dot(K,error_mu)
        #m[:,:,t] = mu_hat
        m[:,t] = mu_hat
        if ndim>1:
            cov = dot((eye(ndim) - dot(K,H)),cov_est)
        else:
            cov = (1-K)*cov_est
        ce[:,t] = cov

    if not f:
        return m
    else:
        figure()
        plot(original_data[:],'ro',ms=1)
        plot(obs[0,:],'g-',ms=1)
        plot(m[0,:],'k-',lw=1)
        plot(m[0,:]+20*ce[0,:],'b--',lw=1)
        plot(m[0,:]-20*ce[0,:],'k--',lw=1)
        #legend(['Noisy Datapoints','Kalman estimate','Covariance'])
        xlabel('Time')


    show()

Kalman(f=True, nsteps=50)

#Kalman(obs=r,f=True)

