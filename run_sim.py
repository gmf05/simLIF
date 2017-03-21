#!/usr/bin/python
#%
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 02:56:23 2015

@author: gmf
"""
import os
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
#import scipy.io as spio
import pp_tools as pp

MGH = os.environ['MGH']

# Time parameters
dt = 1.0e-3
T = 13e3
#T = 13e3
#T = 33e3
Nsec = T*dt
time = np.arange(0,T)*dt

#

def runSim(krate,b0,kint,bint,kspa,bspa1,bspa2):
  
  # Network structure
  N1 = 1 # Number of unconnected spike trains
  N2 = 25 # Number of interior spike trains
  N = N1+N2 # Total number of spike trains
   
  neighbors = np.zeros([N2,4])
  for n in range(N2):
    neighbors[n,0] = N1+n-1
    neighbors[n,1] = N1+n-2
    neighbors[n,2] = N1+n-3
    neighbors[n,3] = N1+n-4

  lint,yint = pp.plotspline(kint,bint)
  lspa,yspa1 = pp.plotspline(kspa,bspa1)
  not_saved,yspa2 = pp.plotspline(kspa,bspa2)
  
  # Set outer edge spike probabilities
  edge_rate = 15
  intensity = np.repeat(edge_rate * dt * np.ones([1,np.round(1/dt)]), Nsec, axis=1)
  intensity_vec = np.repeat(intensity, N1, axis=0)
  edge_dn = scipy.stats.poisson.rvs(intensity_vec)

  p = pp.params()
  p.add_covar('Rate', 0, krate, 'indicator')
  p.add_covar('Intrinsic', 1, kint, 'spline')
  p.add_covar('Spatial1', 2, kspa, 'spline')
  #int_ind = p.index[1]
  int_order = int(kint[-1])
  #spa_ind = p.index[2]
  spa_order = int(kspa[-1])
  burnin = p.get_burnin()
   
  ### BEGIN SIMULATION ###
  sim_dn = np.zeros([N, T])
  # lambda = np.zeros([N, T])
  intensity = np.zeros([N2, T])
  sim_dn[0:N1, :] = edge_dn
  
  print 'Simulation started...'
  for t in range(burnin, int(T)):
    if np.mod(t,1e3)==0: print str(t*dt)
    lambda_t = np.zeros(N2)
    for n in range(N2):
      # baseline rate:
      lambda_t[n] = b0
      # intrinsic effects:
      lambda_t[n] += np.sum(sim_dn[N1+n, t-1-np.arange(0,int_order)] * yint)
      # spatial effects:
      lambda_t[n] += np.sum(sim_dn[neighbors[n,0], t-1-np.arange(0,spa_order)] * yspa1)
      #lambda_t[n] += np.sum(sim_dn[neighbors[n,1], t-1-np.arange(0,spa_order)] * yspa2)
      #lambda_t[n] += np.sum(sim_dn[neighbors[n,2], t-1-np.arange(0,spa_order)] * yspa2)
      #lambda_t[n] += np.sum(sim_dn[neighbors[n,3], t-1-np.arange(0,spa_order)] * yspa2)
      #for k in range(1,6):
      #  if n>k: lambda_t(n) += np.sum(sim_dn[neighbors(n)-k, t-1-np.arange(0,spa_order)] * yspa2)
    intensity[:,t] = np.exp(lambda_t)
    dn_t = ( scipy.stats.poisson.rvs(np.exp(lambda_t)) > 0 )
    sim_dn[N1+np.arange(0,N2),t] = dn_t
  print 'Simulation done!'
  ### END SIMULATION ###
  if np.any(np.isnan(sim_dn)) or np.sum(np.mean(sim_dn,axis=1)>0.1): print 'Blow up??'
  return sim_dn
