/* 
 * Copyright 2015 Jarad Niemi
 * 
 * Licensed under GNU GPLv2
 *
 * MCMC inference for hierarchical normal model
 *   y_i \sim N(\theta_i,1) 
 *   \theta_i \sim N(mu,\tau^2)
 *   p(\mu) \propto 1
 *   \tau ~ Unif(0,10)
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "timer.h"

#define N 1000
#define MEAN 0.0
#define SD 1.0
#define SIGMA2 1.0
#define NREPS 100000

double y[N], theta[N], mu, tau, tau2;

// http://stackoverflow.com/questions/5287009/gaussian-random-number-generator
double rnorm(double mean, double sd) 
{
  double v1,v2,s;

  do {
    v1 = 2.0 * ((double) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((double) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0;
  else
    return mean + sd*(1*sqrt(-2.0 * log(s) / s));
}




int main(int argc, char** argv)
{
  int i,j;
  const unsigned int n=N;
  const unsigned int nreps=NREPS;

  double mean=MEAN, sd=SD, sigma2=SIGMA2, m, V, ybar;

  memset(y,     0, n*sizeof(double));
  memset(theta, 0, n*sizeof(double));
  for (i=0; i<n; i++) 
  {
    y[i] = rnorm(mean,sqrt(1+sd));
  }

  StartTimer();

  // initial values
  mu = 0;
  tau = 1; tau2 = tau*tau;

  for (i=0; i<nreps; i++) {
    // sample theta
    for (j=0; i<n; j++) {
      V = 1.0/(1.0+1.0/tau2);
      m = V*y[j]/1.0;    
      theta[j] = rnorm(m, sqrt(V));
    }
    
    ybar = 0;
    for (j=0; j<n; j++) {
      ybar += theta[j];
    }
    ybar /= n;

    // sample mu and tau
    mu = rnorm(ybar, sqrt(sigma2/n));
  }

  double runtime = GetTimer();
 
  printf(" total: %f s\n", runtime / 1000);
}

