/* 
 * Copyright 2015 Jarad Niemi
 * 
 * Licensed under GNU GPLv2
 *
 * EM for hierarchical normal model
 *   y_i \sim N(\theta_i,\sigma^2) 
 *   \theta_i \sim N(\mu,\tau^2)
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "timer.h"

#define N 1000
#define SIGMA2 1.0
#define MU 0.5
#define TAU2 1.0
#define NREPS 100

// http://stackoverflow.com/questions/5287009/gaussian-random-number-generator
// I don't even use this.
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

  double v, mu=MU, tau2=TAU2, sigma2=SIGMA2, y[n], theta[n];

  memset(y,     0, n*sizeof(double));
  memset(theta, 0, n*sizeof(double));
  
  // Fill with data
  for (i=0; i<n; i++) y[i] = rnorm(mu,sqrt(sigma2+tau2)); // Aren't you calling the function here?

  StartTimer();

  // initial values
  mu = 0.3;
  tau2 = 1.3; 

#pragma acc data copy(mu,tau2) create(v,theta)
#pragma acc kernels
  for (i=0; i<nreps; i++) {
    v = 1.0/(1.0/tau2+1.0/sigma2);
    
    // E-step
    for (j=0; j<n; j++) theta[j] = v*(mu/tau2+y[j]/sigma2);

    // M-step
    mu = 0.0;
    for (j=0;j<n; j++) mu += theta[j];
    mu /= n;
    
    tau2 = 0.0;
    for (j=0;j<n; j++) tau2 += (theta[j]-mu)*(theta[j]-mu);
    tau2 /= n;
    tau2 += v;
  }

  double runtime = GetTimer();
 
  printf(" mu= %f\n", mu);
  printf(" tau2= %f\n", tau2);
  printf(" total: %f s\n", runtime / 1000);
}

