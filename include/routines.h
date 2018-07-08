#ifndef _ROUTINES_H
#define _ROUTINES_H

#include<armadillo>
#include<fstream>
#include<gsl/gsl_rng.h>

using namespace std;


// GIBBS
arma::vec dirichlet(gsl_rng *r, const arma::vec &alpha);

// Useful functions
bool fileExists(const char* file);
double gsl_ran_beta_pdflog(double x, double alpha, double beta);


#endif
