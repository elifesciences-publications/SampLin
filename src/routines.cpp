// Many bugs fixed, everything works.

#include<iostream>
#include <random>
#include <fstream>
#include "../include/routines.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <sys/stat.h>

using namespace std;

// GIBBS
arma::vec dirichlet(gsl_rng *r, const arma::vec &alpha){
    arma::vec theta(alpha.n_elem);
    gsl_ran_dirichlet(r,alpha.n_elem,alpha.memptr(),theta.memptr());
    return(theta);
}

bool fileExists(const char* file) {
    struct stat buf;
    return (stat(file, &buf) == 0);
}

double gsl_ran_beta_pdflog(double x, double alpha, double beta){
    if(x==0 || x==1){
        std::cout<<"beta out of range"<<endl;
        std::abort();
    } else {
        return( (alpha-1)*log(x) + (beta-1)*log(1-x)) - gsl_sf_lnbeta(alpha,beta);
    }
}

double MultiBetaLog(const arma::vec &x,double alpha=0){
    double r=0;
    for(unsigned int i=0;i<x.n_elem;++i){
        r+=gsl_sf_lngamma(x(i)+alpha);
    }
    r+=-gsl_sf_lngamma(arma::accu(x)+x.n_elem*alpha);
    return(r);
}
