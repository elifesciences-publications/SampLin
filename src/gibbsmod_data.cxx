#include<iostream>
#include<cstdio>
#define ARMA_NO_DEBUG
#include<armadillo>
#include "../include/routines.h"
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf.h>
#include<cmath>
#include<omp.h>
#include<sys/stat.h>
#include<ctime>

using namespace std;

int main(int argc, char* argv[]){


    if(argc<8){
        cerr<<"Bayesian analysis"<<endl
            <<"usage:"<<endl
            <<"./gibbs_data <NITER> <BURN_IN> <TRIM> <LINEAGES> <SEED> <INPUT_FILE> <folder>"<<endl;
        return 1;
    }

    // Create output folder
    string proc_folder(argv[7]);
    struct stat sb;
    if (stat(proc_folder.c_str(), &sb) != 0){
        const int dir_err = mkdir(proc_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory!");
            exit(1);
        }
    }

    int NITER=atoi(argv[1]);
    int BURN_IN=atoi(argv[2]);
    int TRIM=atoi(argv[3]);
    int N;
    int M=20;
    int P=atoi(argv[4]); // initial number of lineages
    int RNGSEED=atoi(argv[5]);
    bool PFIX=false;

    arma::fmat s;

    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);

    s.load(argv[6]);

    N=s.n_rows;

    cout<<N<<' '<<M<<endl;

    // streams
    ofstream FreeEnergyOutput(proc_folder+"/F.dat");
    ofstream membershipTraj(proc_folder+"/membership_traj.dat");
    ofstream PTraj(proc_folder+"/P.dat");
    ofstream pmukTraj(proc_folder+"/pmuk.dat");
    ofstream confTraj(proc_folder+"/configurations.dat");

    // output variables

    // INIT

    int ids[N];

    for(unsigned int i=0;i<N;i++) ids[i] = gsl_rng_uniform_int(r,P);

    unsigned int i,j,k,mu,nu;
    double F0;

    gsl_ran_discrete_t* gen=NULL;

    int sample=0;

    // Define stops
    vector<int> stop(N);
    for(i=0;i<N;++i){
        k=4; while(s(i,k)==0) k--;
        stop[i]=k;
    }

    while(sample<NITER){
        sample++;
        //cout<<"sample "<<sample<<'\r'<<flush;

        arma::mat ids_mat(P,N,arma::fill::zeros);
        arma::vec group_sizes(P);
        for(i=0;i<N;++i) ids_mat(ids[i],i)=1;
        group_sizes=arma::sum(ids_mat,1);

        arma::mat tmp_mat(4,P,arma::fill::zeros);
        for(i=0;i<N;++i){
            for(k=0;k<stop[i];k++){
                tmp_mat(k,ids[i])+=s(i,k);
            }
        }

        for(i=0;i<N;++i){
            double tmp_prob_ln[P];
            double tmp_prob[P];

            for(k=0;k<stop[i];k++) tmp_mat(k,ids[i])-=s(i,k);
            group_sizes(ids[i])+=-1;

            for(mu=0;mu<P;++mu){
                tmp_prob_ln[mu]=0;
                for(nu=0;nu<P;nu++){
                    tmp_prob_ln[mu]+=(nu!=mu) ? gsl_sf_lngamma(group_sizes(nu)+1)
                                              : gsl_sf_lngamma(group_sizes(nu)+2);
                    for(k=0;k<stop[i];++k) {
                        //cout<<tmp_mat(k,nu)+s(i,k)+1<<' '<<M*(group_sizes(nu)+1)-tmp_mat(k,nu)-s(i,k)+1<<endl;
                        tmp_prob_ln[mu]+=(nu!=mu) ? gsl_sf_lnbeta(tmp_mat(k,nu)+1,M*(group_sizes(nu))-tmp_mat(k,nu)+1)
                                                  : gsl_sf_lnbeta(tmp_mat(k,nu)+s(i,k)+1,M*(group_sizes(nu)+1)-tmp_mat(k,nu)-s(i,k)+1);
                    }
                }

            }

            for(mu=0;mu<P;++mu) tmp_prob[mu]=exp(tmp_prob_ln[mu]-tmp_prob_ln[0]);

            gen = gsl_ran_discrete_preproc(P, tmp_prob);

            ids[i]=gsl_ran_discrete(r,gen);
            group_sizes(ids[i])++;
            for(k=0;k<stop[i];k++) tmp_mat(k,ids[i])+=s(i,k);
            gsl_ran_discrete_free(gen);
            gen=NULL;

        }

        arma::mat pmuk(4,P,arma::fill::zeros);
        arma::vec conf(N,arma::fill::zeros);
        arma::vec tmp_conf(4);

        for(k=0;k<4;k++){
            for(mu=0;mu<P;mu++){
                pmuk(k,mu)=gsl_ran_beta(r,tmp_mat(k,mu),M*group_sizes(mu)-tmp_mat(k,mu));
            }
        }

        for(i=0;i<N;i++){
            for(k=0;k<4;k++){
                tmp_conf(k) = gsl_ran_binomial(r,pmuk(k,ids[i]),M);
                if(tmp_conf(k)>0) conf(i)+=pow(2,k);
            }
        }


        if(sample%TRIM==0 && sample>BURN_IN){

            double F=0;
            for(nu=0;nu<P;nu++){
                F+=gsl_sf_lngamma(group_sizes(nu)+1);
                for(k=0;k<4;++k) {
                    F+=gsl_sf_lnbeta(tmp_mat(k,nu)+1,M*group_sizes(nu)-tmp_mat(k,nu)+1)
                       -gsl_sf_lngamma(N+P)+gsl_sf_lngamma(P);
                }
            }
            for(i=0;i<N;++i) membershipTraj<<ids[i]<<' '; membershipTraj<<endl;
            FreeEnergyOutput<<F<<endl;

            PTraj<<P<<endl;

			for(mu=0;mu<P;mu++){
				pmukTraj<<pmuk(0,mu)<<' '<<pmuk(1,mu)<<' '<<pmuk(2,mu)<<' '<<pmuk(3,mu)<<' ';
			}
			while(mu<10){
				for(k=0;k<4;k++) pmukTraj<<"NA"<<' '; 
				mu++;
			}
			pmukTraj<<endl;

            for(i=0;i<N;++i) confTraj<<conf(i)<<' '; confTraj<<endl;
        }

        if(!PFIX){
            if(gsl_rng_uniform(r)<0.5){
                if(gsl_rng_uniform(r)<(double)P/(N+P)){
                    P+=1;
                }
            } else {
                for(mu=0;mu<P;mu++){
                    if(group_sizes(mu)==0){
                        P--;
                        for(i=0;i<N;++i){
                            if(ids[i]>mu) ids[i]--;
                        }
                        break;
                    }
                }

            }
        }

    }

    membershipTraj.close();

    return 0;
}
