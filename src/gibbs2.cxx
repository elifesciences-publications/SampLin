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

    // Init message
    if(argc<8){
        cerr<<"Bayesian analysis"<<endl
            <<"usage:"<<endl
            <<"./gibbs_data <NITER> <BURN_IN> <TRIM> <LINEAGES> <SEED> <INPUT_FILE> <folder>"<<endl;
        return 1;
    }
    // ------------------------------------------------------------------------------------------------

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
    // ------------------------------------------------------------------------------------------------
    // Definitions

    unsigned int i,j,k,mu,nu;
    int NITER=atoi(argv[1]);
    int BURN_IN=atoi(argv[2]);
    int TRIM=atoi(argv[3]);
    int N;
    int M=20;
    int sample=0;
    int P=atoi(argv[4]); // initial number of lineages
    vector<int> MxLayer(P);
    int RNGSEED=atoi(argv[5]);
    bool PFIX=false;

    // prior hyper-parameters
    double alpha_n=1;
    double alpha_p=1;
    double beta_p=1;
    double alpha_q=1;
    double beta_q=1;
    double alpha_DP=1;

    arma::fmat s,As;

    // gsl random number generator
    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);
    gsl_ran_discrete_t* gen=NULL;

    // Load count matrix
    s.load(argv[6]);
    As=s;
    arma::vec lineage_size=arma::sum(As,1);

    N=s.n_rows;
    cout<<"count matrix loaded. N=", N<<endl;

    // streams
    ofstream LikelihoodOutput(proc_folder+"/F.dat");
    ofstream membershipTraj(proc_folder+"/membership_traj.dat");
    ofstream PTraj(proc_folder+"/P.dat");
    ofstream pmukTraj(proc_folder+"/pmuk.dat");
    ofstream confTraj(proc_folder+"/configurations.dat");
    ofstream occTraj(proc_folder+"/occupancy.dat");
    ofstream MxLayerTraj(proc_folder+"/MxLayer.dat");

    // ------------------------------------------------------------------------------------------------
    // Initialization

    // memberships
    int ids[N];
    for(unsigned int i=0;i<N;i++) ids[i] = gsl_rng_uniform_int(r,P);

    // Define stops as the first non zero element in s. For instance if the occupancy is [0,0,1,2] then stop=2.
    vector<int> stop(N);
    for(i=0;i<N;++i){
        k=0; while(s(i,k)==0) k++;
        stop[i]=k;
    }

    // ------------------------------------------------------------------------------------------------
    // Gibbs iterations

    sample=0;
    while(sample<NITER){
        sample++;

        // binary matrix of memberships
        arma::mat ids_mat(P,N,arma::fill::zeros);
        for(i=0;i<N;++i) ids_mat(ids[i],i)=1;

        // group sizes
        arma::vec group_sizes(P);
        group_sizes=arma::sum(ids_mat,1);

        // matrix and vector of counts per class per layer
        arma::mat T_mat(4,P,arma::fill::zeros);
        arma::rowvec T_vec(P,arma::fill::zeros);

        for(i=0;i<N;++i){
            for(k=0;k<4;k++){
                T_mat(k,ids[i])+=As(i,k);
            }
            T_vec(i)=sum(T_mat);
        }

        // draw type class for all lineages
        for(i=0;i<N;++i){

            // current class
            int mu0 = ids[i];

            // remove lineage from its class
            group_sizes(mu0)+=-1;

            // remove counts from T_mat and T_vec
            for(k=0;k<4;k++){
                T_mat(k,mu0)-=As(i,k);
                T_vec(mu0)-=As(i,k);
            }

            // draw class
            double tmp_prob[P+1];
            double mu_new;
            for(mu=0;mu<P;mu++) tmp_prob[mu]=group_sizes(mu);
            tmp_prob[P]=alpha_DP;

            gen = gsl_ran_discrete_preproc(P+1, tmp_prob);
            mu_new=gsl_ran_discrete(r,gen);

            gsl_ran_discrete_free(gen);
            gen=NULL;

            // calculate metropolis-hastings alpha
            if(mu_new<P){
                log_alpha_MH=gsl_sf_lnbeta(T_vec(mu_new)+lineage_size(i)+alpha_p,(group_sizes(mu_new)+1)*M-T_vec(mu_new)-lineage_size(i)+beta_p)
                         +gsl_sf_lnbeta(T_vec(mu0)+alpha_p,group_sizes(mu0)*M-T_vec(mu0)+beta_p)
                         -gsl_sf_lnbeta(T_vec(mu_new)+alpha_p,group_sizes(mu_new)*M-T_vec(mu_new)+beta_p)
                         -gsl_sf_lnbeta(T_vec(mu0)+lineage_size(i)+alpha_p,(group_sizes(mu0)+1)*M-T_vec(mu0)-lineage_size(i)+beta_p)

                         +MultiBetaLog(T_mat.col(mu_new)+As.row(i).as_col(),alpha_q)
                         +MultiBetaLog(T_mat.col(mu0),alpha_q)
                         -MultiBetaLog(T_mat.col(mu_new),alpha_q)
                         -MultiBetaLog(T_mat.col(mu0)+As.row(i).as_col(),alpha_q);
            } else {
                log_alpha_MH=gsl_sf_lnbeta(lineage_size(i)+alpha_p,M-lineage_size(i)+beta_p)
                        +gsl_sf_lnbeta(T_vec(mu0)+alpha_p,group_sizes(mu0)*M-T_vec(mu0)+beta_p)
                        -gsl_sf_lnbeta(alpha_p,beta_p)
                        -gsl_sf_lnbeta(T_vec(mu0)+lineage_size(i)+alpha_p,(group_sizes(mu0)+1)*M-T_vec(mu0)-lineage_size(i)+beta_p)

                        +MultiBetaLog(As.row(i).as_col(),alpha_q)
                        +MultiBetaLog(T_mat.col(mu0),alpha_q)
                        -MultiBetaLog(arma::zeros(4),alpha_q)
                        -MultiBetaLog(T_mat.col(mu0)+As.row(i).as_col(),alpha_q);
            }

            if(gsl_rng_uniform(r)<exp(log_alpha_MH)){

                if(mu_new==P){
                    T_mat.insert_cols(P,1);
                    T_vec.insert_cols(P,1);
                    group_sizes.insert_rows(P,1);
                    P=P+1;
                }

                ids[i]=mu_new;
                T_mat(mu_new)+=As(i,k);
                T_vec(mu_new)+=lineage_size(i);
                group_sizes(mu_new)+=1;
            }

        }

        arma::mat qmuk(4,P,arma::fill::zeros);
        arma::vec pmu(P,arma::fill::zeros);
        arma::vec conf(N,arma::fill::zeros);
        arma::vec tmp_conf(4);

        for(mu=0;mu<P;mu++){
            qmuk.col(mu)=dirichlet(r,T_mat.col(mu)+alpha_q);
            pmu(mu)=gsl_ran_beta(r,T_vec(mu)+alpha_p,group_sizes(mu)*M-T_vec(mu)+beta_p);
        }

        // data augmentation here



        if(sample%TRIM==0 && sample>BURN_IN){
            double F=0;
            for(nu=0;nu<P;nu++){
                F+=gsl_sf_lngamma(group_sizes(nu)+1);
                for(k=0;k<4;++k) {
                    F+=gsl_sf_lnbeta(tmp_mat(k,nu)+1,MxLayer(k)*group_sizes(nu)-tmp_mat(k,nu)+1);
                       //-gsl_sf_lngamma(N+P)+gsl_sf_lngamma(P);
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

            for(i=0;i<4;i++) MxLayerTraj<<MxLayer(i)<<' ';
            MxLayerTraj<<endl;

            for(i=0;i<N;i++){
                for(k=0;k<4;k++){
                    if(k>=stop[i]){
                        tmp_conf(k) = gsl_ran_binomial(r,pmuk(k,ids[i]),M);
                        if(tmp_conf(k)>0) conf(i)+=pow(2,k);
                    } else {
                        tmp_conf(k)=0;
                    }
                    occTraj<< tmp_conf(k)<<' ';
                }
            }
            occTraj<<endl;

            for(i=0;i<N;++i) confTraj<<conf(i)<<' '; confTraj<<endl;

        }

        if(!PFIX){
            if(gsl_rng_uniform(r)<0.1){
                if(gsl_rng_uniform(r)<(double)P/(N+P)){
                    P+=1;
                    if(P==11) P=10;
                }
            }

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

    membershipTraj.close();

    return 0;
}
