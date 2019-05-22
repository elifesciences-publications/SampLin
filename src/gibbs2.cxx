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

    unsigned int i,j,k,l,mu,nu;
    int NITER=atoi(argv[1]);
    int BURN_IN=atoi(argv[2]);
    int TRIM=atoi(argv[3]);
    int N;
    int M=20;
    int sample=0;
    int P=atoi(argv[4]); // initial number of lineages
    int RNGSEED=atoi(argv[5]);
    bool PFIX=false;
	bool verb=false;
	double log_alpha_MH;

    // prior hyper-parameters
    double alpha_n=1;
    double alpha_p=1;
    double beta_p=1;
    double alpha_q=1;
    double beta_q=1;
    double alpha_DP=.1;

	// input matrix s and augmented matrix As
    arma::mat s,As;

    // gsl random number generator
    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);
    gsl_ran_discrete_t* gen=NULL;

    // Load count matrix
    s.load(argv[6]);
    As=s;
    arma::vec lineage_size=arma::sum(As,1);
    arma::vec lineage_size_from_s=arma::sum(s,1);

    N=s.n_rows;
    cout<<"count matrix loaded. N="<< N<<endl;

    // streams
    ofstream LikelihoodOutput(proc_folder+"/F.dat");
    ofstream membershipTraj(proc_folder+"/membership_traj.dat");
    ofstream PTraj(proc_folder+"/P.dat");
    ofstream qmukTraj(proc_folder+"/qmuk.dat");
    ofstream pmuTraj(proc_folder+"/pmu.dat");
    ofstream confTraj(proc_folder+"/configurations.dat");
    ofstream occTraj(proc_folder+"/occupancy.dat");

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

	// initial augmentation
	for(i=0;i<N;i++){
		for(k=0;k<stop[i];k++){
			As(i,k) = gsl_ran_binomial(r,0.6,lineage_size(i));
		}
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
        }
        T_vec=sum(T_mat);

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
			arma::vec As_row_i = (As.row(i)).t();

            // calculate metropolis-hastings alpha

            if(mu_new<P){
                log_alpha_MH=gsl_sf_lnbeta(T_vec(mu_new)+lineage_size(i)+alpha_p,(group_sizes(mu_new)+1)*M-T_vec(mu_new)-lineage_size(i)+beta_p)
                         +gsl_sf_lnbeta(T_vec(mu0)+alpha_p,group_sizes(mu0)*M-T_vec(mu0)+beta_p)
                         -gsl_sf_lnbeta(T_vec(mu_new)+alpha_p,group_sizes(mu_new)*M-T_vec(mu_new)+beta_p)
                         -gsl_sf_lnbeta(T_vec(mu0)+lineage_size(i)+alpha_p,(group_sizes(mu0)+1)*M-T_vec(mu0)-lineage_size(i)+beta_p)

                         +MultiBetaLog(T_mat.col(mu_new)+As_row_i,alpha_q)
                         +MultiBetaLog(T_mat.col(mu0),alpha_q)
                         -MultiBetaLog(T_mat.col(mu_new),alpha_q)
                         -MultiBetaLog(T_mat.col(mu0)+As_row_i,alpha_q);
            } else {
                log_alpha_MH=gsl_sf_lnbeta(lineage_size(i)+alpha_p,M-lineage_size(i)+beta_p)
                        +gsl_sf_lnbeta(T_vec(mu0)+alpha_p,group_sizes(mu0)*M-T_vec(mu0)+beta_p)
                        -gsl_sf_lnbeta(alpha_p,beta_p)
                        -gsl_sf_lnbeta(T_vec(mu0)+lineage_size(i)+alpha_p,(group_sizes(mu0)+1)*M-T_vec(mu0)-lineage_size(i)+beta_p)

                        +MultiBetaLog(As_row_i,alpha_q)
                        +MultiBetaLog(T_mat.col(mu0),alpha_q)
                        -MultiBetaLog(arma::zeros(4),alpha_q)
                        -MultiBetaLog(T_mat.col(mu0)+As_row_i,alpha_q);
            }

            if(gsl_rng_uniform(r)<exp(log_alpha_MH)){

                ids[i]=mu_new;
                
				if(mu_new==P){
                    T_mat.insert_cols(P,1);
                    T_vec.insert_cols(P,1);
                    group_sizes.insert_rows(P,1);
                    P=P+1;
                }

            }
            
            for(k=0;k<4;k++){
                T_mat(k,ids[i])+=As(i,k);
                T_vec(ids[i])+=As(i,k);
            }
            group_sizes(ids[i])+=1;
			
			if(group_sizes(mu0)==0){

                P--;
                if(verb) cout<<"delete "<<mu0<<endl;
                for(l=0;l<N;++l){
                    if(ids[l]>mu0) ids[l]--;
                }

                T_mat.shed_col(mu0);
                T_vec.shed_col(mu0);
				
                group_sizes.shed_row(mu0);

            }

        }

        arma::mat qmuk(4,P,arma::fill::zeros);
        arma::vec pmu(P,arma::fill::zeros);
        arma::vec conf(N,arma::fill::zeros);
        arma::Col<unsigned int> tmp_conf(4);
		arma::vec qmuk_mu(4);

        for(mu=0;mu<P;mu++){
            qmuk.col(mu)=dirichlet(r,T_mat.col(mu)+alpha_q);
            pmu(mu)=gsl_ran_beta(r,T_vec(mu)+alpha_p,group_sizes(mu)*M-T_vec(mu)+beta_p);
        }

        // data augmentation here
		for(i=0;i<N;i++){
			// generate clonal size by requiring a binomial number larger than
			// the clonal size from s
			if(stop[i]>0){
			   int clonal_size_tmp=0;
			   int counter=0; 
   			   while(clonal_size_tmp<lineage_size_from_s(i) && counter<1000){ 
				   counter++;
				   clonal_size_tmp = gsl_ran_binomial(r,pmu(ids[i]),M);
			   }
			   
			   if(counter<1000){
				   cout<<clonal_size_tmp<<' '<<lineage_size_from_s(i)<<endl;
				   double* tmp_vec = new double[stop[i]];
				   for(k=0;k<stop[i];k++) tmp_vec[k]=qmuk(k,ids[i]);
				   unsigned int*  n_tmp = new unsigned int[stop[i]];
   				   gsl_ran_multinomial(r, stop[i], clonal_size_tmp-lineage_size_from_s(i),tmp_vec,n_tmp);
				   for(k=0;k<stop[i];k++) As(i,k)=n_tmp[k];
				   delete[] tmp_vec,n_tmp;
			   } else {
				   for(k=0;k<stop[i];k++) As(i,k)=0;
			   }

			}
		}

		lineage_size = sum(As,1);

        if(sample%TRIM==0 && sample>BURN_IN){
            double F=0;
            for(nu=0;nu<P;nu++){
                F+=0;//gsl_sf_lngamma(group_sizes(nu)+1);
                for(k=0;k<4;++k) {
                    //F+=gsl_sf_lnbeta(tmp_mat(k,nu)+1,MxLayer(k)*group_sizes(nu)-tmp_mat(k,nu)+1);
                       //-gsl_sf_lngamma(N+P)+gsl_sf_lngamma(P);
                }
            }
            for(i=0;i<N;++i) membershipTraj<<ids[i]<<' '; membershipTraj<<endl;
            LikelihoodOutput<<F<<endl;

            PTraj<<P<<endl;

            for(mu=0;mu<min(P,10);mu++){
                qmukTraj<<qmuk(0,mu)<<' '<<qmuk(1,mu)<<' '<<qmuk(2,mu)<<' '<<qmuk(3,mu)<<' ';
            }
            while(mu<10){
                for(k=0;k<4;k++) qmukTraj<<"NA"<<' ';
                mu++;
            }
            qmukTraj<<endl;

            for(i=0;i<N;i++){
				qmuk_mu = qmuk.col(ids[i]);
				gsl_ran_multinomial(r,4,lineage_size(i),qmuk_mu.memptr(),tmp_conf.memptr());
                for(k=0;k<4;k++){
                    if(k>=stop[i]){
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

    }

    membershipTraj.close();

    return 0;
}
