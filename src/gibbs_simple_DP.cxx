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

void create_output_folder(char* folder){
    string proc_folder(folder);
    struct stat sb;
    if (stat(proc_folder.c_str(), &sb) != 0){
        const int dir_err = mkdir(proc_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory!");
            exit(1);
        }
    }
}

int main(int argc, char* argv[]){


    if(argc<8){
        cerr<<"Bayesian analysis"<<endl
            <<"usage:"<<endl
            <<"./gibbs_data <NITER> <BURN_IN> <TRIM> <LINEAGES> <SEED> <INPUT_FILE> <folder>"<<endl;
        return 1;
    }

    // Create output folder
	create_output_folder(argv[7]);
    string proc_folder(argv[7]);

    int NITER=atoi(argv[1]);
    int BURN_IN=atoi(argv[2]);
    int TRIM=atoi(argv[3]);
    int N;
    int M=20;
    int P=atoi(argv[4]); // initial number of lineages
    int RNGSEED=atoi(argv[5]);
    bool PFIX=true;
    int verbose=0;
	double alpha_p=1, beta_p=1;
	double alpha_DP=.1;
	double log_alpha_MH;

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
    ofstream occTraj(proc_folder+"/occupancy.dat");

    // output variables

    // INIT

    int ids[N];

    for(unsigned int i=0;i<N;i++) ids[i] = gsl_rng_uniform_int(r,P);

    unsigned int i,j,k,mu,nu;
    gsl_ran_discrete_t* gen=NULL;

    int sample=0;

    // Define stops as the first non zero element in s. For instance if the occupancy is [0,0,1,2] then stop=2.
    vector<int> stop(N);
    for(i=0;i<N;++i){
        k=0; while(s(i,k)==0) k++;
        stop[i]=k;
    }
	// ------------------------------------------------------

    while(sample<NITER){
        sample++;
        //cout<<"sample "<<sample<<'\r'<<flush;

        arma::mat ids_mat(P,N,arma::fill::zeros);
        arma::vec group_sizes(P);
        arma::mat group_sizes_xLayer(P,4,arma::fill::zeros);
        for(i=0;i<N;++i) ids_mat(ids[i],i)=1;

		// filling group sizes per layer
		for(k=0;k<4;k++){
			for(i=0;i<N;++i) if(stop[i]<=k) group_sizes_xLayer(ids[i],k)+=1;
		}
        
		group_sizes=arma::sum(ids_mat,1);

        arma::mat tmp_mat(4,P,arma::fill::zeros);
        for(i=0;i<N;++i){
            for(k=stop[i];k<4;k++){
                tmp_mat(k,ids[i])+=s(i,k);
            }
        }

        for(i=0;i<N;++i){
			// current class;
			int mu0 = ids[i];

			// remove lineage from its class
            group_sizes(ids[i])+=-1;
			for(k=stop[i];k<4;k++){
				tmp_mat(k,ids[i])-=s(i,k);
				group_sizes_xLayer(ids[i],k)-=1;
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

			log_alpha_MH=0;
			if(mu_new<P){
				for(k=stop[i];k<4;k++){
					log_alpha_MH+= gsl_sf_lnbeta(tmp_mat(k,mu_new)+s(i,k)+alpha_p,M*(group_sizes_xLayer(mu_new,k)+1)-tmp_mat(k,mu_new)-s(i,k)+beta_p)
						          +gsl_sf_lnbeta(tmp_mat(k,mu0)+alpha_p,M*group_sizes_xLayer(mu0,k)-tmp_mat(k,mu0)+beta_p)
								  -gsl_sf_lnbeta(tmp_mat(k,mu_new)+alpha_p,M*(group_sizes_xLayer(mu_new,k))-tmp_mat(k,mu_new)+beta_p)
								  -gsl_sf_lnbeta(tmp_mat(k,mu0)+s(i,k)+alpha_p,M*(group_sizes_xLayer(mu0,k)+1)-tmp_mat(k,mu0)-s(i,k)+beta_p);
                } // end for over layers

            } else if(mu_new==P){
				for(k=stop[i];k<4;k++){
                    log_alpha_MH+= gsl_sf_lnbeta(s(i,k)+alpha_p,M-s(i,k)+beta_p)
                                  +gsl_sf_lnbeta(tmp_mat(k,mu0)+alpha_p,M*group_sizes_xLayer(mu0,k)-tmp_mat(k,mu0)+beta_p)
                                  -gsl_sf_lnbeta(alpha_p,beta_p)
                                  -gsl_sf_lnbeta(tmp_mat(k,mu0)+s(i,k)+alpha_p,M*(group_sizes_xLayer(mu0,k)+1)-tmp_mat(k,mu0)-s(i,k)+beta_p);
                } // end for over layers
			}


			if(gsl_rng_uniform(r)<exp(log_alpha_MH)){
				ids[i]=mu_new;
				if(mu_new==P){
					tmp_mat.insert_cols(P,1);
					group_sizes.insert_rows(P,1);
					group_sizes_xLayer.insert_rows(P,1);
					P++;
				}
			}

            group_sizes(ids[i])++;
            for(k=stop[i];k<4;k++){
			   tmp_mat(k,ids[i])+=s(i,k);
			   group_sizes_xLayer(ids[i],k)+=1;
			}

			if(group_sizes(mu0)==0){
				P--;
				if(verbose) cout<<"delete "<<mu0<<endl;
                for(unsigned int l=0;l<N;++l){
                    if(ids[l]>mu0) ids[l]--;
                }

                tmp_mat.shed_col(mu0);
                group_sizes.shed_row(mu0);
                group_sizes_xLayer.shed_row(mu0);

            }

			
        } // end loop over lineages

        arma::mat pmuk(4,P,arma::fill::zeros);
        arma::vec conf(N,arma::fill::zeros);
        arma::vec tmp_conf(4);

        for(k=0;k<4;k++){
            for(mu=0;mu<P;mu++){
                pmuk(k,mu)=gsl_ran_beta(r,tmp_mat(k,mu),M*group_sizes_xLayer(mu,k)-tmp_mat(k,mu));
            }
        }

				
        if(sample%TRIM==0 && sample>BURN_IN){
            double F=0;
            for(nu=0;nu<P;nu++){
                F+=0;
                for(k=0;k<4;++k) {
                    F+=gsl_sf_lnbeta(tmp_mat(k,nu)+alpha_p,M*group_sizes_xLayer(nu,k)-tmp_mat(k,nu)+beta_p);
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

        
    }

    membershipTraj.close();

    return 0;
}
