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
	int Peff;
    int verbose=0;

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
            double tmp_prob_ln[P];
            double tmp_prob[P];

            for(k=stop[i];k<4;k++){
				tmp_mat(k,ids[i])-=s(i,k);
				group_sizes_xLayer(ids[i],k)-=1;
			}

            group_sizes(ids[i])+=-1;
			if(verbose>0)	cout<<i<<" was in "<<ids[i]<<endl;
			if(verbose>0) cout<<"group "<<ids[i]<<" has "<<group_sizes_xLayer.row(ids[i])<<endl;
            for(mu=0;mu<P;++mu){
                tmp_prob_ln[mu]=0;

                for(nu=0;nu<P;nu++){
					// The relevant term in the BoldB function deriving from the integration of the
					// dirichlet distribution is given by the product of gamma. When nu is equal to mu
					// we add 1 to the group_sizes.
                    tmp_prob_ln[mu]+=(nu!=mu) ? gsl_sf_lngamma(group_sizes(nu)+1)
                                              : gsl_sf_lngamma(group_sizes(nu)+2);
				    // here we need to use group sizes that depends on k because of the missing information
					// regarding how many progenitors occupy each layer
					for(k=0;k<4;++k) {
						if(stop[i]<=k){
							if(verbose>1){
								cout<<"mu="<<mu<<' '
									<<"nu="<<nu<<' '
									<<"k ="<<k<<' '
									<<"stop="<<stop[i]<<' '
									<<"s(i,k)="<<s(i,k)<<endl
									<<"gsl_sf_lnbeta("<<tmp_mat(k,nu)+1<<","<<M<<"*"<<(group_sizes_xLayer(nu,k))<<"-"<<tmp_mat(k,nu)<<"+"<<1<<")"<<' '<<endl
								    <<"OR"<<endl
									<<"gsl_sf_lnbeta("<<tmp_mat(k,nu)<<"+"<<s(i,k)+1<<","<<M<<"*"<<(group_sizes_xLayer(nu,k)+1)<<"-"<<tmp_mat(k,nu)<<"-"<<s(i,k)+1<<")"<<endl;
							}

							tmp_prob_ln[mu]+=(nu!=mu) ? gsl_sf_lnbeta(tmp_mat(k,nu)+1,M*(group_sizes_xLayer(nu,k))-tmp_mat(k,nu)+1)
										            	: gsl_sf_lnbeta(tmp_mat(k,nu)+s(i,k)+1,M*(group_sizes_xLayer(nu,k)+1)-tmp_mat(k,nu)-s(i,k)+1);
						} else {
							if(verbose>1){
								cout<<"mu="<<mu<<' '
									<<"nu="<<nu<<' '
									<<"k ="<<k<<' '
									<<"stop="<<stop[i]<<endl
									<<"gsl_sf_lnbeta("<<tmp_mat(k,nu)+1<<","<<M*(group_sizes_xLayer(nu,k))-tmp_mat(k,nu)+1<<")"<<' '<<endl;
							}

							tmp_prob_ln[mu]+=gsl_sf_lnbeta(tmp_mat(k,nu)+1,M*(group_sizes_xLayer(nu,k))-tmp_mat(k,nu)+1);
						}
                    } // end for over layers
                } // end loop over nu for calculating the likelihood

            } // end loop over mu 

            for(mu=0;mu<P;++mu) tmp_prob[mu]=exp(tmp_prob_ln[mu]-tmp_prob_ln[0]);

            gen = gsl_ran_discrete_preproc(P, tmp_prob);

            ids[i]=gsl_ran_discrete(r,gen);
            group_sizes(ids[i])++;
            for(k=stop[i];k<4;k++){
			   tmp_mat(k,ids[i])+=s(i,k);
			   group_sizes_xLayer(ids[i],k)+=1;
			}
            gsl_ran_discrete_free(gen);
            gen=NULL;
			
			if(verbose>0)	cout<<i<<" is now in "<<ids[i]<<endl;
			if(verbose>0) cout<<"group "<<ids[i]<<" has "<<group_sizes_xLayer.row(ids[i])<<endl;
			if(verbose>0) cout<<"sum = "<<arma::sum(group_sizes_xLayer.row(ids[i]),1)<<endl;

        } // end loop over lineages

        arma::mat pmuk(4,P,arma::fill::zeros);
        arma::vec conf(N,arma::fill::zeros);
        arma::vec tmp_conf(4);

        for(k=0;k<4;k++){
            for(mu=0;mu<P;mu++){
                pmuk(k,mu)=gsl_ran_beta(r,tmp_mat(k,mu),M*group_sizes(mu)-tmp_mat(k,mu));
            }
        }

				
        if(sample%TRIM==0 && sample>BURN_IN){
            double F=0;
            for(nu=0;nu<P;nu++){
                F+=gsl_sf_lngamma(group_sizes(nu)+1);
                for(k=0;k<4;++k) {
                    F+=gsl_sf_lnbeta(tmp_mat(k,nu)+1,M*group_sizes(nu)-tmp_mat(k,nu)+1);
                       //-gsl_sf_lngamma(N+P)+gsl_sf_lngamma(P);
                }
            }
            for(i=0;i<N;++i) membershipTraj<<ids[i]<<' '; membershipTraj<<endl;
            FreeEnergyOutput<<F<<endl;

            PTraj<<Peff<<endl;

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

        if(!PFIX){
                if(gsl_rng_uniform(r)<(double)(P+3)/(N+P)){
                    P+=1;
					if(P==11) P=10;
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

			Peff=P;
        } else {
			Peff=P;
			for(mu=0;mu<P;mu++){
				if(group_sizes(mu)<3) Peff--;
			}
		}

    }

    membershipTraj.close();

    return 0;
}
