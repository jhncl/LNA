#include "headers.h"
#include "datain.h"
#include "functions.h"
#include "print.h"
     int
     main (int argc,char *argv[])
     {
	struct_data *data= malloc(sizeof(struct_data));
	struct_para *para= malloc(sizeof(struct_para));
	struct_priors *priors= malloc(sizeof(struct_priors));
	struct_MH *MH = malloc(sizeof(struct_MH));

	int burn,iters,thin, CAPL;
	long seed;
	const gsl_rng_type * T;
	gsl_rng * RNG;


	testargc(argc);

	gsl_rng_env_setup ();
	T = gsl_rng_default;
	RNG = gsl_rng_alloc (T);
	seed = time (NULL) * getpid();    
  	gsl_rng_set (RNG, seed); /*seed*/

	burn=atoi(argv[1]);   /*Burn in*/
	iters=atoi(argv[2]);    /*iterations*/
	thin=atoi(argv[3]);        /*thining*/

	CAPL=atoi(argv[4]);        /*CAP D->L*/

        inzstruct_data(data);
	inzstruct_priors(priors);
	inzstruct_para(para,data,priors);

	inzstruct_MH(MH);
	sprintf(data->latentfilename,"latent_variable_%i_%i_%i_%i.txt",burn,iters,thin,CAPL);
	gibbsandMHloop(burn,1,RNG,data,para,priors,MH,CAPL,0);
	gibbsandMHloop(iters,thin,RNG,data,para,priors,MH,CAPL,1);

       	gsl_rng_free(RNG);
	return 0;
}
