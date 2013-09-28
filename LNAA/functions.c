#include "headers.h"
#include "functions.h"
#include "print.h"


/*Logistic Growth*/
/*r_o_l tuning paramater is a problem*/
/*trun*/
double trun_const_low(double x,double m, double l)
{
  return(1-0.5-(  (pow(3*l,0.5))*(x-m)/(((l*pow(x-m,2))+3)) +atan((pow(l/3,0.5))*(x-m)))/M_PI);
}


double H_func(double x, double xm1, double zeta,double K, double r,double P,double upsilon){         
  double a=exp(r),b=exp(r)/exp(K),output;   
  output=
(a*exp(P)*exp(a*x))/(b*exp(P)*(exp(a*x)-1)+a)
-(a*exp(P)*exp(a*xm1))/(b*exp(P)*(exp(a*xm1)-1)+a)*
(   
exp(a*(x-xm1))*
pow(
(b*exp(P)*(exp(a*xm1)-1)+a)/
(b*exp(P)*(exp(a*x)-1)+a)
,2)
);

  return(output);                                                                                                                      
}

double H_beta_func(double x, double xm1, double zeta,double K, double r,double P,double upsilon){
  double a=exp(r),b=exp(r)/exp(K),output;
  output=
    exp(a*(x-xm1))*
pow(
(b*exp(P)*(exp(a*xm1)-1)+a)/
(b*exp(P)*(exp(a*x)-1)+a)
,2);

  return(output);                                                                
}




double zeta_const(double x, double xm1, double zeta,double K, double r,double P,double upsilon){                                                         
  double a=exp(r),b=exp(r)/exp(K),var;
  var=0.5*exp(-zeta)*a*exp(2*P)*exp(2*a*x)*
    pow(b*exp(P)*(exp(a*x)-1)+a,-4)*
    (
     b*b*exp(2*P)*(exp(2*a*x)-exp(2*a*xm1))+
     4*b*exp(P)*(a-b*exp(P))*(exp(a*x)-exp(a*xm1))+
     2*a*(x-xm1)*pow(a-b*exp(P),2)
     ) ;
  return(var);
}

/*log of normal density (log precision paramater)*/
double log_normdens_E(double yy,double y,double precision)
{double F;
  F=yy-y;
  return(0.5*(precision-F*F*exp(precision)));
}
/*log of normal density*/

double log_normdens(double yy,double y, double precision)
{double F;
  F=yy-y;
  return(0.5*(log(precision)-F*F*precision));
}

double logistic_function(double t,double K, double r, double P){
	double output;
	output=(K*P*exp(r*t))/(K+P*(exp(r*t)-1));
	return(output);
}

double logistic_function_E(double t,double K, double r, double P){
	double output;
	output=logistic_function(t,exp(K), exp(r), exp(P));
	return(output);
}

/*Gibbs*/

double gauss_sample(gsl_rng *RNG, struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
	vec=(tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau)+gsl_ran_gaussian(RNG,1/sqrt(tau_0+Ndouble*tau));
	return(vec);
}

/*MCMC MH*/

double MCMC_base(gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
	double logu,logaprob,can;
	can=para+gsl_ran_gaussian(RNG,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

	logu=log(1-gsl_rng_uniform(RNG));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}


double MCMC_base_truncate_low(double truncate,gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
  double logu,logaprob,can;
    can=para+gsl_ran_gaussian(RNG,*h);
    if(can<=(truncate)){
      can=para;
    }
  logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

  logu=log(1-gsl_rng_uniform(RNG));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);

}


double MCMC_base_truncate_high(double truncate,gsl_rng *RNG, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
  double logu,logaprob,can;
    can=para+gsl_ran_gaussian(RNG,*h);
    if(can>=(truncate)){
      can=para;
    }

  logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

  logu=log(1-gsl_rng_uniform(RNG));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}

double MCMC_K_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para, int l, int m){
  double SUM=0,M=exp(D_para->P),C=0,H,V,Z,Hb;
  int n,nn,mm;
  mm=D->NoSUM[l]+m;
  nn=l*D->M*D->N + m*D->N + 0;
  Z=zeta_const(D->x[nn],0,D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);                                                 
  H=H_func(D->x[nn],0,D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);
  Hb=H_beta_func(D->x[nn],0,D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);
  C=C*Hb*Hb;
  SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
  V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
  M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
  V=V*(C+Z);
  C=C+Z-V;

  for (n=1;n<(D->NoTIME[mm]);n++){
    nn=l*D->M*D->N + m*D->N + n;
    Z=zeta_const(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);
    H=H_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);
    Hb=H_beta_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],para,D_para->r_lm[mm],D_para->P,1);
    C=C*Hb*Hb;
    SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
    V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
    M=H+Hb*M+V*( (D->y[nn])-(H+Hb*M) );
    V=V*(C+Z);
    C=C+Z-V;

  }
    SUM+=-2*log_normdens(para,log(0.1),2);
  return(-0.5*SUM);
}



double MCMC_r_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double SUM=0,M=exp(D_para->P),C=0,H,V,Z,Hb;
  int n,nn,mm;
  mm=D->NoSUM[l]+m;
  nn=l*D->M*D->N + m*D->N + 0;
  Z=zeta_const(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);       
  H=H_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);
  Hb=H_beta_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);
  C=C*Hb*Hb;
  SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
  V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
  M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
  V=V*(C+Z);
  C=C+Z-V;
  for (n=1;n<(D->NoTIME[mm]);n++){
    nn=l*D->M*D->N + m*D->N + n;
    Z=zeta_const(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);       
    H=H_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);
    Hb=H_beta_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],para,D_para->P,1);
    C=C*Hb*Hb;
    SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
    V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
    M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
    V=V*(C+Z);
    C=C+Z-V;
  }
   SUM+=-2*log_normdens(para,log(3),5);
  return(-0.5*SUM);
}


double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para, int l, int m){
  double SUM=0,M=exp(D_para->P),C=0,H,V,Z,Hb;
  int n,nn,mm;
  for (m=0;m<D->NoORF[l];m++){
    M=exp(D_para->P);C=0;
   mm=D->NoSUM[l]+m;
  nn=l*D->M*D->N + m*D->N + 0;
  Z=zeta_const(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);       

  H=H_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);
  Hb=H_beta_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);
  C=C*Hb*Hb;
  SUM+=log(2*M_PI*(exp(-para)+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-para)+C+Z));
  V=(C+Z) / (C+Z+exp(-para));
  M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
  V=V*(C+Z);
  C=C+Z-V;

    for (n=1;n<(D->NoTIME[mm]);n++){
    nn=l*D->M*D->N + m*D->N + n;
    Z=zeta_const(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);       

    H=H_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);
    Hb=H_beta_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,D_para->upsilon_lm[mm]);
    C=C*Hb*Hb;
    SUM+=log(2*M_PI*(exp(-para)+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-para)+C+Z));
    V=(C+Z) / (C+Z+exp(-para));
    M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
    V=V*(C+Z);
    C=C+Z-V;
    }
  }
  SUM+=-2*log_normdens(para,log(1/(0.01*0.01)),0.1);
  return(-0.5*SUM);
}


double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double SUM=0,M=exp(para),C=0,H,V,Z,Hb;
  int n,nn,mm;
  for (l=0;l<D->L;l++){                                                                                                                                     
    for (m=0;m<D->NoORF[l];m++){  
      M=exp(para);C=0;
  mm=D->NoSUM[l]+m;
  nn=l*D->M*D->N + m*D->N + 0;
  Z=zeta_const(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);       

  H=H_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);
  Hb=H_beta_func(D->x[nn],0,D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);
  C=C*Hb*Hb;

  SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));

  V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
  M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
  V=V*(C+Z);
  C=C+Z-V;
  for (n=1;n<(D->NoTIME[mm]);n++){
    nn=l*D->M*D->N + m*D->N + n;
    Z=zeta_const(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);       

    H=H_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);
    Hb=H_beta_func(D->x[nn],D->x[nn-1],D_para->zeta_lm[mm],D_para->K_lm[mm],D_para->r_lm[mm],para,1);
    C=C*Hb*Hb;

    SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
    V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
    M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
    V=V*(C+Z);
    C=C+Z-V;
    }
    }
  }
   SUM+=-2*log_normdens(para,log(0.0001),0.1);

  return(-0.5*SUM);
}




double MCMC_zeta_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double SUM=0,M=exp(D_para->P),C=0,H,V,Z,/*F*/Hb;
  int n,nn,mm;
    mm=D->NoSUM[l]+m;
  nn=l*D->M*D->N + m*D->N + 0;
  Z=zeta_const(D->x[nn],0,para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
  H=H_func(D->x[nn],0,para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
  Hb=H_beta_func(D->x[nn],0,para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
  C=C*Hb*Hb;

  SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
  V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
  M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
  V=V*(C+Z);
  C=C+Z-V;
   for (n=1;n<(D->NoTIME[mm]);n++){
    nn=l*D->M*D->N + m*D->N + n;
    Z=zeta_const(D->x[nn],D->x[nn-1],para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
 
    H=H_func(D->x[nn],D->x[nn-1],para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
    Hb=H_beta_func(D->x[nn],D->x[nn-1],para,D_para->K_lm[mm],D_para->r_lm[mm],D_para->P,1);
    C=C*Hb*Hb;

    SUM+=log(2*M_PI*(exp(-D_para->nu_l[l])+C+Z))+pow((D->y[nn])-(H+Hb*M),2)/((exp(-D_para->nu_l[l])+C+Z));
    V=(C+Z) / (C+Z+exp(-D_para->nu_l[l]));
    M=H+Hb*M+V*((D->y[nn])-(H+Hb*M) );
    V=V*(C+Z);
    C=C+Z-V;
  }
   SUM+=-2*log_normdens(para,log(1/(0.1*0.1)),0.1);
  return(-0.5*SUM);
}


/*Gibbs and MH steps*/



int gibbsandMHloop(int iter,int thin,gsl_rng *RNG,struct_data *D,struct_para *D_para,struct_priors *D_priors ,struct_MH *D_MH,int CAPL,int print){
  int i,j,l,m,mm;
  /*    print=3;*/
    D->NoORF[0]=1;
  D->L=gsl_min(D->L,CAPL);
  if (print==0){/*printheader(D);*/  
  }
  for (i=0;i<iter;i++){
    for (j=0;j<thin;j++){
      D_MH->hK=0.3;
      D_MH->hP=0.5;
      D_para->P=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->P,MCMC_P,-999,-999);
      D_MH->hP=1;

      for (l=0;l<D->L;l++){

	D_para->nu_l[l]=MCMC_base(RNG,D,D_para,D_priors,&D_MH->accept_nu,&D_MH->hnu,D_para->nu_l[l],MCMC_nu_l,l,-999);   

	for (m=0;m<D->NoORF[l];m++){ 
	  mm=D->NoSUM[l]+m;
	  D_MH->hK=0.1;
	  D_para->K_lm[mm]=MCMC_base_truncate_high(0,RNG,D,D_para,D_priors,&D_MH->accept_K,&D_MH->hK,D_para->K_lm[mm],MCMC_K_lm,l,m);
	  D_MH->hK=0.1;
	  D_MH->hr=0.1;
	  D_para->r_lm[mm]=MCMC_base_truncate_high(3.5,RNG,D,D_para,D_priors,&D_MH->accept_r,&D_MH->hr,D_para->r_lm[mm],MCMC_r_lm,l,m);
	  D_MH->hr=0.1;
	  D_MH->hP=0.1;
	  D_para->zeta_lm[mm]=MCMC_base_truncate_low(1.3,RNG,D,D_para,D_priors,&D_MH->accept_P,&D_MH->hP,D_para->zeta_lm[mm],MCMC_zeta_lm,l,m);
	  D_MH->hP=0.2;
	}
      }
    }
    if (print==1){/*printdata(D,D_para,D_MH);*/
      printf("%g %g %g %g %g %g\n",exp(D_para->K_lm[0]),exp(D_para->r_lm[0]),exp(D_para->P),exp(-0.5*D_para->nu_l[0]),exp(-0.5*D_para->zeta_lm[0]),D_para->P);

}
   }
  return 0;
}
