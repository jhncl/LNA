#include "headers.h"
#include "datain.h"


/*TEST*/

int testargc(int argc)
{
  if (argc!=5) {
    perror("argc failed");
    exit(EXIT_FAILURE);
  }
  return 0;
}

int testsame(int a,int b)
{
  if (a!=b) {
    perror("data int failed");
    exit(EXIT_FAILURE);
  }
  return 0;
}



/*READ IN*/


int datadouble(char filename[], double datavec[],int length )
{
  int i;
  char number[20];
  double data;
  FILE *file;
  file=fopen(filename, "r");
  i=0;
  if ( file != NULL ){
    fscanf(file, "%s %lf",number,&data);
    while ( fscanf(file, "%s %lf",number,&data)!=-1){
      datavec[i]=data;   
      i++;
    }
  }
  else{perror(filename);}
  testsame(length,i);
  fclose(file);
  return 0;
}

int dataint(char filename[], int datavec[], int length)
{
  int i;
  char number[20];
  double data;
  FILE *file;
  file = fopen(filename, "r");
  i=0;
  if ( file != NULL ){
    fscanf(file, "%s %lf",number,&data);
    while (fscanf(file, "%s %lf",number,&data)!=-1){
      datavec[i]=data;
      i++;
    }
  }
  else{perror(filename);}
  testsame(length,i);
  fclose(file);
  return 0;
}

int dataLMN(char filename[], int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIME)
{
  char number[20];
  double data;
  FILE *file = fopen(filename, "r");
  if ( file != NULL ){
    fscanf(file, "%s %lf",number,&data);
    fscanf(file, "%s %lf",number,&data);
    *datavecL=data;     
    fscanf(file, "%s %lf",number,&data);
    *datavecM=data;
    fscanf(file, "%s %lf",number,&data);
    *datavecN=data;
    fscanf(file, "%s %lf",number,&data);
    *datavecmaxy=data;
    fscanf(file, "%s %lf",number,&data);
    *datavecmaxTIME=data;
  }
  else{perror(filename);}
  fclose(file);
  return 0;
}

/*INZ*/

int inzstruct_MH(struct_MH *MH)
{
  fillMH(MH);
  return 0;
}

int inzstruct_priors(struct_priors *priors)
{
  fillpriors(priors);
  return 0;
}

int inzstruct_data(struct_data *data)
{
  long size;
  int i;
  dataLMN("LMNmaxdata.txt",&data->L,&data->M,&data->N,&data->maxy,&data->maxNoTIME);  
  testsame(data->L*data->M*data->N,data->maxy);

  size=data->L*data->M*data->N; /*input from file*/ 
  data->y=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
  data->x=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
  data->yy=malloc(size*sizeof(double));
  size=data->L;
  data->NoORF=malloc(size*sizeof(double));    /*Cycle with data->L*/
  data->NoSUM=malloc(size*sizeof(double));    /*Cycle with data->L*/

  /*if (data->y==NULL||data->x==NULL||data->NoORF==NULL||data->NoSUM==NULL||data->NoTIME==NULL) {
    perror("malloc failed");
    exit(EXIT_FAILURE);
    }*/

  datadouble("ydata.txt",data->y,data->L*data->M*data->N);
  datadouble("ydata.txt",data->yy,data->L*data->M*data->N);
  for (i=0;i<(data->L*data->M*data->N);i++){
    /*    if(data->y[i]<=0){ data->y[i]=0.00001;} */ /************/
    /*    if(data->x[i]<0){ data->x[i]=0.00001;}*/
  }

  /*  for (i=0;i<data->L*data->M*data->N;i++){
    if(data->yy[i]<=0){data->yy[i]=0.00001;}
  }*/
  datadouble("xdata.txt",data->x,data->L*data->M*data->N);
  for (i=0;i<data->L*data->M*data->N;i++){
    if(data->x[i]<=0){data->x[i]=(data->x[i+1]+data->x[i])/2+0.00001;}/*Can't have time <=0*/
  }

  dataint("NoORFdata.txt",data->NoORF,data->L);

  filldata(data);
  testsame(data->maxNoTIME,data->SHIFTlmn);

  size=data->SHIFTlmn;/*inputfromfile*/
  data->NoTIME=malloc(size*sizeof(double));   /*Cycle with SHIFTlm*/
  dataint("NoTIMEdata.txt",data->NoTIME,data->SHIFTlmn);


  return 0;
}

int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors)
{
  long size;

  size=data->SHIFTlmn;
  para->K_lm=malloc(size*sizeof(double));
  para->r_lm=malloc(size*sizeof(double));
  para->upsilon_lm=malloc(size*sizeof(double));/*!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  para->zeta_lm=malloc(size*sizeof(double));

  size=data->L;
  para->tau_upsilon_l=malloc(size*sizeof(double));
  para->upsilon_o_l=malloc(size*sizeof(double));
  para->tau_K_l=malloc(size*sizeof(double));
  para->tau_r_l=malloc(size*sizeof(double));
  para->K_o_l=malloc(size*sizeof(double));
  para->r_o_l=malloc(size*sizeof(double));
  para->nu_l=malloc(size*sizeof(double));


  fillpara(para,data,priors);
  return 0;
}

/*FILL*/

int fillMH(struct_MH *MH)
{
  MH->hK=0.1;MH->accept_K=0;
  MH->hr=0.1;MH->accept_r=0;
  MH->hnu=0.1;MH->accept_nu=0;
  MH->hP=0.2;MH->accept_P=0;  /*h sd; accept=0*/
  return 0;
}

int filldata(struct_data *D)
{
  int l;

  D->NoSUM[0]=0;
  for (l=1;l<(D->L);l++){
    D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
  }
  D->SHIFTlmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];/*create mnSHIFT*/
  return 0;
}

int fillpara(struct_para *D_para, struct_data *D,struct_priors *D_priors)
{
  int l,m,mm;
  double SUM=0,SUMa=0;
  /*initials*/
  /*K*/

 for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      if(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]<=0){D_para->K_lm[mm]=D_priors->P_mu;SUM+=D_para->K_lm[mm];}
	else{     
	  D_para->K_lm[mm]=log(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]);SUM+=D_para->K_lm[mm];
	}
    }
    D_para->K_o_l[l]=SUM/D->NoORF[l];
    SUM=0;
    SUMa+=D_para->K_o_l[l];
  }
 D_para->K_p=SUMa/D->L;       /*LMean*/

 for (l=0;l<D->L;l++)          {D_para->tau_K_l[l]=1;/*D_priors->tau_K_mu;*/}/*!!!!!!!!!!!!!!!*/                  /*Precision*/
  
  D_para->sigma_K_o=D_priors->eta_K_o;               /*Precision*/
  /*r*/
  for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      D_para->r_lm[mm]=log(4)/*D_priors->r_mu;*/;/**********************/
      D_para->zeta_lm[mm]=8;
      D_para->upsilon_lm[mm]=D_priors->upsilon_mu;
    }
  }                          /*LMean*/

  for (l=0;l<D->L;l++)          {D_para->tau_r_l[l]=D_priors->tau_r_mu;
    D_para->tau_upsilon_l[l]=D_priors->tau_upsilon_mu;
  }                  /*Precision*/

  for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;
    D_para->upsilon_o_l[l]=D_priors->upsilon_mu;
  }        /*LMean*/

  D_para->sigma_r_o=D_priors->eta_r_o;               /*Precision*/
  D_para->sigma_upsilon_o=D_priors->eta_upsilon_o;               /*Precision*/

  D_para->r_p=D_priors->r_mu;       /*LMean*/
  D_para->upsilon_p=D_priors->upsilon_mu;       /*LMean*/

  /*nu*/
  for (l=0;l<D->L;l++)          {D_para->nu_l[l]=6/*SDE D_priors->nu_mu*/;
}                      /*LMean*/


  D_para->sigma_nu=D_priors->eta_nu;   /*Precision for lMean*/
  D_para->sigma_zeta=10;/*D_priors->eta_zeta;*/   /*Precision for lMean*/

  D_para->nu_p=D_priors->nu_mu;   /*LMean*/
  D_para->zeta_p=15;/*D_priors->zeta_mu;*/   /*LMean*/

  /*P*/
  D_para->P=D_priors->P_mu;      /*LMean*/


  D_para->tau_K_p=D_priors->tau_K_mu;
  D_para->sigma_tau_K=D_priors->eta_tau_K;
  D_para->tau_r_p=D_priors->tau_r_mu;
  D_para->sigma_tau_r=D_priors->eta_tau_r;

  D_para->tau_upsilon_p=D_priors->tau_upsilon_mu;
  D_para->sigma_tau_upsilon=D_priors->eta_tau_upsilon;
  return 0;
}

int fillpriors(struct_priors *D_priors)
{
  /*Priors*/
  char number[20];
  double data;
  FILE *file = fopen("priors.txt", "r");
  if ( file != NULL ){  
    fscanf(file, "%s %lf",number,&data);
    fscanf(file, "%s %lf",number,&data);
    /*K*/
    D_priors->tau_K_mu=data;
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_tau_K_p=data;          
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_K_o=data; 
    fscanf(file, "%s %lf",number,&data);   
    D_priors->psi_K_o=data;            
    /*r*/
    fscanf(file, "%s %lf",number,&data);
    D_priors->tau_r_mu=data;              
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_tau_r_p=data;          
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_r_o=data;              
    fscanf(file, "%s %lf",number,&data);
    D_priors->psi_r_o=data;       
    /*nu*/
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_nu=data;         
    fscanf(file, "%s %lf",number,&data); 
    D_priors->psi_nu=data;      
    /*K*//*r*//*nu*//*P*/
    fscanf(file, "%s %lf",number,&data);
    D_priors->K_mu=exp(data);     
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_K_p=data;      /*Normal  LMean; Precisions */
    fscanf(file, "%s %lf",number,&data);
    D_priors->r_mu=exp(data);   
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_r_p=data;      /*Normal  LMean; Precisions */
    fscanf(file, "%s %lf",number,&data);
    D_priors->nu_mu=data;      
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_nu_p=data;     /*Normal  LMean; Precisions */
    fscanf(file, "%s %lf",number,&data);
    D_priors->P_mu=data;   
    fscanf(file, "%s %lf",number,&data);
    D_priors->eta_P=data;   /*Normal  LMean; Precisions */
    D_priors->df=1;

    D_priors->eta_tau_K=D_priors->eta_tau_K_p;  D_priors->psi_tau_K=D_priors ->eta_tau_K_p;
    D_priors->eta_tau_r=D_priors->eta_tau_r_p;  D_priors->psi_tau_r=D_priors->eta_tau_r_p;

    D_priors->eta_zeta=D_priors->eta_nu;  D_priors->psi_zeta=D_priors->psi_nu;
    D_priors->zeta_mu=D_priors->nu_mu;    D_priors->eta_zeta_p=D_priors->eta_nu_p;

    D_priors->eta_upsilon_o=D_priors->eta_r_o; D_priors->psi_upsilon_o=D_priors->psi_r_o;
    D_priors->upsilon_mu=1;/**SDE*/ D_priors->eta_upsilon_p=D_priors->eta_r_p;
    D_priors->tau_upsilon_mu= D_priors->tau_r_mu; D_priors->eta_tau_upsilon=D_priors->eta_tau_r;
    D_priors->eta_tau_upsilon_p=D_priors->eta_tau_r_p; D_priors->psi_tau_upsilon=D_priors->psi_tau_r;
  }
  else{perror("Priors");}  
  fclose(file);
  return 0;
}

/* eof  */
