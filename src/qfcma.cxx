#include"qfcma.h"

Qfcma::Qfcma(int dimension,
             int data_number,
             int centers_number,
             double fuzzifierEm,
             double fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Sfcma(dimension, data_number, centers_number,fuzzifierEm),
  Efcma(dimension, data_number, centers_number,fuzzifierLambda){
}

void Qfcma::revise_membership(void){
  Tmp_Membership=Membership;
  for(int k=0;k<data_number();k++){
    int numZeroDissimilarities=0;
    Vector indexZeroDissimilarities(centers_number(), 0.0, "all");
    for(int i=0;i<centers_number();i++){
      if(FuzzifierLambda*(1.0-FuzzifierEm)*Dissimilarities[i][k]==1.0){
        numZeroDissimilarities++;
        indexZeroDissimilarities[i]=1.0;
      }
    }
    if(numZeroDissimilarities!=0){
      for(int i=0;i<centers_number();i++){
        Membership[i][k]=indexZeroDissimilarities[i]
          /numZeroDissimilarities;
      }
    }
    else{
      for(int i=0;i<centers_number();i++){
        double denominator=0.0;
        for(int j=0;j<centers_number();j++){
          denominator+=(Clusters_size[j]/Clusters_size[i]*pow(
          (1.0-FuzzifierLambda*(1.0-FuzzifierEm)*Dissimilarities[j][k])
         /(1.0-FuzzifierLambda*(1.0-FuzzifierEm) *Dissimilarities[i][k])
         ,1.0/(1.0-FuzzifierEm)));
        }
        Membership[i][k]=1.0/denominator;
      }
    }//else
  }//k
  return;
}

void Qfcma::revise_centers(void){
  Sfcma::revise_centers();
  return;
}

void Qfcma::revise_clusters_size(void){
  Tmp_Clusters_size=Clusters_size;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    for(int j=0;j<centers_number();j++){
      double tmp1=0.0;
      double tmp2=0.0;
      for(int k=0;k<data_number();k++){
        tmp1+=(pow(Membership[j][k],FuzzifierEm)
               *(1.0-FuzzifierLambda*(1.0-FuzzifierEm)*Dissimilarities[j][k]));
        tmp2+=(pow(Membership[i][k],FuzzifierEm)
               *(1.0-FuzzifierLambda*(1.0-FuzzifierEm)*Dissimilarities[i][k]));
      }
        denominator+=pow(tmp1/tmp2,1.0/FuzzifierEm);
    }
    Clusters_size[i]=1.0/denominator;
  }
  return;
}
