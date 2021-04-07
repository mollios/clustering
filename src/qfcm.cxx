#include"qfcm.h"

Qfcm::Qfcm(int dimension,
           int data_number,
           int centers_number,
           double fuzzifierEm,
           double fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Sfcm(dimension, data_number, centers_number,fuzzifierEm),
  Efcm(dimension, data_number, centers_number,fuzzifierLambda){
}

void Qfcm::revise_membership(void){
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
          denominator+=pow((1.0-FuzzifierLambda*(1.0-FuzzifierEm)
                            *Dissimilarities[j][k])
                           /(1.0-FuzzifierLambda*(1.0-FuzzifierEm)
                             *Dissimilarities[i][k])
                           ,1.0/(1.0-FuzzifierEm));
        }
        Membership[i][k]=1.0/denominator;
      }
    }//else
  }//k
  return;
}

void Qfcm::revise_centers(void){
  Sfcm::revise_centers();
  return;
}
