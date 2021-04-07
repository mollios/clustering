#include"sfcm.h"

Sfcm::Sfcm(const int &dimension,
           const int &data_number,
           const int &centers_number,
           const double &fuzzifierEm) : 
  Hcm(dimension, data_number, centers_number), 
  FuzzifierEm(fuzzifierEm){
  }

double Sfcm::fuzzifierEm(void)const{
  return FuzzifierEm;
}

double &Sfcm::fuzzifierEm(void){
  return FuzzifierEm;
}

void Sfcm::revise_membership(void){
  Tmp_Membership=Membership;
  for(int k=0;k<data_number();k++){
    int numZeroDissimilarities=0;
    Vector indexZeroDissimilarities(centers_number(), 0.0, "all");
    for(int i=0;i<centers_number();i++){
      if(Dissimilarities[i][k]==0.0){
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
          denominator+=pow(Dissimilarities[j][k]/Dissimilarities[i][k]
                           , 1.0/(1.0-FuzzifierEm));
        }
        Membership[i][k]=1.0/denominator;
      }
    }//else
  }//k
  return;
}

void Sfcm::revise_centers(void){
  Tmp_Centers=Centers;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    Vector numerator(dimension(), 0.0, "all");
    for(int k=0;k<data_number();k++){
      denominator+=pow(Membership[i][k], FuzzifierEm);
      numerator+=pow(Membership[i][k], FuzzifierEm)*Data[k];
    }
    Centers[i]=numerator/denominator;
  }
  return;
}
