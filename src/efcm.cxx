#include"efcm.h"
#include <iomanip>
Efcm::Efcm(const int &dimension,
           const int &data_number,
           const int &centers_number,
           const double &fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  FuzzifierLambda(fuzzifierLambda){
  }

double Efcm::fuzzifierLambda(void)const{
  return FuzzifierLambda;
}


double &Efcm::fuzzifierLambda(void){
  return FuzzifierLambda;
}

void Efcm::revise_membership(void){
  Tmp_Membership=Membership;
  for(int k=0;k<data_number();k++){
    int numZeroDissimilarities=0;
    Vector indexZeroDissimilarities(centers_number(), 0.0, "all");
    for(int i=0;i<centers_number();i++){
      double denominator=0.0;
      for(int j=0;j<centers_number();j++){
        denominator+=exp(-FuzzifierLambda*(Dissimilarities[j][k]
                                           -Dissimilarities[i][k]));
        Membership[i][k]=1.0/denominator;
      }
    }
  }//k
  return;
}

void Efcm::revise_centers(void){
  Tmp_Centers=Centers;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    Vector numerator(dimension(), 0.0, "all");
    for(int k=0;k<data_number();k++){
      denominator+=Membership[i][k];
      numerator+=Membership[i][k]*Data[k];
    }
    Centers[i]=numerator/denominator;
  }
  return;
}
