#include"efcma.h"

Efcma::Efcma(const int &dimension,
             const int &data_number,
             const int &centers_number,
             const double &fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Efcm(dimension, data_number, centers_number, fuzzifierLambda){
}


void Efcma::revise_membership(void){
  Tmp_Membership=Membership;
  for(int k=0;k<data_number();k++){
    int numZeroDissimilarities=0;
    Vector indexZeroDissimilarities(centers_number(), 0.0, "all");
    for(int i=0;i<centers_number();i++){
      double denominator=0.0;
      for(int j=0;j<centers_number();j++){
        denominator+=Clusters_size[j]/Clusters_size[i]
          *exp(-FuzzifierLambda*(Dissimilarities[j][k]
                                 -Dissimilarities[i][k]));
      }
      Membership[i][k]=1.0/denominator;
    }
  }
  return;
}

void Efcma::revise_centers(void){
  Efcm::revise_centers();
  return;
}

void Efcma::revise_clusters_size(void){
  Tmp_Clusters_size=Clusters_size;
  for(int i=0;i<centers_number();i++){
    double numerator=0.0;
    for(int k=0;k<data_number();k++){
      numerator+=Membership[i][k];
    }
    Clusters_size[i]=numerator/data_number();
  }
  return;
}


