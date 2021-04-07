#include"hcma.h"

Hcma::Hcma(const int &dimension,
           const int &data_number,
           const int &centers_number):
  Hcm(dimension, data_number, centers_number),
  Clusters_size(centers_number),
  Tmp_Clusters_size(centers_number){
  for(int i=0;i<centers_number;i++)
    Clusters_size[i]=DBL_MAX;
  }

void Hcma::revise_clusters_size(void){
  Tmp_Clusters_size=Clusters_size;
  return;
}

const Vector Hcma::clusters_size(void) const{
  return Clusters_size;
}

const Vector Hcma::tmp_clusters_size(void) const{
  return Tmp_Clusters_size;
}

double &Hcma::clusters_size(const int &index1){
  return Clusters_size[index1];
}
