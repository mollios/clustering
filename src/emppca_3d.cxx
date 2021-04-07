#include"emppca_3d.h"
#include"libLU.h"

Emppca::Emppca(const int &dimension,
               const int &data_number,
               const int &centers_number,
               const int &sub_dimension,
               const double &fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Efcma(dimension, data_number, centers_number, fuzzifierLambda),
  Mppca(dimension, data_number, centers_number, sub_dimension){
}

void Emppca::revise_centers(void){
  Mppca::revise_centers();
  return;
}
