#include"emfa.h"
#include"libLU.h"

Emfa::Emfa(const int &dimension,
           const int &data_number,
           const int &centers_number,
           const int &sub_dimension,
           const double &fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Efcma(dimension, data_number, centers_number, fuzzifierLambda),
  Mppca(dimension, data_number, centers_number, sub_dimension),
  Mfa(dimension, data_number, centers_number, sub_dimension){
}

void Emfa::revise_centers(void){
  Mfa::revise_centers();
  return;
}

