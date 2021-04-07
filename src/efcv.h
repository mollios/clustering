#include"efcma.h"
#include"tensor.h"

#ifndef __EFCV__
#define __EFCV__

class Efcv: virtual public Efcma{
 protected:
  Tensor Basis, Tmp_Basis;
 public:
  Efcv(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierLambda);
  virtual void revise_dissimilarities(void);
  virtual void revise_basis(void);
  double &basis(const int &index1, const int &index2, const int &index3); 
  Matrix &basis(const int &index); 
  Matrix &tmp_basis(const int &index); 
};

#endif
