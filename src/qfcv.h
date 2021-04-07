#include"qfcma.h"
#include"tensor.h"

#ifndef __QFCV__
#define __QFCV__

class Qfcv: virtual public Qfcma{
 protected:
  Tensor Basis, Tmp_Basis;
 public:
  Qfcv(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierEm,
       const double &fuzzifierLambda);
  virtual void revise_dissimilarities(void);
  virtual void revise_basis(void);
  double &basis(const int &index1, const int &index2, const int &index3); 
  Matrix &basis(const int &index); 
  Matrix &tmp_basis(const int &index); 
};

#endif
