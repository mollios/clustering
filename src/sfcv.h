#include"sfcma.h"
#include"tensor.h"

#ifndef __SFCV__
#define __SFCV__

class Sfcv: virtual public Sfcma{
 protected:
  Tensor Basis, Tmp_Basis;
 public:
  Sfcv(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierEm);
  virtual void revise_dissimilarities(void);
  virtual void revise_basis(void);
  double &basis(const int &index1, const int &index2, const int &index3); 
  Matrix &basis(const int &index); 
  Matrix &tmp_basis(const int &index); 
};

#endif
