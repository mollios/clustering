#include"hcma.h"
#include"tensor.h"

#ifndef __MPPCA__
#define __MPPCA__

class Mppca: virtual public Hcma{
 protected:
  Tensor Basis, Tmp_Basis;
  Vector Sigma, Tmp_Sigma;
  Tensor Ez, Ezz;
  Tensor invM;
 public:
  Mppca(const int &dimension,
        const int &data_number,
        const int &centers_number,
        const int &sub_dimension);
  virtual void revise_dissimilarities(void);
  virtual void revise_centers(void);
  virtual void revise_basis(void);
  virtual void revise_sigma(void);
  virtual void revise_ez(void);
  virtual void revise_ezz(void);
  virtual void revise_invm(void);
  double &basis(const int &index1, const int &index2, const int &index3); 
  Matrix &basis(const int &index);
  Matrix &tmp_basis(const int &index);  
  Vector sigma(void);
  Vector tmp_sigma(void);
  double &sigma(const int &index1); 
  double &ez(const int &index1, const int &index2, const int &index3); 
  Matrix &ez(const int &index); 
  double &ezz(const int &index1, const int &index2, const int &index3); 
  Matrix &ezz(const int &index);
  double &invm(const int &index1, const int &index2, const int &index3); 
  Matrix &invm(const int &index); 
};

#endif
