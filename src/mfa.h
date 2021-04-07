#include"mppca.h"
#include"tensor.h"
#include"tensor_4d.h"

#ifndef __MFA__
#define __MFA__

class Mfa: virtual public Mppca{
 protected:
  Matrix Psi, Tmp_Psi, Inv_Psi;
 public:
  Mfa(const int &dimension,
      const int &data_number,
      const int &centers_number,
      const int &sub_dimension);
  virtual void revise_dissimilarities(void);
  virtual void revise_centers(void);
  virtual void revise_basis(void);
  virtual void revise_psi(void);
  virtual void revise_inv_psi(void);
  virtual void revise_ez(void);
  virtual void revise_ezz(void);
  virtual void revise_m(void);
  double &psi(const int &index1, const int &index2); 
  Matrix &psi(void);
  double &inv_psi(const int &index1, const int &index2); 
  Matrix &inv_psi(void);
  Matrix &tmp_psi(void);  
};

#endif
