#include"sfcma.h"
#include"mfa.h"

#ifndef __SMFA__
#define __SMFA__

class Smfa: virtual public Sfcma, public Mfa{
 public:
  Smfa(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierEm);
  virtual void revise_centers(void);
  virtual void revise_basis(void);
  virtual void revise_psi(void);
};

#endif
