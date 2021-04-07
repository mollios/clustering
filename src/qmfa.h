#include"qfcma.h"
#include"mfa.h"

#ifndef __QMFA__
#define __QMFA__

class Qmfa: virtual public Qfcma, public Mfa{
 public:
  Qmfa(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierEm,
       const double &fuzzifierLambda);
  virtual void revise_centers(void);
  virtual void revise_basis(void);
  virtual void revise_psi(void);
};

#endif
