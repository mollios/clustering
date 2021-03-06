#include"qfcma.h"
#include"mppca_eigen.h"

#ifndef __QMPPCA__
#define __QMPPCA__

class Qmppca: virtual public Qfcma, public Mppca{
 public:
  Qmppca(const int &dimension,
         const int &data_number,
         const int &centers_number,
         const int &sub_dimension,
         const double &fuzzifierEm,
         const double &fuzzifierLambda);
  virtual void revise_basis(void);
};

#endif
