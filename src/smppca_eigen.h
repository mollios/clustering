#include"sfcma.h"
#include"mppca_eigen.h"

#ifndef __SMPPCA__
#define __SMPPCA__

class Smppca: virtual public Sfcma, public Mppca{
 public:
  Smppca(const int &dimension,
         const int &data_number,
         const int &centers_number,
         const int &sub_dimension,
         const double &fuzzifierEm);
  virtual void revise_basis(void);
};

#endif
