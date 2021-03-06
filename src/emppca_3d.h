#include"efcma.h"
#include"mppca_3d.h"

#ifndef __EMPPCA__
#define __EMPPCA__

class Emppca: virtual public Efcma, public Mppca{
 public:
  Emppca(const int &dimension,
         const int &data_number,
         const int &centers_number,
         const int &sub_dimension,
         const double &fuzzifierLambda);
  virtual void revise_centers(void);
};

#endif
