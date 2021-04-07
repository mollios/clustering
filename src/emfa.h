#include"efcma.h"
#include"mfa.h"

#ifndef __EMFA__
#define __EMFA__

class Emfa: virtual public Efcma, public Mfa{
 public:
  Emfa(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const int &sub_dimension,
       const double &fuzzifierLambda);
  virtual void revise_centers(void);
};

#endif
