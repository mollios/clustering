#include"hcm.h"

#ifndef __EFCM__
#define __EFCM__

class Efcm: virtual public Hcm{
 protected:
  double FuzzifierLambda;
 public:
  Efcm(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const double &fuzzifierLambda);
  double fuzzifierLambda(void)const;
  double &fuzzifierLambda(void);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
};

#endif

