#include"hcm.h"

#ifndef __SFCM__
#define __SFCM__

class Sfcm: virtual public Hcm{
protected:
  double FuzzifierEm;
public:
  Sfcm(const int &dimension,
       const int &data_number,
       const int &centers_number,
       const double &fuzzifierEm);
  double fuzzifierEm(void)const;
  double &fuzzifierEm(void);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
};

#endif
