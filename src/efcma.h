#include"hcma.h"
#include"efcm.h"

#ifndef __EFCMA__
#define __EFCMA__

class Efcma: virtual public Hcma, public Efcm{
 public:
  Efcma(const int &dimension,
        const int &data_number,
        const int &centers_number,
        const double &fuzzifierLambda);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
  virtual void revise_clusters_size(void);
};

#endif
