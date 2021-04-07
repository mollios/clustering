#include"hcma.h"
#include"sfcm.h"

#ifndef __SFCMA__
#define __SFCMA__

class Sfcma: virtual public Hcma, public Sfcm{
 public:
  Sfcma(const int &dimension,
        const int &data_number,
        const int &centers_number,
        const double &fuzzifierEm);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
  virtual void revise_clusters_size(void);
};
#endif
