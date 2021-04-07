#include"hcm.h"

#ifndef __HCMA__
#define __HCMA__

class Hcma: virtual public Hcm{
 protected:
  Vector Clusters_size, Tmp_Clusters_size;
 public:
  Hcma(const int &dimension,
       const int &data_number,
       const int &centers_number);
  virtual void revise_clusters_size(void);
  const Vector clusters_size(void) const;
  const Vector tmp_clusters_size(void) const;
  double &clusters_size(const int &index1);
};

#endif
