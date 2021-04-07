#include"efcma.h"
#include"sfcma.h"

#ifndef __QFCMA__
#define __QFCMA__

class Qfcma: public Efcma, public Sfcma{

public:
  Qfcma(int dimension,
        int data_number,
        int centers_number,
        double fuzzifierEm,
        double fuzzifierLambda);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
  virtual void revise_clusters_size(void);
};

#endif
