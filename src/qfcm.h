#include"efcm.h"
#include"sfcm.h"

#ifndef __QFCM__
#define __QFCM__

class Qfcm: public Efcm, public Sfcm{

public:
  Qfcm(int dimension,
       int data_number,
       int centers_number,
       double fuzzifierEm,
       double fuzzifierLambda);
  virtual void revise_membership(void);
  virtual void revise_centers(void);
};

#endif
