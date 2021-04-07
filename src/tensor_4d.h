#include<iostream>
#include<cstring>
#include"tensor.h"

#ifndef __TENSOR4D__
#define __TENSOR4D__

class Tensor4d{
 private:
  int Depths;
  Tensor *Element;
 public:
  //Tensor(int rows=0);
  Tensor4d(int depths=0, int heights=0, int rows=0, int cols=0);
  ~Tensor4d(void);
  Tensor4d(const Tensor4d &arg);
  Tensor4d &operator=(const Tensor4d &arg);
  Tensor4d(Tensor4d &&arg);
  Tensor4d &operator=(Tensor4d &&arg);
  int rows(void) const;
  int cols(void) const;
  int heights(void) const;
  int depths(void) const;
  Tensor operator[](int index) const;
  Tensor &operator[](int index);
};


#endif
