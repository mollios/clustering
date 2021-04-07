#include<iostream>
#include<cstring>
#include"matrix.h"

#ifndef __TENSOR__
#define __TENSOR__

class Tensor{
 private:
  int Heights;
  Matrix *Element;
 public:
  //Tensor(int rows=0);
  Tensor(int heights=0, int rows=0, int cols=0);
  ~Tensor(void);
  Tensor(const Tensor &arg);
  Tensor &operator=(const Tensor &arg);
  Tensor(Tensor &&arg);
  Tensor &operator=(Tensor &&arg);
  int rows(void) const;
  int cols(void) const;
  int heights(void) const;
  Matrix operator[](int index) const;
  Matrix &operator[](int index);
};


#endif
