#include<iostream>
#include<cstdlib>
#include<cmath>
#include"vector.h"
#include"matrix.h"
#include"tensor.h"

Tensor::Tensor(int heights, int rows, int cols) try :
  Heights(heights), Element(new Matrix[heights]){
  for(int i=0;i<Heights;i++){
    Element[i]=Matrix(rows, cols);
  }
}
catch(std::bad_alloc){
    std::cerr << "Out of Memory" << std::endl;
    throw;
 }

Tensor::~Tensor(void){
  delete []Element;
}

Tensor::Tensor(const Tensor &arg) try :
  Heights(arg.Heights), Element(new Matrix[Heights]){
  for(int i=0;i<Heights;i++){
    Element[i]=arg.Element[i];
  }
}
catch(std::bad_alloc){
    std::cerr << "Out of Memory" << std::endl;
    throw;
 }

Tensor::Tensor(Tensor &&arg)
  : Heights(arg.Heights), Element(arg.Element){
  arg.Heights=0;
  arg.Element=nullptr;
}

Tensor &Tensor::operator=(Tensor &&arg){
  if(this==&arg){
    return *this;
  }
  else{
    Heights=arg.Heights;
    Element=arg.Element;
    arg.Heights=0;
    arg.Element=nullptr;
    return *this;
  }
}

Tensor &Tensor::operator=(const Tensor &arg){
  if(this==&arg)	return *this;
  //Rows=arg.Rows;ここではRowsを更新してはいけない
  if(this->Heights != arg.Heights ||
     this->cols() != arg.cols()|| this->rows() != arg.rows()){
    Heights=arg.Heights;
    delete []Element;
    try{
      Element=new Matrix[Heights];
    }
    catch(std::bad_alloc){
      std::cerr << "Out of Memory" << std::endl;
      throw;
    }
  }
  for(int i=0;i<Heights;i++){
    Element[i]=arg.Element[i];
  }
  return *this;
}

int Tensor::rows(void) const{
  return Element[0].rows();
}

int Tensor::cols(void) const{
  return Element[0].cols();
}

int Tensor::heights(void) const{
  return Heights;
}

Matrix Tensor::operator[](int index) const{
  return Element[index];
}

Matrix &Tensor::operator[](int index){
  return Element[index];
}




