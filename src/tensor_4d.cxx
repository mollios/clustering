#include<iostream>
#include<cstdlib>
#include<cmath>
#include"tensor_4d.h"

Tensor4d::Tensor4d(int depths, int heights, int rows, int cols) try :
  Depths(depths), Element(new Tensor[depths]){
  for(int i=0;i<Depths;i++){
    Element[i]=Tensor(heights, rows, cols);
  }
}
catch(std::bad_alloc){
    std::cerr << "Out of Memory" << std::endl;
    throw;
 }

Tensor4d::~Tensor4d(void){
  delete []Element;
}

Tensor4d::Tensor4d(const Tensor4d &arg) try :
  Depths(arg.Depths), Element(new Tensor[Depths]){
  for(int i=0;i<Depths;i++){
    Element[i]=arg.Element[i];
  }
}
catch(std::bad_alloc){
    std::cerr << "Out of Memory" << std::endl;
    throw;
 }

Tensor4d::Tensor4d(Tensor4d &&arg)
  : Depths(arg.Depths), Element(arg.Element){
  arg.Depths=0;
  arg.Element=nullptr;
}

Tensor4d &Tensor4d::operator=(Tensor4d &&arg){
  if(this==&arg){
    return *this;
  }
  else{
    Depths=arg.Depths;
    Element=arg.Element;
    arg.Depths=0;
    arg.Element=nullptr;
    return *this;
  }
}

Tensor4d &Tensor4d::operator=(const Tensor4d &arg){
  if(this==&arg)	return *this;
  //Rows=arg.Rows;ここではRowsを更新してはいけない
  if(this->Depths != arg.Depths ||this->heights() != arg.heights() ||
     this->cols() != arg.cols()|| this->rows() != arg.rows()){
    Depths=arg.Depths;
    delete []Element;
    try{
      Element=new Tensor[Depths];
    }
    catch(std::bad_alloc){
      std::cerr << "Out of Memory" << std::endl;
      throw;
    }
  }
  for(int i=0;i<Depths;i++){
    Element[i]=arg.Element[i];
  }
  return *this;
}

int Tensor4d::rows(void) const{
  return Element[0].rows();
}

int Tensor4d::cols(void) const{
  return Element[0].cols();
}

int Tensor4d::heights(void) const{
  return Element[0].heights();
}

int Tensor4d::depths(void) const{
  return Depths;
}
Tensor Tensor4d::operator[](int index) const{
  return Element[index];
}

Tensor &Tensor4d::operator[](int index){
  return Element[index];
}




