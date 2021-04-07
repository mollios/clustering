#include<iostream>
#include<cstdlib>
#include<cmath>
#include"vector.h"

Vector::Vector(int size) try :
  Size(size), Element(new double[Size]){
}
catch(std::bad_alloc){
  std::cerr << "Vector::Vector(int size): Out of Memory!" << std::endl;
  throw;
 }

Vector::~Vector(void){
  delete []Element;
}

Vector::Vector(const Vector &arg) try :
  Size(arg.Size), Element(new double[Size]){
  for(int i=0;i<Size;i++){
    Element[i]=arg.Element[i];
  }
}
catch(std::bad_alloc){
  std::cerr << "Vector::Vector(int size): Out of Memory!" << std::endl;
  throw;
 }

Vector::Vector(Vector &&arg)
  : Size(arg.Size), Element(arg.Element){
  arg.Size=0;
  arg.Element=nullptr;
}

Vector::Vector(int dim, double arg, const char *s) try :
  Size(dim), Element(new double[Size]){
  if(strcmp(s, "all")!=0){
    std::cerr << "Invalid string parameter" << std::endl;
    exit(1);
  }
  for(int i=0;i<Size;i++){
    Element[i]=arg;
  }
}
catch(std::bad_alloc){
  std::cerr << "Vector::Vector(int size): Out of Memory!" << std::endl;
  throw;
 }

Vector &Vector::operator=(const Vector &arg){
  if(this==&arg)	return *this;
  if(this->Size != arg.Size){
    Size=arg.Size;
    delete []Element;
    try{
      Element=new double[Size];
    }
    catch(std::bad_alloc){
      std::cerr << "Out of Memory" << std::endl;
      throw;
    }
  }
  for(int i=0;i<Size;i++){
    Element[i]=arg.Element[i];
  }
  return *this;
}

Vector &Vector::operator=(Vector &&arg){
  if(this!=&arg){
    Size=arg.Size;
	delete []Element;
    Element=arg.Element;
    arg.Size=0;
    arg.Element=nullptr;
  }
  return *this;
}

int Vector::size(void) const{
  return Size;
}

double Vector::operator[](int index) const{
  return Element[index];
}

double &Vector::operator[](int index){
  return Element[index];
}

Vector Vector::operator+(void) const{
  return *this;
}

Vector Vector::operator-(void) const{
  Vector result=*this;
  for(int i=0;i<result.Size;i++){
    result[i]*=-1.0;
  }
  return result;
}

Vector &Vector::operator+=(const Vector &rhs){
  if(rhs.Size==0){
    std::cout << "Vector::operator+=:Size 0" << std::endl;
    exit(1);
  }
  else if(Size!=rhs.Size){
    std::cout << "Vector::operator+=:Size Unmatched" << std::endl;
    exit(1);
  }
  else{
    for(int i=0;i<Size;i++){
      Element[i]+=rhs[i];
    }
  }
  return *this;
}

Vector &Vector::operator*=(double rhs){
  for(int i=0;i<Size;i++){
    Element[i]*=rhs;
  }
  return *this;
}

Vector &Vector::operator/=(double rhs){
  for(int i=0;i<Size;i++){
    Element[i]/=rhs;
  }
  return *this;
}


Vector &Vector::operator-=(const Vector &rhs){
  if(rhs.Size==0){
    std::cout << "Vector::operator-=:Size 0" << std::endl;
    exit(1);
  }
  else if(Size!=rhs.Size){
    std::cout << "Vector::operator-=:Size Unmatched" << std::endl;
    exit(1);
  }
  else{
    for(int i=0;i<Size;i++){
      Element[i]-=rhs[i];
    }
  }
  return *this;
}

Vector Vector::operator+(const Vector &rhs) const{
  Vector result=*this;
  return result+=rhs;
}

Vector Vector::operator-(const Vector &rhs) const{
  Vector result=*this;
  return result-=rhs;
}

Vector operator*(double lhs, const Vector &rhs){
  if(rhs.size()==0){
    std::cout << "Vector operator*:Size 0" << std::endl;
    exit(1);
  }
  Vector result=rhs;
  for(int i=0;i<result.size();i++){
    result[i]*=lhs;
  }
  return result;
}

Vector operator/(const Vector &lhs, double rhs){
  if(lhs.size()==0){
    std::cout << "Vector operator/:Size 0" << std::endl;
    exit(1);
  }
  Vector result=lhs;
  return (result/=rhs);
}


std::ostream &operator<<(std::ostream &os, const Vector &rhs){
  os << "(";
  if(rhs.size()>0){
    for(int i=0;;i++){
      os << rhs[i];
      if(i>=rhs.size()-1) break;
      os << ", ";
    }
  }
  os << ')';
  return os;
}

bool Vector::operator==(const Vector &rhs) const{
  if(Size!=rhs.size())	return false;
  for(int i=0;i<Size;i++){
    if(Element[i]!=rhs[i])	return false;
  }
  return true;
}

double max_norm(const Vector &arg){
  if(arg.size()<1){
    std::cout << "Can't calculate norm for 0-sized vector" << std::endl;
    exit(1);
  }
  double result=fabs(arg[0]);
  for(int i=1;i<arg.size();i++){
    double tmp=fabs(arg[i]);
    if(result<tmp)	result=tmp;
  }
  return result;
}

double squared_norm(const Vector &arg){
  return sqrt(norm_square(arg));
}

double norm_square(const Vector &arg){
  double result=0.0;
  for(int i=0;i<arg.size();i++){
    result+=arg[i]*arg[i];
  }
  return result;
}

double L1norm_square(const Vector &arg){
  double result=0.0;
  for(int i=0;i<arg.size();i++){
    result+=fabs(arg[i]);
  }
  return result;
}

double Vector::operator*(const Vector &rhs) const{
  if(Size<1 || rhs.size()<1 || Size!=rhs.size()){
    std::cout << "Can't calculate innerproduct";
    std::cout << "for 0-sized vector";
    std::cout << "or for different sized vector";
    std::cout << std::endl;
    exit(1);
  }
  double result=Element[0]*rhs[0];
  for(int i=1;i<Size;i++){
    result+=Element[i]*rhs[i];
  }
  return result;
}

bool Vector::operator!=(const Vector &rhs) const{
  if(Size!=rhs.size())	return true;
  for(int i=0;i<Size;i++){
    if(Element[i]!=rhs[i])	return true;
  }
  return false;
}

Vector Vector::sub(int begin, int end) const{
  if(end<begin){
    std::cerr << "Vector::sub:invalid parameter" << std::endl;
    exit(1);
  }
  Vector result(end-begin+1);
  for(int i=0;i<result.size();i++){
    result[i]=Element[begin+i];
  }
  return result;
}

void Vector::set_sub(int begin, int end, const Vector &arg){
  if(end<begin){
    std::cerr << "Vector::sub:invalid parameter" << std::endl;
    exit(1);
  }
  if(end-begin+1!=arg.size()){
    std::cerr << "Vector::sub:invalid parameter" << std::endl;
    exit(1);
  }
  for(int i=0;i<arg.size();i++){
    Element[begin+i]=arg[i];
  }
  return;
}

Vector fraction(const Vector &arg){
  Vector result(arg.size());
  for(int i=0;i<result.size();i++){
    result[i]=1.0/arg[i];
  }
  return result;
}
