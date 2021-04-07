#include<iostream>
#include<cstring>

#ifndef __VECTOR__
#define __VECTOR__

class Matrix;

class Vector{
 private:
  int Size;
  double *Element;
 public:
  Vector(int size=0);
  ~Vector(void);
  explicit Vector(int dim, double arg, const char *s);
  Vector(const Vector &arg);
  Vector &operator=(const Vector &arg);
  Vector(Vector &&arg);
  Vector &operator=(Vector &&arg);
  int size(void) const;
  double operator[](int index) const;
  double &operator[](int index);
  Vector operator+(void) const;
  Vector operator-(void) const;
  Vector &operator+=(const Vector &rhs);
  Vector &operator-=(const Vector &rhs);
  Vector &operator*=(double rhs);
  Vector &operator/=(double rhs);
  Vector operator+(const Vector &rhs) const;
  Vector operator-(const Vector &rhs) const;
  double operator*(const Vector &rhs) const;
  bool operator==(const Vector &rhs) const;
  bool operator!=(const Vector &rhs) const;
  Vector sub(int begin, int end) const;
  void set_sub(int begin, int end, const Vector &arg);
};

Vector operator*(double lhs, const Vector &rhs);
Vector operator/(const Vector &lhs, double rhs);
std::ostream &operator<<(std::ostream &os, const Vector &rhs);
double max_norm(const Vector &arg);
double squared_norm(const Vector &arg);
double norm_square(const Vector &arg);
double L1norm_square(const Vector &arg);
Vector fraction(const Vector &arg);

#endif
