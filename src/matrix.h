#include<iostream>
#include<cstring>
#include"vector.h"

#ifndef __MATRIX__
#define __MATRIX__

class Matrix{
 private:
  int Rows;
  Vector *Element;
 public:
  //Matrix(int rows=0);
  Matrix(int rows=0, int cols=0);
  explicit Matrix(int dim, const char *s);
  explicit Matrix(const Vector &arg, const char *s);
  ~Matrix(void);
  Matrix(const Matrix &arg);
  Matrix &operator=(const Matrix &arg);
  Matrix(Matrix &&arg);
  Matrix &operator=(Matrix &&arg);
  int rows(void) const;
  int cols(void) const;
  Vector operator[](int index) const;
  Vector &operator[](int index);
  Matrix operator+(void) const;
  Matrix operator-(void) const;
  Matrix &operator+=(const Matrix &rhs);
  Matrix &operator-=(const Matrix &rhs);
  Matrix &operator*=(double rhs);
  Matrix &operator/=(double rhs);
  Matrix sub(int row_begin,
		   int row_end,
		   int col_begin,
		   int col_end) const;
  void set_sub(int row_begin,
	       int row_end,
	       int col_begin,
	       int col_end,
	       const Matrix &arg);
};

Matrix operator+(const Matrix &lhs, const Matrix &rhs);
Matrix operator-(const Matrix &lhs, const Matrix &rhs);
Matrix operator*(double lhs, const Matrix &rhs);
Vector operator*(const Matrix &lhs, const Vector &rhs);
Matrix operator*(const Matrix &lhs, const Matrix &rhs);
Matrix operator*(const Matrix &lhs, double rhs);
Matrix operator/(const Matrix &lhs, double rhs);
std::ostream &operator<<(std::ostream &os, const Matrix &rhs);
bool operator==(const Matrix &lhs, const Matrix &rhs);
double max_norm(const Matrix &arg);
double frobenius_norm(const Matrix &arg);
Matrix transpose(const Matrix &arg);
Vector diag(const Matrix &arg);
Matrix pow(const Matrix &arg, double power);
Matrix transpose(const Vector &arg);
Matrix operator*(const Vector &lhs, const Matrix &rhs);

#endif
