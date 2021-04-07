#include "vector.h"
#include "matrix.h"

/***
活用指針
求めたいものが表面的に逆行列であっても、
１．その後に逆行列とベクトルとの積を計算するのであれば、
線形方程式を解く
２．係数行列が同じで右辺のベクトルが異なる線形方程式を複数解くのであれば、
行列方程式を解く。
３．行列式だけを求めるのと特定の連立方程式を併せて解くのは
計算量的にほぼ同等
 ***/

struct vecWithDet{
  Vector solution;
  double det;
};

struct matWithDet{
  Matrix solution;
  double det;
};

struct scalarWithDet{
  double solution;
  double det;
};
  

Matrix make_LU(const Matrix &matrix, int *p);
Matrix make_LU_det(const Matrix &matrix, int *p, double &det);
Vector solve(const Matrix &matrix, const Vector &vector);
Matrix solve(const Matrix &matrix, const Matrix &matrix2);
struct vecWithDet solveWithDet(const Matrix &matrix, const Vector &vector);
struct matWithDet solveWithDet(const Matrix &matrix, const Matrix &matrix2);
struct matWithDet invWithDet(const Matrix &matrix);
Matrix inv(const Matrix &matrix);
double det(const Matrix &matrix);
double inverseWeightedQuadraticForm(const Matrix &matrix, const Vector &vector);
Vector inverseWeightedQuadraticForm(const Matrix &matrix, const Matrix &matrix2);
struct scalarWithDet inverseWeightedQuadraticFormWithDet(const Matrix &matrix, const Vector &vector);
struct vecWithDet inverseWeightedQuadraticFormWithDet(const Matrix &matrix, const Matrix &matrix2);
