#include "vector.h"
#include "matrix.h"
#include <cfloat>
#include <cmath>

double max_abs(const Vector &arg);

#ifndef EIGEN
#define EIGEN

struct Eigen{
	double value;
	Vector vector;
};

struct AllEigen{
	Matrix Lambda;
	Matrix TT;
};

const struct Eigen eigen_power(const Matrix &arg, const double &stop);
void eigen_jacobi0(const Matrix &arg, const double &stop,
                   Matrix &Lambda, Matrix &Q);
void eigen_jacobi(const Matrix &arg, const double &stop,
                  Matrix &Lambda, Matrix &Q);
void eigen_jacobi_oppositetheta(const Matrix &arg, const double &stop,
                                Matrix &Lambda, Matrix &Q);
void eigen_jacobi_cyclic(const Matrix &arg, const double &stop,
                         Matrix &Lambda, Matrix &Q);
void eigen_jacobi_rutishauser(const Matrix &arg, const double &stop,
                              Matrix &Lambda, Matrix &Q);
void eigen_jacobi_cyclic_rutishauser(const Matrix &arg, const double &stop,
                                     Matrix &Lambda, Matrix &Q);
const AllEigen eigen_jacobi_cyclic_rutishauser(const Matrix &arg,
                                               const double &stop);
void sort_alleigen(struct AllEigen &arg);
void sort_abs_alleigen(struct AllEigen &arg);

#endif

