#include <iostream>
#include <cstdlib>
#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "libLU.h"

Matrix make_LU(const Matrix &matrix, int *p){
  int i, j, k, d;
  double max_value, max_value_row, tmp;
  int tmp_p;
  Matrix pmatrix(matrix);

  for(i=0;i<matrix.rows();i++){
    p[i]=i;
  }	

  for(i=0;i<matrix.rows()-1;i++){
    d=i;
    max_value=0.0;
    for(j=i;j<matrix.rows();j++){
      max_value_row=0.0;
      for(k=i;k<matrix.cols();k++){
	tmp=fabs(pmatrix[p[j]][k]);
	if(max_value_row<tmp){
	  max_value_row=tmp;
	}
      }
      tmp=fabs(pmatrix[p[j]][i])/max_value_row;
      if(max_value<tmp){
	max_value=tmp;
	d=j;
      }
    }
    if(d!=i){
      tmp_p=p[i];
      p[i]=p[d];
      p[d]=tmp_p;
    }
    if(pmatrix[p[i]][i]==0.0||pmatrix[p[i]][i]==-0.0){
      std::cerr << "Division by 0(pmatrix["
		<< p[i] << "]["
		<< i << "]=" << pmatrix[i][i] << std::endl;
      exit(1);
    }
    for(j=i+1;j<matrix.rows();j++){
      pmatrix[p[j]][i]/=pmatrix[p[i]][i];
      for(k=i+1;k<matrix.cols();k++){
	pmatrix[p[j]][k]-=pmatrix[p[j]][i]*pmatrix[p[i]][k];
      }
    }
  }
  return pmatrix;
}

Vector solve(const Matrix &matrix, const Vector &vector){
  if(matrix.cols()!=vector.size()){
    std::cerr << "solve:Size Unmatched" << std::endl;
    exit(1);
  }
  if(matrix.rows()!=matrix.cols()){
    std::cerr << "solve:Not Implemented" << std::endl;
    exit(1);
  }
  int *p=new int[vector.size()];
  for(int i=0;i<vector.size();i++){
    p[i]=i;
  }
  Matrix LU=make_LU(matrix, p);
  Vector solution(vector.size());
  for(int j=0;j<solution.size();j++){
    solution[j]=vector[p[j]];
    //    std::cout << "#1:" << solution[j] << std::endl;
    for(int k=0;k<j;k++){
      solution[j]-=LU[p[j]][k]*solution[k];
      //      std::cout << "#2:" << solution[j] << std::endl;
    }
  }
  for(int j=solution.size()-1;j>=0;j--){
    for(int k=j+1;k<solution.size();k++){
      solution[j]-=LU[p[j]][k]*solution[k];
      //      std::cout << "#3:" << solution[j] << std::endl;
    }
    solution[j]/=LU[p[j]][j];
    //    std::cout << "#4:" << solution[j] << std::endl;
  }
  delete []p;
  return solution;
}

Matrix solve(const Matrix &matrix, const Matrix &matrix2){
  if(matrix.cols()!=matrix.rows()){
    std::cerr << "solve:Size Unmatched" << std::endl;
    exit(1);
  }
  if(matrix.rows()!=matrix.cols()){
    std::cerr << "solve:Not Implemented" << std::endl;
    exit(1);
  }
  int *p=new int[matrix2.rows()];
  for(int i=0;i<matrix2.rows();i++){
    p[i]=i;
  }
  Matrix LU=make_LU(matrix, p);
  Matrix solution(matrix2.rows(), matrix2.cols());
  for(int i=0;i<solution.cols();i++){
    for(int j=0;j<solution.rows();j++){
      solution[j][i]=matrix2[p[j]][i];
      //      std::cout << "#1:" << solution[j][i] << std::endl;
      for(int k=0;k<j;k++){
	solution[j][i]-=LU[p[j]][k]*solution[k][i];
	//	std::cout << "#2:" << solution[j][i] << std::endl;
      }
    }
    for(int j=solution.rows()-1;j>=0;j--){
      for(int k=j+1;k<solution.rows();k++){
	solution[j][i]-=LU[p[j]][k]*solution[k][i];
	//	std::cout << "#3:" << solution[j][i] << std::endl;
      }
      solution[j][i]/=LU[p[j]][j];
      //      std::cout << "#4:" << solution[j][i] << std::endl;
    }
  }
  delete []p;
  return solution;
}

Matrix make_LU_det(const Matrix &matrix, int *p, double &det){
  int i, j, k, d;
  double max_value, max_value_row, tmp;
  int tmp_p;
  int replacementCount=0;
  Matrix pmatrix(matrix);

  for(i=0;i<matrix.rows();i++){
    p[i]=i;
  }	

  for(i=0;i<matrix.rows()-1;i++){
    d=i;
    max_value=0.0;
    for(j=i;j<matrix.rows();j++){
      max_value_row=0.0;
      for(k=i;k<matrix.cols();k++){
	tmp=fabs(pmatrix[p[j]][k]);
	if(max_value_row<tmp){
	  max_value_row=tmp;
	}
      }
      tmp=fabs(pmatrix[p[j]][i])/max_value_row;
      if(max_value<tmp){
	max_value=tmp;
	d=j;
      }
    }
    if(d!=i){
      tmp_p=p[i];
      p[i]=p[d];
      p[d]=tmp_p;
      replacementCount++;
    }
    if(pmatrix[p[i]][i]==0.0||pmatrix[p[i]][i]==-0.0){
      std::cerr << "Division by 0(pmatrix["
		<< p[i] << "]["
		<< i << "]=" << pmatrix[i][i] << std::endl;
      exit(1);
    }
    for(j=i+1;j<matrix.rows();j++){
      pmatrix[p[j]][i]/=pmatrix[p[i]][i];
      for(k=i+1;k<matrix.cols();k++){
	pmatrix[p[j]][k]-=pmatrix[p[j]][i]*pmatrix[p[i]][k];
      }
    }
  }
  det=1.0;
  for(int i=0;i<matrix.rows();i++){
    det*=pmatrix[p[i]][i];
  }
  if(replacementCount%2==1){
    det=-det;
  }
  return pmatrix;
}

struct vecWithDet solveWithDet(const Matrix &matrix, const Vector &vector){
  struct vecWithDet result;
  if(matrix.cols()!=vector.size()){
    std::cerr << "solveWithDet:Size Unmatched" << std::endl;
    exit(1);
  }
  if(matrix.rows()!=matrix.cols()){
    std::cerr << "solveWithDet:Not Implemented" << std::endl;
    exit(1);
  }
  int *p=new int[vector.size()];
  for(int i=0;i<vector.size();i++){
    p[i]=i;
  }
  Matrix LU=make_LU_det(matrix, p, result.det);
  result.solution=Vector(vector.size());
  for(int j=0;j<result.solution.size();j++){
    result.solution[j]=vector[p[j]];
    for(int k=0;k<j;k++){
      result.solution[j]-=LU[p[j]][k]*result.solution[k];
    }
  }
  for(int j=result.solution.size()-1;j>=0;j--){
    for(int k=j+1;k<result.solution.size();k++){
      result.solution[j]-=LU[p[j]][k]*result.solution[k];
    }
    result.solution[j]/=LU[p[j]][j];
  }
  delete []p;
  return result;
}

struct matWithDet solveWithDet(const Matrix &matrix, const Matrix &matrix2){
  matWithDet result;
  if(matrix.cols()!=matrix.rows()){
    std::cerr << "solveWithDet(const Matrix &, const Matrix):Size Unmatched" << std::endl;
    exit(1);
  }
  if(matrix.rows()!=matrix.cols()){
    std::cerr << "solveWithDet(const Matrix &, const Matrix):Not Implemented" << std::endl;
    exit(1);
  }
  int *p=new int[matrix2.rows()];
  for(int i=0;i<matrix2.rows();i++){
    p[i]=i;
  }
  Matrix LU=make_LU_det(matrix, p, result.det);
  result.solution=Matrix(matrix2.rows(), matrix2.cols());
  for(int i=0;i<result.solution.cols();i++){
    for(int j=0;j<result.solution.rows();j++){
      result.solution[j][i]=matrix2[p[j]][i];
      for(int k=0;k<j;k++){
	result.solution[j][i]-=LU[p[j]][k]*result.solution[k][i];
      }
    }
    for(int j=result.solution.rows()-1;j>=0;j--){
      for(int k=j+1;k<result.solution.rows();k++){
	result.solution[j][i]-=LU[p[j]][k]*result.solution[k][i];
      }
      result.solution[j][i]/=LU[p[j]][j];
    }
  }
  delete []p;
  return result;
}

struct matWithDet invWithDet(const Matrix &matrix){
  return solveWithDet(matrix, Matrix(matrix.rows(), "I"));
}

Matrix inv(const Matrix &matrix){
  return solve(matrix, Matrix(matrix.rows(), "I"));
}

double det(const Matrix &matrix){
  double result;
  int *p=new int[matrix.rows()];
  for(int i=0;i<matrix.rows();i++){
    p[i]=i;
  }
  make_LU_det(matrix, p, result);
  return result;
}

double inverseWeightedQuadraticForm(const Matrix &matrix, const Vector &vector){
  return vector*solve(matrix, vector);
}

Vector inverseWeightedQuadraticForm(const Matrix &matrix, const Matrix &matrix2){
  /***
これはmatrix2*matrix^{-1}*matrix2を計算するものではない！
共通の逆行列に対して複数のベクトルに対する2次形式を
一括で得るためのものである。
単に呼び出し側で簡潔に書けるだけでなく
1回のLU分解で複数の2次形式を計算する。
   ***/
  Matrix tmp=solve(matrix, matrix2);
  Vector result(matrix2.cols(), 0.0, "all");
  for(int i=0;i<result.size();i++){
    for(int j=0;j<matrix.rows();j++){
      result[i]+=matrix2[j][i]*tmp[j][i];
    }
  }
  return result;
}

struct scalarWithDet inverseWeightedQuadraticFormWithDet(const Matrix &matrix, const Vector &vector){
  struct vecWithDet tmp=solveWithDet(matrix, vector);
  struct scalarWithDet result;
  result.solution=vector*tmp.solution;
  result.det=tmp.det;
  return result;
}

struct vecWithDet inverseWeightedQuadraticFormWithDet(const Matrix &matrix, const Matrix &matrix2){
  struct matWithDet tmp=solveWithDet(matrix, matrix2);
  struct vecWithDet result;
  result.solution=Vector(matrix2.cols(), 0.0, "all");
  for(int i=0;i<result.solution.size();i++){
    for(int j=0;j<matrix.rows();j++){
      result.solution[i]+=matrix2[j][i]*tmp.solution[j][i];
    }
  }
  result.det=tmp.det;
  return result;
}
