#include"efcv.h"
#include"eigen.h"

Efcv::Efcv(const int &dimension,
           const int &data_number,
           const int &centers_number,
           const int &sub_dimension,
           const double &fuzzifierLambda) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Efcma(dimension, data_number, centers_number, fuzzifierLambda),
  Basis(centers_number, sub_dimension, dimension){
}

void Efcv::revise_dissimilarities(void){
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Dissimilarities[i][k]=norm_square(Data[k]-Centers[i]);
      for(int l=0;l<Basis[i].rows();l++){
        Dissimilarities[i][k]-=pow(Basis[i][l]*(Data[k]-Centers[i]),2);
      }
    }
  }
  return;
};

void Efcv::revise_basis(void){
  Tmp_Basis=Basis;
  for(int i=0;i<centers_number();i++){
    Matrix temp=Membership[i][0]*((Data[0]-Centers[i])
                                  *transpose(Data[0]-Centers[i]));
    for(int k=1;k<data_number();k++){
      temp+=Membership[i][k]*(Data[k]-Centers[i])
        *transpose(Data[k]-Centers[i]);
    }
    struct AllEigen result=eigen_jacobi_cyclic_rutishauser(temp, 1.0E-15);
    sort_alleigen(result);
    Matrix transTT=transpose(result.TT);
    for(int l=0;l<Basis[i].rows();l++){
      Basis[i][l]=transTT[l];
    }
  }
  return;
};

double &Efcv::basis(const int &index1, const int &index2,
                    const int &index3){ 
  return Basis[index1][index2][index3];
}

Matrix &Efcv::basis(const int &index){
  return Basis[index];
}

Matrix &Efcv::tmp_basis(const int &index){
  return Tmp_Basis[index];
}
