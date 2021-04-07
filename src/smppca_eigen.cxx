#include"smppca_eigen.h"
#include"libLU.h"
#include"eigen.h"

Smppca::Smppca(const int &dimension,
               const int &data_number,
               const int &centers_number,
               const int &sub_dimension,
               const double &fuzzifierEm) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Sfcma(dimension, data_number, centers_number, fuzzifierEm),
  Mppca(dimension, data_number, centers_number, sub_dimension){
}

void Smppca::revise_basis(void){
  Tmp_Basis=Basis;
  Tmp_Sigma=Sigma;
  for(int i=0;i<centers_number();i++){
    Matrix I(Basis[i].cols(), Basis[i].cols());
    Matrix tmp(Basis[i].cols(), Basis[i].cols());
    for(int l=0;l<Basis[i].cols();l++){
      for(int l2=0;l2<Basis[i].cols();l2++){
        if(l==l2)
          I[l][l2]=1;
        else
          I[l][l2]=0;
      }
    }
    Matrix temp=pow(Membership[i][0],FuzzifierEm)*((Data[0]-Centers[i])
                                  *transpose(Data[0]-Centers[i]));
    for(int k=1;k<data_number();k++){
      temp+=pow(Membership[i][k], FuzzifierEm)*(Data[k]-Centers[i])
        *transpose(Data[k]-Centers[i]);
    }
    temp/=data_number()*clusters_size(i);
    
    struct AllEigen result=eigen_jacobi_cyclic_rutishauser(temp, 1.0E-15);
    sort_alleigen(result);
    Matrix transTT=transpose(result.TT);
#ifdef EIGEN_VECTOR
    std::cout << "Lambda\n" << result.Lambda << std::endl;
    std::cout << "transTT\n" << transTT << std::endl;
#endif
    for(int l=0;l<Basis[i].cols();l++){
      for(int l2=0;l2<Basis[i].cols();l2++){
        tmp[l][l2]=result.Lambda[l][l2];
      }
    }
    tmp-=Sigma[i]*I;
    for(int l=0;l<Basis[i].cols();l++){
      for(int l2=0;l2<Basis[i].cols();l2++){
        tmp[l][l2]=pow(tmp[l][l2],1/2);
      }
    }
    tmp=tmp*I;
    for(int ell=0;ell<dimension();ell++){
      for(int l=0;l<Basis[i].cols();l++){
        Basis[i][ell][l]=transTT[ell][l];
      }
    }
    Basis[i]=Basis[i]*tmp;
    Sigma[i]=0.0;
    for(int l=Basis[i].cols();l<dimension();l++){
      Sigma[i]+=result.Lambda[l][l];
    }
    Sigma[i]/=(dimension()-Basis[i].cols());
  }
  return;
};
