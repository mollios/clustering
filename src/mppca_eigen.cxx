#include"mppca_eigen.h"
#include"libLU.h"
#include"eigen.h"

Mppca::Mppca(const int &dimension,
             const int &data_number,
             const int &centers_number,
             const int &sub_dimension) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Basis(centers_number, dimension, sub_dimension),
  Sigma(centers_number),
  M(centers_number, sub_dimension, sub_dimension){
}

void Mppca::revise_dissimilarities(void){
  for(int i=0;i<centers_number();i++){
    Matrix I(Basis[i].rows(), Basis[i].rows());
    Matrix inv(Basis[i].rows(), Basis[i].rows());
    for(int l=0;l<Basis[i].rows();l++){
      for(int l2=0;l2<Basis[i].rows();l2++){
        if(l==l2)
          I[l][l2]=1.0;
        else
          I[l][l2]=0.0;
        
      }
    }
    for(int k=0;k<data_number();k++){
      double alpha=0.5*dimension()*log(2.0*M_PI);
      Matrix C=Basis[i]*transpose(Basis[i])+Sigma[i]*I;
      struct matWithDet Det=invWithDet(C);
#ifdef DISS     
      std::cout << "det(A):" << Det.det << std::endl;
#endif
      double beta=0.5*log(Det.det);
      inv=(1.0/Sigma[i])*I-(1.0/Sigma[i])*Basis[i]*M[i]*transpose(Basis[i]);
      double gamma=0.5*(Data[k]-Centers[i])*(inv*(Data[k]-Centers[i]));
#ifdef DISS
      std::cout << "alpha\n" << alpha << std::endl;
      std::cout << "beta\n" << beta << std::endl;
      std::cout << "gamma\n" << gamma << std::endl;
#endif
      Dissimilarities[i][k]=alpha+beta+gamma;
    }
  }
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      if(Dissimilarities[i][k]<0){
#ifdef DISS
        std::cout << "d:\n" << Dissimilarities << std::endl;
#endif
        Dissimilarities[i][k]=0.0;
      }
    }
  }
  return;
};

void Mppca::revise_basis(void){
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
    Matrix temp=Membership[i][0]*((Data[0]-Centers[i])
                                  *transpose(Data[0]-Centers[i]));
    for(int k=1;k<data_number();k++){
      temp+=Membership[i][k]*(Data[k]-Centers[i])
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

void Mppca::revise_m(void){
  Tmp_M=M;
  for(int i=0;i<centers_number();i++){
    Matrix I(Basis[i].cols(), Basis[i].cols());
    Matrix inv(Basis[i].cols(), Basis[i].cols());
    for(int l=0;l<Basis[i].cols();l++){
      for(int l2=0;l2<Basis[i].cols();l2++){
        if(l==l2)
          I[l][l2]=1.0;
        else
          I[l][l2]=0.0;
      }
    }
    inv=transpose(Basis[i])*Basis[i]+Sigma[i]*I;
    std::cout << "inv\n" << inv << std::endl;
    for(int l=0;l<Basis[i].cols();l++){
      M[i][l]=solve(inv,I[l]);
    }
  }
  return;
}

double &Mppca::basis(const int &index1, const int &index2, const int &index3){ 
  return Basis[index1][index2][index3];
}

Matrix &Mppca::basis(const int &index){
  return Basis[index];
}

Matrix &Mppca::tmp_basis(const int &index){
  return Tmp_Basis[index];
}

Vector Mppca::sigma(void){
  return Sigma;
}

Vector Mppca::tmp_sigma(void){
  return Tmp_Sigma;
}

double &Mppca::sigma(const int &index1){ 
  return Sigma[index1];
}

double &Mppca::m(const int &index1, const int &index2, const int &index3){ 
  return M[index1][index2][index3];
}

Matrix &Mppca::m(const int &index){
  return M[index];
}

Matrix &Mppca::tmp_m(const int &index){
  return Tmp_M[index];
}
