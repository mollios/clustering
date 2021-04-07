#include"mppca_3d.h"
#include"libLU.h"

Mppca::Mppca(const int &dimension,
             const int &data_number,
             const int &centers_number,
             const int &sub_dimension) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Basis(centers_number, dimension, sub_dimension),
  Sigma(centers_number),
  Ez(centers_number, data_number, sub_dimension),
  Ezz(centers_number, sub_dimension, sub_dimension),
  invM(centers_number, sub_dimension, sub_dimension){
}

void Mppca::revise_dissimilarities(void){
  double min=0.0;
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      
      double alpha=0.5*Basis[i].cols()*log(2.0*M_PI)
        +0.5*dimension()*log(2.0*M_PI*Sigma[i])+0.5*Basis[i].cols();
      double beta=(1.0/(2.0*Sigma[i]))*norm_square(Data[k]-Centers[i]);
      double gamma=(1.0/Sigma[i])*Ez[i][k]
        *(transpose(Basis[i])*(Data[k]-Centers[i]));
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
      for(int l=0;l<Basis[i].cols();l++){
        inv[l]=solve(invM[i],I[l]);
      }
      Matrix delta=inv*Ez[i][k]*transpose(Ez[i][k]);
      double delta2=0.0;
      //#ifdef DEBUG
      std::cout << "alpha\n" << alpha-1.0 << std::endl;
      std::cout << "beta\n" << beta << std::endl;
      std::cout << "gamma\n" << gamma << std::endl;
      std::cout << "delta\n" << delta << std::endl;
      //#endif
      for(int l=0;l<Basis[i].cols();l++){
        delta2+=(1.0/(2.0*Sigma[i]))*delta[l][l];
      }      
      Dissimilarities[i][k]=alpha+beta-gamma+delta2;
      if(Dissimilarities[i][k]<min){
        if(data_number()!=1){
          min=Dissimilarities[i][k];
        }
        else{
          Dissimilarities[i][k]=0.0;
        }
      }
    }
  }
  // for(int i=0;i<centers_number();i++){
  //   for(int k=0;k<data_number();k++){
  //     if(Dissimilarities[i][k]<0){
  //       if(k==2){
  //         std::cout << "d:\n" << min << std::endl;
  //         std::cout << "d:\n" << Dissimilarities << std::endl;
  //         //exit(1);
  //       }
  //     }
  //   }
  // }
  // for(int i=0;i<centers_number();i++){
  //   for(int k=0;k<data_number();k++){
  //     Dissimilarities[i][k]-=min;
  //   }
  // }
  return;
};

void Mppca::revise_centers(void){
  Tmp_Centers=Centers;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    Vector numerator(dimension(), 0.0, "all");
    for(int k=0;k<data_number();k++){
      denominator+=Membership[i][k];
      numerator+=Membership[i][k]*(Data[k]-Basis[i]*Ez[i][k]);
    }
    Centers[i]=numerator/denominator;
  }
  return;
}
void Mppca::revise_basis(void){
  Tmp_Basis=Basis;
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
    for(int l=0;l<Basis[i].cols();l++){
      inv[l]=solve(Ezz[i],I[l]);
    }
    Matrix temp=Membership[i][0]*(Data[0]-Centers[i])*transpose(Ez[i][0]);
    for(int k=1;k<data_number();k++){
      temp+=Membership[i][k]*(Data[k]-Centers[i])*transpose(Ez[i][k]);
    }
    Basis[i]=temp*inv;
    Matrix U(Basis[i].cols(), Basis[i].rows());
    Matrix TW(Basis[i].cols(), Basis[i].rows());
    TW=transpose(Basis[i]);
    for(int l=0;l<Basis[i].cols();l++){
      U[l]=TW[l]/squared_norm(TW[l]);
      for(int l2=0;l2<l-1;l2++){
        U[l]-=(TW[l]*U[l2])*U[l2];
      }
      U[l]=U[l]/squared_norm(U[l]);
    }
    Basis[i]=transpose(U);
#ifdef DEBUG
    std::cout << "temp\n" << temp << std::endl;
#endif
  }
  return;
};

void Mppca::revise_sigma(void){
  Tmp_Sigma=Sigma;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    double numerator=0.0;
    Matrix temp(Basis[i].cols(), Basis[i].cols());
    for(int k=0;k<data_number();k++){
      denominator+=Membership[i][k];
      double alpha=Membership[i][k]*norm_square(Data[k]-Centers[i]);
      double beta=Membership[i][k]*Ez[i][k]
        *(transpose(Basis[i])*(Data[k]-Centers[i]));
      numerator+=(alpha-2.0*beta);
    }
    temp=Ezz[i]*transpose(Basis[i])*Basis[i];
    for(int l=0;l<Basis[i].cols();l++){
      numerator+=temp[l][l];
    }
#ifdef DEBUG
    std::cout << "denominator\n" << denominator << std::endl;
    std::cout << "numerator\n" << numerator << std::endl;
#endif
    Sigma[i]=numerator/(dimension()*denominator);
  }
  return;
}

void Mppca::revise_ez(void){
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Ez[i][k]=invM[i]*transpose(Basis[i])*(Data[k]-Centers[i]);
    }
  }
  return;
}

void Mppca::revise_ezz(void){
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Ezz[i]+=Membership[i][k]*(Sigma[i]*invM[i]+Ez[i][k]*transpose(Ez[i][k]));
    }
  }
  return;
}

void Mppca::revise_invm(void){
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
#ifdef DEBUG
    std::cout <<"M\n"<< inv << std::endl;
#endif
    for(int l=0;l<Basis[i].cols();l++){
      invM[i][l]=solve(inv,I[l]);
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

double &Mppca::ez(const int &index1, const int &index2, const int &index3){ 
  return Ez[index1][index2][index3];
}

Matrix &Mppca::ez(const int &index){
  return Ez[index];
}

double &Mppca::ezz(const int &index1, const int &index2, const int &index3){ 
  return Ezz[index1][index2][index3];
}

Matrix &Mppca::ezz(const int &index){
  return Ezz[index];
}

double &Mppca::invm(const int &index1, const int &index2, const int &index3){ 
  return invM[index1][index2][index3];
}

Matrix &Mppca::invm(const int &index){
  return invM[index];
}
