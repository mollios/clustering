#include"mfa.h"
#include"libLU.h"

Mfa::Mfa(const int &dimension,
         const int &data_number,
         const int &centers_number,
         const int &sub_dimension) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Mppca(dimension, data_number, centers_number, sub_dimension),
  Psi(dimension, dimension),
  Inv_Psi(dimension, dimension){
}

void Mfa::revise_dissimilarities(void){
  double min=0.0;
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      double det_psi=Psi[0][0]*Psi[1][1]-Psi[0][1]*Psi[1][0];      
      double alpha=0.5*log(det_psi);
      double beta=0.5*(Data[k]-Centers[i])
        *(Inv_Psi*(Data[k]-Centers[i]));
      double gamma=Ez[i][k]*(transpose(Basis[i])
        *Inv_Psi*(Data[k]-Centers[i]));
      Matrix delta=transpose(Basis[i])
        *Inv_Psi*Basis[i]*Ezz[i][k];
      double delta2=0.0;
      #ifdef VERBOSE
      std::cout << "alpha\n" << alpha << std::endl;
      std::cout << "beta\n" << beta << std::endl;
      std::cout << "gamma\n" << gamma << std::endl;
      std::cout << "delta\n" << delta << std::endl;
      #endif
      for(int l=0;l<Basis[i].cols();l++){
        delta2+=0.5*delta[l][l];
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
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      if(Dissimilarities[i][k]<0){
        if(k==2){
          //std::cout << "d:\n" << min << std::endl;
          //std::cout << "d:\n" << Dissimilarities << std::endl;
          //exit(1);
        }
      }
    }
  }
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Dissimilarities[i][k]-=min;
      std::cout << "d:\n" << min << std::endl;
    }
  }
  return;
};

void Mfa::revise_centers(void){
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
void Mfa::revise_basis(void){
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
    Matrix temp=Membership[i][0]*(Data[0]-Centers[i])
      *transpose(Ez[i][0]);
    Matrix temp2=Membership[i][0]*Ezz[i][0];
    for(int k=1;k<data_number();k++){
      temp+=Membership[i][k]*(Data[k]-Centers[i])
        *transpose(Ez[i][k]);
      temp2+=Membership[i][k]*Ezz[i][k];
    }
    for(int l=0;l<Basis[i].cols();l++){
      inv[l]=solve(temp2,I[l]);
    }
    Basis[i]=temp*inv;
    // Matrix U(Basis[i].cols(), Basis[i].rows());
    // Matrix TW(Basis[i].cols(), Basis[i].rows());
    // TW=transpose(Basis[i]);
    // for(int l=0;l<Basis[i].cols();l++){
    //   U[l]=TW[l]/squared_norm(TW[l]);
    //   for(int l2=0;l2<l-1;l2++){
    //     U[l]-=(TW[l]*U[l2])*U[l2];
    //   }
    //   U[l]=U[l]/squared_norm(U[l]);
    // }
    //Basis[i]=transpose(U);
#ifdef VERBOSE3
    std::cout << "temp\n" << temp << std::endl;
    std::cout << "temp2\n" << inv << std::endl;
#endif
  }
  return;
};

void Mfa::revise_psi(void){
  Tmp_Psi=Psi;
  for(int ell=0;ell<dimension();ell++){
    for(int ell2=0;ell2<dimension();ell2++){
      Psi[ell][ell2]=0.0;
    }
  }
  for(int i=0;i<centers_number();i++){
    Matrix S=Membership[i][0]*(Data[0]-Centers[i])
      *transpose(Data[0]-Centers[i]);
    Matrix beta=Membership[i][0]*(Basis[i]*Ez[i][0])
      *transpose(Data[0]-Centers[i]);
    for(int k=1;k<data_number();k++){
      S+=Membership[i][k]*(Data[k]-Centers[i])
        *transpose(Data[k]-Centers[i]);
      beta+=Membership[i][k]*(Basis[i]*Ez[i][k])
        *transpose(Data[k]-Centers[i]);
    }
    Psi+=(S-beta)/data_number();
    
    for(int ell=0;ell<dimension();ell++){
      for(int l=0;l<dimension();l++){
        if(l!=ell)
          Psi[ell][l]=0;
      }
    }
  }
  return;
}

void Mfa::revise_inv_psi(void){
  for(int i=0;i<centers_number();i++){
    Matrix I(dimension(), dimension());
    for(int ell=0;ell<dimension();ell++){
      for(int ell2=0;ell2<dimension();ell2++){
        if(ell==ell2)
          I[ell][ell2]=1.0;
        else
          I[ell][ell2]=0.0;
      }
    }
    for(int ell=0;ell<dimension();ell++){
      Inv_Psi[ell]=solve(Psi,I[ell]);
    }
  }
  return;
}

void Mfa::revise_ez(void){
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Ez[i][k]=M[i]*transpose(Basis[i])*Inv_Psi*(Data[k]-Centers[i]);
    }
  }
  return;
}

void Mfa::revise_ezz(void){
  for(int i=0;i<centers_number();i++){
    for(int k=0;k<data_number();k++){
      Ezz[i][k]=M[i]+Ez[i][k]*transpose(Ez[i][k]);
    }
  }
  return;
}

void Mfa::revise_m(void){
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
    inv=transpose(Basis[i])*Inv_Psi*Basis[i]+I;
#ifdef VERBOSE2
    std::cout <<"M\n"<< inv << std::endl;
#endif
    for(int l=0;l<Basis[i].cols();l++){
      M[i][l]=solve(inv,I[l]);
    }
  }
  return;
}

double &Mfa::psi(const int &index1, const int &index2){ 
  return Psi[index1][index2];
}

Matrix &Mfa::psi(void){
  return Psi;
}

double &Mfa::inv_psi(const int &index1, const int &index2){ 
  return Inv_Psi[index1][index2];
}

Matrix &Mfa::inv_psi(void){
  return Inv_Psi;
}

Matrix &Mfa::tmp_psi(void){
  return Tmp_Psi;
}
