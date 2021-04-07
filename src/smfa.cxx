#include"smfa.h"
#include"libLU.h"

Smfa::Smfa(const int &dimension,
           const int &data_number,
           const int &centers_number,
           const int &sub_dimension,
           const double &fuzzifierEm) :
  Hcm(dimension, data_number, centers_number),
  Hcma(dimension, data_number, centers_number),
  Sfcma(dimension, data_number, centers_number, fuzzifierEm),
  Mppca(dimension, data_number, centers_number, sub_dimension),
  Mfa(dimension, data_number, centers_number, sub_dimension){
}

void Smfa::revise_centers(void){
  Tmp_Centers=Centers;
  for(int i=0;i<centers_number();i++){
    double denominator=0.0;
    Vector numerator(dimension(), 0.0, "all");
    for(int k=0;k<data_number();k++){
      denominator+=pow(Membership[i][k], FuzzifierEm);
      numerator+=pow(Membership[i][k], FuzzifierEm)
        *(Data[k]-(Basis[i]*Ez[i][k]));
    }
    Centers[i]=numerator/denominator;
  }
  return;
}

void Smfa::revise_basis(void){
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
    Matrix temp=pow(Membership[i][0], FuzzifierEm)
      *(Data[0]-Centers[i])*transpose(Ez[i][0]);
    Matrix temp2=pow(Membership[i][0], FuzzifierEm)*Ezz[i][0];
    for(int k=1;k<data_number();k++){
      temp+=pow(Membership[i][k], FuzzifierEm)
        *(Data[k]-Centers[i])*transpose(Ez[i][k]);
      temp2+=pow(Membership[i][k], FuzzifierEm)*Ezz[i][k];
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
    // Basis[i]=transpose(U);
#ifdef VERBOSE3
    std::cout << "temp\n" << temp << std::endl;
    std::cout << "temp2\n" << temp2 << std::endl;
#endif
  }
  return;
};

void Smfa::revise_psi(void){
  Tmp_Psi=Psi;
  for(int ell=0;ell<dimension();ell++){
    for(int ell2=0;ell2<dimension();ell2++){
      Psi[ell][ell2]=0.0;
    }
  }
  for(int i=0;i<centers_number();i++){
    Matrix S=pow(Membership[i][0], FuzzifierEm)
      *(Data[0]-Centers[i])*transpose(Data[0]-Centers[i]);
    Matrix beta=pow(Membership[i][0], FuzzifierEm)
      *(Basis[i]*Ez[i][0])*transpose(Data[0]-Centers[i]);
    for(int k=1;k<data_number();k++){
      S+=pow(Membership[i][k], FuzzifierEm)
        *(Data[k]-Centers[i])*transpose(Data[k]-Centers[i]);
      beta+=pow(Membership[i][k], FuzzifierEm)
        *(Basis[i]*Ez[i][k])*transpose(Data[k]-Centers[i]);
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
