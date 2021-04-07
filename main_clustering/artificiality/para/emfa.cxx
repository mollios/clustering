#include<iostream>
#include<fstream>
#include<cstdlib>
#include<random>
#include"../../../src/emfa.h"

#define MAX_ITERATES 100000
#define DIFF_FOR_STOP 1.0E-10

const int centers_number=2;
const int parameter_number=3;

int main(void){
  double Lambda=0.72e+0;
  double para_Lambda[3]={1.0e-1, 1.0e+0, 1.0e+1};
  double stridelambda=2.0;
  std::string filenameData("./data/test/debug3.dat");
#ifdef CHECK_ANSWER
  std::string filenameCorrectCrispMembership
    ("./data/test/debug3.correctCrispMembership");
#endif
  
  std::string::size_type filenameDataSlashPosition=
    filenameData.find_last_of("/")+1;
  std::string::size_type filenameDataDotPosition=
    filenameData.find_last_of(".")-filenameDataSlashPosition; 
  if(filenameDataDotPosition==std::string::npos){
    std::cerr << "File:" << filenameData
              << " needs \".\" and filename-extention." << std::endl;
    exit(1);
  }

#ifdef WHILE
  for(int p=0;p<parameter_number;p++){
      Lambda=para_Lambda[p];
#endif
    
  std::ifstream ifs(filenameData);
  if(!ifs){
    std::cerr << "File:" << filenameData
              << " could not open." << std::endl;
    exit(1);
  }
  int data_number, data_dimension, sub_dimension;
  ifs >> data_number;
  ifs >> data_dimension;
  ifs >> sub_dimension;
  
  Emfa test(data_dimension, data_number, centers_number,
            sub_dimension, Lambda);
#ifdef VERBOSE
  std::cout << test.basis(0).rows() << "\t"
            << test.basis(0).cols() << std::endl;
#endif  
  for(int cnt=0;cnt<data_number;cnt++){
    for(int ell=0;ell<data_dimension;ell++){
      ifs >> test.data(cnt, ell);
    }
  }
#ifdef VERBOSE
  std::cout << test.data()  << std::endl;
#endif

  /***Initial Centers Setting***/
  for(int i=0;i<test.centers_number();i++){
    test.clusters_size(i)=1.0/centers_number;
  }
  // for(int i=0;i<test.centers_number();i++){
  //   for(int ell=0;ell<test.dimension();ell++){
  //       test.centers(i,ell)=0.0;
  //   }
  // }
  // for(int i=0;i<test.centers_number();i++){
  //   for(int ell=0;ell<test.dimension();ell++){
  //     if(ell==2)
  //       test.centers(i,ell)=1.0*pow(-1.0,i);
  //     else
  //       test.centers(i,ell)=0.0;
  //   }
  // }
  for(int i=0;i<test.centers_number();i++){
    for(int ell=0;ell<test.dimension();ell++){
      if(ell==0)
        test.centers(i,ell)=2.0*pow(-1.0,i);
      else
        test.centers(i,ell)=0.0;
    }
  }
  std::ifstream ifs_correctCrispMembership(filenameCorrectCrispMembership);
  if(!ifs_correctCrispMembership){
    std::cerr << "File:" << filenameCorrectCrispMembership
              << " could not open." << std::endl;
    exit(1);
  }
  for(int i=0;i<test.centers_number();i++){
    for(int k=0;k<test.data_number();k++){
      ifs_correctCrispMembership >> test.correctCrispMembership(i,k);
      test.membership(i,k)=test.correctCrispMembership(i,k);
    }
  }
  for(int i=0;i<test.centers_number();i++){
    for(int l=0;l<test.basis(i).cols();l++){
      for(int ell=0;ell<test.dimension();ell++){
        if(l+ell==1.0)
          test.basis(i,ell,l)=1.0;
        else
          test.basis(i,ell,l)=0.0;
      }
    }
  }
    
  for(int l=0;l<test.dimension();l++){
    for(int ell=0;ell<test.dimension();ell++){
      if(l==ell)
        test.psi(ell,l)=1.0;
      else
        test.psi(ell,l)=0.0;
    }
  }
  
#ifdef VERBOSE
  std::cout << "v:\n" << test.centers() << std::endl;
  for(int i=0;i<test.centers_number();i++){
    std::cout << "basis"<< i <<":\n" << test.basis(i) << std::endl;
  }
  std::cout << "psi:\n" << test.psi() << std::endl;
  std::cout << "u:\n" << test.membership() << std::endl;
  std::cout << "a:\n" << test.clusters_size() << std::endl;
#endif
  test.iterates()=0;
  while(1){
    test.revise_inv_psi();
#ifdef VERBOSE
    std::cout << "invPsi:\n" << test.inv_psi() << std::endl;
#endif
    test.revise_m();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      std::cout << "invM" << i <<":\n" << test.m(i) << std::endl;
    }
#endif
    test.revise_ez();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      std::cout << "Ez" << i <<":\n" << test.ez(i) << std::endl;
    }
#endif
    test.revise_ezz();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      Matrix Ezz(sub_dimension, sub_dimension);
      for(int l=0;l<sub_dimension;l++){
        for(int l2=0;l2<sub_dimension;l2++){
          Ezz[l][l2]=0.0;
        }
      }
      for(int k=0;k<test.data_number();k++){
        Ezz+=test.ezz(i)[k];
      }
      std::cout <<"Ezz"<< i <<":\n" << Ezz << std::endl;
    }
#endif
    test.revise_basis();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      std::cout << "basis:\n" << test.basis(i) << std::endl;
    }
#endif

#ifdef VERBOSE3
    for(int i=0;i<test.centers_number();i++){
      Matrix Ezz(sub_dimension, sub_dimension);
      for(int l=0;l<sub_dimension;l++){
        for(int l2=0;l2<sub_dimension;l2++){
          Ezz[l][l2]=0.0;
        }
      }
      for(int k=0;k<test.data_number();k++){
        Ezz+=test.ezz(i)[k];
      }
      std::cout <<"A:\n" << test.basis(i)*Ezz << std::endl;
    }

    for(int i=0;i<test.centers_number();i++){
      Matrix XEz(data_dimension, sub_dimension);
      for(int l=0;l<test.dimension();l++){
        for(int ell=0;ell<test.dimension();ell++){
          XEz[ell][l]=0.0;
        }
      }
      std::cout <<"XEz:\n" << XEz << std::endl;
      for(int k=0;k<test.data_number();k++){
        for(int l=0;l<test.dimension();l++){
          for(int ell=0;ell<test.dimension();ell++){
            XEz[ell][l]+=(test.data(k)[ell]-test.centers(i)[ell])*test.ez(i)[k][l];
          }
        }
      }
      std::cout <<"XEz:\n" << XEz << std::endl;
    }
#endif    
    test.revise_dissimilarities();
#ifdef VERBOSE
    std::cout << "d:\n" << test.dissimilarities() << std::endl;
#endif
    test.revise_membership();
#ifdef VERBOSE
    std::cout << "u:\n" << test.membership() << std::endl;
#endif
    test.revise_centers();
#ifdef VERBOSE
    std::cout << "v:\n" << test.centers() << std::endl;
#endif
    test.revise_psi();
#ifdef VERBOSE
    std::cout << "psi:\n" << test.psi() << std::endl;
#endif
    test.revise_clusters_size();
#ifdef VERBOSE
    std::cout << "a:\n" << test.clusters_size() << std::endl;
#endif
    
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      Matrix temp(test.basis(i).rows(), test.basis(i).cols());
      temp=transpose(test.basis(i));
      for(int l=0;l<test.basis(i).cols();l++){          
        std::cout << "v_" << i << "+t_{" << i <<  l <<"}:\n" <<
          test.centers()[i]+temp[l] << std::endl;
      }
    }
#endif
    
    double diff_u=max_norm(test.tmp_membership()-test.membership());
    double diff_v=max_norm(test.tmp_centers()-test.centers());
    double diff_a=max_norm(test.tmp_clusters_size()-test.clusters_size());
    double diff_W=0;
    double diff_Psi=max_norm(test.tmp_psi()-test.psi());    
    for(int i=0;i<test.centers_number();i++){
      diff_W+=max_norm(test.tmp_basis(i)-test.basis(i));
    }
    double diff=diff_u+diff_v+diff_W+diff_Psi+diff_a;
#ifdef DIFF
    std::cout << "#diff:" << diff << "\t";
    std::cout << "#diff_u:" << diff_u << "\t";
    std::cout << "#diff_v:" << diff_v << "\t";
    std::cout << "#diff_a:" << diff_a << "\t";
    std::cout << "#diff_W:" << diff_W << "\t";
    std::cout << "#diff_Psi:" << diff_Psi << "\n";
    std::cout << "#iterates:" << test.iterates() << "\n";
#endif
    if(diff<DIFF_FOR_STOP)break;
    if(test.iterates()>=MAX_ITERATES)break;
    test.iterates()++;
  }
#ifdef VERBOSE
  std::cout << "v:\n" << test.centers() << std::endl;
#endif
  
#ifdef CHECK_ANSWER
  test.set_crispMembership();
  test.set_contingencyTable();
  std::cout << "Contingency Table:\n"
            << test.contingencyTable() << std::endl;
  std::cout <<Lambda<< "\t"<< test.ARI() << std::endl;
#endif

#ifdef CRFILE
  std::string filenameResultMembership
    =std::string("./debug_data/eMFA/")
    +std::string("eMFA-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_membership");
  std::ofstream ofs_membership(filenameResultMembership);
  if(!ofs_membership){
    std::cerr << "File:" << filenameResultMembership
	      << "could not open." << std::endl;
    exit(1);
  }

  for(int k=0;k<test.data_number();k++){
    for(int ell=0;ell<test.dimension();ell++){
      ofs_membership << test.data()[k][ell] << "\t";
    }
    for(int i=0;i<test.centers_number();i++){
      ofs_membership << test.membership()[i][k] << "\t";
    }
    ofs_membership << std::endl;
  }
  ofs_membership.close();

  std::string filenameResultCenters
    =std::string("./debug_data/eMFA/")
    +std::string("eMFA-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_centers");
  std::ofstream ofs_centers(filenameResultCenters);
  if(!ofs_centers){
    std::cerr << "File:" << filenameResultCenters
	      << "could not open." << std::endl;
    exit(1);
  }
  for(int i=0;i<test.centers_number();i++){
    for(int ell=0;ell<test.dimension();ell++){
      ofs_centers << test.centers()[i][ell] << "\t";
    }
    ofs_centers << std::endl;
  }
  ofs_centers.close();

  std::string filenameResultBasis
    =std::string("./debug_data/eMFA/")
    +std::string("eMFA-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_basis");
  std::ofstream ofs_basis(filenameResultBasis);
  if(!ofs_centers){
    std::cerr << "File:" << filenameResultCenters
	      << "could not open." << std::endl;
    exit(1);
  }
  for(int i=0;i<test.centers_number();i++){
    for(int l=0;l<test.basis(i).cols();l++){
      for(int ell=0;ell<test.dimension();ell++){
        ofs_basis << test.centers()[i][ell] << "\t";
      }
      ofs_basis << std::endl;
      for(int ell=0;ell<test.dimension();ell++){
        ofs_basis << test.centers()[i][ell]+test.basis(i)[ell][l] << "\t";
      }
      ofs_basis << "\n\n"<< std::endl;
    }
  }
  ofs_basis.close();
  
#ifdef CLASSIFICATION_FUNCTION
  //Classification Function
  if(test.dimension()>2){
    std::cerr << "Dimension:" << test.dimension()
              << "is too high for classification function visualization."
              << std::endl;
    exit(1);
  }
  Emfa ClassFunction(test.dimension(), 1, test.centers_number(),
                       1, test.fuzzifierLambda());
  std::string filenameClassificationFunction
    =std::string("./debug_data/eMFA/")
    +std::string("eMFA-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_classificationFunction");
  std::ofstream ofs_classificationFunction(filenameClassificationFunction);
  if(!ofs_classificationFunction){
    std::cerr << "File:" << filenameClassificationFunction
	      << "could not open." << std::endl;
    exit(1);
  }
  for(int i=0;i<test.centers_number();i++){
    ClassFunction.centers(i)=test.centers(i);
  }
  for(int i=0;i<test.centers_number();i++){
    for(int l=0;l<test.basis(i).cols();l++){
      for(int ell=0;ell<test.dimension();ell++){
        ClassFunction.basis(i)[ell][l]=test.basis(i)[ell][l];
      }
    }
  }
  for(int l=0;l<test.dimension();l++){
    for(int ell=0;ell<test.dimension();ell++){
      ClassFunction.psi(ell,l)=test.psi(ell,l);
    }
  }
  for(int i=0;i<test.centers_number();i++){
    ClassFunction.clusters_size(i)=test.clusters_size(i);
  }
  
  Vector Min(test.dimension(), DBL_MAX, "all");
  Vector Max(test.dimension(), -DBL_MAX, "all");
  for(int k=0;k<test.data_number();k++){
    for(int ell=0;ell<test.dimension();ell++){
      if(Min[ell]>test.data(k, ell)){
	Min[ell]=test.data(k, ell);
      }
      if(Max[ell]<test.data(k, ell)){
	Max[ell]=test.data(k, ell);
      }
    }
  }
  Vector Mid=0.5*(Max+Min);
  Vector Width=3.0*(Max-Min);
  Min=Mid-Width;
  Max=Mid+Width;
  int flag=0;
  int flag2=0;
  for(double x0=Min[0];x0<=Max[0];x0+=Width[0]/5.0){
    if(fabs(x0)<0.1&&flag==0){
      x0=-0.05;
      flag++;
    }
    if(x0>=0.8&&flag!=0)
      x0-=3.55;
    for(double x1=Min[1];x1<=Max[1];x1+=Width[1]/5.0){
      if(fabs(x1)<0.1&&flag2==0){
        x1=-0.05;
        flag2++;
      }
      if(x1>=0.8&&flag2!=0)
        x1-=3.55;
#ifdef VERBOSE
      std::cout << "x0:" << x0 << "\t" << "x1:" << x1 << std::endl;
#endif
      ClassFunction.data(0,0)=x0;
      ClassFunction.data(0,1)=x1;
      while(1){
        ClassFunction.revise_inv_psi();
        ClassFunction.revise_m();
        ClassFunction.revise_ez();
        ClassFunction.revise_ezz();
        ClassFunction.revise_dissimilarities();
        ClassFunction.revise_membership();
        double diff_u=frobenius_norm(ClassFunction.tmp_membership()
                                     -ClassFunction.membership());
#ifdef DIFF2
        std::cout << "diff_u:" << diff_u << std::endl;
#endif
        if(diff_u<DIFF_FOR_STOP)break;
      }
      for(int ell=0;ell<ClassFunction.dimension();ell++){
      ofs_classificationFunction << ClassFunction.data()[0][ell] << "\t";
      }
      for(int i=0;i<ClassFunction.centers_number();i++){
      ofs_classificationFunction << ClassFunction.membership()[i][0] <<"\t";
      }
      double max=0.0;
      for(int i=0;i<ClassFunction.centers_number();i++){
        if(max<ClassFunction.membership()[i][0]){
          max=ClassFunction.membership()[i][0];
        }
      }
      ofs_classificationFunction << max << "\t";
      ofs_classificationFunction << std::endl;
      if(x1>=0.01&&flag2!=0){
        x1=0.0;
        flag2=0;
      }
    }
    ofs_classificationFunction << std::endl;
    if(x0>=0.01&&flag!=0){
      x0=0.0;
      flag=0;
    }
  }
#endif
#endif
  
#ifdef WHILE
  //Lambda*=stridelambda;
  }
#endif
  return 0;
}
