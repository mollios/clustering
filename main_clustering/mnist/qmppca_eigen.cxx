#include<iostream>
#include<fstream>
#include<cstdlib>
#include<random>
#include<iomanip>
#include"../../src/qmppca_eigen.h"

using namespace std;
#define MAX_ITERATES 100000
#define DIFF_FOR_STOP 1.0E-10

const int centers_number=10;
const int parameter_number=3;

int main(void){
  double Em=1.0+1.0e-6;
  double Lambda=1.0e-10;
  double para_Em[3]={1.01, 1.5, 2.5};
  double para_Lambda[3]={1.0e-1, 1.0e+0, 1.0e+1};
  double strideEm=1.0e-1;
  double stridelambda=10.0;
  std::string filenameData("./data/image_data/train_0.txt");
#ifdef CHECK_ANSWER
  std::string filenameCorrectCrispMembership
    ("./data/labels/train_0.labels");
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
  std::string filenameResultARI
    =std::string("./mnist/result/")
    +std::string("qMPPCA_eigen-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_ARI");
  std::ofstream ofs_ARI(filenameResultARI,ios::app);
  if(!ofs_ARI){
    std::cerr << "File:" << filenameResultARI
              << "could not open." << std::endl;
    exit(1);
  }
  while(1){
    if(Em>1.6)
      exit(1);
    while(1){
      if(Lambda>1.0e+10)
        break;
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
  
  Qmppca test(data_dimension, data_number, centers_number,
              sub_dimension, Em, Lambda);
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
  
  /***Initial Parameters Setting***/
  for(int i=0;i<test.centers_number();i++){
    test.clusters_size(i)=1.0/centers_number;
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

#ifdef VERBOSE
  std::cout << "v:\n" << test.centers() << std::endl;
  std::cout << "u:\n" << test.membership() << std::endl;
  for(int i=0;i<test.centers_number();i++){
    std::cout << "basis"<< i <<":\n" << test.basis(i) << std::endl;
  }
  std::cout << "sigma:\n" << test.sigma() << std::endl;
  std::cout << "a:\n" << test.clusters_size() << std::endl;
#endif
  test.iterates()=0;
  while(1){
    test.revise_centers();
#ifdef VERBOSE
    std::cout << "v:\n" << test.centers() << std::endl;
#endif
    test.revise_basis();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      std::cout << "basis:\n" << test.basis(i) << std::endl;
    }
    std::cout << "sigma:\n" << test.sigma() << std::endl;
#endif
    test.revise_m();
#ifdef VERBOSE
    for(int i=0;i<test.centers_number();i++){
      std::cout << "invM" << i <<":\n" << test.m(i) << std::endl;
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
    test.revise_clusters_size();
#ifdef VERBOSE
    std::cout << "a:\n" << test.clusters_size() << std::endl;
#endif
    
    double diff_u=max_norm(test.tmp_membership()-test.membership());
    double diff_v=max_norm(test.tmp_centers()-test.centers());
    double diff_a=max_norm(test.tmp_clusters_size()-test.clusters_size());
    double diff_M=0;
    double diff_W=0;
    double diff_sigma=0;
    for(int i=0;i<test.centers_number();i++){
      diff_M+=max_norm(test.tmp_m(i)-test.m(i));
      diff_W+=max_norm(test.tmp_basis(i)-test.basis(i));
      diff_sigma+=max_norm(test.tmp_sigma()-test.sigma());
    }
    double diff=diff_u+diff_v+diff_M+diff_W+diff_sigma+diff_a;
#ifdef DIFF
    std::cout << "#diff:" << diff << "\t";
    std::cout << "#diff_u:" << diff_u << "\t";
    std::cout << "#diff_v:" << diff_v << "\t";
    std::cout << "#diff_a:" << diff_a << "\t";
    std::cout << "#diff_M:" << diff_M << "\t";
    std::cout << "#diff_W:" << diff_W << "\t";
    std::cout << "#diff_sigma:" << diff_sigma << "\n";
    std::cout << "#iterates:" << test.iterates() << "\n";
#endif
    
    if(diff<DIFF_FOR_STOP)break;
    if(test.iterates()>=MAX_ITERATES)break;
    test.iterates()++;
  }
  
#ifdef CHECK_ANSWER
  test.set_crispMembership();
  test.set_contingencyTable();
  // std::cout << "Contingency Table:\n"
  //           << test.contingencyTable() << std::endl;
  std::cout <<  Em <<"\t"<<Lambda<< "\t"<< test.ARI() << std::endl;
#endif
#ifdef WHILE
  ofs_ARI << Em  << "\t" << Lambda << "\t" << test.ARI() << std::endl;
#endif
  
#ifdef CRFILE
  std::string filenameResultMembership
    =std::string("./mnist/qMPPCA_eigen/")
    +std::string("qMPPCA_eigen-")
    +std::to_string(test.fuzzifierEm())+std::string("-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_membership");
  std::ofstream ofs_membership(filenameResultMembership);
  if(!ofs_membership){
    std::cerr << "File:" << filenameResultMembership
	      << " could not open." << std::endl;
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
    =std::string("./mnist/qMPPCA_eigen/")
    +std::string("qMPPCA_eigen-")
    +std::to_string(test.fuzzifierEm())+std::string("-")
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
    =std::string("./mnist/qMPPCA_eigen/")
    +std::string("qMPPCA_eigen-")
    +std::to_string(test.fuzzifierEm())+std::string("-")
    +std::to_string(test.fuzzifierLambda())+std::string("-")
    +filenameData.substr(filenameDataSlashPosition, filenameDataDotPosition)
    +std::string(".result_basis");
  std::ofstream ofs_basis(filenameResultBasis);
  if(!ofs_basis){
    std::cerr << "File:" << filenameResultBasis
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
#endif
  
#ifdef WHILE
  Lambda*=stridelambda;
  }
  Em+=strideEm;
  }
  ofs_ARI.close();
#endif
  return 0;
}
