#include "vector.h"
#include "matrix.h"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include "eigen.h"

double max_abs(const Vector &arg){
	double result=fabs(arg[0]);
	for(int i=1;i<arg.size();i++){
		double tmp=fabs(arg[i]);
		if(result<tmp)	result=tmp;
	}
	return result;
}

const struct Eigen eigen_power(const Matrix &arg, const double &stop){
	if(arg.rows()!=arg.cols()){
		std::cerr << "eigen:Not Square" << std::endl;
		exit(1);
	}
	Vector x(arg.cols());
	for(int i=0;i<x.size();i++){
		x[i]=1.0;
	}
	struct Eigen result;
	double tmp_r=DBL_MAX;

	while(1){
		double ell=max_abs(x);
		result.vector=x/ell;
		x=arg*result.vector;
		result.value=(result.vector*x)/(result.vector*result.vector);
#ifdef EIGEN_DEBUG
		std::cout << "in Eigen:" << fabs(result.value-tmp_r) << std::endl;
#endif
		if(fabs(result.value-tmp_r)<stop)break;
		tmp_r=result.value;
	}

	return result;
}

void eigen_jacobi0(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//Q:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
	//Lambdaの非対角要素のうちで絶対値最大要素の添字を(I,J)とする。
		int I=-1, J=-1; double absmax=-DBL_MAX;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j && absmax<fabs(Lambda[i][j])){
				I=i; J=j; absmax=fabs(Lambda[i][j]);
			}
		}}
#ifdef DEBUG_JACOBI
		std::cout << "I,J:" << I << "\t" << J << std::endl;
#endif
	//Lambda[I][I]==Lambda[J][J] => theta=M_PI_4
	//otherwise => theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]))
		double theta;
		if(Lambda[I][I]==Lambda[J][J]){
			theta=M_PI_4;
		}
		else{
			theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]));
		}
#ifdef DEBUG_JACOBI
		std::cout << "theta:" << theta << std::endl;
#endif
		Matrix T(Lambda.rows(), "I");
		T[I][I]=cos(theta);
		T[I][J]=-sin(theta);
		T[J][I]=sin(theta);
		T[J][J]=cos(theta);
		Lambda=transpose(T)*Lambda*T;
		TT=TT*T;
#ifdef DEBUG_JACOBI
		std::cout << "TT:" << TT << std::endl;
		std::cout << "Lambda:" << Lambda << std::endl;
#endif

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}

void eigen_jacobi(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//Q:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
	//Lambdaの非対角要素のうちで絶対値最大要素の添字を(I,J)とする。
	//対称性から上三角部分(i<j)のみ調べれば良い
		int I=-1, J=-1; double absmax=-DBL_MAX;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=i+1;j<Lambda.cols();j++){
			if(absmax<fabs(Lambda[i][j])){
				I=i; J=j; absmax=fabs(Lambda[i][j]);
			}
		}}
#ifdef DEBUG_JACOBI
		std::cout << "I,J:" << I << "\t" << J << std::endl;
#endif
	//Lambda[I][I]==Lambda[J][J] => theta=M_PI_4
	//otherwise => theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]))
		double theta;
		if(Lambda[I][I]==Lambda[J][J]){
			theta=M_PI_4;
		}
		else{
			theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]));
		}
#ifdef DEBUG_JACOBI
		std::cout << "theta:" << theta << std::endl;
#endif
//TTを更新
		for(int i=0;i<TT.rows();i++){
			double tmpTTI=TT[i][I]*cos(theta)+TT[i][J]*sin(theta);
			double tmpTTJ=-TT[i][I]*sin(theta)+TT[i][J]*cos(theta);
			TT[i][I]=tmpTTI;
			TT[i][J]=tmpTTJ;
		}
#ifdef DEBUG_JACOBI
		std::cout << "TT:" << TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
		for(int j=0;j<Lambda.cols();j++){
			if(j!=I && j!=J){
				double LambdaIj=Lambda[I][j]*cos(theta)+Lambda[J][j]*sin(theta);
				double LambdaJj=-Lambda[I][j]*sin(theta)+Lambda[J][j]*cos(theta);
				Lambda[I][j]=LambdaIj;
				Lambda[J][j]=LambdaJj;
			}
		}
//I,J列目
		for(int i=0;i<Lambda.rows();i++){
			if(i!=I && i!=J){
				double LambdaiI=Lambda[i][I]*cos(theta)+Lambda[i][J]*sin(theta);
				double LambdaiJ=-Lambda[i][I]*sin(theta)+Lambda[i][J]*cos(theta);
				Lambda[i][I]=LambdaiI;
				Lambda[i][J]=LambdaiJ;
			}
		}
//(I,I),(J,J)成分
		double LambdaII=0.5*(Lambda[I][I]+Lambda[J][J])
			+0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			+Lambda[I][J]*sin(2.0*theta);
		double LambdaJJ=0.5*(Lambda[I][I]+Lambda[J][J])
			-0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			-Lambda[I][J]*sin(2.0*theta);
		Lambda[I][I]=LambdaII;
		Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
		Lambda[I][J]=Lambda[J][I]=0.0;
//それ以外は更新なし
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << Lambda << std::endl;
#endif

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}

void eigen_jacobi_oppositetheta(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//Q:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
	//Lambdaの非対角要素のうちで絶対値最大要素の添字を(I,J)とする。
	//対称性から上三角部分(i<j)のみ調べれば良い
		int I=-1, J=-1; double absmax=-DBL_MAX;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=i+1;j<Lambda.cols();j++){
			if(absmax<fabs(Lambda[i][j])){
				I=i; J=j; absmax=fabs(Lambda[i][j]);
			}
		}}
#ifdef DEBUG_JACOBI
		std::cout << "I,J:" << I << "\t" << J << std::endl;
#endif
	//Lambda[I][I]==Lambda[J][J] => theta=M_PI_4
	//otherwise => theta=-0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]))
		double theta;
		if(Lambda[I][I]==Lambda[J][J]){
			theta=M_PI_4;
		}
		else{
			theta=-0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]));
		}
#ifdef DEBUG_JACOBI
		std::cout << "theta:" << theta << std::endl;
#endif
//TTを更新
		for(int i=0;i<TT.rows();i++){
			double tmpTTI=TT[i][I]*cos(theta)-TT[i][J]*sin(theta);
			double tmpTTJ=TT[i][I]*sin(theta)+TT[i][J]*cos(theta);
			TT[i][I]=tmpTTI;
			TT[i][J]=tmpTTJ;
		}
#ifdef DEBUG_JACOBI
		std::cout << "TT:" << TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
		for(int j=0;j<Lambda.cols();j++){
			if(j!=I && j!=J){
				double LambdaIj=Lambda[I][j]*cos(theta)-Lambda[J][j]*sin(theta);
				double LambdaJj=Lambda[I][j]*sin(theta)+Lambda[J][j]*cos(theta);
				Lambda[I][j]=LambdaIj;
				Lambda[J][j]=LambdaJj;
			}
		}
//I,J列目
		for(int i=0;i<Lambda.rows();i++){
			if(i!=I && i!=J){
				double LambdaiI=Lambda[i][I]*cos(theta)-Lambda[i][J]*sin(theta);
				double LambdaiJ=Lambda[i][I]*sin(theta)+Lambda[i][J]*cos(theta);
				Lambda[i][I]=LambdaiI;
				Lambda[i][J]=LambdaiJ;
			}
		}
//(I,I),(J,J)成分
		double LambdaII=0.5*(Lambda[I][I]+Lambda[J][J])
			+0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			-Lambda[I][J]*sin(2.0*theta);
		double LambdaJJ=0.5*(Lambda[I][I]+Lambda[J][J])
			-0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			+Lambda[I][J]*sin(2.0*theta);
		Lambda[I][I]=LambdaII;
		Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
		Lambda[I][J]=Lambda[J][I]=0.0;
//それ以外は更新なし
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << Lambda << std::endl;
#endif

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}

void eigen_jacobi_cyclic(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//Q:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
	//枢軸(I,J)を巡回
	//対称性から上三角部分(i<j)のみ調べれば良い
		for(int I=0;I<Lambda.rows()-1;I++){
		for(int J=I+1;J<Lambda.cols();J++){
		if(Lambda[I][J]==0.0)continue;
	//Lambda[I][I]==Lambda[J][J] => theta=M_PI_4
	//otherwise => theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]))
		double theta;
		if(Lambda[I][I]==Lambda[J][J]){
			theta=M_PI_4;
		}
		else{
			theta=0.5*atan(2.0*Lambda[I][J]/(Lambda[I][I]-Lambda[J][J]));
		}
#ifdef DEBUG_JACOBI
		std::cout << "theta:" << theta << std::endl;
#endif
//TTを更新
		for(int i=0;i<TT.rows();i++){
			double tmpTTI=TT[i][I]*cos(theta)+TT[i][J]*sin(theta);
			double tmpTTJ=-TT[i][I]*sin(theta)+TT[i][J]*cos(theta);
			TT[i][I]=tmpTTI;
			TT[i][J]=tmpTTJ;
		}
#ifdef DEBUG_JACOBI
		std::cout << "TT:" << TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
		for(int j=0;j<Lambda.cols();j++){
			if(j!=I && j!=J){
				double LambdaIj=Lambda[I][j]*cos(theta)+Lambda[J][j]*sin(theta);
				double LambdaJj=-Lambda[I][j]*sin(theta)+Lambda[J][j]*cos(theta);
				Lambda[I][j]=LambdaIj;
				Lambda[J][j]=LambdaJj;
			}
		}
//I,J列目
		for(int i=0;i<Lambda.rows();i++){
			if(i!=I && i!=J){
				double LambdaiI=Lambda[i][I]*cos(theta)+Lambda[i][J]*sin(theta);
				double LambdaiJ=-Lambda[i][I]*sin(theta)+Lambda[i][J]*cos(theta);
				Lambda[i][I]=LambdaiI;
				Lambda[i][J]=LambdaiJ;
			}
		}
//(I,I),(J,J)成分
		double LambdaII=0.5*(Lambda[I][I]+Lambda[J][J])
			+0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			+Lambda[I][J]*sin(2.0*theta);
		double LambdaJJ=0.5*(Lambda[I][I]+Lambda[J][J])
			-0.5*(Lambda[I][I]-Lambda[J][J])*cos(2.0*theta)
			-Lambda[I][J]*sin(2.0*theta);
		Lambda[I][I]=LambdaII;
		Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
		Lambda[I][J]=Lambda[J][I]=0.0;
//それ以外は更新なし
		}}//巡回の終り
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << Lambda << std::endl;
#endif

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}

void eigen_jacobi_rutishauser(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//TT:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
//Lambdaの非対角要素のうちで絶対値最大要素の添字を(I,J)とする。
//対称性から上三角部分(i<j)のみ調べれば良い
		int I=-1, J=-1; double absmax=-DBL_MAX;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=i+1;j<Lambda.cols();j++){
			if(absmax<fabs(Lambda[i][j])){
				I=i; J=j; absmax=fabs(Lambda[i][j]);
			}
		}}
//Rutishauserの計算式
			if(Lambda[I][J]==0.0)continue;
			double z=(Lambda[J][J]-Lambda[I][I])
				/(2.0*Lambda[I][J]);
			int signz=((z>0.0) ? 1: -1);
			double t=signz/(fabs(z)+sqrt(1+z*z));
			double c=1.0/sqrt(1.0+t*t);
			double s=c*t;
//TTを更新
			for(int i=0;i<TT.rows();i++){
				double tmpTTI=TT[i][I]*c-TT[i][J]*s;
				double tmpTTJ=TT[i][I]*s+TT[i][J]*c;
				TT[i][I]=tmpTTI;
				TT[i][J]=tmpTTJ;
			}
#ifdef DEBUG_JACOBI
			std::cout << "TT:" << TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
//I,J列目
			for(int j=0;j<Lambda.cols();j++){
				if(j!=I && j!=J){
					double LambdaIj=Lambda[I][j]*c-Lambda[J][j]*s;
					double LambdaJj=Lambda[I][j]*s+Lambda[J][j]*c;
					Lambda[I][j]=Lambda[j][I]=LambdaIj;
					Lambda[J][j]=Lambda[j][J]=LambdaJj;
//					Lambda[I][j]=LambdaIj;
//					Lambda[J][j]=LambdaJj;
				}
			}
//(I,I),(J,J)成分
			double LambdaII=Lambda[I][I]*c*c
				-2.0*Lambda[I][J]*s*c
				+Lambda[J][J]*s*s;
			double LambdaJJ=Lambda[I][I]
				+Lambda[J][J]
				-LambdaII;
			Lambda[I][I]=LambdaII;
			Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
			Lambda[I][J]=Lambda[J][I]=0.0;
//それ以外は更新なし
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << Lambda << std::endl;
#endif

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}
void eigen_jacobi_cyclic_rutishauser(const Matrix &A, const double &stop, Matrix &Lambda, Matrix &TT){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//TT:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	Lambda=A;
	for(int i=0;i<TT.rows();i++){
	for(int j=0;j<TT.cols();j++){
		TT[i][j]=0.0;
	}}
	for(int i=0;i<TT.rows();i++){
		TT[i][i]=1.0;
	}
	int cnt=0;
	while(1){
//枢軸(I,J)を巡回させる
//対称性から上三角部分(i<j)のみ調べれば良い
		for(int I=0;I<Lambda.rows()-1;I++){
		for(int J=I+1;J<Lambda.cols();J++){
//Rutishauserの計算式
			if(Lambda[I][J]==0.0)continue;
			double z=(Lambda[J][J]-Lambda[I][I])
				/(2.0*Lambda[I][J]);
			int signz=((z>0.0) ? 1: -1);
			double t=signz/(fabs(z)+sqrt(1+z*z));
			double c=1.0/sqrt(1.0+t*t);
			double s=c*t;
//TTを更新
			for(int i=0;i<TT.rows();i++){
				double tmpTTI=TT[i][I]*c-TT[i][J]*s;
				double tmpTTJ=TT[i][I]*s+TT[i][J]*c;
				TT[i][I]=tmpTTI;
				TT[i][J]=tmpTTJ;
			}
#ifdef DEBUG_JACOBI
			std::cout << "TT:" << TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
//I,J列目
			for(int j=0;j<Lambda.cols();j++){
				if(j!=I && j!=J){
					double LambdaIj=Lambda[I][j]*c-Lambda[J][j]*s;
					double LambdaJj=Lambda[I][j]*s+Lambda[J][j]*c;
					Lambda[I][j]=Lambda[j][I]=LambdaIj;
					Lambda[J][j]=Lambda[j][J]=LambdaJj;
//					Lambda[I][j]=LambdaIj;
//					Lambda[J][j]=LambdaJj;
				}
			}
//(I,I),(J,J)成分
			double LambdaII=Lambda[I][I]*c*c
				-2.0*Lambda[I][J]*s*c
				+Lambda[J][J]*s*s;
			double LambdaJJ=Lambda[I][I]
				+Lambda[J][J]
				-LambdaII;
			Lambda[I][I]=LambdaII;
			Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
			Lambda[I][J]=Lambda[J][I]=0.0;
//それ以外は更新なし
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << Lambda << std::endl;
#endif
		}}//巡回の終り

		double residual=0.0;
		for(int i=0;i<Lambda.rows();i++){
		for(int j=0;j<Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DEBUG_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
	return;
}

const struct AllEigen eigen_jacobi_cyclic_rutishauser(const Matrix &A, const double &stop){
//Lambda:Aで初期化しておいて固有値を対角に持つ行列にする。
//TT:Eで初期化しておいて固有ベクトルを列成分に持つ行列にする。
	struct AllEigen result;
	result.Lambda=A;
	result.TT=Matrix(A.rows(),"I");
	int cnt=0;
	while(1){
//枢軸(I,J)を巡回させる
//対称性から上三角部分(i<j)のみ調べれば良い
		for(int I=0;I<result.Lambda.rows()-1;I++){
		for(int J=I+1;J<result.Lambda.cols();J++){
//Rutishauserの計算式
			if(result.Lambda[I][J]==0.0)continue;
			double z=(result.Lambda[J][J]-result.Lambda[I][I])
				/(2.0*result.Lambda[I][J]);
			int signz=((z>0.0) ? 1: -1);
			double t=signz/(fabs(z)+sqrt(1+z*z));
			double c=1.0/sqrt(1.0+t*t);
			double s=c*t;
//TTを更新
			for(int i=0;i<result.TT.rows();i++){
				double tmpTTI=result.TT[i][I]*c-result.TT[i][J]*s;
				double tmpTTJ=result.TT[i][I]*s+result.TT[i][J]*c;
				result.TT[i][I]=tmpTTI;
				result.TT[i][J]=tmpTTJ;
			}
#ifdef DEBUG_JACOBI
			std::cout << "TT:" << result.TT << std::endl;
#endif

//Lambdaの更新
//2行(列)づつ、中間変数を用いて更新
//I,J行目
//I,J列目
			for(int j=0;j<result.Lambda.cols();j++){
				if(j!=I && j!=J){
					double LambdaIj=result.Lambda[I][j]*c-result.Lambda[J][j]*s;
					double LambdaJj=result.Lambda[I][j]*s+result.Lambda[J][j]*c;
					result.Lambda[I][j]=result.Lambda[j][I]=LambdaIj;
					result.Lambda[J][j]=result.Lambda[j][J]=LambdaJj;
//					result.Lambda[I][j]=LambdaIj;
//					result.Lambda[J][j]=LambdaJj;
				}
			}
//(I,I),(J,J)成分
			double LambdaII=result.Lambda[I][I]*c*c
				-2.0*result.Lambda[I][J]*s*c
				+result.Lambda[J][J]*s*s;
			double LambdaJJ=result.Lambda[I][I]
				+result.Lambda[J][J]
				-LambdaII;
			result.Lambda[I][I]=LambdaII;
			result.Lambda[J][J]=LambdaJJ;
//(I,J),(J,I)成分
			result.Lambda[I][J]=result.Lambda[J][I]=0.0;
//それ以外は更新なし
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:" << result.Lambda << std::endl;
#endif
		}}//巡回の終り

		double residual=0.0;
		for(int i=0;i<result.Lambda.rows();i++){
		for(int j=0;j<result.Lambda.cols();j++){
			if(i!=j){
				residual+=fabs(result.Lambda[i][j]);
			}
		}}
		if(residual<stop)break;
		cnt++;
#ifdef DIFF_JACOBI
		std::cout << "cnt:" << cnt << "\t" << residual << std::endl;
#endif
	}
#ifdef DEBUG_JACOBI
		std::cout << "Lambda:\n" << result.Lambda << "\nTT:\n" << result.TT << std::endl;
#endif
	return result;
}

void sort_abs_alleigen(struct AllEigen &arg){
	for(int k=0;k<arg.Lambda.rows()-1;k++){
		double absmax=-1.0; int absmax_index=-1;
		for(int i=k;i<arg.Lambda.rows();i++){
			double tmp=fabs(arg.Lambda[i][i]);
			if(absmax<tmp){
				absmax=tmp;
				absmax_index=i;
			}
		}
#ifdef DEBUG_SORT_ALLEIGEN
		std::cout << "absmax:" << absmax << std::endl;
		std::cout << "absmax_index:" << absmax_index << std::endl;
#endif
		double tmp=arg.Lambda[absmax_index][absmax_index];
		arg.Lambda[absmax_index][absmax_index]=arg.Lambda[k][k];
		arg.Lambda[k][k]=tmp;
		Vector tmpVector(arg.TT.rows());
		for(int j=0;j<tmpVector.size();j++){
			tmpVector[j]=arg.TT[j][absmax_index];
		}
		for(int j=0;j<tmpVector.size();j++){
			arg.TT[j][absmax_index]=arg.TT[j][k];
		}
		for(int j=0;j<tmpVector.size();j++){
			arg.TT[j][k]=tmpVector[j];
		}
#ifdef DEBUG_SORT_ALLEIGEN
		std::cout << "Lambda[" << k << "][" << k << "]:" << arg.Lambda[k][k] << std::endl;
		std::cout << "Lambda[" << absmax_index << "][" << absmax_index << "]:" << arg.Lambda[absmax_index][absmax_index] << std::endl;
#endif
	}
	return;
}

void sort_alleigen(struct AllEigen &arg){
	for(int k=0;k<arg.Lambda.rows()-1;k++){
		double max=-DBL_MAX; int max_index=-1;
		for(int i=k;i<arg.Lambda.rows();i++){
			double tmp=arg.Lambda[i][i];
			if(max<tmp){
				max=tmp;
				max_index=i;
			}
		}
#ifdef DEBUG_SORT_ALLEIGEN
		std::cout << "max:" << max << std::endl;
		std::cout << "max_index:" << max_index << std::endl;
#endif
		double tmp=arg.Lambda[max_index][max_index];
		arg.Lambda[max_index][max_index]=arg.Lambda[k][k];
		arg.Lambda[k][k]=tmp;
		Vector tmpVector(arg.TT.rows());
		for(int j=0;j<tmpVector.size();j++){
			tmpVector[j]=arg.TT[j][max_index];
		}
		for(int j=0;j<tmpVector.size();j++){
			arg.TT[j][max_index]=arg.TT[j][k];
		}
		for(int j=0;j<tmpVector.size();j++){
			arg.TT[j][k]=tmpVector[j];
		}
#ifdef DEBUG_SORT_ALLEIGEN
		std::cout << "Lambda[" << k << "][" << k << "]:" << arg.Lambda[k][k] << std::endl;
		std::cout << "Lambda[" << max_index << "][" << max_index << "]:" << arg.Lambda[max_index][max_index] << std::endl;
#endif
	}
	return;
}
