# https://github.com/mollios/some-datasets
# $ make (vectorクラス等の中間ファイル作成、-Bで強制コンパイル) 
# $ make ターゲット(.out)で実行ファイルコンパイル
# 実データの場合:
# $ make ターゲット dmakeata=マクロ名
# 例(klfcsをBookCrossingに適用する場合):w
# $ make rklfcs.out data=BOOK
# クラスの呼び出し等デバグしたいとき
# $ make ターゲット data=マクロ名 class=1
# 人工データの場合:
# $ make ターゲット
CXX = g++
CXXFLAGS = -std=c++17 -DDIFF -DCRFILE -DCHECK_ANSWER
FS = -lstdc++fs
objects = .o/vector.o .o/matrix.o .o/tensor.o \
.o/tensor4d.o .o/eigen.o .o/libLU.o
efcv = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/efcv.o
sfcv = $(objects) .o/hcm.o .o/hcma.o \
.o/sfcm.o .o/sfcma.o .o/sfcv.o
qfcv = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/sfcm.o .o/sfcma.o \
.o/qfcm.o .o/qfcma.o .o/qfcv.o
emppca_eigen = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/mppca_eigen.o .o/emppca_eigen.o
smppca_eigen = $(objects) .o/hcm.o .o/hcma.o \
.o/sfcm.o .o/sfcma.o .o/mppca_eigen.o .o/smppca_eigen.o
qmppca_eigen = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/sfcm.o .o/sfcma.o \
.o/qfcm.o .o/qfcma.o .o/mppca_eigen.o .o/qmppca_eigen.o
emppca = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/mppca.o .o/emppca.o
emppca_3d = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/mppca_3d.o .o/emppca_3d.o
smppca = $(objects) .o/hcm.o .o/hcma.o \
.o/sfcm.o .o/sfcma.o .o/mppca.o .o/smppca.o
qmppca = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/sfcm.o .o/sfcma.o \
.o/qfcm.o .o/qfcma.o .o/mppca.o .o/qmppca.o
emfa = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/mppca.o .o/mfa.o .o/emfa.o
smfa = $(objects) .o/hcm.o .o/hcma.o \
.o/sfcm.o .o/sfcma.o .o/mppca.o .o/mfa.o .o/smfa.o
qmfa = $(objects) .o/hcm.o .o/hcma.o \
.o/efcm.o .o/efcma.o .o/sfcm.o .o/sfcma.o \
.o/qfcm.o .o/qfcma.o .o/mppca.o .o/mfa.o .o/qmfa.o

method_all = $(all) \
$(efcv) $(sfcv) $(qfcv) \
$(emppca_eigen) $(smppca_eigen) $(qmppca_eigen) \
$(emppca) $(smppca) $(qmppca) $(emppca_3d)\
$(emfa) $(smfa) $(qmfa) \
cross_efcv.out \
cross_sfcv.out \
cross_qfcv.out \
cross_emppca.out \
cross_emppca_3d.out \
cross_smppca.out \
cross_qmppca.out \
cross_emppca_eigen.out \
cross_smppca_eigen.out \
cross_qmppca_eigen.out \
cross_emfa.out \
cross_smfa.out \
cross_qmfa.out \
para_efcv.out \
para_sfcv.out \
para_qfcv.out \
para_emppca.out \
para_smppca.out \
para_qmppca.out \
para_emppca_eigen.out \
para_smppca_eigen.out \
para_qmppca_eigen.out \
para_emfa.out \
para_smfa.out \
para_qmfa.out \
mnist_efcv.out \
mnist_sfcv.out \
mnist_qfcv.out \
mnist_emppca_eigen.out \
mnist_smppca_eigen.out \
mnist_qmppca_eigen.out \

ifdef data
	DATASET=-D$(data) 
endif
ifdef class
	MACRO=-DCHECK_CLASS 
endif
ifdef log
	LOG=-DVERBOSE 
endif
ifdef while
	WHILE=-DWHILE 
endif
ifdef a
	A=-D$(a) 
endif

all : $(objects) 

method_all : $(method_all)

.o/vector.o : src/vector.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/matrix.o : src/matrix.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/tensor.o : src/tensor.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/tensor4d.o : src/tensor_4d.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/eigen.o : src/eigen.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/libLU.o : src/libLU.cxx 
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/hcm.o : src/hcm.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/hcma.o : src/hcma.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/efcm.o : src/efcm.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/efcma.o : src/efcma.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/sfcm.o : src/sfcm.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/sfcma.o : src/sfcma.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/qfcm.o : src/qfcm.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/qfcma.o : src/qfcma.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/efcv.o : src/efcv.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/sfcv.o : src/sfcv.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/qfcv.o : src/qfcv.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) -c $^ -o $@
.o/mppca_eigen.o : src/mppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/emppca_eigen.o : src/emppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/smppca_eigen.o : src/smppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/qmppca_eigen.o : src/qmppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/mppca.o : src/mppca.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/mppca_3d.o : src/mppca_3d.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/emppca.o : src/emppca.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/emppca_3d.o : src/emppca_3d.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/smppca.o : src/smppca.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/qmppca.o : src/qmppca.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/mfa.o : src/mfa.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/emfa.o : src/emfa.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/smfa.o : src/smfa.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@
.o/qmfa.o : src/qmfa.cxx
	$(CXX) $(CXXFLAGS) $(MACRO) $(LOG) $(WHILE)-c $^ -o $@


#クラスタリング人工データ(クロス)
cross_efcv.out : $(efcv) \
main_clustering/artificiality/cross/efcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_sfcv.out : $(sfcv) \
main_clustering/artificiality/cross/sfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_qfcv.out : $(qfcv) \
main_clustering/artificiality/cross/qfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_emppca_eigen.out : $(emppca_eigen) \
main_clustering/artificiality/cross/emppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_smppca_eigen.out : $(smppca_eigen) \
main_clustering/artificiality/cross/smppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_qmppca_eigen.out : $(qmppca_eigen) \
main_clustering/artificiality/cross/qmppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_emppca.out : $(emppca) \
main_clustering/artificiality/cross/emppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_emppca_3d.out : $(emppca_3d) \
main_clustering/artificiality/cross/emppca_3d.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_smppca.out : $(smppca) \
main_clustering/artificiality/cross/smppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_qmppca.out : $(qmppca) \
main_clustering/artificiality/cross/qmppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_emfa.out : $(emfa) \
main_clustering/artificiality/cross/emfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_smfa.out : $(smfa) \
main_clustering/artificiality/cross/smfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
cross_qmfa.out : $(qmfa) \
main_clustering/artificiality/cross/qmfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@


#クラスタリング人工データ（平行）
para_efcv.out : $(efcv) \
main_clustering/artificiality/para/efcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_sfcv.out : $(sfcv) \
main_clustering/artificiality/para/sfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_qfcv.out : $(qfcv) \
main_clustering/artificiality/para/qfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_emppca_eigen.out : $(emppca_eigen) \
main_clustering/artificiality/para/emppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_smppca_eigen.out : $(smppca_eigen) \
main_clustering/artificiality/para/smppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_qmppca_eigen.out : $(qmppca_eigen) \
main_clustering/artificiality/para/qmppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_emppca.out : $(emppca) \
main_clustering/artificiality/para/emppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_smppca.out : $(smppca) \
main_clustering/artificiality/para/smppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_qmppca.out : $(qmppca) \
main_clustering/artificiality/para/qmppca.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_emfa.out : $(emfa) \
main_clustering/artificiality/para/emfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_smfa.out : $(smfa) \
main_clustering/artificiality/para/smfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@
para_qmfa.out : $(qmfa) \
main_clustering/artificiality/para/qmfa.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	-DCLASSIFICATION_FUNCTION $(FS) -o $@


#クラスタリングmnist
mnist_efcv.out : $(efcv) \
main_clustering/mnist/efcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@
mnist_sfcv.out : $(sfcv) \
main_clustering/mnist/sfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@
mnist_qfcv.out : $(qfcv) \
main_clustering/mnist/qfcv.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@
mnist_emppca_eigen.out : $(emppca_eigen) \
main_clustering/mnist/emppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@
mnist_smppca_eigen.out : $(smppca_eigen) \
main_clustering/mnist/smppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@
mnist_qmppca_eigen.out : $(qmppca_eigen) \
main_clustering/mnist/qmppca_eigen.cxx
	$(CXX) $(CXXFLAGS) $(LOG) $(WHILE) $^ \
	$(DATASET) $(FS) -o $@

#クラスタリングoil

clean :
	rm -f *.out
clean.o :
	rm -f .o/*.o
