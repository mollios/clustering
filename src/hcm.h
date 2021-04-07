#include<cmath>
#include<cfloat>
#include"matrix.h"

#ifndef __HCM__
#define __HCM__

class Hcm{
 protected:
  Matrix Data, Centers, Tmp_Centers;
  Matrix Membership, Tmp_Membership, Dissimilarities;
  Matrix CrispMembership, CorrectCrispMembership, ContingencyTable;
  int Iterates;
  double Objective;
 public:
  Hcm(const int &dimension,
      const int &data_number,
      const int &centers_number);
  virtual void revise_membership(void);
  virtual void revise_dissimilarities(void);
  virtual void revise_centers(void);
  int dimension(void) const;
  int data_number(void) const;
  int centers_number(void) const;
  const Matrix centers(void) const;
  const Matrix tmp_centers(void) const;
  const Matrix data(void) const;
  const Matrix membership(void) const;
  const Matrix tmp_membership(void) const;
  int &iterates(void);
  const Matrix dissimilarities(void) const;
  double &data(const int &index1, const int &index2); 
  Vector &data(const int &index1); 
  double &centers(const int &index1, const int &index2); 
  Vector &centers(const int &index1); 
  double &membership(const int &index1, const int &index2); 
  double &dissimilarities(const int &index1, const int &index2);
  virtual void set_objective(void);
  const double objective(void) const;
  void set_crispMembership(void);
  const Matrix crispMembership(void) const;
  double &crispMembership(const int &index1, const int &index2);
  const Matrix correctCrispMembership(void) const;
  double &correctCrispMembership(const int &index1, const int &index2);
  void set_contingencyTable(void);
  const Matrix contingencyTable(void) const;
  const double ARI(void) const;
};

#endif
