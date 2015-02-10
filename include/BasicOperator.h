#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include "global_data.h"
#include "Perambulator.h"
#include "typedefs.h"
#include "VdaggerV.h"

// struct for Look-up table in create_gamma and get_operator. To read as
// "in column i the row[i]-element is non-zero and its value is value[i]"
// As Gamma matrices are 4x4 matrices, row and value are 4-vectors

struct lookup {
  int row[4];
  std::complex<double> value[4];
};

class BasicOperator {

public:
  BasicOperator();
  ~BasicOperator () {};

  void init_operator(const char dilution, 
                     const LapH::VdaggerV& vdaggerv,
                     const LapH::Perambulator& peram);

  void init_operator_u(const char dilution, 
                     const LapH::VdaggerV& vdaggerv,
                     const LapH::Perambulator& peram);

  void init_operator_d(const char dilution, 
                     const LapH::VdaggerV& vdaggerv,
                     const LapH::Perambulator& peram);


  // returns (P^(b) rho_i)^dagger D_d^-1 V^dagger V Gamma D_u^-1 (P^(b) rho_j)
  inline const Eigen::MatrixXcd& get_operator(const int t1, const int t2,
                                 const int t3, const size_t index,
                                 const size_t rnd_i, const size_t rnd_j) const {
    return Q2[t1][t2][t3][index][rnd_i][rnd_j];
  }

  // returns (P^(b) rho_i)^dagger V^dagger V Gamma D_u^-1 (P^(b) rho_j)
  inline const Eigen::MatrixXcd& get_operator_u(const int t1, const int t2,
                                 const size_t index,
                                 const size_t rnd_i, const size_t rnd_j) const {
    return Q1_u[t1][t2][index][rnd_i][rnd_j];
  }

   // returns (P^(b) rho_i)^dagger V^dagger V Gamma D_u^-1 (P^(b) rho_j)
  inline const Eigen::MatrixXcd& get_operator_d(const int t1, const int t2,
                                 const size_t index,
                                 const size_t rnd_i, const size_t rnd_j) const {
    return Q1_d[t1][t2][index][rnd_i][rnd_j];
  }


  void mult_dirac(const Eigen::MatrixXcd& matrix, Eigen::MatrixXcd& reordered,
                  const size_t index) const;

  size_t order_dirac(const size_t index, const size_t block) const;
  void value_dirac(const size_t index, const size_t block, 
                   cmplx& value) const;

  std::vector<struct lookup>  gamma;

private:
  array_Xcd_d6_eigen Q2;
  array_Xcd_d5_eigen Q1_u;
  array_Xcd_d5_eigen Q1_d;

};

#endif /* BASICOPERATOR_H_ */
