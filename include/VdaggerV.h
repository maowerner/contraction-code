#ifndef VDAGGERV_H_
#define VDAGGERV_H_

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "boost/multi_array.hpp"
#include "Eigen/Dense"

#include "EigenVector.h"
#include "global_data.h"
#include "RandomVector.h"
#include "typedefs.h"

namespace LapH {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class VdaggerV {

private:
  array_Xcd_d2_eigen vdaggerv;
  array_Xcd_d3_eigen rvdaggerv;
  array_Xcd_d3_eigen vdaggervr;
  array_Xcd_d4_eigen rvdaggervr;
  array_cd_d2  momentum;
  bool is_vdaggerv_set;
  void create_momenta();

public:
  VdaggerV ();
  ~VdaggerV () {};

  void build_vdaggerv(const int config_i);
  void build_rvdaggervr(const int config_i, 
                        const std::vector<LapH::RandomVector>& rnd_vec);

  // return reference on vdaggerv
  inline const Eigen::MatrixXcd& return_vdaggerv(const size_t index,
                                                 const size_t t) const {
    return vdaggerv[index][t];
  }
  // return reference on rvdaggerv
  inline const Eigen::MatrixXcd& return_rvdaggerv(const size_t index,
                                                 const size_t t,
                                                 const size_t rnd1) const {
    return rvdaggerv[index][t][rnd1];
  }
  // return reference on vdaggervr
  inline const Eigen::MatrixXcd& return_vdaggervr(const size_t index,
                                                 const size_t t,
                                                 const size_t rnd1) const {
    return vdaggervr[index][t][rnd1];
  }
  // return reference on rvdaggervr
  inline const Eigen::MatrixXcd& return_rvdaggervr(const size_t index, 
                                                   const size_t t, 
                                                   const size_t rnd1,
                                                   const size_t rnd2) const {
    return rvdaggervr[index][t][rnd1][rnd2];
  }

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif // VDAGGERV_H_ 


