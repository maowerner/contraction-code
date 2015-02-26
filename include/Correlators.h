#ifndef CORRELATORS_H_
#define CORRELATORS_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "CorrelatorIo2pt.h"
#include "CrossOperator.h" 
#include "global_data.h"
#include "IoHelpers.h"
#include "BasicOperator.h"
#include "Perambulator.h"
#include "typedefs.h"

#include "omp.h"

namespace LapH {

class Correlators{

public:
  Correlators();
  ~Correlators() {};

  void compute_correlators(const size_t config_i);

private:
  BasicOperator basic;
  LapH::Perambulator peram;
  std::vector<LapH::RandomVector> rnd_vec;
  LapH::VdaggerV vdaggerv;
  array_cd_d2 C4_mes;
  array_cd_d2 C4_Dd_11_mes;
  array_cd_d2 C4_Dd_12_mes;
  array_cd_d2 C4_Dd_21_mes;
  array_cd_d2 C4_Dd_22_mes;
  array_cd_d2 C4_Du_11_mes;
  array_cd_d2 C4_Du_12_mes;
  array_cd_d2 C4_Du_21_mes;
  array_cd_d2 C4_Du_22_mes;
  array_cd_d2 C3_A_mes;
  array_cd_d2 C3_B_mes;
  array_cd_d2 C2_mes;
  array_cd_d3 Q1_u_trace;
  array_cd_d3 Q1_d_trace;
  array_cd_d6 Q2_trace;
  array_cd_d6 Q2_trace_uncharged;

  void set_corr(const size_t config){
    read_rnd_vectors_from_file(config);
    vdaggerv.build_vdaggerv(config);
    vdaggerv.build_rvdaggervr(config, rnd_vec);
    peram.read_perambulators_from_file(config);
  }
  void read_rnd_vectors_from_file (const int config_i);
  void compute_meson_small_traces(const size_t id_si, 
                                  const Eigen::MatrixXcd& Q2,
                                  const Eigen::MatrixXcd& rVdaggerVr, 
                                  cmplx& Q2_trace);

  void build_Q1_trace();
  void build_Q2_trace();
  void build_Q2_trace_uncharged();

  void build_and_write_2pt(const size_t config_i);
  void build_and_write_C4_1(const size_t config_i);
  void build_and_write_C4_2(const size_t config_i);

  void compute_meson_3pt_cross_trace(LapH::CrossOperator& X);
  void compute_meson_3pt_cross_trace_verbose(LapH::CrossOperator& X);
  void compute_meson_4pt_cross_trace(LapH::CrossOperator& X);
  void compute_meson_4pt_box_trace(LapH::CrossOperator& X);

  void write_C3(const size_t config_i);
  void write_C3_verbose(const size_t config_i);
  void write_C4_cross(const size_t config_i);
  void write_C4_box(const size_t config_i);

};

}// end of namespace

#endif // CORRELATORS_H_
