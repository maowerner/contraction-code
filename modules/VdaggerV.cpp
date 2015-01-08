#include "VdaggerV.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
namespace {

inline void read_eigenvectors_from_file(LapH::EigenVector& V, 
                                        const int config_i, const int t) {
  char name[200];
  std::string filename = global_data->get_path_eigenvectors() + "/"
      + global_data->get_name_eigenvectors();
  sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t);

  V.read_eigen_vector(name, 0, 0);
}

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
LapH::VdaggerV::VdaggerV() : vdaggerv(), rvdaggervr(), momentum(), nb_mom(1),
                             is_vdaggerv_set(false) {

  // is needed in the whole class
  nb_mom = global_data->get_number_of_momenta();

  const size_t Lt = global_data->get_Lt();
  const size_t Vs = global_data->get_Lx() * global_data->get_Ly() * 
                 global_data->get_Lz();         
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;

  const size_t nb_VdaggerV = global_data->get_number_of_VdaggerV();
  const size_t nb_rVdaggerVr = global_data->get_number_of_rVdaggerVr();

  // only half of the array is stored to save memory. But be carefull, it 
  // must be mapped correctly from outside by addressing the memomentum
  // correctly and daggering
  vdaggerv.resize(boost::extents[nb_VdaggerV][Lt]);
  rvdaggervr.resize(boost::extents[nb_rVdaggerVr][Lt][nb_rnd][nb_rnd]);

  momentum.resize(boost::extents[nb_mom][Vs]);
  create_momenta();

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::create_momenta () {

  const int Lx = global_data->get_Lx();
  const int Ly = global_data->get_Ly();
  const int Lz = global_data->get_Lz();

  static const std::complex<double> I(0.0, 1.0);

  //const int number_of_max_mom = global_data->get_number_of_max_mom();
  const int max_mom_in_one_dir = global_data->get_max_mom_in_one_dir();

  // helper variables for momenta
  const double px = 2. * M_PI / (double) Lx;
  const double py = 2. * M_PI / (double) Ly;
  const double pz = 2. * M_PI / (double) Lz;
  int p = 0;
  int max_mom_squared = global_data->get_number_of_max_mom();

  // running over all momentum components
  for(int ipx = -max_mom_in_one_dir; ipx <= max_mom_in_one_dir; ++ipx){
    for(int ipy = -max_mom_in_one_dir; ipy <= max_mom_in_one_dir; ++ipy){
      for(int ipz = -max_mom_in_one_dir; ipz <= max_mom_in_one_dir; ++ipz){
        if((ipx * ipx + ipy * ipy + ipz * ipz) > max_mom_squared) {
          continue;
        }
        //TODO: for Lx == Ly == Lz ipxH and ipxHipyH may be integers and px, 
        //py get multiplied in the exponential
        // running over all lattice points
        for(int x = 0; x < Lx; ++x){
          const int xH = x * Ly * Lz; // helper variable
          const double ipxH = ipx * px * x; // helper variable
          for(int y = 0; y < Ly; ++y){
            const int xHyH = xH + y * Lz; // helper variable
            const double ipxHipyH = ipxH + ipy * py * y; // helper variable
            for(int z = 0; z < Lz; ++z){
              momentum[p][xHyH + z] = exp(-I * (ipxHipyH + ipz * pz * z));
            }
          }
        }
        ++p;
      }
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_vdaggerv (const int config_i) {

  clock_t t2 = clock();
  std::cout << "\tbuild vdaggerv:";

  const size_t Lt = global_data->get_Lt();
  const size_t dim_row = global_data->get_dim_row();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const size_t nb_VdaggerV = global_data->get_number_of_VdaggerV();
  const size_t id_unity = global_data->get_index_of_unity();

  const size_t nb_dis = global_data->get_number_of_displ();

  std::fill(vdaggerv.origin(), vdaggerv.origin() + vdaggerv.num_elements(), 
            Eigen::MatrixXcd::Zero(nb_ev, nb_ev));

#pragma omp parallel
{
  Eigen::VectorXcd mom = Eigen::VectorXcd::Zero(dim_row);
  LapH::EigenVector V_t(1, dim_row, nb_ev);// each thread needs its own copy
  #pragma omp for schedule(dynamic)
  for(size_t t = 0; t < Lt; ++t){

    read_eigenvectors_from_file(V_t, config_i, t);

    for(size_t op = 0; op < nb_VdaggerV; op++){
      if(op != id_unity){
        // momentum vector contains exp(-i p x). Divisor 3 for colour index. 
        // All three colours on same lattice site get the same momentum.
        for(size_t x = 0; x < dim_row; ++x) {
          mom(x) = momentum[op/nb_dis][x/3];
        }
        vdaggerv[op][t] = V_t[0].adjoint() * mom.asDiagonal() * V_t[0];

      }
      else{
        // zero momentum 
        (vdaggerv[op][t]) = Eigen::MatrixXcd::Identity(nb_ev, nb_ev);
      }
    }

//    }} // loop over momentum and displacement
  } // loop over time
}// pragma omp parallel ends here

  // set flag that vdaggerv is set
  is_vdaggerv_set = true;

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::VdaggerV::build_rvdaggervr(const int config_i,
                             const std::vector<LapH::RandomVector>& rnd_vec) {
  // check of vdaggerv is already build
  if(not is_vdaggerv_set){
    std::cout << "\n\n\tCaution: vdaggerv is not set and rvdaggervr cannot be" 
              << " computed\n\n" << std::endl;
    exit(0);
  }

  clock_t t2 = clock();
  std::cout << "\tbuild rvdaggervr:";

  const size_t Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t nb_VdaggerV = global_data->get_number_of_VdaggerV();
  const std::list<std::pair<size_t, size_t> > op_rVdaggerVr = 
    global_data->get_op_rVdaggerVr();

  std::fill(rvdaggervr.data(), rvdaggervr.data() + rvdaggervr.num_elements(), 
            Eigen::MatrixXcd::Zero(dilE, 4*dilE));

  // TODO: just a workaround
  // can be changed to op by running over p = op/nb_dg, but dis currently
  // not supported.

  #pragma omp parallel for schedule(dynamic)
  for(size_t t = 0; t < Lt; t++){

  for(size_t op = 0; op < nb_VdaggerV; op++){

    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(nb_ev, 4*dilE);
      // dilution from left
      for(size_t block= 0; block < 4; block++){
      for(size_t vec_i = 0; vec_i < nb_ev; ++vec_i) {
        size_t blk_i =  block + vec_i * 4 + 4 * nb_ev * t;
        
        M.block(0, vec_i%dilE + dilE*block, nb_ev, 1) += 
             vdaggerv[op][t].col(vec_i) * 
             rnd_vec[rnd_i][blk_i];
      }}// end of dilution
      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
      if(rnd_i != rnd_j){
        // dilution from right
        for(size_t block = 0; block < 4; block++){
        for(size_t vec_j = 0; vec_j < nb_ev; ++vec_j) {
          size_t blk_j =  block + vec_j * 4 + 4 * nb_ev * t;
          rvdaggervr[op][t][rnd_j][rnd_i]
                        .block(vec_j%dilE, dilE*block , 1, dilE) +=
              M.block(vec_j, dilE*block, 1, dilE) * 
              std::conj(rnd_vec[rnd_j][blk_j]);
        }}// end of dilution
      }}// rnd_j loop ends here
    }// rnd_i loop ends here
  }

  // building the other half of momenta
  for(const auto& op : op_rVdaggerVr){
    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
    for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j){
    if(rnd_i != rnd_j){
      //is .adjoint().transpose() faster?
      for(size_t block = 0; block < 4; block++){
        rvdaggervr[op.second][t][rnd_j][rnd_i]
                            .block(0, block*dilE, dilE, dilE) =
          (rvdaggervr[op.first][t][rnd_i][rnd_j]
                            .block(0, block*dilE, dilE, dilE)).adjoint();
      }
    }}}// loops over rnd vecs
  }

  }// time, momemtum and displacement loop ends here

  t2 = clock() - t2;
  std::cout << std::setprecision(1) << "\t\tSUCCESS - " << std::fixed 
    << ((float) t2)/CLOCKS_PER_SEC << " seconds" << std::endl;

}





