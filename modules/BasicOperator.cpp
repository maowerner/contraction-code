/*
 * BasicOperator.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#include "BasicOperator.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

namespace { // some internal namespace

static const std::complex<double> I(0.0, 1.0);

// Look-up table for gamma matrices. For every Gamma structure (currently 0-15)
// the four non-zero values are specified.

static void create_gamma (std::vector<struct lookup>& gamma, const int i) {
  try {
    switch(i) {
    case 0: // gamma_0
      gamma[0].row[0] = 2;
      gamma[0].value[0] = -1;
      gamma[0].row[1] = 3;
      gamma[0].value[1] = -1;
      gamma[0].row[2] = 0;
      gamma[0].value[2] = -1;
      gamma[0].row[3] = 1;
      gamma[0].value[3] = -1;
      break;

    case 1: // gamma_1
      gamma[1].row[0] = 3;
      gamma[1].value[0] = I;
      gamma[1].row[1] = 2;
      gamma[1].value[1] = I;
      gamma[1].row[2] = 1;
      gamma[1].value[2] = -I;
      gamma[1].row[3] = 0;
      gamma[1].value[3] = -I;
      break;

    case 2: // gamma_2
      gamma[2].row[0] = 3;
      gamma[2].value[0] = -1;
      gamma[2].row[1] = 2;
      gamma[2].value[1] = 1;
      gamma[2].row[2] = 1;
      gamma[2].value[2] = 1;
      gamma[2].row[3] = 0;
      gamma[2].value[3] = -1;
      break;

    case 3: // gamma_3
      gamma[3].row[0] = 2;
      gamma[3].value[0] = I;
      gamma[3].row[1] = 3;
      gamma[3].value[1] = -I;
      gamma[3].row[2] = 0;
      gamma[3].value[2] = -I;
      gamma[3].row[3] = 1;
      gamma[3].value[3] = I;
      break;

    case 4: // unity
      gamma[4].row[0] = 0;
      gamma[4].value[0] = 1;
      gamma[4].row[1] = 1;
      gamma[4].value[1] = 1;
      gamma[4].row[2] = 2;
      gamma[4].value[2] = 1;
      gamma[4].row[3] = 3;
      gamma[4].value[3] = 1;
      break;

    case 5: // gamma_5
      gamma[5].row[0] = 0;
      gamma[5].value[0] = 1;
      gamma[5].row[1] = 1;
      gamma[5].value[1] = 1;
      gamma[5].row[2] = 2;
      gamma[5].value[2] = -1;
      gamma[5].row[3] = 3;
      gamma[5].value[3] = -1;
      break;

    case 6: // gamma_0 * gamma_5
      gamma[6].row[0] = 2;
      gamma[6].value[0] = -1;
      gamma[6].row[1] = 3;
      gamma[6].value[1] = -1;
      gamma[6].row[2] = 0;
      gamma[6].value[2] = 1;
      gamma[6].row[3] = 1;
      gamma[6].value[3] = 1;
      break;

    case 7: // gamma_1 * gamma_5
      gamma[7].row[0] = 3;
      gamma[7].value[0] = I;
      gamma[7].row[1] = 2;
      gamma[7].value[1] = I;
      gamma[7].row[2] = 1;
      gamma[7].value[2] = I;
      gamma[7].row[3] = 0;
      gamma[7].value[3] = I;
      break;

    case 8: // gamma_2 * gamma_5
      gamma[8].row[0] = 3;
      gamma[8].value[0] = -1;
      gamma[8].row[1] = 2;
      gamma[8].value[1] = 1;
      gamma[8].row[2] = 1;
      gamma[8].value[2] = -1;
      gamma[8].row[3] = 0;
      gamma[8].value[3] = 1;
      break;

    case 9: // gamma_3 * gamma_5
      gamma[9].row[0] = 2;
      gamma[9].value[0] = I;
      gamma[9].row[1] = 3;
      gamma[9].value[1] = -I;
      gamma[9].row[2] = 0;
      gamma[9].value[2] = I;
      gamma[9].row[3] = 1;
      gamma[9].value[3] = -I;
      break;

    case 10: // gamma_0 * gamma_1
      gamma[10].row[0] = 1;
      gamma[10].value[0] = -I;
      gamma[10].row[1] = 0;
      gamma[10].value[1] = -I;
      gamma[10].row[2] = 3;
      gamma[10].value[2] = I;
      gamma[10].row[3] = 2;
      gamma[10].value[3] = I;
      break;

    case 11: // gamma_0 * gamma_2
      gamma[11].row[0] = 1;
      gamma[11].value[0] = 1;
      gamma[11].row[1] = 0;
      gamma[11].value[1] = -1;
      gamma[11].row[2] = 3;
      gamma[11].value[2] = -1;
      gamma[11].row[3] = 2;
      gamma[11].value[3] = 1;
      break;

    case 12: // gamma_0 * gamma_3
      gamma[12].row[0] = 0;
      gamma[12].value[0] = -I;
      gamma[12].row[1] = 1;
      gamma[12].value[1] = I;
      gamma[12].row[2] = 2;
      gamma[12].value[2] = I;
      gamma[12].row[3] = 3;
      gamma[12].value[3] = -I;
      break;

    case 13: // gamma_1 * gamma_2
      gamma[13].row[0] = 0;
      gamma[13].value[0] = I;
      gamma[13].row[1] = 1;
      gamma[13].value[1] = -I;
      gamma[13].row[2] = 2;
      gamma[13].value[2] = I;
      gamma[13].row[3] = 3;
      gamma[13].value[3] = -I;
      break;

    case 14: // gamma_1 * gamma_3
      gamma[14].row[0] = 1;
      gamma[14].value[0] = 1;
      gamma[14].row[1] = 0;
      gamma[14].value[1] = -1;
      gamma[14].row[2] = 3;
      gamma[14].value[2] = 1;
      gamma[14].row[3] = 2;
      gamma[14].value[3] = -1;
      break;

    case 15: // gamma_2 * gamma_3
      gamma[15].row[0] = 1;
      gamma[15].value[0] = I;
      gamma[15].row[1] = 0;
      gamma[15].value[1] = I;
      gamma[15].row[2] = 3;
      gamma[15].value[2] = I;
      gamma[15].row[3] = 2;
      gamma[15].value[3] = I;
      break;

    case 16: // gamma_2 * gamma_0 * gamma_5
      gamma[16].row[0] = 1;
      gamma[16].value[0] = -1;
      gamma[16].row[1] = 0;
      gamma[16].value[1] = 1;
      gamma[16].row[2] = 3;
      gamma[16].value[2] = -1;
      gamma[16].row[3] = 2;
      gamma[16].value[3] = 1;
      break;
    default:
      printf("Dirac component %d not found in BasicOperator::create_gamma\n",i);
      exit(0);
    }
  return;
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: BasicOperator::create_gamma\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
} // internal namespace ends here

/******************************************************************************/
// constructor ****************************************************************/
/******************************************************************************/
BasicOperator::BasicOperator() : gamma(), Q2() {

  // creating gamma matrices
  gamma.resize(16);
  for(int i = 0; i < 16; ++i)
    create_gamma(gamma, i);

  std::cout << "\tallocated memory in BasicOperator" << std::endl;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::reset_operator() {

  //TODO Q2.clear(); ??
  Q2.resize(boost::extents[0][0][0][0][0][1]);
  Q1_u.resize(boost::extents[0][0][0][0][1]);
  Q1_d.resize(boost::extents[0][0][0][0][1]);

}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::init_operator_verbose(
                              const std::map<size_t, size_t> map_required_Q2,
                              const std::map<size_t, size_t> map_required_times,
                              const LapH::VdaggerV& vdaggerv, 
                              const LapH::Perambulator& peram){

  const int Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const int dilT = quarks[0].number_of_dilution_T;
  const size_t Q2_size = 4 * dilE;
  
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  const size_t nb_op = map_required_Q2.size();
  const size_t nb_ti = map_required_times.size();

  std::cout << "\n" << std::endl;
  clock_t time = clock();

  Q2.resize(boost::extents[Lt][Lt/dilT][nb_ti][nb_op][nb_rnd][nb_rnd]);
  std::fill(Q2.data(), Q2.data() + Q2.num_elements(), 
            Eigen::MatrixXcd::Zero(Q2_size, Q2_size));

  #pragma omp parallel 
  {
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(Q2_size, 4 * nb_ev);
  #pragma omp for schedule(dynamic)
  for(int t_0 = 0; t_0 < Lt; t_0++){

    if(omp_get_thread_num() == 0)
      std::cout << "\tcomputing double quarkline: " 
                << std::setprecision(2) << (float) t_0/Lt*100 << "%\r" 
                << std::flush;

    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {

      // ti denotes times of Quarkline:
      // 5 - t_0 -> t_0-1 -> t  I=1   3pt,4pt
      // 6 - t_0 -> t_0+1 -> t  I=1   3pt,4pt
      for(const auto& ti : map_required_times){
      // getting the neighbour blocks
        int tend;
        tend = ((Lt + t_0 + 2*(5-ti.first) + 1)%Lt)/dilT;

        for(const auto& q2 : map_required_Q2){
          // new momentum -> recalculate M[0]
          // M only depends on momentum and displacement. first_vdv 
          // prevents repeated calculation for different gamma structures
          if(op_Corr[q2.first].first_vdv == true){

            for(size_t col = 0; col < 4; ++col) {
            for(size_t row = 0; row < 4; ++row) {
              if(op_Corr[q2.first].negative_momentum == false){
                M.block(dilE * col, nb_ev * row, dilE, nb_ev) = 
                  (peram[rnd_i].block(nb_ev * (4 * t_0 + row), 
                                      dilE * (4 * tend + col), 
                                      nb_ev, dilE)).adjoint() *
                  vdaggerv.return_vdaggerv(op_Corr[q2.first].id_vdv, t_0);
              }
              else {
                M.block(dilE * col, nb_ev * row, dilE, nb_ev) = 
                  (peram[rnd_i].block(nb_ev*(4*t_0 + row), dilE*(4*tend + col), 
                                      nb_ev, dilE)).adjoint() *
                  // TODO: calculate V^daggerV Omega from op.negative_momentum 
                  // == false and multiply Omega from the left
                  // and then (V^daggerV Omega)^dagger * Omega
                  (vdaggerv.return_vdaggerv(op_Corr[q2.first].id_vdv, t_0)).adjoint();
              }  

              // gamma_5 trick. It changes the sign of the two upper right and two
              // lower left blocks in dirac space
              if( ((row + col) == 3) || (abs(row - col) > 1) )
                M.block(dilE * col, row * nb_ev, dilE, nb_ev) *= -1.;
            }}// loops over row and col end here
          }//if over same gamma structure ends here

          for(int t = 0; t < Lt/dilT; t++){
            for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) {
              if(rnd_i != rnd_j){

                //dilution of d-quark from left
                for(size_t block_dil = 0; block_dil < 4; block_dil++){
                  cmplx value = 1.;
                  value_dirac(q2.first, block_dil, value);

                  for(size_t row = 0; row < 4; row++){
                  for(size_t col = 0; col < 4; col++){

                    Q2[t_0][t][ti.second][q2.second][rnd_i][rnd_j]
                        .block(row*dilE, col*dilE, dilE, dilE) += value *
                      M.block(row*dilE, order_dirac(q2.first, block_dil)* nb_ev, dilE, nb_ev) * 
                      peram[rnd_j]
                        .block(4*nb_ev*t_0 + block_dil*nb_ev, 
                          Q2_size*t + col*dilE, nb_ev, dilE);

              }}}//dilution ends here

          }}}// loops over rnd_j and t block end here 
        }//loop operators
      }// loop over ti ends here
    }// loop over rnd_i ends here
  }// loops over t_0 ends here
  }// pragma omp ends

  std::cout << "\tcomputing double quarkline: 100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::init_operator(
                              const std::map<size_t, size_t> map_required_Q2,
                              const std::map<size_t, size_t> map_required_times,
                              const LapH::VdaggerV& vdaggerv, 
                              const LapH::Perambulator& peram){

  const int Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const int dilT = quarks[0].number_of_dilution_T;
  const size_t Q2_size = 4 * dilE;
  
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  const size_t nb_op = map_required_Q2.size();
  const size_t nb_ti = map_required_times.size();

  std::cout << "\n" << std::endl;
  clock_t time = clock();

  Q2.resize(boost::extents[Lt][Lt/dilT][nb_ti][nb_op][nb_rnd][nb_rnd]);
  std::fill(Q2.data(), Q2.data() + Q2.num_elements(), 
            Eigen::MatrixXcd::Zero(Q2_size, Q2_size));

  #pragma omp parallel 
  {
  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(Q2_size, 4 * nb_ev);
  #pragma omp for schedule(dynamic)
  for(int t_0 = 0; t_0 < Lt; t_0++){

    if(omp_get_thread_num() == 0)
      std::cout << "\tcomputing double quarkline: " 
                << std::setprecision(2) << (float) t_0/Lt*100 << "%\r" 
                << std::flush;

    for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(int t = 0; t < Lt/dilT; t++){

        for(const auto& q2 : map_required_Q2){
          // new momentum -> recalculate M[0]
          // M only depends on momentum and displacement. first_vdv 
          // prevents repeated calculation for different gamma structures
          if(op_Corr[q2.first].first_vdv == true){

            for(size_t col = 0; col < 4; ++col) {
            for(size_t row = 0; row < 4; ++row) {
              if(op_Corr[q2.first].negative_momentum == false){
                M.block(dilE * col, nb_ev * row, dilE, nb_ev) = 
                  (peram[rnd_i].block(nb_ev * (4 * t_0 + row), 
                                      dilE * (4 * t + col), 
                                      nb_ev, dilE)).adjoint() *
                  vdaggerv.return_vdaggerv(op_Corr[q2.first].id_vdv, t_0);
              }
              else {
                M.block(dilE * col, nb_ev * row, dilE, nb_ev) = 
                  (peram[rnd_i].block(nb_ev * (4 * t_0 + row), 
                                      dilE * (4 * t + col), 
                                      nb_ev, dilE)).adjoint() *
                  // TODO: calculate V^daggerV Omega from op.negative_momentum 
                  // == false and multiply Omega from the left
                  // and then (V^daggerV Omega)^dagger * Omega
                  (vdaggerv.return_vdaggerv(op_Corr[q2.first].id_vdv, t_0)).adjoint();
              }  

              // gamma_5 trick. It changes the sign of the two upper right and two
              // lower left blocks in dirac space
              if( ((row + col) == 3) || (abs(row - col) > 1) )
                M.block(dilE * col, row * nb_ev, dilE, nb_ev) *= -1.;
            }}// loops over row and col end here
          }//if over same gamma structure ends here

          // ti denotes times of Quarkline:
          // 0 - t -> t_0 -> t-1    I=2   4pt
          // 1 - t -> t_0 -> t      I=1,2 2pt
          // 2 - t -> t_0 -> t+1    I=2   4pt
          // 3 - t -> t_0 -> t_0-1  I=1   3pt,4pt
          // 4 - t -> t_0 -> t_0+1  I=1   3pt,4pt
          for(const auto& ti : map_required_times){
            // getting the neighbour blocks
            int tend;
            if(ti.first < 3)
              tend = (Lt/dilT+t + ti.first - 1)%(Lt/dilT);  
            else
              tend = ((Lt + t_0 + 2*(ti.first-3) - 1)%Lt)/dilT;

            for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) {
              if(rnd_i != rnd_j){

                //dilution of d-quark from left
                for(size_t block_dil = 0; block_dil < 4; block_dil++){
                  cmplx value = 1.;
                  value_dirac(q2.first, block_dil, value);

                  for(size_t row = 0; row < 4; row++){
                  for(size_t col = 0; col < 4; col++){

                    Q2[t_0][t][ti.second][q2.second][rnd_i][rnd_j]
                        .block(row*dilE, col*dilE, dilE, dilE) += value *
                      M.block(row*dilE, order_dirac(q2.first, block_dil)* nb_ev, 
                              dilE, nb_ev) * 
                      peram[rnd_j]
                        .block(4*nb_ev*t_0 + block_dil*nb_ev, 
                               Q2_size*tend + col*dilE, nb_ev, dilE);

              }}}}//dilution ends here

          }}// loops over rnd_j and ti block end here 
        }//loop operators
      }// loop over t ends here
    }// loop over rnd_i ends here
  }// loops over t_0 ends here
  }// pragma omp ends

  std::cout << "\tcomputing double quarkline: 100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::init_operator_u(const std::map<size_t, size_t> map_required_u,
                                    const LapH::VdaggerV& vdaggerv, 
                                    const LapH::Perambulator& peram){

  const int Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const int dilT = quarks[0].number_of_dilution_T;
  const size_t Q2_size = 4 * dilE;
  
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  const size_t nb_op = map_required_u.size();

  std::cout << "\n" << std::endl;
  clock_t time = clock();

  Q1_u.resize(boost::extents[Lt][Lt/dilT][nb_op][nb_rnd][nb_rnd]);
  std::fill(Q1_u.data(), Q1_u.data() + Q1_u.num_elements(), 
            Eigen::MatrixXcd::Zero(Q2_size, Q2_size));

  #pragma omp parallel 
  {
//  Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(Q2_size, 4 * nb_ev);
  #pragma omp for schedule(dynamic)
  for(int t_0 = 0; t_0 < Lt; t_0++){

    if(omp_get_thread_num() == 0)
      std::cout << "\tcomputing single quarkline: " 
                << std::setprecision(2) << (float) t_0/Lt*100 << "%\r" 
                << std::flush;

    for(const auto& u : map_required_u){

      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) {
        for(int t = 0; t < Lt/dilT; t++){
          // new momentum -> recalculate M[0]
          // M only depends on momentum and displacement. first_vdv 
          // prevents repeated calculation for different gamma structures
          
          //dilution of u-quark from left
          for(size_t block_dil = 0; block_dil < 4; block_dil++){
            cmplx value = 1.;
            value_dirac(u.first, block_dil, value);

              Q1_u[t_0][t][u.second][rnd_i][rnd_j].block(order_dirac(u.first, 
                  block_dil)*dilE, 0, dilE, 4*dilE) += value * 
                vdaggerv.return_rvdaggerv(op_Corr[u.first].id_rvdvr, t_0, rnd_i)
                  .block(0, order_dirac(u.first, block_dil)*nb_ev, dilE, nb_ev) * 
                peram[rnd_j].block(4*nb_ev*t_0 + block_dil*nb_ev, Q2_size*t, 
                                   nb_ev, 4*dilE);

          }//dilution ends here
        }// loop over t ends here
     }}// loop over rnd ends here
    }//loop operators
  }// loops over t_0 ends here
  }// pragma omp ends

  std::cout << "\tcomputing single quarkline: 100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::init_operator_d(const std::map<size_t, size_t> map_required_d,
                                    const LapH::VdaggerV& vdaggerv, 
                                    const LapH::Perambulator& peram){

  const int Lt = global_data->get_Lt();
  const size_t nb_ev = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;
  const int dilT = quarks[0].number_of_dilution_T;
  const size_t Q2_size = 4 * dilE;
  
  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  const size_t nb_op = map_required_d.size();

  std::cout << "\n" << std::endl;
  clock_t time = clock();

  Q1_d.resize(boost::extents[Lt/dilT][Lt][nb_op][nb_rnd][nb_rnd]);
  std::fill(Q1_d.data(), Q1_d.data() + Q1_d.num_elements(), 
            Eigen::MatrixXcd::Zero(Q2_size, Q2_size));

  #pragma omp parallel 
  {
  #pragma omp for schedule(dynamic)
  for(int t_0 = 0; t_0 < Lt/dilT; t_0++){

    if(omp_get_thread_num() == 0)
      std::cout << "\tcomputing single quarkline: " 
                << std::setprecision(2) << (float) t_0/Lt*100 << "%\r" 
                << std::flush;

    for(const auto& d : map_required_d){

      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
      for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) {
        for(int t = 0; t < Lt; t++){
          // new momentum -> recalculate M[0]
          // M only depends on momentum and displacement. first_vdv 
          // prevents repeated calculation for different gamma structures
          
          //dilution of u-quark from left
          for(size_t block_dil = 0; block_dil < 4; block_dil++){
            cmplx value = 1.;
            value_dirac(d.first, block_dil, value);

              Q1_d[t_0][t][d.second][rnd_i][rnd_j]
                  .block(0, block_dil*dilE, 4*dilE, dilE) += value * 
                (peram[rnd_i].block(4*nb_ev*t + order_dirac(d.first, block_dil) *
                  nb_ev, Q2_size*t_0 , 
                  nb_ev, 4*dilE)).adjoint() *
                vdaggerv.return_vdaggervr(op_Corr[d.first].id_rvdvr, t, rnd_j).block(0, 
                  block_dil*dilE, nb_ev, dilE);

          }//dilution ends here

          // gamma_5 trick
          for(size_t row = 0; row < 4; row++)
          for(size_t col = 0; col < 4; col++)
            if( ((row + col) == 3) || (abs(row - col) > 1) )
              Q1_d[t_0][t][d.second][rnd_i][rnd_j]
                .block(row*dilE, col*dilE, dilE, dilE) *= -1.;

        }// loop over t ends here
     }}// loop over rnd ends here
    }//loop operators
  }// loops over t_0 ends here
  }// pragma omp ends

  std::cout << "\tcomputing single quarkline: 100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time) / CLOCKS_PER_SEC 
            << " seconds" << std::endl;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//void BasicOperator::swap_operators(){
//
//  const std::vector<quark> quarks = global_data->get_quarks();
//  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
//  const size_t nb_mom = global_data->get_number_of_momenta();
//  #pragma omp parallel for schedule(dynamic)
//  for(size_t p = 0; p < nb_mom; p++){
//    for(size_t dir = 0; dir < 1; ++dir) {// TODO: dirac still hard coded
//      for(size_t rnd_i = 0; rnd_i < nb_rnd; ++rnd_i) {
//        for(size_t rnd_j = 0; rnd_j < nb_rnd; ++rnd_j) { 
//          Q2[0][p][dir][rnd_i][rnd_j].swap(
//          Q2[1][p][dir][rnd_i][rnd_j]);
//        }
//      }
//    }
//  }
//}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
size_t BasicOperator::order_dirac(const size_t index, size_t block) const {

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  size_t block_reordered = block;

  for(const auto& dirac : op_Corr[index].gamma){
    if(dirac != 4){
      block_reordered = gamma[dirac].row[block_reordered];
    }
  }

  return block_reordered;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::value_dirac(const size_t index, const size_t block, 
                                cmplx& value) const{

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  size_t block_reordered = block;

  for(const auto& dirac : op_Corr[index].gamma){
    if(dirac != 4){
      value = value * gamma[dirac].value[block_reordered];
      block_reordered = gamma[dirac].row[block_reordered];
    }
  }

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void BasicOperator::mult_dirac(const Eigen::MatrixXcd& matrix,
                               Eigen::MatrixXcd& reordered, 
                               const size_t index) const {

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();

  const size_t rows = matrix.rows();
  const size_t cols = matrix.cols();
  const size_t col = cols/4;
  if( (rows != reordered.rows()) || (cols != reordered.cols()) ){
    std::cout << "Error! In BasicOperator::mult_dirac: size of matrix and "
                 "reordered must be equal" << std::endl;
    exit(0);
  }

  for(const auto& dirac : op_Corr[index].gamma){
    if(dirac != 4){
      for(size_t block = 0; block < 4; block++){
        reordered.block(0, block * col, rows, col) = 
          gamma[dirac].value[block] *
          matrix.block(0, gamma[dirac].row[block]*col, rows, col);
      }
    }
  }

}


