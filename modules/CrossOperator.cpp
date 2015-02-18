#include "CrossOperator.h"  

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

LapH::CrossOperator::CrossOperator(const size_t number) : X(number) {

  const vec_index_IO_1 op_C4 = global_data->get_lookup_4pt_3_IO();
  const size_t nb_op = op_C4.size();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t nb_rnd = quarks[0].number_of_rnd_vec;
  const size_t dilE = quarks[0].number_of_dilution_E;

  for(auto& xx : X){
    xx.resize(boost::extents[nb_op][nb_rnd][nb_rnd][nb_rnd]);

    std::fill(xx.data(), xx.data() + xx.num_elements(), 
                Eigen::MatrixXcd(4 * dilE, 4 * dilE));
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::construct(const BasicOperator& basic, 
                                    const VdaggerV& vdaggerv, const size_t nb,
                                    const int t_source, const int t_sink,
                                    const size_t type){

  const int Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;
  const size_t dilT = quarks[0].number_of_dilution_T;

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();

  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();
  
  if(not ((nb == 0) || (nb == 2) || (nb == 3) || (nb == 4))){
    std::cout << "\tIn LapH::CrossOperator::construct: nb must be 0 or 1" 
              << std::endl;
    exit(0);
  }

  //TODO implement cases 0, 2, 3, 4 in correspondence to type of Q2
  size_t tu, td, t2;
  // case second visited sink index is larger than first one
  if(type == 0){
    tu = (t_sink/dilT);
    if (tu == (((t_sink+1)%Lt)/dilT))
      td = 1;
    else
      td = 2;
    t2 = (t_sink + 1)%Lt;
  }
  // case second visited sink index is smaller than first one
  else if(type == 2){
    t2 = (t_sink + 1)%Lt;
    tu = (t2/dilT);
    if (tu == t_sink/dilT)
      td = 1;
    else
      td = 0;
    t2 = t_sink;
  }
  
  else if(type == 3){
    tu = t_sink/dilT;
    td = 3;
    t2 = (t_source + Lt - 1)%Lt;
  }
  else {
    tu = t_sink/dilT;
    td = 4;
    t2 = (t_source + 1)%Lt;
  }

#pragma omp parallel
#pragma omp single
{     
  for(const auto& op : op_C4_IO){
  for(const auto& i : op.index_pt){

//    size_t id_Q2 = op_C4[i].index_Q2[nb];
    size_t id_Q2 = op_C4[i].index_Q2[1];
//    size_t id_Corr = op_C4[i].index_Corr[nb];
    size_t id_Corr = op_C4[i].index_Q2[0];

    #pragma omp task shared (op)
    for(auto& rnd_it : rnd_vec_index) {
      compute_X(basic, id_Corr, 
                basic.get_operator(t_source, tu, td, id_Q2, rnd_it[0], rnd_it[1]),
                vdaggerv.return_rvdaggervr(op_Corr[id_Corr].id_rvdvr, t2, 
                    rnd_it[1], rnd_it[2]),
                X[nb][op.id][rnd_it[0]][rnd_it[1]][rnd_it[2]]);

    } // loop over random vectors
  }}// loops operators
}// end of parallel region

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::compute_X(const BasicOperator& basic, 
                                    const size_t id_Corr, 
                                    const Eigen::MatrixXcd& Q2, 
                                    const Eigen::MatrixXcd& rVdaggerVr,
                                    Eigen::MatrixXcd& X) {

  const std::vector<quark> quarks = global_data->get_quarks();
  const size_t dilE = quarks[0].number_of_dilution_E;

  X = Eigen::MatrixXcd::Zero(4*dilE, 4*dilE);

  for(size_t block = 0; block < 4; block++){

    cmplx value = 1;
    basic.value_dirac(id_Corr, block, value);

//    for(size_t row = 0; row < 4; row++){
//    for(size_t col = 0; col < 4; col++){
//
//      X.block(row*dilE, col*dilE, dilE, dilE) +=
//        value * Q2.block(row*dilE, basic.order_dirac(id_Corr, block)*dilE, dilE, dilE) * 
//        rVdaggerVr.block(basic.order_dirac(id_Corr, block)*dilE, 
////        rVdaggerVr.block(block*dilE, 
//            basic.order_dirac(id_Corr, col)*dilE, dilE, dilE);
//
//    }}

    X.block(0, block*dilE, 4*dilE, dilE) = 
      value * Q2.block(0, block*dilE, 4*dilE, dilE) * 
      rVdaggerVr.block(block*dilE,
basic.order_dirac(id_Corr, block)*dilE, 
          dilE, dilE);

    }// loop block ends here

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void LapH::CrossOperator::swap(const size_t nb1, const size_t nb2){
  
  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();
  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();
  // TODO: }
  // TODO: Think about for each loop

  // omp parallel for cannot handle autoloops,
  // therefore this workaround is implemented
  #pragma omp parallel
  #pragma omp single
  {
    for(auto& op : op_C4_IO){
      #pragma omp task shared(op)
      for(auto& rnd_it : rnd_vec_index) {
        X[nb1][op.id][rnd_it[0]][rnd_it[1]][rnd_it[2]].swap(
        X[nb2][op.id][rnd_it[0]][rnd_it[1]][rnd_it[2]]);

      } // random vector loop
    }//loops operator
  }
}
