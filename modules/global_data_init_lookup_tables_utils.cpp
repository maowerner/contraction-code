#include "global_data.h"
#include "global_data_utils.h"

namespace {

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// helper functions
static void copy_quantum_numbers(const pdg& in, std::array<int, 6>& out){

  out[0] = in.dis3[0];
  out[1] = in.dis3[1];
  out[2] = in.dis3[2];
  out[3] = in.p3[0];
  out[4] = in.p3[1]; 
  out[5] = in.p3[2];
}

} // end of unnamed namespace

namespace global_data_utils {

// helper functions to execute vector addition with arrays. 
// TODO: use std::vector?
std::array<int, 3> add_p3(const pdg& in1, const pdg& in2){

  std::array<int, 3> result;
  for(size_t i = 0; i < 3; i++){
    result[i] = in1.p3[i] + in2.p3[i];
  }

  return result;
}

// helper functions to obtain absolute values of the momentum of a certain
// lookup_corr number.
int abs_p3(const pdg& in){

  int result = 0;
  for(size_t i = 0; i < 3; i++){
    result += in.p3[i] * in.p3[i];
  }

  return result;
}

// helper function to obtain absolute value of array like for vectors
// TODO: use std::vector?
int abs_p3(const std::array<int, 3> in){

  int result = 0;
  for(size_t i = 0; i < 3; i++){
    result += in[i] * in[i];
  }

  return result;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// function that compares two pdg structs and checks if the corresponding 
// entries of lookup_corr coincide
bool compare_quantum_numbers_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3) && 
      (in1.gamma == in2.gamma))
    return true;
  else
    return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// function that compares two pdg structs and checks if the corresponding 
// momenta and displacements coincide
bool compare_mom_dis_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3))
    return true;
  else
    return false;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// function that compares a pdg struct with an operator as defined by the input
// file and checkes if the quantum numbers of pdg are contained in the physical
// situation described by the operator
bool compare_quantum_numbers_of_pdg(const pdg& in1, const Operators& in2){

  if( in1.gamma == in2.gammas &&
      in1.dis3 == in2.dil_vec){
    for(const auto& mom_vec : in2.mom_vec){
      for(auto mom : mom_vec){
        if(in1.p3 == mom){
          return true;
          //TODO: is that safer or just spam?
          break;
        }
      }
    }
  }

  return false;
}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// comparision function for lists of size_t
bool compare_index_list(index_IO_1& in1, index_IO_1& in2) {

  if(in1.index_pt.size() != in2.index_pt.size())
    return false;

  auto it1 = in1.index_pt.begin();
  auto it2 = in2.index_pt.begin();
  while(it1 != in1.index_pt.end()){
    if(*it1 != *it2){
      return false;
      break;
    }
    it1++;
  }

  return true;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// comparision function for lists of two size_t's
bool compare_index_list(index_IO_2& in1, index_IO_2& in2) {

  if(in1.index_pt.size() != in2.index_pt.size())
    return false;

  auto it1 = in1.index_pt.begin();
  auto it2 = in2.index_pt.begin();
  while(it1 != in1.index_pt.end()){
    if(*it1 != *it2){
      return false;
      break;
    }
    it1++;
  }

  return true;

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// function that sets the indices of lookup_vdv and lookup_rvdvr. These are 
// different from lookup_corr as V^dagger*V is independent of the gamma 
// structure it thus is not necessary to calculate this quantity for every 
// operator. Additionally momenta related by a factor (-1) can be obtained by
// daggering V^dagger*V. This function truncates lookup_corr to only contain
// the necessary indices.
void set_index_corr(vec_pdg_Corr& lookup_corr, vec_pd_VdaggerV& lookup_vdv,
                    vec_pd_rVdaggerVr& lookup_rvdvr) {

  // first number is the operator id, the next three is the displacement vector
  // and the last three the momentum vector
  std::vector<std::array<int, 6> > rvdaggervr_qu_nb;
  std::vector<std::array<int, 6> > vdaggerv_qu_nb;
  size_t counter_rvdvr = 0;
  size_t counter_vdv = 0;
  for(auto& op : lookup_corr){
    std::array<int, 6> write;
    if(op.id != 0){
      copy_quantum_numbers(op, write);
      // ######################################################################
      // check if quantum numbers are already stored in rvdaggervr_qu_nb
      bool is_known_rvdvr = false;
      size_t fast_counter_rvdvr = 0;// this gives the Op id if QN are duplicate
      for(const auto& rvdvr : rvdaggervr_qu_nb){
        if(rvdvr == write){
          is_known_rvdvr = true;
          break;
        }
        fast_counter_rvdvr++;
      }
      if(!is_known_rvdvr){ // setting the unknown quantum numbers
        op.id_rvdvr = counter_rvdvr;
        counter_rvdvr++;
        rvdaggervr_qu_nb.push_back(write);
      }
      else
        op.id_rvdvr = fast_counter_rvdvr;

      // ######################################################################
      // check if quantum numbers are already stored in vdaggerv_qu_nb
      bool is_known_vdv = false;
      op.first_vdv = true;
      op.negative_momentum = false;
      size_t fast_counter_vdv = 0;// this gives the Op id if QN are duplicate

      // first check for duplicate quantum numbers
      for(const auto& vdv : vdaggerv_qu_nb){
        if(vdv == write){
          is_known_vdv = true;
          break;
        }
        fast_counter_vdv++;
      }

      // second check for complex conjugate momenta
      if(!is_known_vdv){ 
        fast_counter_vdv = 0;
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;
        for(const auto& vdv : vdaggerv_qu_nb){
          if(vdv == write){
            is_known_vdv = true;
            break;
          }
          fast_counter_vdv++;
        }
        // undo the signchange of momentum
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;

        // case quantum numbers are unknown -> create new entry for 
        // vdaggerv_qu_nb
        if(!is_known_vdv){
          op.id_vdv = counter_vdv;
          vdaggerv_qu_nb.push_back(write);
          counter_vdv++;
//          op.first_vdv = true;
        }
        // case quantum numbers are unknown, but opposite momente exist. 
        // VdaggerV can be obtained from the negative momentum.
        else{
          for(const auto& op2 : lookup_corr){
            if( (op2.negative_momentum == true) && 
                (op2.id_vdv == fast_counter_vdv) ){
              op.first_vdv = false;
              break;
            }
          }

          op.negative_momentum = true;
          op.id_vdv = fast_counter_vdv;
        }
      }
      // case same quantum numbers already exist. 
      // TODO: I don't think that works for displacements
      else{
        op.negative_momentum = lookup_corr[op.id-1].negative_momentum;
        op.id_vdv = fast_counter_vdv;
        op.first_vdv = false;
      }
    }
    else{ // setting the very first entry
      copy_quantum_numbers(op, write);
      rvdaggervr_qu_nb.push_back(write);
      vdaggerv_qu_nb.push_back(write);
      op.id_vdv = counter_vdv;
      op.id_rvdvr = counter_rvdvr;
      op.first_vdv = true;
      counter_vdv++;
      counter_rvdvr++;
    }
  }

  // setting the lookuptables to be able to reconstruct the quantum numbers
  // when computing VdaggerV and rVdaggerVr
  lookup_vdv.resize(vdaggerv_qu_nb.size());
  lookup_rvdvr.resize(rvdaggervr_qu_nb.size());

  // set identification number of lookup_vdv and corresponding index of 
  // lookup_corr to be able to relate to the quantum numbers
  size_t counter = 0;
  for(auto& op_vdv : lookup_vdv){
    op_vdv.id = counter;
    for(const auto& op : lookup_corr){
      if( (counter == op.id_vdv) && (op.first_vdv == true) && 
          (op.negative_momentum == false) )
        op_vdv.index = op.id;
    }
    counter++;
  }

  // test output
  std::cout << "lookup_vdv" << std::endl;
  for(const auto& op_vdv : lookup_vdv)
    std::cout << op_vdv.id << "\t" << op_vdv.index << std::endl;

  // set identification number of lookup_rvdvr and corresponding index of 
  // lookup_corr to be able to relate to the quantum numbers. Also set adjoint
  // flag and id in case the opposite momentum is also built and can be 
  // obtained from adjoining
  counter = 0;
  for(auto& op_rvdvr : lookup_rvdvr){
    op_rvdvr.id = counter;
    for(const auto& op : lookup_corr){
      if(counter == op.id_rvdvr){
        op_rvdvr.index = op.id;
        // if the momentum was only build for the negative in VdaggerV, the 
        // adjoint has been taken
        if(op.negative_momentum == true){
          op_rvdvr.adjoint = true;
          op_rvdvr.id_adjoint = lookup_corr[lookup_vdv[op.id_vdv].index].id_rvdvr;
        }
        else{
          op_rvdvr.adjoint = false;
        }
      }
    }
    counter++;
  }

  // test output
  std::cout << "lookup_rvdvr" << std::endl;
  for(const auto& op_rvdvr : lookup_rvdvr)
    std::cout << op_rvdvr.id << "\t" << op_rvdvr.adjoint << "\t" 
              << op_rvdvr.id_adjoint << "\t" << op_rvdvr.index << std::endl;

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// for the 2pt function a list of all quantum numbers specified in the infile 
// for any 2pt correlator is writen in lookup_2pt. This function takes two 
// operators for source and sink and adds all quantum numbers they contain to
// lookup_2pt
void set_index_2pt(const Operators& in1, const Operators& in2, 
                   const vec_pdg_Corr& lookup_corr, vec_index_2pt& lookup_2pt) {

  index_2pt write;
  std::array<int, 3> zero = {{0, 0, 0}};

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
        write.index_Q2[0] = op1.id;
        write.index_Corr[0] = op2.id;
  
        lookup_2pt.push_back(write);
    }} //loops over sink end here
  }} //loops over source end here

}

// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// for the 3pt function a list of all quantum numbers specified in the infile 
// for any 3pt correlator is writen in lookup_3pt. This function takes two 
// operators for source and one for sink and adds all quantum numbers they 
// contain to lookup_3pt
void set_index_3pt(const Operators& in1, const Operators& in2, 
                   const Operators& in3, const vec_pdg_Corr& lookup_corr, 
                   vec_index_3pt& lookup_3pt) {

  std::array<int, 3> zero = {{0, 0, 0}};

  index_3pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      for(const auto& op3 : lookup_corr){
      if(compare_quantum_numbers_of_pdg(op3, in3)){

// HERE THE REFERENCE FRAME IS STILL HARDCODED! ONLY COMMENT IN ONE BLOCK OF IF 
// STATEMENTS FOR THE DESIRED MOVING FRAME AND CHOOSE THE SAME IN 
// set_index_4pt() AND THE INFILE


        // only include quantum numbers, if momenta correspond to CMS system
//        if( (add_p3(op1, op3) == zero) && (abs_p3(op2) == 0) ){
          
        // only include quantum numbers, if momenta correspond to first moving 
        // frame
//        std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//        if( (abs_p3(add_p3(op1, op3)) == 1) && 
//            (op2_p3 == add_p3(op1, op3)) ){

        // only include quantum numbers, if momenta correspond to second moving
        // frame
//        std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//        if( (((abs_p3(op1) == 2) && (abs_p3(op3) == 0)) || 
//            ((abs_p3(op1) == 0) && (abs_p3(op3) == 2))) &&
//            (op2_p3 == add_p3(op1, op3)) ){

        // only include quantum numbers, if momenta correspond to third moving
        // frame
//        std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//        if( (((abs_p3(op1) == 3) && (abs_p3(op3) == 0)) || 
//            ((abs_p3(op1) == 0) && (abs_p3(op3) == 3))) &&
//            (op2_p3 == add_p3(op1, op3)) ){

        // only include quantum numbers, if momenta correspond to fourth
        // moving frame
        std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
        if( (((abs_p3(op1) == 4) && (abs_p3(op3) == 0)) || 
            ((abs_p3(op1) == 0) && (abs_p3(op3) == 4))) &&
            (op2_p3 == add_p3(op1, op3)) ){

          write.index_Q2[0] = op1.id;
          write.index_Corr[0]  = op2.id;
          write.index_Q2[1] = op3.id;

          lookup_3pt.push_back(write);
        }

      }} //loops over sink end here
    }}
  }} //loops over source end here

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************

// for the 4pt function a list of all quantum numbers specified in the infile 
// for any 4pt correlator is writen in lookup_4pt. This function takes four
// operators for source and sink and adds all quantum numbers they contain to
// lookup_4pt
void set_index_4pt(const Operators& in1, const Operators& in2, 
                   const Operators& in3, const Operators& in4,
                   const vec_pdg_Corr& lookup_corr, vec_index_4pt& lookup_4pt) {

  std::array<int, 3> zero = {{0, 0, 0}};

  index_4pt write;

  for(const auto& op1 : lookup_corr){
  if(compare_quantum_numbers_of_pdg(op1, in1)){
    for(const auto& op2 : lookup_corr){
    if(compare_quantum_numbers_of_pdg(op2, in2)){
      for(const auto& op3 : lookup_corr){
      if(compare_quantum_numbers_of_pdg(op3, in3)){
        for(const auto& op4 : lookup_corr){
        if(compare_quantum_numbers_of_pdg(op4, in4)){

// HERE THE REFERENCE FRAME IS STILL HARDCODED! ONLY COMMENT IN ONE BLOCK OF IF 
// STATEMENTS FOR THE DESIRED MOVING FRAME AND CHOOSE THE SAME IN 
// set_index_3pt() AND THE INFILE

          // only include quantum numbers, if momenta correspond to CMS system
//          if( (add_p3(op1, op3) == zero) && (add_p3(op2, op4) == zero) ){

          // only include quantum numbers, if momenta correspond to first moving 
          // frame
//          std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//          std::array<int, 3> op4_p3 = {{(-1) *op4.p3[0], (-1)*op4.p3[1], (-1)*op4.p3[2]}};
//          if( (abs_p3(add_p3(op1, op3)) == 1) && (abs_p3(add_p3(op2, op4)) == 1)
//              && ( (add_p3(op1, op3) == op2_p3) || (add_p3(op1, op3) == op4_p3)) ){

          // only include quantum numbers, if momenta correspond to second moving 
          // frame
//          std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//          std::array<int, 3> op4_p3 = {{(-1) *op4.p3[0], (-1)*op4.p3[1], (-1)*op4.p3[2]}};
//          if( (((abs_p3(op1) == 2) && (abs_p3(op3) == 0)) || 
//              ((abs_p3(op1) == 0) && (abs_p3(op3) == 2))) &&
//              (((abs_p3(op2) == 2) && (abs_p3(op4) == 0)) || 
//              ((abs_p3(op2) == 0) && (abs_p3(op4) == 2)))
//              && ( (add_p3(op1, op3) == op2_p3) || (add_p3(op1, op3) == op4_p3)) ){

          // only include quantum numbers, if momenta correspond to third moving 
          // frame
//          std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
//          std::array<int, 3> op4_p3 = {{(-1) *op4.p3[0], (-1)*op4.p3[1], (-1)*op4.p3[2]}};
//          if( (((abs_p3(op1) == 3) && (abs_p3(op3) == 0)) || 
//              ((abs_p3(op1) == 0) && (abs_p3(op3) == 3))) &&
//              (((abs_p3(op2) == 3) && (abs_p3(op4) == 0)) || 
//              ((abs_p3(op2) == 0) && (abs_p3(op4) == 3)))
//              && ( (add_p3(op1, op3) == op2_p3) || (add_p3(op1, op3) == op4_p3)) ){

          // only include quantum numbers, if momenta correspond to fourth 
          // moving frame
          std::array<int, 3> op2_p3 = {{(-1) *op2.p3[0], (-1)*op2.p3[1], (-1)*op2.p3[2]}};
          std::array<int, 3> op4_p3 = {{(-1) *op4.p3[0], (-1)*op4.p3[1], (-1)*op4.p3[2]}};
          if( (((abs_p3(op1) == 4) && (abs_p3(op3) == 0)) || 
              ((abs_p3(op1) == 0) && (abs_p3(op3) == 4))) &&
              (((abs_p3(op2) == 4) && (abs_p3(op4) == 0)) || 
              ((abs_p3(op2) == 0) && (abs_p3(op4) == 4)))
              && ( (add_p3(op1, op3) == op2_p3) || (add_p3(op1, op3) == op4_p3)) ){

            write.index_Q2[0]   = op1.id;
            write.index_Corr[0] = op2.id;
            write.index_Q2[1]   = op3.id;
            write.index_Corr[1] = op4.id;

            lookup_4pt.push_back(write);
          }

        }}
      }} //loops over sink end here
    }}
  }} //loops over source end here

}


} // end of namespace global_data_utils
