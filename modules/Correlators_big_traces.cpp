#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C3() ?
void LapH::Correlators::compute_meson_3pt_trace_verbose(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int dilT = quarks[0].number_of_dilution_T;

  const vec_pdg_Corr op_Corr = global_data->get_lookup_corr();
  const vec_index_3pt op_C3 = global_data->get_lookup_3pt_trace();
  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();

  size_t counter_1 = 0;
  size_t counter_2 = 0;
  std::pair<std::map<size_t, size_t>::iterator, bool> ret;
  map_required_Q2.clear();
  map_required_u.clear();
  map_required_times.clear();

  const vec_index_IO_1 op_C3_IO = global_data->get_lookup_3pt_IO();

  std::cout << "\n\tcomputing the 3pt function:\r";
  clock_t time = clock();

  for(const auto& op : op_C3_IO){
    for(const auto& i : op.index_pt){
      for(const auto q2 : op_C3[i].index_Q2){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(q2, counter_1));
        if(ret.second == true)
        // case element did not exist yet
          counter_1++;
      }
      for(const auto corr : op_C3[i].index_Corr){
        ret = map_required_u.insert(std::pair<size_t, size_t>(corr, counter_2));
        if(ret.second == true)
        // case element did not exist yet
          counter_2++;
      }
    }
  }

  // only need quarklines without fiertz rearangement
  map_required_times.insert(std::pair<size_t, size_t>(5, 0));
  map_required_times.insert(std::pair<size_t, size_t>(6, 1));

  basic.reset_operator();
  basic.init_operator_verbose(map_required_Q2, map_required_times, vdaggerv, peram);
  basic.init_operator_u(map_required_u, vdaggerv, peram);

  std::fill(C3_A_mes.data(), C3_A_mes.data() + C3_A_mes.num_elements(), 
            cmplx(0.0, 0.0));
  std::fill(C3_B_mes.data(), C3_B_mes.data() + C3_B_mes.num_elements(), 
            cmplx(0.0, 0.0));

  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the 3pt function: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t = 0; t < Lt; t++){
      const int t_source = (t_sink + 1 + t)%Lt;
      const int t_source_1 = (t_source + 1) % Lt;

      X.construct_3pt(map_required_Q2, map_required_times, basic, vdaggerv, 0, 
                      t_source_1, t_sink, 5);
      X.construct_3pt(map_required_Q2, map_required_times, basic, vdaggerv, 1, 
                      t_source, t_sink, 6);

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C3_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
        for(const auto& i : op.index_pt){

          size_t op_Corr = op_C3[i].index_Corr[0];

          #pragma omp task shared(op, i)
          {
          cmplx priv_C3_A(0.0,0.0);
          cmplx priv_C3_B(0.0,0.0);

          for(const auto& rnd_it : rnd_vec_index) {

            priv_C3_A += (X(0, op.id, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          basic.get_operator_u(t_sink, t_source_1/dilT,
                                               map_required_u[op_Corr], 
                                               rnd_it[2], rnd_it[0])).trace();

            priv_C3_B += (X(1, op.id, rnd_it[0], rnd_it[1], rnd_it[2]) *
                         basic.get_operator_u(t_sink, t_source/dilT, 
                                              map_required_u[op_Corr], 
                                              rnd_it[2], rnd_it[0])).trace();

          }
          #pragma omp critical
          {
            C3_A_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C3_A;
            C3_B_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C3_B;
          }
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the 3pt function: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C3() ?
void LapH::Correlators::compute_meson_3pt_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int dilT = quarks[0].number_of_dilution_T;

  const vec_index_3pt op_C3 = global_data->get_lookup_3pt_trace();
  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();

  size_t counter_1 = 0;
  size_t counter_2 = 0;
  std::pair<std::map<size_t, size_t>::iterator, bool> ret;
  map_required_Q2.clear();
  map_required_d.clear();
  map_required_times.clear();

  const vec_index_IO_1 op_C3_IO = global_data->get_lookup_3pt_IO();

  std::cout << "\n\tcomputing the 3pt function:\r";
  clock_t time = clock();

  for(const auto& op : op_C3_IO){
    for(const auto& i : op.index_pt){
      for(const auto q2 : op_C3[i].index_Q2){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(q2, counter_1));
        if(ret.second == true)
        // case element did not exist yet
          counter_1++;
      }
      for(const auto corr : op_C3[i].index_Corr){
        ret = map_required_d.insert(std::pair<size_t, size_t>(corr, counter_2));
        if(ret.second == true)
        // case element did not exist yet
          counter_2++;
      }
    }
  }

  // only need quarklines without fiertz rearangement
  map_required_times.insert(std::pair<size_t, size_t>(3, 0));
  map_required_times.insert(std::pair<size_t, size_t>(4, 1));

  basic.reset_operator();
  basic.init_operator(map_required_Q2, map_required_times, vdaggerv, peram);
  basic.init_operator_d(map_required_d, vdaggerv, peram);

  // setting the correlation functions to zero
  std::fill(C3_A_mes.data(), C3_A_mes.data() + C3_A_mes.num_elements(), 
            cmplx(0.0, 0.0));
  std::fill(C3_B_mes.data(), C3_B_mes.data() + C3_B_mes.num_elements(), 
            cmplx(0.0, 0.0));

  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the 3pt function: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t = 0; t < Lt; t++){
      const int t_source = (t_sink + 1 + t)%Lt;
      const int t_source_1 = (t_source + 1) % Lt;

      X.construct_3pt(map_required_Q2, map_required_times, basic, vdaggerv, 0, 
                      t_source_1, t_sink, 3);
      X.construct_3pt(map_required_Q2, map_required_times, basic, vdaggerv, 1, 
                      t_source, t_sink, 4);

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C3_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
        for(const auto& i : op.index_pt){

          size_t id_Corr = op_C3[i].index_Corr[0];
          #pragma omp task shared(op, i)
          {
          cmplx priv_C3_A(0.0,0.0);
          cmplx priv_C3_B(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C3_A += (X(0, op.id, rnd_it[0], rnd_it[1], rnd_it[2]) *
                          basic.get_operator_d(t_source/dilT, t_sink,
                                map_required_d[id_Corr], rnd_it[2], rnd_it[0]))
                            .trace();

            priv_C3_B += (X(1, op.id, rnd_it[0], rnd_it[1], rnd_it[2]) *
                         (basic.get_operator_d(t_source_1/dilT, t_sink, 
                                               map_required_d[id_Corr], 
                                               rnd_it[2], rnd_it[0]))).trace();

          }
          #pragma omp critical
          {
            C3_A_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C3_A;
            C3_B_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C3_B;
          }
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the 3pt function: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C4_3() ?
void LapH::Correlators::compute_meson_4pt_box_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();

  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  // TODO: must be changed in GlobalData {
  // TODO: }

  size_t counter = 0;
  std::pair<std::map<size_t, size_t>::iterator, bool> ret;
  map_required_Q2.clear();
  map_required_times.clear();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_c4i10_IO();

  std::cout << "\n\tcomputing the 4pt box-diagram:\r";
  clock_t time = clock();

  for(const auto& op : op_C4_IO){
    for(const auto& i : op.index_pt){
      for(const auto q2 : op_C4[i].index_Q2){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(q2, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
      for(const auto corr : op_C4[i].index_Corr){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(corr, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
    }
  }

  // only need quarklines without fiertz rearangement
  map_required_times.insert(std::pair<size_t, size_t>(3, 0));
  map_required_times.insert(std::pair<size_t, size_t>(4, 1));

  basic.reset_operator();
  basic.init_operator(map_required_Q2, map_required_times, vdaggerv, peram);

  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the 4pt box-diagram: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;

    for(int t = 0; t < Lt; t++){
      const int t_source = (t_sink + 1 + t)%Lt;
      const int t_source_1 = (t_source + 1) % Lt;

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source_1, t_sink_1, 3);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink, t_source, 4);

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
        for(const auto& i : op.index_pt){

//          size_t i_0 = op_C4[i].index_Q2[0];
//          size_t i_1 = op_C4[i].index_Corr[0];
//          size_t i_2 = op_C4[i].index_Q2[1];
//          size_t i_3 = op_C4[i].index_Corr[1];
          // complete diagramm. combine X and Y to four-trace
          // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
          //             D_d^-1(t_sink       | t_source + 1) Gamma 
          //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
          //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Dd_11_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region


      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source, t_sink, 4);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink_1, t_source_1, 3);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Dd_12_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source_1, t_sink, 3);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink_1, t_source, 3);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Dd_21_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source, t_sink_1, 4);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink, t_source_1, 4);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Dd_22_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the 4pt box-diagram: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
void LapH::Correlators::compute_meson_4pt_box_trace_verbose(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();

  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  // TODO: must be changed in GlobalData {
  // TODO: }

  size_t counter = 0;
  std::pair<std::map<size_t, size_t>::iterator, bool> ret;
  map_required_Q2.clear();
  map_required_times.clear();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_c4i10_IO();

  std::cout << "\n\tcomputing the 4pt box-diagram:\r";
  clock_t time = clock();

  for(const auto& op : op_C4_IO){
    for(const auto& i : op.index_pt){
      for(const auto q2 : op_C4[i].index_Q2){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(q2, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
      for(const auto corr : op_C4[i].index_Corr){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(corr, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
    }
  }

  // only need quarklines without fiertz rearangement
  map_required_times.insert(std::pair<size_t, size_t>(5, 0));
  map_required_times.insert(std::pair<size_t, size_t>(6, 1));

  basic.reset_operator();
  basic.init_operator_verbose(map_required_Q2, map_required_times, vdaggerv, peram);

  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the 4pt box-diagram: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;

    for(int t = 0; t < Lt; t++){
      const int t_source = (t_sink + 1 + t)%Lt;
      const int t_source_1 = (t_source + 1) % Lt;

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source, t_sink_1, 6);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink_1, t_source, 5);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Du_11_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source_1, t_sink, 5);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink, t_source_1, 6);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Du_12_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source, t_sink, 6);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink, t_source, 6);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Du_21_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 0, t_source_1, t_sink_1, 5);
      X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                  op_C4_IO, 1, t_sink_1, t_source_1, 5);

      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        for(const auto& i : op.index_pt){
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

            priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                        X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
          }
          #pragma omp critical
          {
            C4_Du_22_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region

    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the 4pt box-diagram: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C4_3() ?
void LapH::Correlators::compute_meson_4pt_cross_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();

  const vec_index_4pt op_C4 = global_data->get_lookup_4pt_trace();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();

  std::fill(C4_mes.data(), C4_mes.data() + C4_mes.num_elements(), 
            cmplx(0.0, 0.0));
  size_t counter = 0;
  std::pair<std::map<size_t, size_t>::iterator, bool> ret;
  map_required_Q2.clear();
  map_required_times.clear();

  std::cout << "\n\tcomputing the traces of 2 pi_+/-:\r";
  clock_t time = clock();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();

  for(const auto& op : op_C4_IO){
    for(const auto& i : op.index_pt){
      for(const auto q2 : op_C4[i].index_Q2){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(q2, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
      for(const auto corr : op_C4[i].index_Corr){
        ret = map_required_Q2.insert(std::pair<size_t, size_t>(corr, counter));
        if(ret.second == true)
        // case element did not exist yet
          counter++;
      }
    }
  }

  // only need quarklines without fiertz rearangement
  map_required_times.insert(std::pair<size_t, size_t>(0, 0));
  map_required_times.insert(std::pair<size_t, size_t>(1, 1));
  map_required_times.insert(std::pair<size_t, size_t>(2, 2));

  basic.reset_operator();
  basic.init_operator(map_required_Q2, map_required_times, vdaggerv, peram);

  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of 2 pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t = 0; t < Lt; t++){
      const int t_source = (t_sink + 1 + t)%Lt;
      const int t_source_1 = (t_source + 1) % Lt;

//      if(t != 0){
//        if(t%2 == 0){
//          X.swap(1, 0);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 0);
//        }
//        else{
//          X.swap(0, 1);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
//        }
//      }
//      else{
//        if(t_source%2 == 0){
          X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                      op_C4_IO, 0, t_source, t_sink, 2);
          X.construct(map_required_Q2, map_required_times, basic, vdaggerv, 
                      op_C4_IO, 1, t_source_1, t_sink, 0);
//        }
//        else{
//          X.construct(basic, vdaggerv, 0, t_source,   t_sink, 0);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
//        }
//      }

//      if(t_source != (t_sink+1)%Lt){
//      if(t != 0){
//        if(t_source%2 == 0){
//          X.swap(1, 0);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 0);
//        }
//        else{
//          X.swap(0, 1);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
//        }
//      }
//      else{
//        if(t_source%2 == 0){
//          X.construct(basic, vdaggerv, 0, t_source,   t_sink, 1);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 0);
//        }
//        else{
//          X.construct(basic, vdaggerv, 0, t_source,   t_sink, 0);
//          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
//        }
//      }

      // The parallelisation is not done with #pragma omp for because it is 
      // incompatible with auto loops
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4_IO){
        // different quantum number combinations denoted by the index i may be 
        // handled by different tasks. Their results are summed up in the end.
        // TODO: Further speedup by different parallelisation might be possible
        for(const auto& i : op.index_pt){

//          size_t i_0 = op_C4[i].index_Q2[0];
//          size_t i_1 = op_C4[i].index_Corr[0];
//          size_t i_2 = op_C4[i].index_Q2[1];
//          size_t i_3 = op_C4[i].index_Corr[1];
          // complete diagramm. combine X and Y to four-trace
          // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
          //             D_d^-1(t_sink       | t_source + 1) Gamma 
          //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
          //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          #pragma omp task shared(op, i)
          {
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {

//            if(t%2 == 0)

              priv_C4 += (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
                          X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace();
//            else
//              priv_C4 += std::conj(
//                         (X(1, op.id, rnd_it[3], rnd_it[0], rnd_it[1]) *
//                          X(0, op.id, rnd_it[1], rnd_it[2], rnd_it[3])).trace());
//
          }
          #pragma omp critical
          {
            C4_mes[op.id][abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
          }
          }
      }}//loops operators
      } // end parallel region
    }// loop t_source
  }// loop t_sink

  std::cout << "\tcomputing the traces of 2 pi_+/-: " << "100.00%" << std::endl;
  time = clock() - time;
  std::cout << "\t\tSUCCESS - " << ((float) time)/CLOCKS_PER_SEC 
            << " seconds" << std::endl;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C3_verbose(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C3_IO = global_data->get_lookup_3pt_IO();

  if(op_C3_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C3_A_mes.data(); i < (C3_A_mes.data()+C3_A_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C3_B_mes.data(); i < (C3_B_mes.data()+C3_B_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C3_mes  - boost structure containing all correlators
  sprintf(outfile, "%s/C3_Au_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C3_IO, "C3I10", C3_A_mes);

  sprintf(outfile, "%s/C3_Bu_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C3_IO, "C3I10", C3_B_mes);

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C3(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_3 rnd_vec_index = global_data->get_rnd_vec_3pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C3_IO = global_data->get_lookup_3pt_IO();

  if(op_C3_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C3_A_mes.data(); i < (C3_A_mes.data()+C3_A_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C3_B_mes.data(); i < (C3_B_mes.data()+C3_B_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C3_mes  - boost structure containing all correlators
  sprintf(outfile, "%s/C3_Ad_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C3_IO, "C3I10", C3_A_mes);

  sprintf(outfile, "%s/C3_Bd_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C3_IO, "C3I10", C3_B_mes);

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C4_box(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_c4i10_IO();

  if(op_C4_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C4_Dd_11_mes.data(); i < (C4_Dd_11_mes.data()+C4_Dd_11_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Dd_12_mes.data(); i < (C4_Dd_12_mes.data()+C4_Dd_12_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Dd_21_mes.data(); i < (C4_Dd_21_mes.data()+C4_Dd_21_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Dd_22_mes.data(); i < (C4_Dd_22_mes.data()+C4_Dd_22_mes.num_elements()); i++)
    *i /= norm1;

  for(auto i = C4_Du_11_mes.data(); i < (C4_Du_11_mes.data()+C4_Du_11_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Du_12_mes.data(); i < (C4_Du_12_mes.data()+C4_Du_12_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Du_21_mes.data(); i < (C4_Du_21_mes.data()+C4_Du_21_mes.num_elements()); i++)
    *i /= norm1;
  for(auto i = C4_Du_22_mes.data(); i < (C4_Du_22_mes.data()+C4_Du_22_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4I10_Dd_11_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Dd_11_mes);
  sprintf(outfile, "%s/C4I10_Dd_12_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Dd_12_mes);
  sprintf(outfile, "%s/C4I10_Dd_21_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Dd_21_mes);
  sprintf(outfile, "%s/C4I10_Dd_22_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Dd_22_mes);

  sprintf(outfile, "%s/C4I10_Du_11_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Du_11_mes);
  sprintf(outfile, "%s/C4I10_Du_12_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Du_12_mes);
  sprintf(outfile, "%s/C4I10_Du_21_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Du_21_mes);
  sprintf(outfile, "%s/C4I10_Du_22_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I10", C4_Du_22_mes);

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

//TODO: is that still necessary?
void LapH::Correlators::write_C4_cross(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_4pt();
  const size_t norm1 = Lt*rnd_vec_index.size();

  const vec_index_IO_1 op_C4_IO = global_data->get_lookup_4pt_3_IO();

  if(op_C4_IO.size() == 0)
    return;

  // normalisation
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to lime file
  // outfile - filename
  // C4_mes  - boost structure containing all correlators

  sprintf(outfile, "%s/C4I2+_3_conf%04d.dat", outpath.c_str(), (int)config_i);
  export_corr_IO(outfile, op_C4_IO, "C4I2+", C4_mes);

}


