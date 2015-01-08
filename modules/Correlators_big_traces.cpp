#include "Correlators.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//TODO: Call that build_C4_3() ?
void LapH::Correlators::compute_meson_4pt_cross_trace(LapH::CrossOperator& X) {

  const int Lt = global_data->get_Lt();
  const vec_pdg_C4 op_C4 = global_data->get_op_C4();
  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();
  // TODO: must be changed in GlobalData {
  // TODO: }

  std::cout << "\n\tcomputing the traces of 2 pi_+/-:\r";
  clock_t time = clock();
  for(int t_sink = 0; t_sink < Lt; ++t_sink){
    std::cout << "\tcomputing the traces of 2 pi_+/-: " 
        << std::setprecision(2) << (float) t_sink/Lt*100 << "%\r" 
        << std::flush;
    int t_sink_1 = (t_sink + 1) % Lt;
    for(int t_source = 0; t_source < Lt; ++t_source){
      const int t_source_1 = (t_source + 1) % Lt;

      if(t_source != 0){
        if(t_source%2 == 0){
          X.swap(1, 0);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
        }
        else{
          X.swap(0, 1);
          X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 0);
        }
      }
      else{
        X.construct(basic, vdaggerv, 0, t_source,   t_sink, 0);
        X.construct(basic, vdaggerv, 1, t_source_1, t_sink, 1);
      }
      if(t_source == t_sink)
        continue;
    
      #pragma omp parallel
      #pragma omp single
      {
      for(const auto& op : op_C4){
      for(const auto& i : op.index){
        // complete diagramm. combine X and Y to four-trace
        // C4_mes = tr(D_u^-1(t_source     | t_sink      ) Gamma 
        //             D_d^-1(t_sink       | t_source + 1) Gamma 
        //             D_u^-1(t_source + 1 | t_sink + 1  ) Gamma 
        //             D_d^-1(t_sink + 1   | t_source    ) Gamma)
          cmplx priv_C4(0.0,0.0);
          for(const auto& rnd_it : rnd_vec_index) {
            #pragma omp task shared(rnd_it, i)
            if(t_source%2 == 0)
              priv_C4 += (X(0, i[2], i[1], rnd_it[2], rnd_it[1], rnd_it[3]) *
                          X(1, i[3], i[0], rnd_it[3], rnd_it[0], rnd_it[2])).trace();
            else
              priv_C4 += std::conj(
                         (X(0, i[2], i[1], rnd_it[2], rnd_it[1], rnd_it[3]) *
                          X(1, i[3], i[0], rnd_it[3], rnd_it[0], rnd_it[2])).trace());
            }
          #pragma omp critical
          {
            C4_mes[op.p_sq_cm][op.p_sq_so_1][op.p_sq_si_1][op.dg_so][op.dg_si]
                [abs((t_sink - t_source) - Lt) % Lt] += priv_C4;
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
void LapH::Correlators::write_C4_3(const size_t config_i){

  char outfile[400];
  FILE *fp = NULL;
  std::string outpath = global_data->get_output_path() + "/" + 
      global_data->get_name_lattice();

  const int Lt = global_data->get_Lt();
  const size_t nb_mom_sq = global_data->get_number_of_momentum_squared();
  const std::vector<int> dirac_ind {5};
  const size_t nb_dir = dirac_ind.size();

  const indexlist_4 rnd_vec_index = global_data->get_rnd_vec_C4();
  const size_t norm1 = Lt*rnd_vec_index.size();

  // normalisation
  for(auto i = C4_mes.data(); i < (C4_mes.data()+C4_mes.num_elements()); i++)
    *i /= norm1;

  // output to binary file
  for(size_t dirac_1 = 0; dirac_1 < nb_dir; ++dirac_1){     
  for(size_t dirac_2 = 0; dirac_2 < nb_dir; ++dirac_2){
    for(size_t p = 0; p < nb_mom_sq; p++){
      sprintf(outfile, 
          "%s/dirac_%02d_%02d_p_%01d_%01d_displ_%01d_%01d/"
          "C4_3_conf%04d.dat", 
          outpath.c_str(), dirac_ind.at(dirac_1), dirac_ind.at(dirac_2), 
          (int)p, (int)p, 0, 0, (int)config_i);
//      std::cout << outfile << std::endl;
      if((fp = fopen(outfile, "wb")) == NULL)
        std::cout << "fail to open outputfile" << std::endl;
      fwrite((double*) &(C4_mes[0][p][p][dirac_1][dirac_2][0]), 
             sizeof(double), 2 * Lt, fp);
      fclose(fp);
    }// loop p
  }}// loops dirac_1 dirac_2

}

