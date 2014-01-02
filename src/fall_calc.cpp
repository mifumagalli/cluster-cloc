

#include "def.h"

using namespace std;
int main(int argc, char* argv[]){
  // Temporarily hardcode these values
  double t0 = 0.;
  double t1 = 0.1;
  double mmin = 100;
  double mmax = 1e7;
  double cmfslope = -2;
  double gamma_min = 5.05e18;
  double gamma_max = 1.198e20;
  double eta = -0.688;
  double age_slope = -0.9;
  // initialize the model
  double start = 15, stop = 35, step=0.0125; //these are in terms of
						//log_10 lum
  double sfr_err=0.5;
  double obs_err=log(sqrt(2.512));
  //double sfr=0;
  double f_c=.01;
  double sfr_start = -4, sfr_stop = 4, sfr_step=.0125;//sfr_step=.0125;
  double length = 1e9; //length of time in yr for SFR event
			//important for nhat calcuation from SFR
  if(argc > 1) t0 =  atof(argv[1]);
//cout<<argc<<"\t"<<t0*2<<endl;
  if(argc > 2) t1 = atof(argv[2]);
  if(argc > 3) mmin = atof(argv[3]);
  if(argc > 4) mmax = atof(argv[4]);
  if(argc > 5) cmfslope = atof(argv[5]);
  if(argc > 6) gamma_min = atof(argv[6]);
  if(argc > 7) gamma_max = atof(argv[7]);
  if(argc > 8) eta = atof(argv[8]);
  if(argc > 9) age_slope = atof(argv[9]);
  if(argc > 10) start = atof(argv[10]);
  if(argc > 11) stop = atof(argv[11]);
  if(argc > 12) step = atof(argv[12]);
  if(argc > 13) length = atof(argv[13]);
  if(argc > 14) f_c = atof(argv[14]);
  if(argc > 15) obs_err = atof(argv[15]);
  if(argc > 16) sfr_start = atof(argv[16]);
  if(argc > 17) sfr_stop = atof(argv[17]);
  if(argc > 18) sfr_step = atof(argv[18]);
  if(argc > 19) sfr_err = atof(argv[19]); 
  string grid_out="grid_def2.dat";
  if(argc > 20) grid_out = (string) argv[20];
//  if(argc > 20) t0 = atof(*argv[20]);
  //---



printtime();
  //double start = 15, stop = 30, step=0.05;
  model mod(t0, t1, mmin, mmax, cmfslope, gamma_min, gamma_max, eta,
             age_slope, start, stop, step, length, f_c);
  //Now I need to mix the PDF by the observational errora
  mod.obs_pdf(obs_err);

  

//need to make a grid of SFRs and then convolve on that dimension to get
// the final including uncertain SFR
  mod.mk_sfr_grid(sfr_start, sfr_stop, sfr_step);
  //cout<<"starting to mix sfr error"<<endl;
 // printtime();
  mod.mix_sfr_grid(sfr_err);

  mod.printstate();
  mod.print();
  mod.print_sfr_grid(grid_out);
  printtime();
  exit(0);
}
