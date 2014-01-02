
//include standard libraries
#include <sstream>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h> 
#include <fstream>
#include <vector>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <algorithm>

//include gsl stuff
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
using namespace std;

class model;

double tsum(vector<double>, vector<double>);

class model {
public:
  double gamma_min, gamma_max, mmin, mmax, eta, omega, t0, t1, cmf_slope;
  double age_slope, cmf_slope1, mean_m;
  vector<double> pdf;
  vector<double> cdf;
  vector<double> obs_err_pdf;
  vector<double> obs_err_cdf;
  vector<double> log_obs_err_cdf;
  vector<double> x;
  vector< vector<double> > sfr_grid;
  vector< vector<double> > mixed_sfr_grid;
  double sfr_start, sfr_stop, sfr_step, sfr_err;
  vector<double> sfr_gridx;
  double obs_err, start, stop, step, f_c, length;
  double gamt(double t);
  double eval(double  z0, double z1, double t0, double t1,
                    double L, double C, double p);
  double pdf_eval(double L, double t0, double t1);
  double eval2(double z0, double z1, double t0, double t1,  
                     double L, double C, double p);
  double cdf_eval1(double  L, double t0, double t1);
  model(double, double, double, double, double, double, 
        double,  double, double, double, double, double, double, double);
  ~model(void);
  void printstate(void);
  void get_pdf(void);
  void get_cdf(void);
  vector<double> sfr_pdf(double);
  void obs_pdf(double);
  void print(void);
  void print_pdf(vector<double>);
  void mk_sfr_grid(double, double, double);
  void mix_sfr_grid(double);
  void print_sfr_grid(string);
};

void printtime(void);
void printvec(vector<double>, string);
#include "tsum.cpp"
#include "model.cpp"
#include "eval_step.cpp"
#include "printtime.cpp"
#include "printvec.cpp"
