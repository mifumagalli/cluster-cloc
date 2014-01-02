;/*z0, z1: the range of the integration 
	;e.g. log mmin*gamma_min
;t0, t1: the range of the dust step function
	; in log lum... ie.  should be ~50 not exp(50)
;L : the log luminosity being evaluated
;C: constant in front of expression being integrated
;p: the power that L is raised to*/



double model::gamt(double t){
return ( 7.782e24*pow(t,eta));
}


double model::eval(double  z0, double z1, double t0, double t1, 
                    double L, double C, double p){ 
double out=0.;
if(p == 0){
  cout<<"Power p is equal to zero"<< p<<endl;
  exit(1);
}
if (L > z1-t0 || L < z0-t1) return( 0);


if(z0-t0 > z1-t1){
  if (L > z0-t1 && L < z1-t1) out=(C/p)*(exp(p*(t1+L)) - exp(p*z0));
  if (L > z1-t1 && L < z0-t0) out=(C/p)*(exp(p*z1) - exp(p*z0));
  if (L > z0-t0 && L < z1-t0) out=(C/p)*(exp(p*z1) - exp(p*(t0+L)));
} else {
  if( L > z0-t1 && L < z0-t0) out=(C/p)*(exp(p*(t1+L)) - exp(p*z0));
  if( L > z0-t0 && L < z1-t1) out=(C/p)*(exp(p*(t1+L)) - exp(p*(t0+L)));
  if( L > z1-t1 && L < z1-t0) out=(C/p)*(exp(p*z1) - exp(p*(t0+L)));
}

return (out);
}
//;-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
double model::pdf_eval(double L, double t0, double t1){
//common steps, mmin, mmax, gamma_min, gamma_max, cmf_slope1, omega
double outpdf=0.;
//branch 1
outpdf+=eval(log(mmax*gamma_min), log(mmax*gamma_max), 
		t0, t1, L, pow(gamma_max, omega), cmf_slope1);
outpdf+=eval(log(mmax*gamma_min),log( mmax*gamma_max), t0,
	 t1, L, -pow(mmax, -omega), cmf_slope1+omega);

//cout<<outpdf<<endl;
//branch 2
outpdf+=eval(log(mmin*gamma_max), log(mmax*gamma_min),
                t0, t1, L, pow(gamma_max, omega)-
		pow(gamma_min, omega), cmf_slope1);

//cout<<outpdf<<endl;
//branch 3
outpdf += eval(log(mmin*gamma_min), log(mmin*gamma_max), t0, 
            t1, L, pow(mmin, -omega), cmf_slope1+omega);
//cout<<outpdf<<endl;
outpdf += eval(log(mmin*gamma_min), log(mmin*gamma_max),
            t0, t1, L, -pow(gamma_min, omega), cmf_slope1);

//cout<<outpdf<<endl;
//cout<<L<<endl;

//printstate();
//exit(1);
return (outpdf);
}
//;-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

double model::eval2(double z0, double z1, double t0, double t1, 
                     double L, double C, double p){
double outcdf=0.;
double LL=0.;
if (p == 0){
  cout<<"Power p is equal to zero\t"<<endl;
  exit(1);
}
if(z0-t0 > z1-t1) {
	LL=min(max(L,z0-t1), z1-t1);
	outcdf+=exp(p*t1)*(exp(p*LL)-exp(p*(z0-t1)))/p ;
	outcdf-= exp(p*z0)*(LL-(z0-t1));

	LL=min( max(L,z1-t1) , z0-t0);
	outcdf+=(exp(p*z1)-exp(p*z0))*(LL-(z1-t1));

	LL=min(max(L, z0-t0), z1-t0);
	outcdf+=exp(p*z1)*(LL-(z0-t0)) ;
	outcdf-= exp(p*t0)*(exp(p*LL)-exp(p*(z0-t0)))/p;
} else {
	LL=min(max(L, z0-t1), z0-t0);
	outcdf+=exp(p*t1)*(exp(p*LL)-exp(p*(z0-t1)))/p; 
	outcdf-= exp(p*z0)*(LL-(z0-t1));

	LL=min(max(L, z0-t0), z1-t1);
	outcdf+=(exp(p*t1)-exp(p*t0))*(exp(p*LL)-exp(p*(z0-t0))) /p;

	LL=min(max(L, z1-t1), z1-t0);
	outcdf+=exp(p*z1)*(LL-(z1-t1)); 
	outcdf-= exp(p*t0)*(exp(p*LL)-exp(p*(z1-t1)))/p;
}

return (outcdf*(C/p));
}

//;-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
double model::cdf_eval1(double  L, double t0, double t1){
double outcdf2=0.;
//branch 1
outcdf2+=eval2(log(mmax*gamma_min), log(mmax*gamma_max), t0,
	 t1, L, pow(gamma_max, omega), cmf_slope1);
outcdf2+=eval2(log(mmax*gamma_min),log( mmax*gamma_max), t0, 
	t1, L, -pow(mmax, -omega), cmf_slope1+omega);

//branch 2
outcdf2+=eval2(log(mmin*gamma_max), log(mmax*gamma_min), t0, 
                t1, L, pow(gamma_max, omega)-
                pow(gamma_min, omega), cmf_slope1);

//branch 3
outcdf2+=eval2(log(mmin*gamma_min), log(mmin*gamma_max), t0, 
	t1, L, pow(mmin, -omega), cmf_slope1+omega);
outcdf2+=eval2(log(mmin*gamma_min),log( mmin*gamma_max), t0,
	 t1, L, -pow(gamma_min, omega), cmf_slope1);

//cout<<L<<"\t"<<outcdf2<<endl;
//exit(1);
return (outcdf2);
}

void model::get_pdf(void){
//pdf=dblarr(n_elements(l))
  for(int i=0;i<(long)x.size();i++)
     pdf.push_back(pdf_eval(x[i], t0, t1));
  //long last = (long)x.size()-1;

  double sum = tsum(x, pdf);
  for(int i=0;i<(long)x.size();i++)
     pdf[i] /= sum;
  
  //for(int i=0;i<(long)x.size();i++)
  //   pdf[i] /= cdf[last];
  
}

void model::get_cdf(void){
  for(int i=0;i<(long)x.size();i++)
     cdf.push_back(cdf_eval1(x[i], t0, t1));
  long last = (long)x.size()-1;

  for(int i=0;i<(long)x.size();i++)
     cdf[i] /= cdf[last];
}

