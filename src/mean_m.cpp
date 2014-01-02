
double mean_m(double mmin, double mmax, double cmf_slope) {
  double mean_m;
  if (cmf_slope == -1){
    mean_m=(pow(mmax,cmf_slope+2)-pow(mmin,cmf_slope+2))/(cmf_slope+2);
    mean_m/=log(mmax/mmin);
  } else if( cmf_slope == -2){
    mean_m=log(mmax/mmin);
    mean_m/=(pow(mmax,cmf_slope+1)-pow(mmin,cmf_slope+1))/(cmf_slope+1);
  }else{
    mean_m=(pow(mmax,cmf_slope+2)-pow(mmin,cmf_slope+2))/(cmf_slope+2);
    mean_m/=(pow(mmax,cmf_slope+1)-pow(mmin,cmf_slope+1))/(cmf_slope+1);
  }
   
  return(mean_m);
}
