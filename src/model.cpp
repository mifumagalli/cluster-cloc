model::model(double in_t0, double in_t1, double in_mmin, 
                    double in_mmax, double in_cmfslope,
                    double in_gamma_min, double in_gamma_max, 
                    double in_eta, double in_ageslope,
                    double in_start, double in_stop, double in_step, 
		    double in_length, double in_fc){
    f_c = in_fc;
    length = in_length;
    t0 = in_t0;
    t1 = in_t1;
    mmin = in_mmin;
    mmax = in_mmax;
    cmf_slope = in_cmfslope;
    cmf_slope1 = cmf_slope + 1.;
    eta = in_eta;
    gamma_min = in_gamma_min;
    gamma_max = in_gamma_max;
    age_slope = in_ageslope;
    start = in_start;
    stop = in_stop;
    step = in_step;
//make array of luminosity convertying to ln L
    for(double v =log(pow(10,start)); v<log(pow(10,stop));v+=log(pow(10,step))){
          x.push_back(v);
    }
//cout<<x.size()<<endl;
//exit(1);
    omega = (age_slope + 1.)/eta - cmf_slope1;
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
//for(int i=0;i<(long)x.size();i++)cout<<x[i]<<endl;
    get_cdf();   //CDF needs to be first to get normalization right
    get_pdf();

/*    sfr_grid.push_back(x);
    sfr_grid.push_back(x);
    sfr_grid.push_back(x);
    sfr_grid.push_back(x);
    cout<<x[0]<<"\t"<<x[1]<<endl;
    cout<<sfr_grid[0][0]<<"\t"<<sfr_grid[0][1]<<endl;
    cout<<sfr_grid[1][0]<<"\t"<<sfr_grid[1][1]<<endl;
    cout<<sfr_grid[2][0]<<"\t"<<sfr_grid[2][1]<<endl;
    cout<<sfr_grid[3][0]<<"\t"<<sfr_grid[3][1]<<endl;
exit(1);*/
}
//-----------------------------------------------------------------
void model::obs_pdf(double err){
   obs_err_pdf.clear();
   obs_err_cdf.clear();
   log_obs_err_cdf.clear();
   obs_err_pdf.resize(x.size());
   obs_err_cdf.resize(x.size());
   log_obs_err_cdf.resize(x.size());
   obs_err = err;
   vector<double> wgts_y;

   double fac = 1./(obs_err*sqrt(2*3.14159265359));
   long zero_ind=99999999999, ct=0;
   long num_ind=0;
   double diff_val=0;
   for(int i =0;i<(long)x.size();i++){
//cout<<x[i]<<endl;
      if(x[i]-x[0] > 5*obs_err){
         num_ind = i;
         diff_val = x[i]-x[0];
         break;
      }
   }
   for(double xx= -diff_val;xx<diff_val+log(pow(10,step))*.98; xx+=log(pow(10,step))){  
          if(fabs(xx) < 1e-8) zero_ind=ct;
          ct++;
          wgts_y.push_back(fac*exp(-0.5*(xx*xx)/(obs_err*obs_err)));
   }
//exit(10);
   double pdf_tmp, cdf_tmp;
   long xsz = (long) x.size();
   long wgts_ysz = (long) wgts_y.size();
   long wgts_wing = wgts_ysz - zero_ind-1;
   if(wgts_ysz != wgts_wing*2+1){
     cout<<"model::obs_pdf wgts is not the right dimension.. logic bug"<<endl;
     cout<<wgts_ysz<<"\t"<<wgts_wing<<"\t"<<zero_ind<<endl;
     exit(1);
   }
   double sum=0;
   for(int i=0; i<(long)x.size();i++){
       sum=0;
       pdf_tmp=0;
       cdf_tmp=0;
       for(int j= max((long)0,i-wgts_wing); j< min(i+wgts_wing, xsz-1); j++){
           //cout<<j-i+zero_ind<<"\t"<<j<<"\t"<<i<<"\t"<<zero_ind<<endl;
           sum += wgts_y[j-i+zero_ind];
           pdf_tmp+=wgts_y[j-i+zero_ind]*pdf[j];
           cdf_tmp+=wgts_y[j-i+zero_ind]*cdf[j];
       }
       obs_err_pdf[i]=(pdf_tmp/sum);
       obs_err_cdf[i]=(cdf_tmp/sum);
   }
   long last = (long)x.size()-1;
   double sum2 =tsum(x, obs_err_pdf);
   for(int i=0; i<(long)x.size();i++){
       obs_err_pdf[i] /= sum2;
//have to do pdf first
       obs_err_cdf[i] /= obs_err_cdf[last];
       log_obs_err_cdf[i] = log(obs_err_cdf[i]);
   }
 
   
}

//-----------------------------------------------------------------
vector<double> model::sfr_pdf(double sfr){
    vector<double> out;
    out.resize(x.size());
    double expr =0;
    //double buffer =0;
    double nhat = f_c*length*pow(10., sfr)/mean_m;
      double tmp_log_nhat;
//cout<<sfr<<"\t"<<nhat<<"\t"<<f_c<<"\t"<<length<<"\t"<<pow(10., sfr)<<"\t"<<mean_m<<endl;
    if(nhat < 100){
    //if nhat is less than 100 then we can just brute force 300 iterations
    //to get the full distribution
    //I assume at least 1 cluser in each galaxy looked at because if there
    // wasn't a cluster it probably wasn't included in this study
    // THis is why the for loop starts at 1.... I will thus have to renormalize
      tmp_log_nhat = log(nhat);
      double jm1, tmp_fac;
      for(double j=1; j<300;j+=1){
//           j=floor(jj + 0.5); //round
           jm1=j-1;
           tmp_fac= j*exp( tmp_log_nhat*j-  lgamma(j+1)-nhat);
           for(int i=0;i<(long)x.size(); i++){
              expr = exp( log_obs_err_cdf[i]*(jm1))*obs_err_pdf[i]*
                          tmp_fac;
              if(expr != expr) expr=0;
              out[i]+=expr;
           }
       }
       //case of no clusters
cout<<exp(-nhat)<<endl;
       out[0] += exp(-nhat);


     /*   for(int i=0;i<(long)x.size(); i++){
           buffer = 0;
//           tmp_log_obs_err_cdf = log(obs_err_cdf[i]);
           tmp_log_nhat = log(nhat);
           for(int j=1;j < 300; j++){
               expr = j*exp( log_obs_err_cdf[i]*(j-1))*obs_err_pdf[i]*
                     exp( tmp_log_nhat*j-  lgamma(j+1)-nhat);
              if(expr != expr) expr=0;
              buffer += expr;
           }
           out[i]=buffer;
        }*/

    } else {
//majority of the code lives in this loop!!!!
      double sigma=sqrt(nhat);
      double sig_step=sigma/10.;
      double j=0;
      double tmp_log_nhat, tmp_fac, jm1;
      tmp_log_nhat = log(nhat);
      for(double jj=nhat-5*sigma; jj<( nhat+5*sigma);jj+=sig_step){
           j=floor(jj + 0.5); //round
           jm1=j-1;
           tmp_fac= j*exp( tmp_log_nhat*j-  lgamma(j+1)-nhat);
           for(int i=0;i<(long)x.size(); i++){
              expr = exp( log_obs_err_cdf[i]*(jm1))*obs_err_pdf[i]*
			  tmp_fac;
              if(expr != expr) expr=0;
              out[i]+=expr;
           }
       }
    }

 double sum = tsum(x, out);
 for(int i = 0;i<(long)x.size();i++) out[i] /= sum;
//cout<<out.size()<<endl;
 return(out);
}

//-----------------------------------------------------------------
void model::mk_sfr_grid(double in_sfr_start, double in_sfr_stop, 
                        double in_sfr_step){
    sfr_start = in_sfr_start;
    sfr_stop = in_sfr_stop;
    sfr_step = in_sfr_step;
 //   sfr_gridx.resize(x.size());
//    long ct = 0;
    for(double s = sfr_start; s<sfr_stop; s+=sfr_step){
       sfr_gridx.push_back(s);
  //     ct++;
    }
    sfr_grid.resize(sfr_gridx.size());
    long ct2=0;
    for(double s = sfr_start; s<sfr_stop; s+=sfr_step){
//       cout<<ct2<<endl;
       sfr_grid[ct2] = (sfr_pdf(s));
//if( ct2 == 554) cout<<"adsfads"<<sfr_grid[ct2].size()<<endl;
       ct2++;
    }
}

//-----------------------------------------------------------------

void model::mix_sfr_grid(double err){
   sfr_err = err;
   vector<double> wgts_y;
   //calculate the vector the determines the weights for the mixing
   double fac = 1./(sfr_err*sqrt(2*3.14159265359));
   long zero_ind=99999999999, ct=0;
   long num_ind=0;
   double diff_val=0;
   for(int i =0;i<(long)sfr_gridx.size();i++){
      if(sfr_gridx[i]-sfr_gridx[0] > 5*sfr_err){
         num_ind = i;
         diff_val = sfr_gridx[i]-sfr_gridx[0];
         break;
      }
   }
//   cout<<diff_val<<endl;
   for(double xx= -diff_val;xx<diff_val + sfr_step*.98; xx+=sfr_step){  
          if(fabs(xx) < 1e-8) zero_ind=ct;
          ct++;
          wgts_y.push_back(fac*exp(-0.5*(xx*xx)/(sfr_err*sfr_err)));
   }
   printvec(wgts_y, "wgts_y.txt");
   //calculate some stuff to do the mixing
   long xsz = (long) sfr_gridx.size();
   long wgts_ysz = (long) wgts_y.size();
   long wgts_wing = wgts_ysz - zero_ind-1; // we go +- wgts_wing of central point
   //define the mixed_sfr_grid and size it appropriately
   mixed_sfr_grid.resize(sfr_gridx.size());
   for(int i=0; i<(long)sfr_gridx.size();i++){
         mixed_sfr_grid[i].resize(x.size());
   } //things get automatically intiailized to 0
   if(wgts_ysz != wgts_wing*2+1){
     cout<<"model::mix_sfr_grid wgts is not the right dimension.. logic bug"<<endl;
     cout<<wgts_ysz<<"\t"<<wgts_wing<<"\t"<<zero_ind<<endl;
     exit(1);
   }
   double sum=0;
   double buffer;
   for(int i=0; i<(long)sfr_gridx.size();i++){ //loop over SFR
         sum=0;
         //find the sum of the wgts for renormalization later
         for(int j= max((long)0,i-wgts_wing); j< min(i+wgts_wing, xsz-1); j++){
             sum += wgts_y[j-i+zero_ind];
         }
//if(sum == 0){
//cout<<"WTF"<<endl;
//cout<<sum<<endl;
//exit(1);
//}
         for(int j= max((long)0,i-wgts_wing); j< min(i+wgts_wing, xsz-1); j++){
             //loop over smoothing kernel
             for(int k =0;k<(long)x.size();k++){
                buffer = wgts_y[j-i+zero_ind]*sfr_grid[j][k];
                if(buffer != buffer) buffer =0;
                mixed_sfr_grid[i][k] += buffer;
             }
         }
         //divide by the sum from above
         for(int k =0;k<(long)x.size();k++){
            mixed_sfr_grid[i][k] /= sum;
         }
   }

    
}

//-----------------------------------------------------------------
/*vector<double> model::sfr_pdf(double sfr){
    vector<double> out;
    out.resize(x.size());
   vector<double> wgts_y;
   double fac = 1./(sfr_err*sqrt(2*3.14159265359));
 //  for(int i =0; i<(long)wgts_x.size();i++)
   //        wgts_y.push_back(fac*exp(-0.5*(wgts_x[i]*wgts_x[i])/obs_err));
   long zero_ind=99999999999, ct=0;
   for(double xx= -5;xx<5; xx+=step){
          if(fabs(xx) < 1e-8) zero_ind=ct;
          ct++;
          wgts_y.push_back(fac*exp(-0.5*(xx*xx)/(obs_err*sfr_err)));
   }


}*/
//-----------------------------------------------------------------
void model::print_pdf(vector<double> p){
  ofstream fout;
  fout.open("test2.txt");
  for(int i =0; i<(long)x.size();i++){
     fout<<x[i]<<"\t"<<p[i]<<endl;
  }
  fout.close();
}

//-----------------------------------------------------------------
void model::printstate(void){
  cout<<"t0:\t"<<t0<<endl;
  cout<<"t1:\t"<<t1<<endl;
  cout<<"mmin:\t"<<mmin<<endl;
  cout<<"mmax:\t"<<mmax<<endl;
  cout<<"cmf_slope:\t"<<cmf_slope<<endl;
  cout<<"cmf_slope1:\t"<<cmf_slope1<<endl;
  cout<<"mean_m:\t"<<mean_m<<endl;
  cout<<"eta:\t"<<eta<<endl;
  cout<<"gamma_min:\t"<<gamma_min<<endl;
  cout<<"gamma_max:\t"<<gamma_max<<endl;
  cout<<"age_slope:\t"<<age_slope<<endl;
  cout<<"omega:\t"<<omega<<endl;
  cout<<"sfr_start:\t"<<sfr_start<<endl;
  cout<<"sfr_stop:\t"<<sfr_stop<<endl;
  cout<<"sfr_step:\t"<<sfr_step<<endl;
  cout<<"start:\t"<<start<<endl;
  cout<<"stop:\t"<<stop<<endl;
  cout<<"step:\t"<<step<<endl;
  cout<<"f_c:\t"<<f_c<<endl;
  cout<<"length:\t"<<length<<endl;
  cout<<"obs_err:\t"<<obs_err<<endl;
  cout<<"sfr_err:\t"<<sfr_err<<endl;
   
}


//-----------------------------------------------------------------

void model::print(void){
    ofstream fout;
//cout<<x.size()<<endl;
//cout<<pdf.size()<<endl;
//cout<<cdf.size()<<endl;
    if(obs_err_pdf.size() < 10){
       cout<<"MODEL::PRINT obs_err_pdf not yet set"<<endl;
       cout<<"you must run model::obs_pdf before model::print"<<endl;
       exit(1);
    }
    fout.open("test.txt");
long half = (long)floor(sfr_grid.size()/2.+0.5);
//cout<<"SFR:"<<sfr_gridx[half]<<endl;
//cout<<sfr_grid.size()<<"\t"<<mixed_sfr_grid.size()<<endl;

//exit(1);
    for(int i=0;i<(long)x.size();i++){
      fout<<x[i]<<"\t"<<pdf[i]<<"\t"<<cdf[i]<<"\t"<<obs_err_pdf[i]<<"\t"<<obs_err_cdf[i]<<"\t"<<sfr_grid[half][i]<<"\t"<<mixed_sfr_grid[half][i]<<endl;;
   //   cout<<x[i]<<"\t"<<pdf[i]<<"\t"<<cdf[i]<<endl;
    }
    fout.close();
}

//-----------------------------------------------------------------
void model::print_sfr_grid(string name){
   cout<<"Writing Mixed SFR grid to "<<name<<endl;
   ofstream fout;
   double buffer;
   fout.open(name.data(), ofstream::binary);
   //write info on SFR axis
   buffer = (double) sfr_gridx.size();
   //the number of entries
   fout.write((char *)(&(buffer)),sizeof(buffer));

   //the values of each
   for(int i = 0; i<(long)sfr_gridx.size();i++){
      buffer = sfr_gridx[i];
      fout.write((char *)(&(buffer)),sizeof(buffer));
   }
   //the number of entrix on the lum axis
   buffer = (double) x.size();
   fout.write((char *)(&(buffer)),sizeof(buffer));
   //the entries
   for(int i = 0; i<(long)x.size();i++){
      buffer = x[i];
      fout.write((char *)(&(buffer)),sizeof(buffer));
   }

   for(int i = 0; i<(long)x.size();i++){
      buffer = obs_err_pdf[i];
      fout.write((char *)(&(buffer)),sizeof(buffer));
   }
   for(int i = 0; i<(long)x.size();i++){
      buffer = pdf[i];
      fout.write((char *)(&(buffer)),sizeof(buffer));
   }

   //now the grid itself
   for(int i = 0; i<(long)sfr_gridx.size();i++){
      for(int j = 0; j<(long)x.size();j++){
         buffer = mixed_sfr_grid[i][j];
         fout.write((char *)(&(buffer)),sizeof(buffer));
      }
   }
   fout.close();

}

//-----------------------------------------------------------------


model::~model(void){
}
