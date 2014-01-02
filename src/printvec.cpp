void printvec(vector<double> in, string name){
  ofstream fout;
  fout.open(name.data());
  for(int i =0;i<(long)in.size();i++){
    fout<<in[i]<<endl;
  }
  fout.close();

}
