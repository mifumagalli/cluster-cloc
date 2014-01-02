
double tsum(vector<double> x, vector<double>y) {
    if (x.size() != y.size()) {
      cout<<("x and y must be the same size")<<endl;
      exit(1);
    }
    double sum = 0.0;
    for (int i = 1; i < (long) x.size(); i++) {
        sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);
    }
    return sum * 0.5;
}
