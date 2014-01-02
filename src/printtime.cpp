
void printtime(void){
time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo=localtime(&rawtime);
  cout<<asctime(timeinfo)<<endl;
}
