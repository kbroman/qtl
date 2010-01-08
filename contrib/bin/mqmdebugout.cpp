#include "mqm.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <limits>
#include <getopt.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

using namespace std;

FILE* redirect_info;  // Redirect output for testing

bool checkfileexists(const char *filename) {
  ifstream myfile;
  bool exists;
  myfile.open(filename);
  exists = myfile.is_open();
  myfile.close();
  return exists;
}

double mytruncate(double n, double p = 3){
    int sign = 0;
	if(n >= 0){
        sign = 1;
    }else{
        sign = -1;
    }
	double val = fabs((pow(10,p)) * n);
	val = floor(val);
    val /= pow(10,p);
	return (double) sign * val;
}

void unittest_pbeta(void){
  for(int df1=1;df1<50;df1++){
    for(int df2=1;df2<df1;df2++){
      for(double halfway = 5.0; halfway < 100;halfway+=0.5){
        double prob = pbeta(df2/(df2+df1*halfway), df2/2.0, df1/2.0, 1, 0);
        if(prob > 0.0001) fprintf(redirect_info,"df1:%d df2:%d hw:%f prob:%.5f\n",df1,df2,halfway,mytruncate(prob,5));
      }
    }
  }
}


void unittest_dnorm(void){
  for(double variance=1.0;variance < 100.0;variance += 0.5){
    for(double residual=1.0;residual < variance;residual += 0.5){
      double prob = dnorm(residual,0,sqrt(variance),0);
      if(prob > 0.00001) fprintf(redirect_info,"%f %f %f\n",variance,residual,mytruncate(prob,5));
    }
  }
}

static struct option long_options[] = {
  {0, 0, 0, 0}
};

int main(int argc, char** argv){
  printf("testing external functions used in mqm\n");
  int option_index = 0;
  char c;
  char* outputfile = NULL;
  bool pbetaflag= false;
  bool dnormflag= false;
  //Parsing of arguments
  while ((c = getopt_long(argc, argv, "dpo:",long_options, &option_index)) != -1)
    switch (c) {
    case 'p':
      pbetaflag = true;
      break;
    case 'd':
      dnormflag = true;
      break;
    case 'o':
      outputfile = optarg;
      printf("outputfile=%s\n",outputfile);
      break;      
    default:
      break;
    }
    //Check the output file
    if (outputfile != NULL){
    // Open outputstream if specified - using C type for redirection
    FILE *fout = stdout;
    if (outputfile){
      fout = fopen(outputfile,"w");
      redirect_info = fout;
    }    
    if(pbetaflag) unittest_pbeta();
    if(dnormflag) unittest_dnorm();
  
  
    printf("done\n");
    return 0;
  }else{
    printf("No output file\n");
    return 1;
  }
}
