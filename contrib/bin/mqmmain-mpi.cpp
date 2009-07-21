/**********************************************************************
 *
 * mqmmain-mpi.cpp - standalone MQM MPI version
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * last modified Apr,2009
 * first written Feb, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <string>

#define MAXCOMMANDLENGTH   1000

using namespace std;

void ourerror(char *s) {
  printf("\n** Error: %s **\n",s);
  exit(-1);
}

int threader(int num,int batchsize,int start,int stop,int thread_id,char *command,int verbose) {
  if (verbose) {
    printf("Executing # %d on %d on thread: %d with command: %s\n",start,num,thread_id, command);
  }
  int e;
  char MQMcommand[MAXCOMMANDLENGTH];

  for (int i=start; i<stop; i++) {
    sprintf(MQMcommand,"%s -T=%d",command,i);
    printf("%s\n",MQMcommand);
    e = system(MQMcommand);
    if (e !=1) {
      e = i;
      break;
    }
  }
  return e;
}

int main(int argc, char *argv[]) {
  char command[MAXCOMMANDLENGTH];
  int n=1000;
  int nthreads = 10;
  int nbatch = 50;
  int verbose = 0;
  int log = 0;
  int itemstodo = 10000;
  int lbatch = 0;
  int lcores = 0;
  int itemspbatch = 0;
  int nroftodos;
  char c;
  int thread_id;
  strcpy(command,"N");
  printf("C++ Multiprocessor for sMQM V0.1 (c) Danny Arends\n");
  for (int i=1; i<argc; i++) {
    if (!strcmp(argv[i],argv[0])) continue;
    if (argv[i][0] != '-') ourerror("Dash needed at argument, argument structure for booleans: -V or for strings/numbers -C=<command>\n");

    c = toupper(argv[i][1]);

    // On-Off flags----------------
    if (c == 'L') {
      log=1;
      continue;
    }
    if (c == 'V') {
      verbose =1;
      continue;
    }

    // -argum=value
    if (argv[i][2]!='=') {
      ourerror("Equal symbol needed at argument, argument structure for booleans: -V or for strings/numbers -C=<command>");
    }

    switch (c) {
    case 'C':
      strcpy(command,&argv[i][3]);
      break;
    case 'P':
      nthreads = atoi(&argv[i][3]);
      break;
    case 'B':
      nbatch = atoi(&argv[i][3]);
      break;
    case 'N':
      itemstodo = atoi(&argv[i][3]);
      break;
    default:
      ourerror("Unknown parameter");
    }
  }
  if (command[0]=='N') {
    ourerror("Please supply a command using the -C=<command> switch");
  }
  printf("Requesting %d threads\n", nthreads);
  printf("Batchsize = %d\n", nbatch);
  printf("itemstodo = %d\n", itemstodo);
  itemspbatch = nbatch*nthreads;
  printf("# items per run = %d\n", itemspbatch);
  nroftodos = (int)ceil((float)itemstodo / (nbatch*nthreads));
  printf("# runs = %d\n", nroftodos);
  lcores = ((itemstodo % (nbatch*nthreads))/nbatch);
  printf("# threads in last run = %d\n", lcores);
  if (lcores > 0 ) {
    lbatch = itemstodo-(itemspbatch*(nroftodos-1)+(lcores-1)*nbatch);
  } else {
    lbatch = itemstodo-(itemspbatch*(nroftodos-1));
  }
  printf("# items in last run = %d\n", lbatch);
  printf("Command = %s\n", command);
  omp_set_num_threads(nthreads);
  for (int x=0;x<nroftodos;x++) {
    if (x==(nroftodos-1) && lbatch != 0) {
      nthreads=(lcores+1);
      if (verbose) {
        printf("l-cores set\n");
      }
    }
#pragma omp parallel shared(n,command)
    {
      thread_id = omp_get_thread_num();
#pragma omp for
      for (int i=0; i<nthreads; i++) {
        if (lbatch > 0 && x==(nroftodos-1)) {
          nbatch=lbatch;
          if (verbose) {
            printf("l-batch set\n");
          }
        }
        int ret = threader(i,nbatch,i*nbatch+x*itemspbatch,(i+1)*nbatch+x*itemspbatch,thread_id,command,verbose);
        if (ret != 1) {
          printf("Threader %d on thread %d doing jobs [%d...%d] produced an error at %d\n",x,i,i*nbatch+x*itemspbatch,(i+1)*nbatch+x*itemspbatch,ret);
        } else {
          printf("Threader %d on thread %d doing jobs [%d...%d] completed\n",x,i,i*nbatch+x*itemspbatch,(i+1)*nbatch+x*itemspbatch);
        }
      }
    }
  }
  return 1;
}
