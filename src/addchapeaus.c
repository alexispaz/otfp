#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapeau_obj.h"
#include "chapeau.h"
// Este programa lee los .bsp.retart y saca un .bsp que corresponde al restart
// todo, que lea tambien el .bsp y saque un .bsp en el tiempo 

//#define _PARANOIA_ 1
//


int main ( int argc, char * argv[] ) {
  chapeau ** ch;
  int i,j;
  
  
  ch=(chapeau**)malloc((argc-1)*sizeof(chapeau*));
  
  for (i=1;i<argc;i++) {
    fprintf(stdout,"reading %s\n",argv[i]); 
    ch[i-1]=chapeau_allocloadstate(argv[i]);
  
    for (j=0;j<ch[i-1]->dm;j++) {
      fprintf(stderr,"--- ch dm %d/%d\n",j,ch[i-1]->dm); 
      fprintf(stderr,"--- ch N %d\n",ch[i-1]->N[j]); 
      fprintf(stderr,"--- ch rmin %.5f\n",ch[i-1]->rmin[j]); 
      fprintf(stderr,"--- ch rmax %.5f\n",ch[i-1]->rmax[j]); 
      fprintf(stderr,"--- ch dr %.5f\n",ch[i-1]->dr[j]); 
      fprintf(stderr,"--- ch idr %.5f\n",ch[i-1]->idr[j]); 
    }
     
  }
  
  // Output of the first chapeau before add in it
  ch[0]->nupdate=1;
  chapeau_solve(ch[0]);
  chapeau_output(ch[0],1);

  //Prepare the output of the rest to the same file
  for (i=1;i<argc-1;i++) {
    ch[i]->outputFreq =ch[0]->outputFreq;
    ch[i]->outputLevel=ch[0]->outputLevel;
    ch[i]->ofp=ch[0]->ofp;
    ch[i]->nupdate=1;
  }

  // start adding and output of each chapeau
  for (i=1;i<argc-1;i++) {
    chapeau_solve(ch[i]);
    chapeau_output(ch[i],1);
    chapeau_sum(ch[0],ch[i]);
  }

  // Output of the sum
  chapeau_setupoutput(ch[0],"added.bsp","added",1,1);
  chapeau_solve(ch[0]);
  chapeau_output(ch[0],1);
  chapeau_savestate (ch[0],"added.ch");
                   
  fflush(ch[0]->ofp);
                     
  return 0;

}
