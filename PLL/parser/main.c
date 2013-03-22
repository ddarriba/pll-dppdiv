/*
   Parser usage example file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"

int 
main (int argc, char * argv[])
{
  struct pllPhylip * phylip;

  if (argc != 2)
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }
  
  phylip = pllPhylipParse (argv[1]);
  if (!phylip) 
   {
     printf ("Error while parsing\n");
     return (EXIT_FAILURE);
   }

  printf ("Taxa: %d SeqLen: %d\n", phylip->nTaxa, phylip->seqLen);
  
  
  pllPhylipDump (phylip);
  pllPhylipRemoveDuplicate (phylip);
  
  pllPhylipDestroy (phylip);


  return (EXIT_SUCCESS);
}
