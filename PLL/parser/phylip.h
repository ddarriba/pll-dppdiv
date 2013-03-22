#ifndef PHYLIP_H
#define PHYLIP_H

#define PHYLIP_KEEP_UNIQUE      1 << 0
#define PHYLIP_LEX_SORT         1 << 1
#define PHYLIP_SITE_WEIGHTS     1 << 2
#define PHYLIP_DNA_DATA         1 << 3
#define PHYLIP_PROT_DATA        1 << 4


struct pllPhylip
 {
   int              nTaxa;
   int              seqLen;
   char          ** label;
   unsigned char ** seq;
   int            * weights;
 };

struct pllPhylip * pllPhylipParse (const char *);
void pllPhylipRemoveDuplicate (struct pllPhylip *);
void pllPhylipDestroy (struct pllPhylip *);
void usage (const char * cmd_name);
void pllPhylipDump (struct pllPhylip *);

#endif
