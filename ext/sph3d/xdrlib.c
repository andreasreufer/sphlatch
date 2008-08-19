/******************************************************************************
 * File:      xdrlib.cc                                                       *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <rpc/types.h>
#include <rpc/xdr.h>


/* store file pointers for xdr files */
FILE *outfile;
FILE *infile;
XDR   xdrs[1];
int   openforwrite = 0;
int   openforread  = 0;

/* openwritexdrfile_ opens the requested xdrfile for output */
/* particles will be written to this file by calling writeparticle_() */
void openwritexdrfile_(char *outname, char **varnam, int *numvars,
			int *numpart, float *time,  char *title) {
  int  i;
  
  /* make sure file handles arent currently being used */
  if (openforwrite || openforread ) {
    printf("ERROR opening %s - another file is already open\n",outname);
    exit(1);
  }
  
  outfile=fopen(outname,"w");
  if (outfile == NULL) {
    printf("ERROR opening  file %s.\n",outname);
    exit(1);
  }
  printf("Opening file |%s|\n",outname);
  openforwrite = TRUE;
  
  /* dump header */
  fprintf(outfile, "%s\n", title);
  fprintf(outfile, "%d, %d\n", *numpart, -(*numvars));
  fprintf(outfile, "3\n");
  fprintf(outfile, "TIME\n");
  fprintf(outfile, "%f\n", *time);
  fprintf(outfile, "XDRD\n");
  fprintf(outfile, "0\n");
  // If you need double precision change the next line to
  // fprintf(outfile, "DOUBLE\n");
  fprintf(outfile, "FLOAT\n");
  fprintf(outfile, "0\n");
  
  for (i=0; i< *numvars; i++) fprintf(outfile, "%s\n", &varnam[i*7]);
  
  fprintf(outfile, "^^^\n");
  fflush(outfile);
  
  /* now setup for XDR stuff */
  xdrstdio_create(xdrs, outfile, XDR_ENCODE);
}

void writeparticle_(float *store, int *numvars) {
  if (!openforwrite) {
    printf("ERROR - no file open for writing.\n");
    exit(1);
  }
  xdr_vector(xdrs, (char *)store, *numvars, sizeof(float), xdr_float); 
}


void closexdrfile_( void ) {
  if (openforwrite || openforread) fclose(outfile);
  else printf("WARNING - closexdrfile() called, but no files open.\n");

  openforread = openforwrite = FALSE;
}


void openreadxdrfile_(char *inname,  char **varnam, int *numvars,
		       int *numpart, float *time,  char *title) {
  int  i, l, m;
  char tempstr[80];

  if (openforwrite || openforread) {
    printf("ERROR opening %s - another file is already open\n",inname);
    exit(1);
  }
  
  infile=fopen(inname,"r");
  if (infile == NULL) {
    printf("ERROR opening  file %s.\n", inname);
    exit(1);
  }
  
  printf("Opening file |%s|\n", inname);
  openforread = TRUE;
  
  /* read header */
  fscanf(infile, "%s\n", title);
  
  /* get number of particles and number of variables stored in file */
  i=fscanf(infile,"%d, %d", &l, &m);
  *numpart =  l;
  *numvars = -m;
  printf("numpart, numvars = %d  %d\n", *numpart, *numvars);
  
  /* skipover some lines */
  fscanf(infile,"%s\n", tempstr);
  fscanf(infile,"%s\n", tempstr);

  /* scan in time of dump */
  fscanf(infile,"%f\n",  time);
  printf("Time = %f\n", *time);

  /* skipover some lines */
  fscanf(infile,"%s\n", tempstr);
  fscanf(infile,"%s\n", tempstr);  

  /* skip over end of header */
  fscanf(infile,"%s\n", tempstr);
  fscanf(infile,"%s\n", tempstr);
      
  /* grab names of variables */
  for (i = 0; i < *numvars; i++) fscanf(infile, "%s\n", &varnam[i*7]);

  /* grab ^^^*/
   fscanf(infile,"%s\n", tempstr); 
  fflush(infile);

  /* now do XDR stuff */
  xdrstdio_create(xdrs, infile, XDR_DECODE);
}

void readparticle_(float *store, int *numvars) {
  if (!openforread) {
    printf("ERROR - no file open for writing.\n");
    exit(1);
  }
  xdr_vector(xdrs, (char *)store, *numvars, sizeof(float), xdr_float); 
}

void writecracks_(float *epsmin, float *xm, int *nflaws, float *acoef, float *young,
			      float *grav) {
	 static FILE * fd0=NULL;
     if (fd0==NULL) 
	 	fd0=fopen("in.cracks","w");

	 fwrite(epsmin, sizeof(*epsmin),1,fd0);
	 fwrite(xm, sizeof(*xm),1,fd0);
	 fwrite(nflaws, sizeof(*nflaws),1,fd0);
	 fwrite(acoef, sizeof(*acoef),1,fd0);
	 fwrite(young, sizeof(*young),1,fd0);
	 fwrite(grav, sizeof(*grav),1,fd0);
/* ascii output
	 fprintf(fd0,"%g\t%g\t%d\t%g\t%g\t%g\n",*epsmin,*xm,*nflaws,*acoef,*young,*grav);
*/
}	 

