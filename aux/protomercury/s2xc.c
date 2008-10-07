#include <stdio.h>
#include <string.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

/* store file pointers for xdr files */
FILE *outfile;
FILE *infile;
XDR  xdrs[1];
int  openforwrite=0;
int  openforread=0;

/* openwritexdrfile_ opens the requested xdrfile for output */
/* particles will be written to this file by calling writeparticle_() */
void openwritexdrfile_( char *outname,  char **varnam, int *numvars,
	        int *numpart, double *dtime, char *title)
{
  int  i;

  float  ft   = (float)(*dtime);
  float *time = &ft;

  printf("dtime = %f\n", *dtime);

  printf("time = %f\n", *time);

  /* make sure file handles arent currently being used */
  if (openforwrite || openforread ) 
    {
      printf("ERROR opening %s - another file is already open\n",outname);
      exit;
    }
  
  outfile=fopen(outname,"w");
  if (outfile==NULL)
    {
      printf("ERROR opening  file %s.\n",outname);
      exit;
    }
  printf("Opening file |%s|\n",outname);
  openforwrite=TRUE;
  
  /* dump header */
  fprintf(outfile,"%s\n",title); /* title */
  fprintf(outfile,"%d, %d\n", *numpart, -(*numvars));
  fprintf(outfile,"2\n");
  fprintf(outfile,"TIME\n");
  fprintf(outfile,"%f\n", *time);
  fprintf(outfile,"XDRD\n");
  fprintf(outfile,"0\n");
  
  for (i=0; i< *numvars; i++)
    {
      fprintf(outfile,"%s\n",&varnam[i*7]);
    }
  
  fprintf(outfile,"^^^\n");
  fflush(outfile);
  
  /* now setup for XDR stuff */
  xdrstdio_create(xdrs, outfile, XDR_ENCODE);

} /*openwritefile*/

void writeparticle_( float *store, int *numvars)
{
  /* make sure file is open */
  if (!openforwrite) 
    {
      printf("ERROR - no file open for writing.\n");
      exit;
    }
  xdr_vector(xdrs, (char *)store, *numvars, sizeof(float),
	     xdr_float);
  
} /*writeparticle_*/


void closexdrfile_( void )
{
  /* all done, so close file */
  if (openforwrite || openforread)
    fclose(outfile);
  else
    printf("WARNING - closexdrfile() called, but no files open.\n");

  openforread=openforwrite=FALSE;
} /* closexdrfile_ */


void openreadxdrfile_( char *inname,  char **varnam, int *numvars,
	        int *numpart, double *dtime, char *title)
{
  int  i, l, m;
  char tempstr[80];
  float *time;

  /* make sure file handles arent currently being used */
  if (openforwrite || openforread) 
    {
      printf("ERROR opening %s - another file is already open\n",inname);
      exit;
    }
  
  infile=fopen(inname,"r");
  if (infile==NULL)
    {
      printf("ERROR opening  file %s.\n",inname);
      exit;
    }
  
  printf("Opening file |%s|\n",inname);
  openforread=TRUE;
  
  /* read header */
  fscanf(infile,"%s\n",title);
  
  /* get number of particles and number of variables stored in file */
  i=fscanf(infile,"%d, %d", &l, &m);
  *numpart = l;
  *numvars = -m;
  printf("numpart, numvars = %d  %d\n",*numpart,*numvars);
  
  /* skipover some lines */
  fscanf(infile,"%s\n",tempstr);
  fscanf(infile,"%s\n",tempstr);

  /* scan in time of dump */
  fscanf(infile,"%f\n", time);
  printf("Time = %f\n",*time);
  
  *dtime = *time;

  /* skip over end of header */
  fscanf(infile,"%s\n",tempstr);
  fscanf(infile,"%s\n",tempstr);
      
  /* grab names of variables */
  for (i=0; i< *numvars; i++)
    {
      fscanf(infile,"%s\n",&varnam[i*7]);
    }
  fflush(infile);

  /* now do XDR stuff */
  xdrstdio_create(xdrs, infile, XDR_DECODE);
} /* openreadxdrfile_ */


void readparticle_( float *store, int *numvars)
{
  /* make sure file is open */
  if (!openforread) 
    {
      printf("ERROR - no file open for writing.\n");
      exit;
    }
  xdr_vector(xdrs, (char *)store, *numvars, sizeof(float),
	     xdr_float);
  
} /*readparticle_*/













