/*
 ** companion.c --- Version 1.0, 07.27.04 (ZML and DCR)
 ** searches for binaries 
 */

/*
 ** This code is Copyright &copy; 2004 by Z. M. Leinhardt and D. C. 
 ** Richardson.
 ** Under the terms of the GNU Public License, you are free to	
 ** redistribute, modify, or even sell this code, but we ask that all
 ** headers identifying the original authors of this code be left intact.
 **
 ** This software comes with NO WARRANTY.  The authors cannot be held
 **	responsible for any undesirable consequences of using this code.
 */

/*
 ** Memory storage for particles will be inflated by the following
 ** factor.  This is needed because we cannot realloc() the particle
 ** storage (in order to add com particles) after the final data read;
 ** otherwise pointers in the tree and binary list to the particle data
 ** won't work.  A better method is to store particle array indices in
 ** the tree and binary list data fields and pass a pointer to the
 ** particle array to whatever routines need the particle data.  This is
 ** annoying because the functions affected include recursive tree
 ** functions and the sorting comparison functions, so for the moment we
 ** stick with this crude buffer inflation and hope it's enough storage
 ** (an assert() is used to make sure).
 */
#define EXTRA_STORE 2

/*
 ** Fraction of particles allowed to remain since last tree build before
 ** rebuilding the tree.  Used for hierarchical system search.  A large
 ** value forces more frequent tree rebuilds.  A smaller value relies
 ** on older tree data for longer, with larger inefficiencies.  More
 ** tests are needed to optimize this value, though it may vary
 ** depending on the specific problem.
 */
#define TREE_REBUILD_FRAC 0.9

/*
 ** If mass ratio between two components of a com particle is extreme
 ** force a tree rebuild.
 */
#define REBUILD_MASS_RATIO 1.0e6

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() */
#include <string.h>
#include <math.h>
#include <limits.h> /* may or may not contain DBL_MAX (values.h obsolete?) */
#include <assert.h>

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ANA_EXT ".ana"
#define PR_EXT ".pr"
#define HIER_EXT ".hier"
#define EXT_EXT ".ext"
#define HIER_EXT_EXT ".hext"
#define VEC_EXT ".vec"

#define BUF_SIZE_INIT 256
#define BUF_SIZE_MULT 2
#define BUF_NPART_INIT 10000

#define CHILD_PER_NODE 8 /* don't change this! (oct-tree) */

#define FILE_TYPE_STR_MAX_LEN 4 /* 3 chars plus null char */

#ifdef SS_CORE

#include <rpu.h>
#include <ss.h>
#include <vector.h>

enum {FileTypeTxt,FileTypeBin,FileTypeSS};
const char FileTypeStr[][FILE_TYPE_STR_MAX_LEN] = {"txt","bin","ss"};

#define DFLT_FILE_TYPE FileTypeSS
#define DFLT_IN_SYS TRUE
#define DFLT_IN_CGS FALSE
#define DFLT_OUT_SYS TRUE
#define DFLT_OUT_CGS FALSE
#define DFLT_LENGTH_CONV	(1.0)		/* default in sys units if ss_core defined */
#define DFLT_MASS_CONV		(1.0)
#define DFLT_TIME_CONV		(1.0)
#else /* define macros and types from ss_core */

#define N_DIM 3
#define TWO_PI (2*M_PI)

#define BOOLEAN int
#define FALSE 0
#define TRUE 1

/* Fundamental constants from 1998 Astronomical Almanac, in mks */

#define AU		1.49597870e11	/* One A.U. in metres [1.4959787066e11] */
#define M_SUN	1.9891e30		/* Solar mass in kilograms */
#define SID_YR	3.15581497632e7	/* One sidereal year in seconds (1998.0) */

typedef double VECTOR[N_DIM];

#define X 0
#define Y 1
#define Z 2

/* Vector function prototypes */

void SET_VEC(VECTOR,double,double,double);
void ZERO_VEC(VECTOR);
void COPY_VEC(VECTOR,VECTOR);
void ADD_VEC(VECTOR,VECTOR,VECTOR);
void SUB_VEC(VECTOR,VECTOR,VECTOR);
void SCALE_VEC(VECTOR,double);
void NORM_VEC(VECTOR,double);
double DOT(VECTOR,VECTOR);
void CROSS(VECTOR,VECTOR,VECTOR);
double MAG_SQ(VECTOR);
double MAG(VECTOR);

/* Assigns a value (x,y,z) to vector v */

#define SET_VEC(v,x,y,z) {\
  (v)[X] = (x);\
  (v)[Y] = (y);\
  (v)[Z] = (z);\
}

/* Assigns zero to vector v */

#define ZERO_VEC(v) SET_VEC((v),0,0,0)

/* Copies vector v1 to vector v2 */

#define COPY_VEC(v1,v2) {\
  (v2)[X] = (v1)[X];\
  (v2)[Y] = (v1)[Y];\
  (v2)[Z] = (v1)[Z];\
}

/* Adds vectors v1 & v2 and puts the result in vector v */

#define ADD_VEC(v1,v2,v) {\
  (v)[X] = (v1)[X] + (v2)[X];\
  (v)[Y] = (v1)[Y] + (v2)[Y];\
  (v)[Z] = (v1)[Z] + (v2)[Z];\
}
/* Subtracts vector v2 from vector v1 and puts the result in vector v */

#define SUB_VEC(v1,v2,v) {\
  (v)[X] = (v1)[X] - (v2)[X];\
  (v)[Y] = (v1)[Y] - (v2)[Y];\
  (v)[Z] = (v1)[Z] - (v2)[Z];\
}

/* Multiplies vector v by scalar a */

#define SCALE_VEC(v,a) {\
  double _scalar = (a);\
  (v)[X] *= _scalar;\
  (v)[Y] *= _scalar;\
  (v)[Z] *= _scalar;\
}

/* Divides vector v by scalar a */

#define NORM_VEC(v,a) SCALE_VEC((v),1.0/(a))

/* Returns dot product of vectors v1 & v2 */

#define DOT(v1,v2) ((v1)[X]*(v2)[X] + (v1)[Y]*(v2)[Y] + (v1)[Z]*(v2)[Z])

/* Returns cross product of vectors v1 & v2 in vector v */

#define CROSS(v1,v2,v) {\
  (v)[X] = (v1)[Y]*(v2)[Z] - (v1)[Z]*(v2)[Y];\
  (v)[Y] = (v1)[Z]*(v2)[X] - (v1)[X]*(v2)[Z];\
  (v)[Z] = (v1)[X]*(v2)[Y] - (v1)[Y]*(v2)[X];\
}

/* Returns square magnitude of vector v */

#define MAG_SQ(v) (DOT((v),(v)))

/* Returns magnitude of vector v */

#define MAG(v) (sqrt(MAG_SQ(v)))

/* Struct for storing particle data */

typedef struct {
	double mass;
	double radius;
	double pos[N_DIM];
	double vel[N_DIM];
	double spin[N_DIM];
	int color;
	int org_idx;  
	} SSDATA;

enum {FileTypeTxt,FileTypeBin};
const char FileTypeStr[][FILE_TYPE_STR_MAX_LEN] = {"txt","bin"};

#define DFLT_FILE_TYPE FileTypeTxt
#define DFLT_IN_SYS FALSE
#define DFLT_IN_CGS TRUE
#define DFLT_OUT_SYS FALSE
#define DFLT_OUT_CGS TRUE
#define DFLT_LENGTH_CONV	(AU*1.0e2)		/* default in cgs units */
#define DFLT_MASS_CONV		(M_SUN*1.0e3)
#define DFLT_TIME_CONV		(SID_YR/TWO_PI)
#define SQ(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#endif /* !SS_CORE */

#define SUPPORT_OLD_FORMAT

/* defaults */

#define DFLT_HIER FALSE
#define DFLT_TIPSY_FILE FALSE
#define DFLT_IN_MKS FALSE
#define DFLT_OUT_MKS FALSE
#define DFLT_DO_PERI_CUT FALSE
#define DFLT_DO_ORBIT_CUT FALSE
#define DFLT_APPLY_HIER_CUT FALSE
#define DFLT_ECC_CUT 1.0
#define DFLT_ENG_CUT 0.0
#define DFLT_HILL_CUT 0.0
#define DFLT_PERI_CUT 0.0
#define DFLT_ORBIT_CUT 0.0
#define DFLT_OPEN_ANGLE 0.5
#define DFLT_HIER_EXTRACT_INDEX (-1)
#define DFLT_EXTRACT_INDEX (-1)

typedef struct {
	/* args */
	BOOLEAN Hier,TipsyFile,InCgsUnits,InMksUnits,InSysUnits,OutCgsUnits,OutMksUnits,OutSysUnits;
	int FileType;
	long ExtIdx,HierExtIdx;
	double EccCut,EngCut,HillCut,PeriCut,OrbitCut,OpenAng;
	/* prompted */
	double avga,starmass;
	/* derived */
	BOOLEAN DoPeriCut,DoOrbitCut,ApplyHierCut;
	double InLengthConv,InMassConv,InTimeConv,InEnergyConv/*NOT USED*/;
	double OutLengthConv,OutMassConv,OutTimeConv,OutEnergyConv;
	} PARAMS;

struct compdata { /* contains particle info and a pointer to parent and children */
	struct compdata *prim;
	struct compdata *sat;
	struct compdata *com;
	SSDATA data;
	long index;
	};

typedef struct compdata COMPDATA;

typedef struct { /* contains the primary and satellite of a bound pair */
	COMPDATA *prim;
	COMPDATA *sat;
	double period; /* only used for hierarchical grouping */
	} BINARY;

typedef struct { /* contains summary data for each hierarchical system */
	const COMPDATA *com;
	double sys_mass;
	double max_a;
	double bind_E;
	} HIER_OUTPUT; /* only used in write_hier_output */

struct node {
	VECTOR pos;
	VECTOR vel;
	double size,mass;
	double half_size, eff_half_size,eff_size_sq;
	struct node *child[CHILD_PER_NODE];
	COMPDATA *leaf[CHILD_PER_NODE];
	int n_part;
	};

typedef struct node NODE;

#define BIN_ENERGY		(1 << 0) /* do not change these! */
#define BIN_ANGMOM		(1 << 1)
#define BIN_SEMI		(1 << 2 | BIN_ENERGY)
#define BIN_PERIOD		(1 << 3 | BIN_SEMI)
#define	BIN_ECC			(1 << 4 | BIN_ANGMOM | BIN_SEMI)
#define BIN_PERIAPSE	(1 << 5 | BIN_SEMI | BIN_ECC)
#define BIN_INCL		(1 << 6 | BIN_ANGMOM)

#define BIN_CUT (BIN_ENERGY|BIN_SEMI|BIN_ECC|BIN_PERIAPSE) /* (BIN_SEMI for Hill cut) */
#define BIN_ALL (BIN_ENERGY|BIN_ANGMOM|BIN_SEMI|BIN_PERIOD|BIN_ECC|BIN_PERIAPSE|BIN_INCL)

int BIT_ON(int,int);
#define BIT_ON(flag,mask) (((flag) & (mask)) == (mask))

typedef struct {
	double E; /* total energy */
	double a; /* semimajor axis */
	double P; /* period */
	double e; /* eccentricity */
	double q; /* periapse */
	double i; /* inclination */
	} BIN_STATS;

/* handy macros for hierarchical system search */

BOOLEAN IS_COM_MBR(const COMPDATA *);
#define IS_COM_MBR(p) ((p)->com != NULL)

BOOLEAN HAS_COM_MBR(const COMPDATA *);
#define HAS_COM_MBR(p) ((p)->prim != NULL || (p)->sat != NULL)

/*BOOLEAN IS_COM_PART(const COMPDATA *);
#define IS_COM_PART(p) ((HAS_COM_MBR(p) && !IS_COM_MBR(p)))*/

BOOLEAN IS_SAME_BINARY(const BINARY *,const BINARY *);
#define IS_SAME_BINARY(b1,b2) ((b1)->prim == (b2)->prim && (b1)->sat == (b2)->sat)

void calc_bin_stats(const SSDATA *p1,const SSDATA *p2,int flag,BIN_STATS *bs)
{
	VECTOR r,v;
	double M_inv;

	assert(p1 && p2 && bs);
	assert(flag > 0);

	SUB_VEC(p1->pos,p2->pos,r);
	SUB_VEC(p1->vel,p2->vel,v);

	M_inv = p1->mass + p2->mass; /* not needed if just BIN_INCL */
	assert(M_inv > 0.0);
	M_inv = 1.0/M_inv;

	if (BIT_ON(flag,BIN_ENERGY)) {
		double mag_r_inv,v2overM;

		mag_r_inv = MAG(r);
		assert(mag_r_inv > 0.0);
		mag_r_inv = 1.0/mag_r_inv;
		v2overM = MAG_SQ(v)*M_inv;
		bs->E = p1->mass*p2->mass*(0.5*v2overM - mag_r_inv);
		assert(bs->E < 0.0); /* must be bound to be a binary! */
		if (BIT_ON(flag,BIN_SEMI)) {
			bs->a = 1.0/(2.0*mag_r_inv - v2overM);
			assert(bs->a > 0.0);
			}
		if (BIT_ON(flag,BIN_PERIOD)) {
			bs->P = TWO_PI*sqrt(CUBE(bs->a)*M_inv);
			assert(bs->P > 0.0);
			}
		}

	if (BIT_ON(flag,BIN_ANGMOM)) {
		VECTOR h;
		double h2;

		CROSS(r,v,h); /* ang mom per unit reduced mass */
		h2 = MAG_SQ(h);
		if (BIT_ON(flag,BIN_ECC)) {
			double x = h2*M_inv/bs->a;
			assert(x <= 1.0);
			bs->e = sqrt(1 - x);
			assert(bs->e >= 0.0 && bs->e <= 1.0); /* allow e=1 for now, will be cut */
			}
		if (BIT_ON(flag,BIN_INCL)) {
			double mag_h = sqrt(h2);
			assert(mag_h > 0.0);
			bs->i = acos(h[Z]/mag_h);
			}
		}

	if (BIT_ON(flag,BIN_PERIAPSE)) {
		bs->q = bs->a*(1 - bs->e);
		assert(bs->q >= 0.0); /* allow q=0 for now, will be cut */
		}
	}

int sort_bin(const void *a,const void *b)
{
	/* sort function for "normal" output (not hierarchical) */

	const BINARY *b1,*b2;
	const SSDATA *p1,*p2;

	b1 = (const BINARY *) a; /* pointers to binaries */
	b2 = (const BINARY *) b;

	p1 = &(b1->prim->data); /* pointers to primaries */
	p2 = &(b2->prim->data);

	if (p1->mass < p2->mass) return 1; /* largest primary masses first */
	if (p1->mass > p2->mass) return -1;

	if (b1->prim->index < b2->prim->index) return -1; /* smallest indices first */
	if (b1->prim->index > b2->prim->index) return 1;

	/* within each system, sort by binding energy */

	{
	const SSDATA *s1 = &(b1->sat->data),*s2 = &(b2->sat->data); /* satellite pointers */
	BIN_STATS bs1,bs2;

	calc_bin_stats(p1,s1,BIN_ENERGY,&bs1);
	calc_bin_stats(p2,s2,BIN_ENERGY,&bs2);

	if (bs1.E < bs2.E) return -1; /* more bound first */
	if (bs1.E > bs2.E) return 1;

	if (b1->sat->index < b2->sat->index) return -1; /* smallest indices first */
	if (b1->sat->index > b2->sat->index) return 1;
	}

	assert(0); /* shouldn't be here (duplicates not allowed) */

	return 0;
	}

double hill(const PARAMS *p,double mass_i,double mass_c)
{
	/* Calculate and return value of Hill sphere */

	return pow((mass_i + mass_c)/(3.0*(p->starmass + mass_i + mass_c)),
			   (1.0/3.0))*p->avga;  
	}

const char *myBasename(const char *path)
{
	char *p;
  
	assert(path);
	p = strrchr(path,'/');
	if (p) return p + 1;
	else return path;
	}


int myNewExt(const char *infile,const char *inext,
			 char *outfile,const char *outext)
{
	/* adds (or replaces) extension to filename */

	const char *basename;
	char *c;
	size_t n;
	
	assert(infile && inext && outfile && outext);
	basename = myBasename(infile);
	if ((c = strrchr(basename,'.')) && strstr(c,inext))
		n = c - basename;
	else
		n = strlen(basename);
	if (n + strlen(outext) >= (size_t) MAXPATHLEN)
		return 1;
	(void) strncpy(outfile,basename,n); /* not null terminated */
	(void) strcpy(outfile + n,outext);
	return 0;
	}

void add_to_list(COMPDATA *prim,COMPDATA *sat,BINARY **list,long *list_size,long *list_posn) 
{
	/* Add bound system to list */
	COMPDATA *temp;

	if (*list_posn >= *list_size) { 
		(void) printf("Growing list space\n");
		if (*list_size == 0) {
			assert(*list == NULL);
			*list_size = BUF_SIZE_INIT;
			}
		else {
			assert(*list != NULL);
			*list_size *= BUF_SIZE_MULT;
			}
		*list = (BINARY *) realloc((void *) (*list),(*list_size)*sizeof(BINARY));
		assert(*list != NULL);
		(void) printf("New list size = %li\n",*list_size);
		}
/*	if ((IS_COM_PART(prim) || IS_COM_PART(sat)) && (prim->data.mass<sat->data.mass)) {*/
	if ((HAS_COM_MBR(prim) || HAS_COM_MBR(sat)) && prim->data.mass < sat->data.mass) {
		temp = prim;
		prim = sat;
		sat = temp;
		} /*DEBUG*/
	else if (!HAS_COM_MBR(prim) && !HAS_COM_MBR(sat)) /* make sure that normal particle */
		assert(prim->data.mass >= sat->data.mass);     /* binaries have prim > sat mass */

	(*list)[*list_posn].prim = prim;
	(*list)[*list_posn].sat = sat;
	++(*list_posn);
	}

void find_companion(const PARAMS *p,const NODE *node,COMPDATA *part,
					BINARY **list,long *list_size,long *list_posn)
{
	/* Identifies particles with speeds less than the escape speed */

	const SSDATA *pd,*ld;
	VECTOR r,v;
	double r2,v2;
	int i;

	pd = &part->data;

	for (i=0;i<CHILD_PER_NODE;i++)
		if (node->child[i] != NULL) {
			SUB_VEC(node->child[i]->pos,pd->pos,r);
			r2 = MAG_SQ(r);
			assert(r2 > 0.0);
			if (node->child[i]->eff_size_sq/r2 > p->OpenAng)
				find_companion(p,node->child[i],part,list,list_size,list_posn);
			else {
				SUB_VEC(node->child[i]->vel,pd->vel,v); /*v=rel vel*/
				v2 = MAG_SQ(v);
				/* do it this way to avoid sqrt()s... */
				if (v2*v2 < 4.0*SQ(pd->mass + node->child[i]->mass)/r2)
					find_companion(p,node->child[i],part,list,list_size,list_posn);
				}
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != part && !IS_COM_MBR(node->leaf[i])) {
			/* (note: member check in previous line only necessary for hierarchical search) */
			ld = &node->leaf[i]->data;
/*			if (!IS_COM_PART(part) && (ld->mass > pd->mass ||*/ /*DCR please check*/
			if (!HAS_COM_MBR(part) && (ld->mass > pd->mass ||
				(ld->mass == pd->mass && node->leaf[i]->index < part->index)))
				continue; /* to prevent double counting, but only for non-hierarchical searching case  */
			SUB_VEC(pd->pos,ld->pos,r);
			SUB_VEC(pd->vel,ld->vel,v);
			r2 = MAG_SQ(r);
			assert(r2 > 0.0);
			v2 = MAG_SQ(v);
			if (v2*v2 < 4.0*SQ(pd->mass + ld->mass)/r2)
				add_to_list(part,node->leaf[i],list,list_size,list_posn);
			}
	}

void get_com_vel(NODE *node)
{
	/*gets vel moments, com vel of node, and total mass in node*/

	VECTOR v;
	double m;
	int i;

	node->mass = 0.0;
	ZERO_VEC(node->vel);
	node->n_part = 0;
	for (i=0;i<CHILD_PER_NODE;i++){
		if (node->child[i]) {
			get_com_vel(node->child[i]);
			m = node->child[i]->mass;
			node->mass += m;
			COPY_VEC(node->child[i]->vel,v);
			SCALE_VEC(v,m);
			ADD_VEC(node->vel,v,node->vel);
			node->n_part += node->child[i]->n_part;
			}
		else if (node->leaf[i]) {
			m = node->leaf[i]->data.mass;
			node->mass += m;
			COPY_VEC(node->leaf[i]->data.vel,v);
			SCALE_VEC(v,m);
			ADD_VEC(node->vel,v,node->vel);
			++node->n_part;
			}
		}
	assert(node->mass > 0.0);
	NORM_VEC(node->vel,node->mass);
	assert(node->n_part > 0);
	}

void make_node(const VECTOR pos,double size,NODE **node)
{
	/* creates new nodes for the tree */

	int i;

	assert(size > 0.0);

	*node = (NODE *) malloc(sizeof(NODE)); /*make space for a node*/
	assert(*node != NULL);

	COPY_VEC(pos,(*node)->pos);
	ZERO_VEC((*node)->vel);
	(*node)->mass = 0.0;
	(*node)->size = size;
	(*node)->half_size = (*node)->eff_half_size = 0.5*size;
	assert((*node)->size > 0.0); /* check for underflow */
	(*node)->eff_size_sq = SQ(size);
	assert((*node)->eff_size_sq > 0.0); /* ditto */

	for (i=0;i<CHILD_PER_NODE;i++) {
		(*node)->child[i] = NULL; /*set children and leaf pointers to null*/ 
		(*node)->leaf[i] = NULL;
		}
	}

void add_to_tree(NODE *node,COMPDATA *p)
{
	/* adds particles to tree */

	int i,idx,idy,idz;

	idx = (p->data.pos[X] < node->pos[X] ? -1 : 1); /*locates the particle in one */
	idy = (p->data.pos[Y] < node->pos[Y] ? -1 : 1); /* of eight quadrants */
	idz = (p->data.pos[Z] < node->pos[Z] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1)); /*sets i=0-7 depending on quadrant*/
	if (node->child[i])  /*if node contains children open the node*/
		add_to_tree(node->child[i],p);
	else if (node->leaf[i]) {   /*if node already contains a particle*/
		VECTOR v;                 /*create children */
		SET_VEC(v,idx,idy,idz);
		SCALE_VEC(v,0.5*node->half_size);
		ADD_VEC(v,node->pos,v);
		make_node(v,node->half_size,&(node->child[i]));
		add_to_tree(node->child[i],node->leaf[i]);
		add_to_tree(node->child[i],p);
		node->leaf[i] = NULL;
		}
	else {
		node->leaf[i] = p; /*if particle is in an empty node make it a leaf*/
		if (p->data.radius > node->eff_half_size) { /*if particle is large make cell*/
			node->eff_half_size = p->data.radius;    /*large - scales with mass*/
			node->eff_size_sq = 4.0*SQ(p->data.radius); 
			}
		}
	}

void kill_node(NODE *node)
{
	/* nodes are no longer needed: release memory used for nodes */

	int i;

	assert(node != NULL);

	for (i=0;i<CHILD_PER_NODE;i++)
		if (node->child[i])
			kill_node(node->child[i]);

	free((void *) node);
	}

int read_data(const PARAMS *p,const char *file_in,COMPDATA **part,long *n,
			  double *m_tot,VECTOR root_center,double *root_size)
{
	/* read data from ss file */

	SSDATA *d;
	FILE *fp; /* for txt & bin file types only */
	double xmin,ymin,zmin,xmax,ymax,zmax;
	long i;

#ifdef SS_CORE
	SSIO ssio; /* for ss file type only */
	SSHEAD h;
#endif

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		if ((fp = fopen(file_in,"r"))  == NULL) {
			(void) fprintf(stderr,"Problem opening file %s\n",file_in);
			return 1;
			}
		*n = BUF_NPART_INIT;
		break;
#ifdef SS_CORE
	case FileTypeSS:
		if (ssioOpen(file_in,&ssio,SSIO_READ)) {
			(void) fprintf(stderr,"Unable to open %s for reading\n",file_in);
			return 1;
			}
		if (ssioHead(&ssio,&h)){
			(void) fprintf(stderr,"Corrupt header\n");
			(void) ssioClose(&ssio);
			return 1;
			}
		if (h.n_data <= 0) {
			(void) fprintf(stderr,"Invalid file size\n");
			(void) ssioClose(&ssio);
			return 1;
			}
		*n = h.n_data;
		break;
#endif
	default:
		(void) fprintf(stderr,"file type %i is invalid\n",p->FileType);
		return 1;
		}

	*part = (COMPDATA *) malloc((*n)*sizeof(COMPDATA)); /*allocate space for part*/
	assert(*part != NULL); /*make sure space allocation worked*/

	*m_tot = 0.0;

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		i=0;
		while (feof(fp) == 0) {
			if (i>=(*n)) {
				(void) printf("Growing particle space\n");
				*n *= BUF_SIZE_MULT;
				*part = (COMPDATA *) realloc((void *) (*part),(*n)*sizeof(COMPDATA));
				assert(*part != NULL);
				(void) printf("New data space = %li\n",*n);
				}
			d = &((*part)[i].data);
			switch (p->FileType) {
			case FileTypeTxt:
				if (fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf\n", 
						   &d->mass,&d->radius,
						   &d->pos[X],&d->pos[Y],&d->pos[Z],
						   &d->vel[X],&d->vel[Y],&d->vel[Z]) != 8) goto error;
				break;
			case FileTypeBin:
				if (fread(&d->mass,sizeof(double),1,fp) != 1) goto error;
				if (fread(&d->radius,sizeof(double),1,fp) != 1) goto error;
				if (fread(d->pos,sizeof(double),3,fp) != 3) goto error;
				if (fread(d->vel,sizeof(double),3,fp) != 3) goto error;
				break;
			default:
				assert(0); /* should never get here */
				}
			d->mass = d->mass/p->InMassConv;
			d->radius = d->radius/p->InLengthConv;
			NORM_VEC(d->pos,p->InLengthConv);
			NORM_VEC(d->vel,p->InLengthConv);
			d->org_idx = i;
			*m_tot += d->mass;
			++i;
			}
		*n=i;
		/* release unused buffer space */
		*part = (COMPDATA *) realloc((void *) (*part),(*n)*sizeof(COMPDATA));
		assert(*part != NULL);
		(void) printf("Number of particles = %li\n",*n);
		break;
#ifdef SS_CORE
	case FileTypeSS:
		for (i=0;i<*n;i++) {
			d = &((*part)[i].data);
			if (ssioData(&ssio,d) != 0) {  
				(void) fprintf(stderr,"Corrupt data\n");
				(void) ssioClose(&ssio);
				return 1;
				}
#ifdef SUPPORT_OLD_FORMAT
			if (d->org_idx == -1) /* assign index for each particle */
				d->org_idx = i;
#endif
			*m_tot += d->mass;
			}
		break;
#endif /* SS_CORE */
	default:
		assert(0); /* invalid file type */
		}

	/* pad storage -- see comment at top of file */
	*part = (COMPDATA *) realloc((void *) (*part),(*n)*EXTRA_STORE*sizeof(COMPDATA));

	xmin = ymin = zmin = DBL_MAX;
	xmax = ymax = zmax = - DBL_MAX;
	for (i=0;i<*n;i++) {
		(*part)[i].prim = (*part)[i].sat = NULL;
		(*part)[i].com = NULL;
		(*part)[i].index = i;
		d = &((*part)[i].data);
		if (d->pos[X] < xmin) xmin = d->pos[X]; /*find max extent of particles*/
		if (d->pos[Y] < ymin) ymin = d->pos[Y]; /*for size of root cell*/ 
		if (d->pos[Z] < zmin) zmin = d->pos[Z];
		if (d->pos[X] > xmax) xmax = d->pos[X];
		if (d->pos[Y] > ymax) ymax = d->pos[Y];
		if (d->pos[Z] > zmax) zmax = d->pos[Z];
		}

	*root_size = xmax - xmin;
	if (ymax - ymin > *root_size) *root_size = ymax - ymin;
	if (zmax - zmin > *root_size) *root_size = zmax - zmin;
	SET_VEC(root_center,(xmin + xmax)/2,(ymin + ymax)/2,(zmin + zmax)/2);
	(void) printf("root center = (%g,%g,%g), size = %g\n",root_center[X],
				  root_center[Y],root_center[Z],*root_size);

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		fclose(fp);
		break;
#ifdef SS_CORE
	case FileTypeSS:
		(void) ssioClose(&ssio);
		break;
#endif
	default:
		assert(0);
		}

	(void) printf("%li particle%s read\n",*n,*n==1 ? "" : "s");
  
	return 0;

 error:
	(void) fprintf(stderr,"Error during read.\n");
	(void) fclose(fp);
	return 1;
	}

void make_com_part(BINARY *tightest,long *n_part,long *part_buf_size,COMPDATA **part)
{
	/* 
	 ** Creates com particle, adds it to the particle list, modifies
	 ** primary and satellite structures so they will no longer be
	 ** used in the search tree.
	 */

	const SSDATA *ptr_prim,*ptr_sat;
	SSDATA *ptr_com;
	COMPDATA *comp_prim,*comp_sat,*comp_com;
	BIN_STATS bs;
	VECTOR v1,v2;

	/* check if buffer needs to grow */

	if (*n_part == *part_buf_size) {
		(void) printf("Growing particle list space\n");
		/* for now, particle list not allowed to grow -- see comment at top */
		assert(0); /* particle realloc() forbidden */
		*part_buf_size *= BUF_SIZE_MULT;
		*part = (COMPDATA *) realloc((void *) (*part),(*part_buf_size)*sizeof(COMPDATA));
		assert(*part != NULL);
		(void) printf("New part list size = %li\n",*part_buf_size);
		}

	/* abbreviation for COMPDATA */

	comp_prim = tightest->prim;
	comp_sat = tightest->sat;
	comp_com = &((*part)[*n_part]); /* new com particle */

	/* abbreviation for SSDATA */

	ptr_prim = &comp_prim->data;
	ptr_sat = &comp_sat->data;
	ptr_com = &comp_com->data;

	/* store pointers */

	comp_com->com = NULL;
	comp_com->prim = comp_prim;
	comp_com->sat = comp_sat;
	comp_com->index = *n_part;

	comp_prim->com = comp_com;
	comp_sat->com = comp_com;

	++(*n_part); /* particle list grows by one for com particle */

	ptr_com->mass = ptr_prim->mass + ptr_sat->mass;
	calc_bin_stats(ptr_prim,ptr_sat,BIN_SEMI,&bs);
	ptr_com->radius = bs.a; /* com "radius" is semimajor axis of binary */

	/* com position */

	COPY_VEC(ptr_prim->pos,v1);
	SCALE_VEC(v1,ptr_prim->mass);
	COPY_VEC(ptr_sat->pos,v2);
	SCALE_VEC(v2,ptr_sat->mass);
	ADD_VEC(v1,v2,ptr_com->pos);
	NORM_VEC(ptr_com->pos,ptr_com->mass);

	/* com velocity */

	COPY_VEC(ptr_prim->vel,v1);
	SCALE_VEC(v1,ptr_prim->mass);
	COPY_VEC(ptr_sat->vel,v2);
	SCALE_VEC(v2,ptr_sat->mass);
	ADD_VEC(v1,v2,ptr_com->vel);
	NORM_VEC(ptr_com->vel,ptr_com->mass);

	/* com spin zeroed (could store orbital ang vel instead) */

	ZERO_VEC(ptr_com->spin);

	/* cycle colors */

	ptr_com->color = 2 + (ptr_prim->color + ptr_sat->color - 2)%14;

	/* assign index of primary */

	ptr_com->org_idx = ptr_prim->org_idx;
	}

/*BOOLEAN ok_to_cut(const PARAMS *p,const SSDATA *pd,const SSDATA *sd,const BIN_STATS *bs)*/
BOOLEAN ok_to_cut(const PARAMS *p,const COMPDATA *prim,const COMPDATA *sat,const BIN_STATS *bs)
{
	/* returns 1 if any list cut criteria met, 0 otherwise */

	if ((p->DoPeriCut && !HAS_COM_MBR(prim) && 
		 bs->q < (p->PeriCut < 0.0 ? - p->PeriCut*prim->data.radius : 
				  p->PeriCut == 0.0 ? prim->data.radius + sat->data.radius :
				  p->PeriCut)) || (p->DoOrbitCut && HAS_COM_MBR(prim) &&
		bs->q < (p->OrbitCut < 0.0 ? - p->OrbitCut*prim->data.radius :
				 p->OrbitCut == 0.0 ? prim->data.radius + 
				 (HAS_COM_MBR(sat) ? sat->data.radius : 0.0) : 
				 p->OrbitCut*prim->data.radius)) ||
		bs->e >= p->EccCut || bs->E >= p->EngCut ||
		(p->HillCut > 0.0 && bs->a >= p->HillCut*hill(p,prim->data.mass,sat->data.mass)))
		return 1;
	else
		return 0;
	}

void cut_list(const PARAMS *p,BINARY *list,long *list_length)
{
	/*cuts out all systems with q<, e>, a>, eng> params recursively*/

	SSDATA *pd,*sd;
	BIN_STATS bs;
	long i,n;
	BOOLEAN *flag;

	flag = (BOOLEAN *) malloc((*list_length)*sizeof(BOOLEAN));

	/*removes binaries from list*/
	for (i=0;i<(*list_length);i++) {
		pd = &list[i].prim->data;
		sd = &list[i].sat->data;
		calc_bin_stats(pd,sd,BIN_CUT,&bs);
		flag[i] = ok_to_cut(p,list[i].prim,list[i].sat,&bs);
		}

	/*shrinks list and counts up total number of surviving binaries*/
	n = 0;
	for (i=0;i<(*list_length);i++) 
		if (flag[i] == 0) {
			if (n < i) list[n] = list[i];
			n++;
			}
  
	(*list_length) = n;

	free((void *) flag);
	}

int extract(const PARAMS *p,const char *filename,const BINARY *list,long nbin)
{
	/* extracts system with index = primary index and creates new data file */

	SSDATA *part;
	FILE *fp; /* for txt & bin file types only */
#ifdef SS_CORE
	SSHEAD head; /* for ss file type only */
	SSIO ssio_out;
#endif
	long n,i;
	int primary_done;
	
	(void) printf("Extracting system...\n");

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		if ((fp = fopen(filename,"w")) == NULL) {
			(void) fprintf(stderr,"Can't open %s\n",filename);
			return 1;
			}
		break;
#ifdef SS_CORE
	case FileTypeSS:
		if (ssioOpen(filename,&ssio_out,SSIO_WRITE)) {
			(void) fprintf(stderr,"Unable to open %s for writing\n",filename);
			return 1;
			}
		break;
#endif
	default:
		(void) fprintf(stderr,"file type %i is invalid\n",p->FileType);
		return 1;
		}

	for (n=i=0;i<nbin;i++)
		if (list[i].prim->index == p->ExtIdx)
			++n;
  
	if (n == 0) {
		(void) fprintf(stderr,"Index number %li not found\n",p->ExtIdx);
		return 1;
		}

#ifdef SS_CORE
	if (p->FileType == FileTypeSS) {
		head.time = 0.0;
		head.n_data = n + 1; /* one primary plus n satellites */
		head.pad = -1;
 
		if (ssioHead(&ssio_out,&head)) {
			(void) fprintf(stderr,"Unable to write header.\n");
			return 1;
			}
		}
#endif

	primary_done = 0;
	for (i=0;i<nbin;i++)
		if (list[i].prim->index == p->ExtIdx) {
			if (!primary_done) {
				part = &list[i].prim->data;
				primary_done = 1;
				--i; /* go back and do satellite */
				}
			else
				part = &list[i].sat->data;
			switch (p->FileType) {
			case FileTypeTxt:
				if (fprintf(fp,"%e %e %e %e %e %e %e %e\n",
							part->mass,part->radius, 
							part->pos[0],part->pos[1],part->pos[2],
							part->vel[0],part->vel[1],part->vel[2]) < 1) goto error;
				break;
			case FileTypeBin:
				if (fwrite(&part->mass,sizeof(double),1,fp) != 1) goto error;
				if (fwrite(&part->radius,sizeof(double),1,fp) != 1) goto error;
				if (fwrite(part->pos,sizeof(double),3,fp) != 3) goto error;
				if (fwrite(part->vel,sizeof(double),3,fp) != 3) goto error;
				break;
#ifdef SS_CORE
			case FileTypeSS:
				if (ssioData(&ssio_out,part) != 0) {
					(void) fprintf(stderr,"Error writing particle %li.\n",list[i].prim->index);
					return 1;
					}
				break;
#endif
			default:
				(void) fprintf(stderr,"file type %i is invalid\n",p->FileType);
				return 1;
				}
			}

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		fclose(fp);
		break;
#ifdef SS_CORE
	case FileTypeSS:
		(void) ssioClose(&ssio_out);
		break;
#endif
	default:
		assert(0);
		}

	return 0;

 error:
	(void) fprintf(stderr,"Error during write.\n");
	(void) fclose(fp);
	return 1;
	}


#ifdef SS_CORE
int find_real_part(const PARAMS *p,FILE *fp,SSIO *ssio_out,COMPDATA *part,long *n)
#else
int find_real_part(const PARAMS *p,FILE *fp,COMPDATA *part,long *n)
#endif
{
	if (part->prim == NULL)	{ /*reached bottom of tree*/
		++*n;
		switch (p->FileType) {
		case FileTypeTxt:
			if (fprintf(fp,"%e %e %e %e %e %e %e %e\n",
						part->data.mass,part->data.radius, 
						part->data.pos[0],part->data.pos[1],part->data.pos[2],
						part->data.vel[0],part->data.vel[1],part->data.vel[2]) < 1) goto error;
			break;
		case FileTypeBin:
			if (fwrite(&part->data.mass,sizeof(double),1,fp) != 1) goto error;
			if (fwrite(&part->data.radius,sizeof(double),1,fp) != 1) goto error;
			if (fwrite(part->data.pos,sizeof(double),3,fp) != 3) goto error;
			if (fwrite(part->data.vel,sizeof(double),3,fp) != 3) goto error;
			break;
#ifdef SS_CORE
		case FileTypeSS:
			if (ssioData(ssio_out,&part->data) != 0) {
				(void) fprintf(stderr,"Error writing particle %li.\n",part->index);
				return 1;
				}
			break;
#endif
		default:
			(void) fprintf(stderr,"file type %i is invalid\n",p->FileType);
			return 1;
			}

		return 0;
		}

	if (p->ApplyHierCut) {
		BIN_STATS bs;
		const SSDATA *pd,*sd;

		pd = &part->prim->data;
		sd = &part->sat->data;
	
		calc_bin_stats(pd,sd,BIN_CUT|BIN_INCL|BIN_PERIOD,&bs);
		if (ok_to_cut(p,part->prim,part->sat,&bs))
			return 0;
		}
#ifdef SS_CORE
	if (find_real_part(p,fp,ssio_out,part->prim,n) != 0) return 1;
	if (find_real_part(p,fp,ssio_out,part->sat,n) != 0) return 1;
#else
	if (find_real_part(p,fp,part->prim,n) != 0) return 1;
	if (find_real_part(p,fp,part->sat,n) != 0) return 1;
#endif
	return 0;

 error:
	(void) fprintf(stderr,"Error during write.\n");
	(void) fclose(fp);
	return 1;
	}

int hier_extract(const PARAMS *p,const char *filename,COMPDATA *part)
{
	/* extracts hierarchy system with top com particle index specified by user */
	long n;
	FILE *fp; /* for txt & bin file types only */
#ifdef SS_CORE
	SSHEAD head; /* for ss file type only */
	SSIO ssio_out;
#endif
	
	/* check that particle past is top com part */
	if (!HAS_COM_MBR(part) || IS_COM_MBR(part)) {
		(void) fprintf(stderr,"Particle %ld is not a top center of mass particle\n", part->index);
		return 1;
		}

	(void) printf("Extracting hierarchy system...\n");

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		if ((fp = fopen(filename,"w")) == NULL) {
			(void) fprintf(stderr,"Can't open %s\n",filename);
			return 1;
			}
		break;
#ifdef SS_CORE
	case FileTypeSS:
		if (ssioOpen(filename,&ssio_out,SSIO_WRITE)) {
			(void) fprintf(stderr,"Unable to open %s for writing\n",filename);
			return 1;
			}
		break;
#endif
	default:
		(void) fprintf(stderr,"file type %i is invalid\n",p->FileType);
		return 1;
		}

#ifdef SS_CORE /* write a dummy header */
	if (p->FileType == FileTypeSS) {
		head.time = 0.0;
		head.n_data = 1; /* one primary plus n satellites */
		head.pad = -1; 
 
		if (ssioHead(&ssio_out,&head)) {
			(void) fprintf(stderr,"Unable to write dummy header.\n");
			return 1;
			}
		}
#endif
	n = 0; /* number of real particles */

#ifndef SS_CORE
	find_real_part(p,fp,part,&n);
#else
	find_real_part(p,fp,&ssio_out,part,&n);
	if (p->FileType == FileTypeSS) {
		ssioRewind(&ssio_out);
		head.time = 0.0;
		head.n_data = n;
		head.pad = -1;

		if (ssioHead(&ssio_out,&head)) {
			(void) fprintf(stderr,"Unable to write true header.\n");
			return 1;
			}
		}
#endif
	(void) printf("Found %ld real particles in extracted system\n",n);

	switch (p->FileType) {
	case FileTypeTxt:
	case FileTypeBin:
		fclose(fp);
		break;
#ifdef SS_CORE
	case FileTypeSS:
		(void) ssioClose(&ssio_out);
		break;
#endif
	default:
		assert(0);
		}

	return 0;
	}

int write_tipsy(const char *filename,long npart,const BINARY *list,long nbin)
{
	/* creates a tipsy vector file of binary energy to use with tipsy software */

	/* NOTE: this function assumes list has been sorted properly! */

	BIN_STATS bs;
	FILE *fp;
	long i,j;

	if ((fp = fopen(filename,"w")) == NULL) {
		(void) fprintf(stderr,"Can't open %s\n",filename);
		return 1;
		}

	(void) printf("Creating tipsy file\n");
	if (fprintf(fp,"%li\n",npart) < 1) goto error;
	for (i=0;i<npart;i++) {
		bs.E = 0.0;
		for (j=0;j<nbin;j++)
			if (list[j].prim->index == i) { /* should be most bound if list sorted */
				calc_bin_stats(&list[j].prim->data,&list[j].sat->data,BIN_ENERGY,&bs);
				break;
				}
		if (fprintf(fp,"%e\n",bs.E) < 1) goto error;
		}

	(void) fclose(fp);

	return 0;

 error:
	(void) fprintf(stderr,"Error during write.\n");
	(void) fclose(fp);
	return 1;
	}

void usage(const char *progname)
{
	/* explains usage of companion and flag options */

	(void) printf("Usage: %s [-H [-z cutoff] [-g index[-a]]|-t] [-c|-m|-s] [-C|-M|-S] [ -f filetype ] [ -e cutoff ] [ -E cutoff ]\n",progname);
	(void) printf("       [ -h cutoff ] [ -q cutoff ] [ -o angle ] [ -x index ] file [ file ... ]\n");
	(void) printf("\n");
	(void) printf("Options: -H = search for hierarchies\n");
	(void) printf("         -z = close approach cutoff for center of mass particles (0 to eliminate orbit crossers,\n");
	(void) printf("              < 0 for semimajor axis\n");
	(void) printf("         -g = extract hierarchy system\n");
	(void) printf("         -a = apply cuts to extraction of hierachy system\n");
	(void) printf("         -t = create Tipsy vector file of binding energy\n");
	(void) printf("         -c | -m | -s = input in cgs, mks or system units (default \"%s\")\n",
#ifdef SS_CORE
				  "system");
#else
	              "cgs");
#endif
	(void) printf("         -C | -M | -S = output in cgs, mks or system units (default \"%s\")\n",
#ifdef SS_CORE
				  "system");
	(void) printf("         -f = file type: plain text (\"%s\"), binary (\"%s\"), or \"%s\" (default \"%s\")\n",
				  FileTypeStr[FileTypeTxt],FileTypeStr[FileTypeBin],FileTypeStr[FileTypeSS],FileTypeStr[DFLT_FILE_TYPE]);
#else
				  "cgs");
	(void) printf("         -f = file type: plain text (\"%s\") or binary (\"%s\") (default \"%s\")\n",
				  FileTypeStr[FileTypeTxt],FileTypeStr[FileTypeBin],FileTypeStr[DFLT_FILE_TYPE]);
#endif
	(void) printf("         -e = eccentricity cutoff\n");
	(void) printf("         -E = binding energy cutoff\n");
	(void) printf("         -h = Hill sphere cutoff, in Hill radii (prompts for semimajor axis and star mass)\n");
    (void) printf("         -q = close approach cutoff for normal particles (0 to eliminate colliders, < 0 for\n");
    (void) printf("              primary radii)\n");
	(void) printf("         -o = opening angle (default %g rad)\n",DFLT_OPEN_ANGLE);
	(void) printf("         -x = extracts system of given primary index\n");
	(void) printf("\n");
	(void) printf("NOTE: cutoff limits taken to be in output units where applicable.\n");

	exit(1);
	}

void set_defaults(PARAMS *params) 
{	
	/* parameters defaults */

	params->Hier = DFLT_HIER;
	params->TipsyFile = DFLT_TIPSY_FILE;
	params->InCgsUnits = DFLT_IN_CGS;
	params->InMksUnits = DFLT_IN_MKS;
	params->InSysUnits = DFLT_IN_SYS;
	params->OutCgsUnits = DFLT_OUT_CGS;
	params->OutMksUnits = DFLT_OUT_MKS;
	params->OutSysUnits = DFLT_OUT_SYS;
	params->FileType = DFLT_FILE_TYPE;
	params->EccCut = DFLT_ECC_CUT;
	params->EngCut = DFLT_ENG_CUT;
	params->HillCut = DFLT_HILL_CUT;
	params->PeriCut = DFLT_PERI_CUT;
	params->OrbitCut = DFLT_ORBIT_CUT;
	params->OpenAng = DFLT_OPEN_ANGLE;
	params->HierExtIdx = DFLT_HIER_EXTRACT_INDEX;
	params->ExtIdx = DFLT_EXTRACT_INDEX;

	/* derived parameters that need to be preset */

	params->DoPeriCut = DFLT_DO_PERI_CUT;
	params->DoOrbitCut = DFLT_DO_ORBIT_CUT;

	/* default units -- conversions between sys units and default units */

	params->InLengthConv = params->OutLengthConv = DFLT_LENGTH_CONV;
	params->InMassConv = params->OutMassConv = DFLT_MASS_CONV;
	params->InTimeConv = params->OutTimeConv = DFLT_TIME_CONV;
	}

void parse_in(int argc,char *argv[],PARAMS *p)
{
	extern int getopt(int,char *const *,const char *); /* in case unistd.h unavailable */

	extern int optind;
	extern char *optarg;

	char file_ext[FILE_TYPE_STR_MAX_LEN];
	int c;

	(void) strncpy(file_ext,FileTypeStr[DFLT_FILE_TYPE],FILE_TYPE_STR_MAX_LEN);

	while ((c = getopt(argc,argv,"HatcmsCMSz:g:f:e:E:h:q:o:x:")) != EOF)
		switch (c) {
		case 'H':
			p->Hier = TRUE;
			break;
		case 'a':
			p->ApplyHierCut = TRUE;
			break;
		case 't':
			p->TipsyFile = TRUE;
			break;
		case 'c':
			p->InCgsUnits = TRUE;
#ifdef SS_CORE
			if (p->InSysUnits)
				p->InSysUnits = FALSE;
#endif
			break;
		case 'm':
			p->InMksUnits = TRUE;
#ifndef SS_CORE
			if (p->InCgsUnits)
				p->InCgsUnits = FALSE;
#else
			if (p->InSysUnits)
				p->InSysUnits = FALSE;
#endif
			break;
		case 's':
			p->InSysUnits = TRUE;
#ifndef SS_CORE
			if (p->InCgsUnits)
				p->InCgsUnits = FALSE;
#endif
			break;
		case 'C':
			p->OutCgsUnits = TRUE;
#ifdef SS_CORE
			if (p->OutSysUnits)
				p->OutSysUnits = FALSE;
#endif
			break;
		case 'M':
			p->OutMksUnits = TRUE;
#ifndef SS_CORE
			if (p->OutCgsUnits)
				p->OutCgsUnits = FALSE;
#else
			if (p->OutSysUnits)
				p->OutSysUnits = FALSE;
#endif
			break;
		case 'S':
			p->OutSysUnits = TRUE;
#ifndef SS_CORE
			if (p->OutCgsUnits)
				p->OutCgsUnits = FALSE;
#endif
			break;
		case 'z':
			p->OrbitCut = atof(optarg);
			p->DoOrbitCut = TRUE;
			break;
		case 'g':
			p->HierExtIdx = atoi(optarg);
			break;
		case 'f':
			strcpy(file_ext,optarg);
			break;
		case 'e':
			p->EccCut = atof(optarg);
			break;
		case 'E':
			p->EngCut = atof(optarg);
			break;
		case 'h':
			p->HillCut = atof(optarg);
			break;
		case 'q':
			p->PeriCut = atof(optarg);
			p->DoPeriCut = TRUE;
			break;
		case 'o':
			p->OpenAng = atof(optarg);
			break;
		case 'x':
			p->ExtIdx = atoi(optarg);
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind >= argc) 
		usage(argv[0]);

	if (strcmp(file_ext,FileTypeStr[FileTypeTxt]) == 0)
		p->FileType = FileTypeTxt;
	else if (strcmp(file_ext,FileTypeStr[FileTypeBin]) == 0)
		p->FileType = FileTypeBin;
#ifdef SS_CORE
	else if (strcmp(file_ext,FileTypeStr[FileTypeSS]) == 0) 
		p->FileType = FileTypeSS;
#endif
	else 
		usage(argv[0]);

	/* sanity checks */

	if (p->Hier && p->TipsyFile)
		usage(argv[0]);

	if ((p->InCgsUnits == TRUE && p->InMksUnits == TRUE) ||
		(p->InCgsUnits == TRUE && p->InSysUnits == TRUE) ||
		(p->InMksUnits == TRUE && p->InSysUnits == TRUE) ||
		(p->OutCgsUnits == TRUE && p->OutMksUnits == TRUE) ||
		(p->OutCgsUnits == TRUE && p->OutSysUnits == TRUE) ||
		(p->OutMksUnits == TRUE && p->OutSysUnits == TRUE))
		usage(argv[0]);

	if (p->EccCut < 0.0 || p->EccCut > 1.0) {
		(void) fprintf(stderr,"Eccentricity cut must be between 0 and 1.\n");
		exit(1);
		}

	if (p->HillCut < 0.0) {
		(void) fprintf(stderr,"Hill sphere cut must be positive.\n");
		exit(1);
		}

	if (p->EngCut > 0.0) {
		(void) fprintf(stderr,"Energy cut must be negative.\n");
		exit(1);
		}

	if (p->DoOrbitCut == TRUE && p->Hier == FALSE) {
		(void) fprintf(stderr,"Close approach cut for center of mass particles must be used with hierarchy\n");
		exit(1);
		}

	if (p->OpenAng < 0.0) {
		(void) fprintf(stderr,"Opening angle must be positive or zero.\n");
		exit(1);
		}

	if (p->HierExtIdx >= 0 && p->Hier == FALSE) {
		(void) fprintf(stderr,"Hierarchy extraction must be used with hierarchy option\n");
		exit(1);
		}

	if (p->ApplyHierCut == TRUE && p->HierExtIdx >= 0) {
		(void) fprintf(stderr,"ApplyHierCut option must be used with hierarchy extraction\n");
		exit(1);
		}

	if (p->HierExtIdx < -1) { /* -1 is default, i.e., no extraction */
		(void) fprintf(stderr,"Hierarchy extraction index must be non-negative.\n");
		exit(1);
		}

	if (p->ExtIdx < -1) { /* -1 is default, i.e., no extraction */
		(void) fprintf(stderr,"Extraction index must be non-negative.\n");
		exit(1);
		}

	/* Default I/O in cgs units (sys if SS_CORE defined); data stored internally in system units */

	if (p->InCgsUnits == TRUE) {
		p->InLengthConv = AU*1.0e2;
		p->InMassConv = M_SUN*1.0e3;
		p->InTimeConv = SID_YR/TWO_PI;
		(void) printf("Input in cgs units.\n");
	} else if (p->InMksUnits == TRUE) {
		p->InLengthConv = AU;
		p->InMassConv = M_SUN;
		p->InTimeConv = SID_YR/TWO_PI;
		(void) printf("Input in mks units.\n");
	} else if (p->InSysUnits == TRUE){
		p->InLengthConv = 1.0;
		p->InMassConv = 1.0;
		p->InTimeConv = 1.0;
		(void) printf("Input in system units.\n");
		}

	p->InEnergyConv = p->InMassConv*SQ(p->InLengthConv/p->InTimeConv); /* NOT USED */

	if (p->OutCgsUnits == TRUE) {
		p->OutLengthConv = AU*1.0e2;
		p->OutMassConv = M_SUN*1.0e3;
		p->OutTimeConv = SID_YR/TWO_PI;
		(void) printf("Output in cgs units.\n");
	} else if (p->OutMksUnits == TRUE) {
		p->OutLengthConv = AU;
		p->OutMassConv = M_SUN;
		p->OutTimeConv = SID_YR/TWO_PI;
		(void) printf("Output in mks units.\n");
	} else if (p->OutSysUnits == TRUE) {
		p->OutLengthConv = 1.0;
		p->OutMassConv = 1.0;
		p->OutTimeConv = 1.0;
		(void) printf("Output in sys units.\n");
		}

	p->OutEnergyConv = p->OutMassConv*SQ(p->OutLengthConv/p->OutTimeConv);

	p->OpenAng = SQ(p->OpenAng); /* store square of opening angle */

	/* convert cuts to system units as needed */

	p->EngCut /= p->OutEnergyConv;

	if (p->HillCut > 0.0){
		(void) printf("What is the average semimajor axis (in AU)?\n");
		(void) scanf("%lf",&(p->avga));
		if (p->avga <= 0.0) {
			(void) fprintf(stderr,"Semimajor axis must be positive.\n");
			exit(1);
			}
		(void) printf("What is the mass of the star (in M_Sun)?\n");
		(void) scanf("%lf",&(p->starmass));
		if (p->starmass <= 0.0) {
			(void) fprintf(stderr,"Star mass must be positive.\n");
			exit(1);
			}
		}

	if (p->PeriCut > 0.0)
		p->PeriCut /= p->OutLengthConv;

	if (p->OrbitCut > 0.0)
		p->OrbitCut /= p->OutLengthConv;
	}

int write_output(const PARAMS *p,const char *filename_in,const BINARY *list,long nbin,double m_tot)
{
	BIN_STATS bs;
	const SSDATA *pd,*sd;
	FILE *fp_pr,*fp_ana;
	char pr_outfile[MAXPATHLEN],ana_outfile[MAXPATHLEN]; 
	long nsys,i,last_index;

	if (myNewExt(filename_in,FileTypeStr[p->FileType],pr_outfile,PR_EXT)) {
		(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
					   filename_in);
		return 1;	
		}
	
	if (myNewExt(filename_in,FileTypeStr[p->FileType],ana_outfile,ANA_EXT)) {
		(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
					   filename_in);
		return 1;
		}

	if ((fp_pr = fopen(pr_outfile,"w")) == NULL) {
		(void) fprintf(stderr,"Can't open %s\n",pr_outfile);
		return 1;
		}

	if ((fp_ana = fopen(ana_outfile,"w")) == NULL) {
		(void) fprintf(stderr,"Can't open %s\n",ana_outfile);
		return 1;
		}

	(void) fprintf(fp_pr," M_p/M_t     p_ind    p_rad  M_s/M_p     s_ind    s_rad  bind_eng        a    e    i      per\n");
	(void) fprintf(fp_pr,"-------- --------- -------- -------- --------- -------- --------- -------- ---- ---- --------\n");
	last_index = -1;
	for (nsys=i=0;i<nbin;i++) {
		pd = &list[i].prim->data;
		sd = &list[i].sat->data;
		calc_bin_stats(pd,sd,BIN_ENERGY|BIN_SEMI|BIN_ECC|BIN_INCL|BIN_PERIOD,&bs);
		/* "pretty" output */
		if (list[i].prim->index != last_index) {
			if (fprintf(fp_pr,"%8.2e %9li %8.2e ",pd->mass/m_tot,list[i].prim->index,
						pd->radius*p->OutLengthConv) < 1) goto error;
			last_index = list[i].prim->index;
			++nsys;
			}
		else
			(void) fprintf(fp_pr,"%28s",""); /* pad 28 spaces */
		if (fprintf(fp_pr,"%8.2e %9li %8.2e %9.2e %8.2e %4.2f %4.2f %8.2e\n",
					sd->mass/pd->mass,list[i].sat->index,sd->radius*p->OutLengthConv,
					bs.E*p->OutEnergyConv,bs.a*p->OutLengthConv,bs.e,bs.i,bs.P*p->OutTimeConv) < 1) goto error;
		/* machine output */
		if (fprintf(fp_ana,"%e %li %e %e %li %e %e %e %e %e %e\n",pd->mass/m_tot,
					list[i].prim->index,pd->radius*p->OutLengthConv,sd->mass/pd->mass,
					list[i].sat->index,sd->radius*p->OutLengthConv,bs.E*p->OutEnergyConv,
					bs.a*p->OutLengthConv,bs.e,bs.i,bs.P*p->OutTimeConv) < 1) goto error;
		}

	(void) fprintf(fp_pr,"Summary: %li system%s, %li binar%s, total mass considered = %e\n",
				   nsys,nsys==1?"":"s",nbin,nbin==1?"y":"ies",m_tot*p->OutMassConv);
	
	(void) fclose(fp_ana);
	(void) fclose(fp_pr);

	return 0;

 error:
	(void) fprintf(stderr,"Error during write.\n");
	(void) fclose(fp_ana);
	(void) fclose(fp_pr);
	return 1;
	}

int walk_system(const PARAMS *p,FILE *fp,const COMPDATA *part,double m_tot,BOOLEAN write_header,int *n_layer)
{
	BIN_STATS bs;
	const SSDATA *pd,*sd;
	const COMPDATA *primary,*satellite;

	if (part->prim == NULL)
		return 0;

	primary = part->prim;
	satellite = part->sat;
	pd = &primary->data;
	sd = &satellite->data;
	calc_bin_stats(pd,sd,BIN_CUT|BIN_INCL|BIN_PERIOD,&bs);
	if (ok_to_cut(p,primary,satellite,&bs))
		return 0;

	if (write_header) {
		*n_layer = 0;
		(void) fprintf(fp,"    c_ind  M_p/M_t     p_ind    p_rad  M_s/M_p     s_ind    s_rad  bind_eng        a    e    i      per\n");
		(void) fprintf(fp,"--------- -------- --------- -------- -------- --------- -------- --------- -------- ---- ---- --------\n");
		}

	if (fprintf(fp,"%9li %8.2e %9li %8.2e %8.2e %9li %8.2e %9.2e %8.2e %4.2f %4.2f %8.2e\n",
				part->index,pd->mass/m_tot,part->prim->index,pd->radius*p->OutLengthConv,
				sd->mass/pd->mass,part->sat->index,sd->radius*p->OutLengthConv,bs.E*p->OutEnergyConv,
				bs.a*p->OutLengthConv,bs.e,bs.i,bs.P*p->OutTimeConv) < 1) return 1;
	++(*n_layer);
	if (walk_system(p,fp,part->prim,m_tot,FALSE,n_layer) != 0) return 1;
	if (walk_system(p,fp,part->sat,m_tot,FALSE,n_layer) != 0) return 1;

	return 0;
	}

void create_summary(const PARAMS *p,const COMPDATA *part,HIER_OUTPUT *sum,double *max_a)
{
	BIN_STATS bs;
	const SSDATA *pd,*sd;

	if (part->prim == NULL)
		return;

	pd = &part->prim->data;
	sd = &part->sat->data;
	calc_bin_stats(pd,sd,BIN_CUT,&bs);
	if (ok_to_cut(p,part->prim,part->sat,&bs))
		return;

	if (bs.a > *max_a) {
		sum->max_a = bs.a;
		*max_a = bs.a;
		}

	sum->bind_E += bs.E;
	
	create_summary(p,part->prim,sum,max_a);
	create_summary(p,part->sat,sum,max_a);
	}

int sort_mass(const void *a, const void *b)
{

	const HIER_OUTPUT *s1,*s2;

	s1 = (const HIER_OUTPUT *) a; /* pointers to systems that made the cut */
	s2 = (const HIER_OUTPUT *) b;

	if (s1->sys_mass < s2->sys_mass) return 1; /* largest system first */
	if (s1->sys_mass > s2->sys_mass) return -1;

	if (s1->com->index < s2->com->index) return -1; /* if total system masses are equal */
	if (s1->com->index > s2->com->index) return 1;  /* smallest com indices first */

	assert(0); /* can't get here, top com particles can not have same indices */

	return 0;
	}

int write_hier_output(const PARAMS *p,const char *filename,const COMPDATA *part,long n_part,double m_tot,long n_orig)
{
	BIN_STATS bs;
	const SSDATA *pd,*sd;
	HIER_OUTPUT *sum;
	FILE *fp;
	char outfile[MAXPATHLEN];
	int n_layer;
	long i,n_sys,n_true_hier,sum_buf_size;

	if (myNewExt(filename,FileTypeStr[p->FileType],outfile,HIER_EXT)) {
		(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
					   filename);
		return 1;	
		}
	
	if ((fp = fopen(outfile,"w")) == NULL) {
		(void) fprintf(stderr,"Can't open %s\n",outfile);
		return 1;
		}

	n_sys = 0;
	sum_buf_size = BUF_SIZE_INIT;
	sum = (HIER_OUTPUT *) malloc(sum_buf_size*sizeof(HIER_OUTPUT));
	assert(sum != NULL);

	for (i=0;i<n_part;i++) {
		pd = &(part[i].prim->data);
		sd = &(part[i].sat->data);
		if (part[i].com == NULL && part[i].prim != NULL) { /* top of system tree */
			calc_bin_stats(pd,sd,BIN_CUT,&bs);
			if (ok_to_cut(p,part[i].prim,part[i].sat,&bs))
				continue; /* if first layer of system does not make cut ignore system altogether */

			sum[n_sys].sys_mass = pd->mass + sd->mass;
			sum[n_sys].max_a = bs.a;
			sum[n_sys].bind_E = bs.E;
			sum[n_sys].com = &part[i];
			create_summary(p,part[i].prim,&sum[n_sys],&(bs.a));
			create_summary(p,part[i].sat,&sum[n_sys],&(bs.a));

			if (++n_sys == sum_buf_size) {
				sum_buf_size *= BUF_SIZE_MULT;
				sum = (HIER_OUTPUT *) realloc((void *) sum,sum_buf_size*sizeof(HIER_OUTPUT));
				assert(sum != NULL);
				} 
			}
		}

	qsort((void *) sum,n_sys,sizeof(HIER_OUTPUT),sort_mass); /* sort systems with most massive first */
	n_true_hier = 0;
	for (i=0;i<n_sys;i++) {
		if (walk_system(p,fp,sum[i].com,m_tot,TRUE,&n_layer) != 0) {
			(void) fprintf(stderr,"Error during write.\n");
			(void) fclose(fp);
			return 1;
			}
		if (fprintf(fp,"System summary: mass = %8.2e, max semimajor axis = %8.2e, total binding energy = %9.2e\n",
					sum[i].sys_mass*p->OutMassConv,sum[i].max_a*p->OutLengthConv,sum[i].bind_E*p->OutEnergyConv) < 1) return 1;
		fprintf(fp,"\n");
		if (n_layer > 1)
			++n_true_hier;
		}	
	if (fprintf(fp,"%ld system%s found: %ld 2-particle system%s and %ld multi-particle system%s\n", 
				n_sys, n_sys==1?"":"s",n_sys-n_true_hier,n_sys-n_true_hier==1?"":"s", 
				n_true_hier, n_true_hier==1?"":"s") < 1) return 1;
	if (fprintf(fp,"Total number of original particle%s: %ld\n", n_orig==1?"":"s",n_orig) < 1) return 1;
	if (fprintf(fp,"Total mass in original particle%s: %8.2e\n", 
				n_orig==1?"":"s",m_tot*p->OutMassConv) < 1) return 1;

	(void) fclose(fp);
	free((void *) sum); 
	return 0;
	}

int sort_per(const void *a,const void *b)
{
	const BINARY *b1,*b2;
	const SSDATA *p1,*p2;

	b1 = (const BINARY *) a; /* pointers to binaries */
	b2 = (const BINARY *) b;

	p1 = &(b1->prim->data); /* pointers to primaries */
	p2 = &(b2->prim->data);

	if (b1->period < b2->period) return 1; /* largest periods first */
	if (b1->period > b2->period) return -1;

	if (p1->mass < p2->mass) return 1; /* largest primary masses first */
	if (p1->mass > p2->mass) return -1;

	if (b1->prim->index < b2->prim->index) return -1; /* smallest indices first */
	if (b1->prim->index > b2->prim->index) return 1;

	/* within each system, sort by binding energy, but check for duplicate first */

	{
	    const SSDATA *s1 = &(b1->sat->data),*s2 = &(b2->sat->data); /* satellite pointers */

		if (b1->sat->index == b2->sat->index) return 0; /* oops! */

		{
		    BIN_STATS bs1,bs2;

			calc_bin_stats(p1,s1,BIN_ENERGY,&bs1);
			calc_bin_stats(p2,s2,BIN_ENERGY,&bs2);

			if (bs1.E < bs2.E) return -1; /* more bound first */
			if (bs1.E > bs2.E) return 1;

			if (b1->sat->index < b2->sat->index) return -1; /* smallest indices first */
			if (b1->sat->index > b2->sat->index) return 1;
			}
		}

	assert(0); /* can't get here */

	return 0;
	}

void calc_per(BINARY *list,long n)
{
	BIN_STATS bs;
	int i;

	for (i=0;i<n;i++) {
		calc_bin_stats(&list[i].prim->data,&list[i].sat->data,BIN_PERIOD,&bs);
		list[i].period = bs.P;
		}
	}

void find_systems(const PARAMS *p,COMPDATA **part,long *npart,
				  VECTOR root_center,double root_size,NODE **root,
				  BINARY **bin,long *bin_buf_size,long *nbin)
{
	BINARY *binary,*add,*new,*ptr;
	COMPDATA *com;
	long *del;
	long part_buf_size,del_buf_size,add_buf_size,new_buf_size,buf;
	long nmax,nloops,ntree,ncmp,ncom,ndel,nadd,nnew,idel,ibin,iadd,n,i;

	(void) printf("Starting hierarchical search for systems...\n");

	/*
	 ** Do the following just once: compute periods for each existing
	 ** binary and sort the binaries in decreasing order of period.
	 ** Subsequently binaries will be deleted and possibly added to
	 ** the list while preserving the sort order.
	 */

	calc_per(*bin,*nbin);
	qsort((void *) *bin,*nbin,sizeof(BINARY),sort_per);

	/* initialize */
	ntree = ncmp = *npart; /* used to monitor tree state */
	ncom = 0; /* ditto */
	part_buf_size = (*npart)*EXTRA_STORE; /* see comment at top of file */

	/* allocate space for maintenance lists */
	del_buf_size = BUF_SIZE_INIT;
	del = (long *) malloc(del_buf_size*sizeof(long));
	assert(del != NULL);
	add_buf_size = BUF_SIZE_INIT;
	add = (BINARY *) malloc(add_buf_size*sizeof(BINARY));
	assert(add != NULL);
	new_buf_size = *bin_buf_size;
	new = (BINARY *) malloc(new_buf_size*sizeof(BINARY));
	assert(new != NULL);

	/*
	 ** Now loop, finding "tightest" (shortest period) binary each
	 ** time, and updating the binary list as required, until no
	 ** binaries remain.  Periodically rebuild the tree to improve
	 ** efficiency.
	 */

	nmax = (*npart > 0xffff ? INT_MAX : (*npart - 1)*(*npart)/2); /* worst-case scenario */
	nloops = 0;

	while (*nbin > 0) {
		--(*nbin); /* truncate list */
		binary = &((*bin)[*nbin]); /* last binary in list had shortest period */
		make_com_part(binary,npart,&part_buf_size,part);
		--ncmp; /* 2 particles replaced by 1 com particle */
		assert(ncmp > 0); /* can't run out of particles! */
		com = &((*part)[*npart - 1]); /* last particle in list is new com particle */

		/* 
		 ** Check to see if tree should be rebuilt:
		 ** 1) if the ratio of the number of original particles left 
		 **    to the number of particles in the tree since the last
		 **    tree build is less than TREE_REBUILD_FRAC.
		 ** 2) if the mass ratio of the two components of the new com
		 **    particle exceed REBUILD_MASS_RATIO.
		 */

		if ((double) ncmp/ntree < TREE_REBUILD_FRAC || 
			binary->prim->data.mass/binary->sat->data.mass > REBUILD_MASS_RATIO) {
			(void) printf("Rebuilding tree... (N = %li)\n",ncmp);
			kill_node(*root);
			make_node(root_center,root_size,root); /* note: recomputing center & size would improve efficiency */
			for (i=0;i<*npart;i++)
				if (!IS_COM_MBR(&((*part)[i]))) /* no child particles */
					add_to_tree(*root,&((*part)[i]));
			get_com_vel(*root);
			assert((*root)->n_part == ncmp); /* particle conservation check */
			ntree = ncmp; /* number of particles in tree for this rebuild */
			ncom = 0;
			}
		else {
			add_to_tree(*root,com); /* add com particle to tree */
			++ncom; /* number of com particles added since last tree rebuild */
			get_com_vel(*root);
			assert((*root)->n_part == ntree + ncom);
			}

		/* create list of binaries to remove */

		for (ndel=ibin=0;ibin<*nbin;ibin++) {
			binary = &((*bin)[ibin]);
			/*
			 ** Record any binary whose primary or satellite was either
			 ** of the children of the new com particle, and mark the
			 ** members of that binary to be resent to find_companion().
			 */
			if (IS_COM_MBR(binary->prim) || IS_COM_MBR(binary->sat)) {
				del[ndel] = ibin;
				if (++ndel == del_buf_size) {
					del_buf_size *= BUF_SIZE_MULT;
					del = (long *) realloc((void *) del,del_buf_size*sizeof(long));
					assert(del != NULL);
					}
				}
			}

		/*
		 ** Now call find_companion() for the new com particle, storing
		 ** results in new list.
		 */

		nadd = 0; /* (reuse existing storage) */
		find_companion(p,*root,com,&add,&add_buf_size,&nadd);

		/* compute periods of binaries to add, then sort */

		calc_per(add,nadd);
		qsort((void *) add,nadd,sizeof(BINARY),sort_per);

		/*
		 ** Update binary list by deleting old binaries and adding
		 ** new binaries all in a single pass, being careful to
		 ** reject any duplicated entries in the new list.
		 */

		idel = ibin = iadd = nnew = 0;
		while (ibin < *nbin || iadd < nadd) {
			/* omit current binary from new list? */
			if (ibin < *nbin) {
				binary = &((*bin)[ibin]);
				if (idel < ndel && del[idel] == ibin) {
					++idel; /* increment and don't copy */
					++ibin;
					continue;
					}
				}
			/*
			 ** Following "while" cascade does not consider binding energy,
			 ** unlike sort_per() -- we'd rather avoid two calc_bin_stats()
			 ** calls here.  Normally only an artificial test should lead
			 ** to this being a problem.
			 */
			while (iadd < nadd &&
				   (ibin == *nbin ||
					(add[iadd].period > binary->period ||
					 (add[iadd].period == binary->period &&
					  (add[iadd].prim->data.mass > binary->prim->data.mass ||
					   (add[iadd].prim->data.mass == binary->prim->data.mass &&
						(add[iadd].prim->index < binary->prim->index ||
						 (add[iadd].prim->index == binary->prim->index &&
						  add[iadd].sat->index < binary->sat->index)))))))) {
				new[nnew] = add[iadd];
				/* check for duplicate add -- they will always be paired together */
				if (++iadd < nadd && IS_SAME_BINARY(&add[iadd],&new[nnew]))
					++iadd; /* skip it */
				/* increment and check for possible buffer overflow */
				if (++nnew == new_buf_size) {
					new_buf_size *= BUF_SIZE_MULT;
					new = (BINARY *) realloc((void *) new,new_buf_size*sizeof(BINARY));
					assert(new != NULL);
					}
				}
			/* copy current binary to new list? */
			if (ibin < *nbin) {
				new[nnew] = *binary;
				++ibin;
				/* increment and check for possible buffer overflow */
				if (++nnew == new_buf_size) {
					new_buf_size *= BUF_SIZE_MULT;
					new = (BINARY *) realloc((void *) new,new_buf_size*sizeof(BINARY));
					assert(new != NULL);
					}
				}
			}

		/* swap new binary list with original list via pointers to save time */

		ptr = *bin;
		buf = *bin_buf_size;
		n = *nbin;
		*bin = new;
		*bin_buf_size = new_buf_size;
		*nbin = nnew;
		new = ptr;
		new_buf_size = buf;
		nnew = n;

		if (*nbin % 100 == 0) (void) printf("%li binaries remaining (%li/%li active/total particles)\n",*nbin,ncmp,*npart);

		if (++nloops == nmax) {
			(void) fprintf(stderr,"Maximum number of loops exceeded\n");
			break;
			}
		}

	free((void *) new);
	free((void *) add);
	free((void *) del);
	}

int main(int argc,char *argv[])
{
	extern int optind;

	COMPDATA *part;
	NODE *root;
	BINARY *list;
	PARAMS params;
	VECTOR root_center;
	double m_tot,root_size;
	long list_size,list_posn;
	long n_part,n_orig,i;

	/*Defaults*/
	set_defaults(&params);

	/* Parse command line arguments */
	parse_in(argc,argv,&params);

	for(;optind<argc;optind++) {

		(void) printf("Reading data...\n");
		if (read_data(&params,argv[optind],&part,&n_part,&m_tot,root_center,&root_size) != 0)
			return 1;

		(void) printf("Building tree...\n");
		make_node(root_center,root_size,&root);
		for (i=0;i<n_part;i++) 
			add_to_tree(root,&part[i]);

		(void) printf("Computing center of mass...\n");
		get_com_vel(root);
		assert(root->n_part == n_part); /* particle conservation check */

		(void) printf("Beginning satellite search...\n");
		list = NULL;
		list_size = list_posn = 0;
		for (i=0;i<n_part;i++)
			find_companion(&params,root,&part[i],&list,&list_size,&list_posn);
		(void) printf("%li binar%s found.\n",list_posn,list_posn==1?"y":"ies");

		if (list_posn == 0) goto done; /* no point in continuing */

		if (params.ExtIdx >= 0 && params.ExtIdx < n_part) {
			char ext_outfile[MAXPATHLEN];

			if (myNewExt(argv[optind],FileTypeStr[params.FileType],ext_outfile,EXT_EXT)) {
				(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
							   argv[optind]);
				return 1;
				}

			(void) extract(&params,ext_outfile,list,list_posn);
			}

		if (params.Hier == TRUE) {
			n_orig = n_part;
			find_systems(&params,&part,&n_part,root_center,root_size,&root,&list,&list_size,&list_posn);
			/*add hierarchy extraction here*/
			if (params.HierExtIdx >= 0) {
				if (params.HierExtIdx >= n_orig && params.HierExtIdx < n_part) {
					char hier_ext_outfile[MAXPATHLEN];
				
					if (myNewExt(argv[optind],FileTypeStr[params.FileType],hier_ext_outfile,HIER_EXT_EXT)) {
						(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
									   argv[optind]);
						return 1;
						}

					if (hier_extract(&params,hier_ext_outfile,&part[params.HierExtIdx])) {
						(void) fprintf(stderr,"Unable to extract hierarchy system of particle %ld\n",
									   params.HierExtIdx);
						return 1;
						}
					}
				else {
					(void) fprintf(stderr,"Hierarchy extraction index must be a center of mass particle\n");
					(void) fprintf(stderr,"Index >= %ld\n",n_orig);
					return 1;
					}
				}
			(void) write_hier_output(&params,argv[optind],part,n_part,m_tot,n_orig); /* cuts applied */
			}
		else { /* do a normal cull */

			(void) printf("Sorting...\n");
			qsort((void *) list,list_posn,sizeof(BINARY),sort_bin);

			cut_list(&params,list,&list_posn); /*applies any cuts stored in params*/
			(void) printf("%li binar%s survived the cut.\n",list_posn,list_posn==1?"y":"ies");

			(void) write_output(&params,argv[optind],list,list_posn,m_tot); /* output normally */

			if (params.TipsyFile == TRUE) {
				char tip_outfile[MAXPATHLEN];

				if (myNewExt(argv[optind],FileTypeStr[params.FileType],tip_outfile,VEC_EXT)) {
					(void) fprintf(stderr,"Unable to generate output filename for %s.\n",
								   argv[optind]);
					return 1;	
					}

				(void) write_tipsy(tip_outfile,n_part,list,list_posn);
				}

			}

	done:

		/* all done */

		free((void *) list);
		kill_node(root);
		free((void *) part);
		}

	return 0;
	}
