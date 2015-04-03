// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FoldConstraints.hh
/// @brief Abinitio-Folding under (distance-)constraints
/// @details
///	similar to classic Foldconstraints Protocol
///
///
/// @author Oliver Lange


#ifndef INCLUDED_devel_simpleoptions_filenm_HH
#define INCLUDED_devel_simpleoptions_filenm_HH

#include <core/types.hh>
#include <vector>
/* This macro calculates the size of a array */
#define asize(a) (sizeof(a)/sizeof((a)[0]))


/* this enum should correspond to the array deffile in gmxlib/filenm.c */
enum {
  efMDP, efGCT,
  efTRX, efTRN, efTRR, efTRJ, efXTC, efG87,
  efENX, efEDR, efENE,
  efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT, efESP, efPQR,
  efLOG, efXVG, efOUT,
  efNDX,
  efTOP, efITP,
  efTPX, efTPS, efTPR, efTPA, efTPB,
  efTEX, efRTP, efATP, efHDB,
  efDAT, efDLG,
  efMAP, efEPS, efMAT, efM2P,
  efMTX,
  efEDI, efEDO,
  efPPA, efPDO,
  efHAT,
  efXPM,
  efNR
};

typedef struct {
  int  ftp;		/* File type (see enum above)		*/
  char *opt;		/* Command line option			*/
  char *fn;		/* File name (as set in source code)	*/
  unsigned long flag;	/* Flag for all kinds of info (see defs)*/
  int  nfiles;		/* number of files			*/
  char **fns;		/* File names				*/
} t_filenm;

#define ffSET	1<<0
#define ffREAD	1<<1
#define ffWRITE	1<<2
#define ffOPT	1<<3
#define ffLIB	1<<4
#define ffMULT	1<<5
#define ffRW	(ffREAD	| ffWRITE)
#define ffOPTRD	(ffREAD	| ffOPT)
#define ffOPTWR	(ffWRITE| ffOPT)
#define ffOPTRW	(ffRW	| ffOPT)
#define ffLIBRD	(ffREAD	| ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
#define ffRDMULT   (ffREAD  | ffMULT)
#define ffOPTRDMULT   (ffRDMULT | ffOPT)
#define ffWRMULT   (ffWRITE  | ffMULT)
#define ffOPTWRMULT   (ffWRMULT | ffOPT)


/* This structure is used for parsing arguments off the comand line */
enum {
  etINT, etREAL, etTIME, etSTR,    etBOOL, etRVEC,   etENUM, etNR
};
/* names to print in help info */
static char *argtp[etNR] = {
  "int", "real", "time", "string", "bool", "vector", "enum"
};

typedef struct {
  char *option;
  bool bSet;
  int  type;
  union {
    void *v;   /* This is a nasty workaround, to be able to use initialized */
    int  *i;   /* arrays */
    core::Real *r;
    char **c;  /* Must be pointer to string (when type == etSTR)         */
               /* or null terminated list of enums (when type == etENUM) */
    bool *b;
  } u;
  char *desc;
} t_pargs;

typedef  std::vector< t_pargs > CmdLinePargs;
void fill_cmd_vector ( t_pargs* source, CmdLinePargs target ) {
  int NPA = asize(source);
  for (int i=0; i<NPA; i++ ) {
    target.push_back( source[i] );
  }
}


void write_man(FILE *out,char *mantp,
	       char *program,
	       int nldesc,char **desc,
	       int nfile,t_filenm *fnm,
	       int npargs,t_pargs *pa,
  int nbug,char **bugs,bool bHidden);


#define PCA_CAN_VIEW       (1<<5)
/* add option -w to view output files (must be implemented in program) */
#define PCA_CAN_BEGIN      (1<<6)
#define PCA_CAN_END        (1<<7)
#define PCA_CAN_DT         (1<<14)
#define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
/* adds options -b and -e for begin and end time for reading trajectories */
#define PCA_TIME_UNIT      (1<<15)
/* set time unit for output */
#define PCA_KEEP_ARGS      (1<<8)
/* keep parsed args in argv (doesn't make sense without NOEXIT_ON_ARGS) */
#define PCA_SILENT         (1<<9)
/* don't print options by default */
#define PCA_CAN_SET_DEFFNM (1<<10)
/* does something for non-master mdrun nodes */
#define PCA_NOEXIT_ON_ARGS (1<<11)
/* no fatal_error when invalid options are encountered */
#define PCA_QUIET          (1<<12)
/* does something for non-master mdrun nodes */
#define PCA_BE_NICE        (1<<13)
/* Default to low priority, unless configured with --disable-nice */


void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],book bKeep );


#endif /* #ifdef ... _HH */
