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

#include "devel/simple_options/filenm.h"
#include <stdio.h>
#include <string>
#include "utility/exit.hh"
#include "core/types.hh"

extern char *mk_desc(t_pargs *pa, char *time_unit_str);
extern void get_pargs(int *argc, char *argv[], CmdLinePargs, bool bKeepArgs);
static char  *cmdline      = NULL;

void parse_common_args(int *argc,char *argv[],unsigned long Flags,
					 int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
					 int ndesc,char **desc,int nbugs,char **bugs)
{
	static bool bHelp=false,bHidden=false,bQuiet=false;
	static char *manstr[]      = { NULL, "no", "html", "tex", "nroff", "ascii", "completion", "py", "xml", NULL };
	static char *not_nicestr[] = { NULL, "0", "4", "10", "19", NULL };
	static char *nicestr[]     = { NULL, "19", "10", "4", "0", NULL };
	static char *not_npristr[] = { NULL, "0", "128", "100", "200", "250", NULL };
	static char *npristr[]     = { NULL, "128", "250", "200", "100", "0", NULL };
	static int  nicelevel=0;//,mantp=0,npri=0;
	static bool /*bGUI=false,*/bDebug=false;
	//static char *deffnm=NULL;


	CmdLinePargs all_pa;

	t_pargs nice_paX  = { "-nice", false, etENUM,  {not_nicestr},
					 "Set the nicelevel" };
	t_pargs nice_pa   = { "-nice", false, etINT,   {&nicelevel},
					 "Set the nicelevel" };

#define EXTRA_PA 16

	t_pargs pca_pa[] = {
		{ "-h",    false, etBOOL, {&bHelp},
			"Print help info and quit" },
		{ "-hidden", false, etBOOL, {&bHidden},
			"HIDDENPrint hidden options" },
		{ "-quiet",false, etBOOL, {&bQuiet},
			"HIDDENDo not print help info" },
		{ "-man",  false, etENUM,  {manstr},
			"HIDDENWrite manual and quit" },
		{ "-debug",false, etBOOL, {&bDebug},
			"HIDDENWrite file with debug information" },
	};

#define NPCA_PA asize(pca_pa)

	//FILE *fp;
	bool bPrint,bExit;//,bXvgr;
	int  i,j,/*k,npall,*/max_pa,cmdlength;
	//char *ptr,*newdesc;
	//char *envstr;

	std::string program;
#define FF(arg) ((Flags & arg)==arg)

	cmdlength = strlen(argv[0]);
	/* Check for double arguments */
	for (i=1; (i<*argc); i++) {
		cmdlength += strlen(argv[i]);
		if (argv[i] && (strlen(argv[i]) > 1) && (!isdigit(argv[i][1]))) {
			for (j=i+1; (j<*argc); j++) {
	if ( (argv[i][0]=='-') && (argv[j][0]=='-') &&
			 (strcmp(argv[i],argv[j])==0) ) {
		if (FF(PCA_NOEXIT_ON_ARGS)) {
			fprintf(stderr,"Double command line argument %s\n",argv[i]);
		} else {
			utility_exit_with_message( "Double command line argument "+std::string(argv[i])+ "\n" );
		}
	}
			}
		}
	}

	program = std::string(argv[0]);

	/* Fill the cmdline string */
	cmdline = new char[cmdlength+*argc+1];
	for (i=0; (i<*argc); i++) {
		strcat(cmdline,argv[i]);
		strcat(cmdline," ");
	}

	/* Handle the flags argument, which is a bit field
	 * The FF macro returns whether or not the bit is set
	 */
	bPrint        = !FF(PCA_SILENT);

	/* Check ALL the flags ... */
	//max_pa = NPCA_PA + EXTRA_PA + npargs;
	all_pa.clear();
	// add common flags
	fill_cmd_vector( pca_pa, all_pa );

	if (FF(PCA_BE_NICE))
		nicelevel=19;
	all_pa.push_back( nice_pa );

	/* Now append the program specific arguments */
	fill_cmd_vector( pa, all_pa );

	/* set etENUM options to default */
	for(uint i=0; (i<all_pa.size()); i++) {
		if (all_pa[i].type==etENUM) {
			all_pa[i].u.c[0]=all_pa[i].u.c[1];
		}
	}

	/* Now parse all the command-line options */
	get_pargs(argc,argv,all_pa,false);

	/* Parse the file args */
	parse_file_args(argc,argv,nfile,fnm,false);

	/*  for(uint i=0; (i<all_pa.size() ); i++)
		all_pa[i].desc = mk_desc(&(all_pa[i]), time_unit() );
	*/
	bExit = bHelp || (strcmp(manstr[0],"no") != 0);

#ifdef HAVE_UNISTD_H
	if (nicelevel != 0 && !bExit)
		nice(nicelevel);
#endif

	if (!(FF(PCA_QUIET) || bQuiet )) {
		int npall;
		if (bHelp)
			write_man(stderr,"help",program,ndesc,desc,nfile,fnm,npall,all_pa,
		nbugs,bugs,bHidden);
		else if (bPrint) {
			pr_fns(stderr,nfile,fnm);
			print_pargs(stderr,npall,all_pa);
		}
	}

	if ( bExit ) {
		utility_exit();
	};
#undef FF
}

