// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#define CPLUSPLUS
#include <devel/simple_options/option.hh>
// AUTO-REMOVED #include <stdio.h>
// AUTO-REMOVED #include <stdlib.h>

#ifndef WIN32
	#include <pwd.h>
#endif

#include <cstring>
#include <iostream>
#include <queue>

#include <utility/assert.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#else
	#include <unistd.h>
#endif




using namespace std;
#ifdef USE_OLLI_HELP
START_OPT(opt_help,"help","Options to manipulate splash screen behaviour")
REG(BOPT("h",bHelp,false, "just print this help and exit"));
BEGIN_VAR_LIST
bool bHelp;
END_OPT(opt_help)
#endif
namespace devel {
namespace option {




OptionFile::OptionFile(string _in_file,string _out_file) : in_file(_in_file),out_file(_out_file) {
	read_options();
	bFileWritten=false;
}

OptionFile::~OptionFile() {
	if (!bFileWritten)
		write_options();
}

/* base function of virtual function, will just register
	 the module_description and check consistency */
void OptionModule::register_block(OptionBlock& block) {
	if (block.name()!=id) {
		cerr << "try to register block " << block.name().c_str() << " with "
	 << "option module "<< id << ". That should not have happend. Bailing out" << endl;
		exit(1);
	};
	block.set_descr(module_description);
	block.register_module( this );
	//... in derived classes this is the place where all the individual options are registered.
}

void OptionBackend::register_opt( OptionModule const& opt_module ) {
//get module identifier (this is also the blockname in the options file '[ block ]'
	string opt_module_id=opt_module.get_id();
	//module already registered ?
	if (modules.find(opt_module_id)==modules.end()) {
		//no? than create a new module entry

		//get a pointer to an instance of the module, new_block() uses polymorphy
		// to return a pointer to the correct subclass of OptionModule
		OptionModuleOP popt_module =  opt_module.clone();
		modules[opt_module_id]=popt_module;

		//do we have a corresponding block in the option file?
		if (blocks.find(opt_module_id)==blocks.end()) {//no? than create a new entry
			blocks[opt_module_id]= new_block( opt_module_id );
		}

		//connect the module_variables with the entries in the option_file
		popt_module->register_block(*(blocks[opt_module_id]));
	};
	//  opt_module=*(modules[opt_module_id]); //just get a copy
}

void OptionBackend::get_options(OptionModule &opt_module) {
 //get module identifier (this is also the blockname in the options file '[ block ]'
	string opt_module_id=opt_module.get_id();
	opt_module=*(modules[opt_module_id]); //just get a copy
}

void OptionFile::get_options(OptionModule &opt_module) {
	register_opt( opt_module );
	OptionBackend::get_options( opt_module );
}

void CommandLineOptions::get_options(OptionModule &opt_module) {
	if ( processed_ ) {
		OptionBackend::get_options( opt_module );
	}
}

void OptionFile::read_options() {
	FILE *infile=NULL;
	infile=fopen(in_file.c_str(),"r");

	bool bEOF=false;
	string block_name="test_options";
	//get first block_name in file:
	bEOF=OptionFileBlock("").read_entries(infile,block_name);//get first block name
	while(!bEOF) {
		if (blocks.find(block_name)==blocks.end()) {
			//block_name did not exist in list, create new block
			blocks[block_name]= new_block(block_name);
			//bEOF=blocks[block_name]->read_entries(infile,block_name);
		} else {
			cerr << "block : "<< block_name << " has already been read in option file" << endl;
			exit(1);
		}
	}
}

void nice_header (FILE *out,const char *fn)
/* called by write_options()
 it produces a header file for the output-option file */
{
	#ifndef WIN32
	const char *unk = "unknown";
	time_t clock;
	char   *user=NULL;
	int    gh;
	uid_t  uid;
	char   buf[256] = "Unknown";

	struct passwd *pw;
	/* Print a nice header above the file */
	clock = time (0);
	fprintf (out,"%c\n",COMMENTSIGN);
	fprintf (out,"%c\tFile '%s' was generated\n",COMMENTSIGN,fn ? fn : unk);

	#ifndef __CYGWIN__
	uid = getuid();
	pw  = getpwuid(uid);
	#ifndef __native_client__
	gh  = gethostname(buf,255);
	#endif

	user= pw->pw_name;

	fprintf (out,"%c\tBy user: %s (%d)\n",COMMENTSIGN,
		 user ? user : unk,(int) uid);
	fprintf(out,"%c\tOn host: %s\n",COMMENTSIGN,(gh == 0) ? buf : unk);
	#endif

	fprintf (out,"%c\tAt date: %s",COMMENTSIGN,ctime(&clock));
	fprintf (out,"%c\n",COMMENTSIGN);
	#endif
}

void OptionBackend::write_options( std::ostream& out ) {
	char buf[1000];
	for (t_map_blocks::iterator it=blocks.begin();it!=blocks.end();++it) {
		OptionBlock &ablock=*(it->second);
		if (ablock.isKnown()) {
			sprintf(buf,"\n\n---------------------------------------------------------------------\n");out << buf;
			sprintf(buf," %40s \n;%20s\n",ablock.name().c_str(),ablock.descr().c_str()); out << buf;
			sprintf(buf,"---------------------------------------------------------------------\n");out << buf;
			ablock.write_entries(out);
			out << std::endl;
		} else {
			cerr << "unknown block [ " << ablock.name() << " ] is ignored" << endl;
		}
	};
}

void OptionFile::write_options() {
	FILE *out;
	out=fopen(out_file.c_str(),"w");
	nice_header(out,out_file.c_str());
	for (t_map_blocks::iterator it=blocks.begin();it!=blocks.end();++it) {
		OptionBlock &ablock=*(it->second);
		if (ablock.isKnown()) {
			fprintf(out,"[ %s ]\n;%s\n",ablock.name().c_str(),ablock.descr().c_str());
			ablock.write_lines(out);
			fprintf(out,"\n\n");
		} else
			cerr << "unknown block [ " << ablock.name() << " ] is ignored" << endl;
	};
	fclose(out);
	bFileWritten=true;
}

void OptionBlock::write_entries( std::ostream &out) {
	bool bHaltOnUnknown = false;
	char buf[1000];
	/* we want to print the lines in the same sequence as they have been
		 registered. All additional (unknown) lines are sorted behind that and
		 suppressed in the output
	*/
	/* we use a priority queue but it might be fine to use another container
		 or sort the map directly (how)
	*/
	priority_queue< OptionEntry, vector< OptionEntry>, greater< OptionEntry> > sorted_lines;

	/* give the unknown lines a unique entry_number that is bigger than any
		 entry_number of the known lines */
	int mm=-1;
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		mm=max(mm,it->second.entry_number);
	};
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		OptionEntry &line=it->second;
		if (line.entry_number<0)
			line.entry_number=mm++;
		sorted_lines.push(line);
	};
	sprintf(buf,"%12s %6s %6s  %s\n","Option","Type","Value","Description"); out<<buf;
	sprintf(buf,"------------------------------------------------------\n"); out<<buf;

	/* print lines in sorted sequence */
	while(!sorted_lines.empty()) {
		const OptionEntry &line=sorted_lines.top();
		if (line.theOption) {
			base_optOP opt = line.theOption;
			if(line.name[0]==';' || (line.name.length()>2 && line.name[1]==';')) {
	sprintf(buf,"%-24s\n",line.theOption->full_name().c_str()); out << buf;
			} else {
	sprintf(buf,"%12s %6s %6s  %s\n",("-"+opt->full_name()).c_str(),opt->type().c_str(),opt->value().c_str(),opt->descr().c_str()); out << buf;
		//	sprintf(buf,"%-24s = %s\n",line.theOption->full_name().c_str(),line.theOption->value().c_str()); out << buf;
			}
		} else {
			fprintf(stderr,"WARNING: Unknown left-hand '%s' in parameter file\n",
				line.name.c_str());
			if (bHaltOnUnknown) {
	exit(1);
			};
		};
		sorted_lines.pop();
	};
}


void OptionBlock::write_lines(FILE *out) {
	bool bHaltOnUnknown = false;

	/* we want to print the lines in the same sequence as they have been
		 registered. All additional (unknown) lines are sorted behind that and
		 suppressed in the output
	*/
	/* we use a priority queue but it might be fine to use another container
		 or sort the map directly (how)
	*/
	priority_queue<OptionEntry,vector<OptionEntry>,greater<OptionEntry> > sorted_lines;

	/* give the unknown lines a unique entry_number that is bigger than any
		 entry_number of the known lines */
	int mm=-1;
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		mm=max(mm,it->second.entry_number);
	};
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		OptionEntry &line=it->second;
		if (line.entry_number<0)
			line.entry_number=mm++;
		sorted_lines.push(line);
	};

	/* print lines in sorted sequence */
	while(!sorted_lines.empty()) {
		const OptionEntry &line=sorted_lines.top();
		if (line.bKnown) {
			if(line.name[0]==';' || (line.name.length()>2 && line.name[1]==';'))
	fprintf(out,"%-24s\n",line.name.c_str());
			else {
	//int bla=1;
	//fprintf(out,"%-24s = %s\n",line.name.c_str(),line.value.c_str());
			}
		} else {
			fprintf(stderr,"WARNING: Unknown left-hand '%s' in parameter file\n",
				line.name.c_str());
			if (bHaltOnUnknown) {
	exit(1);
			};
		};
		sorted_lines.pop();
	}
}

void CommandLineOptions::process_options_(int& argc, char**argv) {
	processed_ = true;
	for (t_map_blocks::iterator it=blocks.begin();it!=blocks.end();++it) {
		OptionBlock &ablock=*(it->second);
		if (ablock.isKnown()) {
			ablock.process_entries( argc, argv);
		} else
			cerr << "unknown block [ " << ablock.name() << " ] is ignored" << endl;
	};
	write_options( );
	#ifdef USE_OLLI_HELP
	opt_help opts;
	if (opts.bHelp) exit(1);
	#endif
}

void OptionBlock::process_entries( int &argc, char* argv[] ) {
	// int pos = 1;
	std::string name;
	for ( int pos = 1; pos < argc; ++pos ) {
		if ( argv[pos][0]=='-' ) {
			string opt = string ( argv[pos]+1 );
			if ( opt.find("no",0) <=1 ) {
	name = opt.substr(2);
			} else {
	name = opt;
			};
			if ( lines.find( name ) != lines.end() ) {
	pos = lines[ name ].process( opt, pos, argc, argv );
			}
		}
	}
	//my_module->register_block( *this );
}

int OptionEntry::process( std::string const& opt, int pos, int &argc, char* argv[] ) {
	int start_pos = pos;
	pos ++; //skip "-opt"
	if (theOption) theOption->set_name( opt ); // this sets the booleans e.g., nohelp / help

	for ( ; pos < argc; ++pos) {
		if ( argv[pos][0] == '-' ) break;
		if (theOption) theOption->set_value( argv[ pos ]);
	};
	int j,i;
	for (i=pos, j=start_pos; i<argc; i++) {
		argv[j++]=argv[i];
	}
	argc = j ;
	return start_pos - 1;
}

/* ################## block and line implementation ##################*/

/* first some forward declarations of helper functions for string handling */

char *fgets2(char *line, int n, FILE *stream);
/* This routine reads a string from stream of max length n
 * and zero terminated, without newlines
 * line should be long enough (>= n)
 */
void ltrim (char *str);
void rtrim (char *str);
void trim (char *str);

/* this is the actual parser routine, could be changed to base it
	 on the iostream library at some later point

	 brief description:
	 *   lines are read and stripped from pre- and post-trailing whitespace
	 *   comments are ignored.
	 *   the lines are parsed for '[' ']' characters, which constitutes
				 the beginning of a new block.
				 The new block_name is given back in 'string &next_block' and the routine
				 terminated.
	 *  if no new block starts any line containing '=' will
				 be separated into <left of '='> and <right of '='>
				 and stored as an instance of  OptionEntry in lines[]
*/
bool OptionFileBlock::read_entries(FILE *in,string &next_block) {
	char buf[STRLEN], lbuf[STRLEN], rbuf[STRLEN];
	char *ptr, *cptr;
	int /*nin,*/ lc, i, j, k;

	/* leftovers from the GROMACS code */
	const char *fn="option_file"; //used in fprintf(debug,...
	FILE *debug = NULL;

	//nin = 0;  // set but never used ~Labonte
	lc  = 0;
	do {
		ptr=fgets2(buf,STRLEN-1,in);
		lc++;
		if (ptr) {
			/* Strip comment */
			if ((cptr=strchr(buf,COMMENTSIGN)) != NULL)
	*cptr='\0';
			/* Strip spaces */
			trim(buf);

			for(j=0; (buf[j] != '=') && (buf[j] != '\0'); j++)
	;
			if (buf[j] == '\0') {
	if (j > 0) {
		if (debug)
			fprintf(debug,"No = on line %d in file %s, ignored\n",lc,fn);
		// try whether there is a '[' sign
		int kk;
		for(kk=0; (buf[kk] != '[') && (buf[kk] != '\0'); kk++)
			;
		if (buf[kk] != '\0') {
			int rr;
			for(rr=0; (buf[rr] != ']') && (buf[rr] != '\0'); rr++)
				;
			if (buf[rr] != '\0') {
				/* okay, there is '[ name ]' entry */
				// remove the brackets
				buf[rr]='\0';
				buf[kk]=' ';
				trim(buf);
				next_block=string(buf);
				return false; //no EOF
			};
		}
	}
			}
			else { // this line has an '=' sign
	for(i=0; (i<j); i++)
		lbuf[i]=buf[i];
	lbuf[i]='\0';
	trim(lbuf);
	if (lbuf[0] == '\0') {
		if (debug)
			fprintf(debug,"Empty left hand side on line %d in file %s, ignored\n",lc,fn);
	}
	else {
		for(i=j+1,k=0; (buf[i] != '\0'); i++,k++)
			rbuf[k]=buf[i];
		rbuf[k]='\0';
		trim(rbuf);
		if (rbuf[0] == '\0') {
			if (debug)
				fprintf(debug,"Empty right hand side on line %d in file %s, ignored\n",lc,fn);
		}
		else {
			/* Now finally something sensible */
			string name(lbuf);
			lines[name]=OptionEntry(name,string(rbuf));
		}
	}
			}
		}
	} while (ptr);
	next_block="ENDE";
	return true; //EOF
}

/* reg_opt is called (from the REG macro) to
	 register options
	 it uses the polymorphy of set_value()
	 to connect the string line.value with
	 an user defined object depending on the option
*/

void OptionBlock::reg_opt(const base_optOP opt) {
	bKnown=true; // any option makes this block known
	OptionEntry &line=lines[ opt->name() ];
	line.theOption=opt;
	line.bKnown=true;
	// line.name=opt.name();
	// line.comment=opt.comment();
	line.entry_number=number_of_known_entries++;

	//connect option string with variable and set default value if empty string
	//  opt.set_value(line.values,line.truefalse);
}


/* definitions of the c-string helper functions */

char *fgets2(char *line, int n, FILE *stream)
/* This routine reads a string from stream of max length n
 * and zero terminated, without newlines
 * line should be long enough (>= n)
 */
{
	char *c;
	if (fgets(line,n,stream)==NULL) return NULL;
	if ((c=strchr(line,'\n'))!=NULL) *c=0;
	return line;
}


void ltrim (char *str)
{
	char *tr;
	int c;

	if (!str)
		return;

	tr = strdup (str);
	c  = 0;
	while ((tr[c] == ' ') || (tr[c] == '\t'))
		c++;

	strcpy (str,tr+c);
	free (tr);
}

void rtrim (char *str)
{
	int nul;

	if (!str)
		return;

	nul = strlen(str)-1;
	while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
		str[nul] = '\0';
		nul--;
	}
}

void trim (char *str)
{
	ltrim (str);
	rtrim (str);
}


}
}
