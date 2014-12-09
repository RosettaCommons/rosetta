// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#define CPLUSPLUS
#include "option.h"
#include <stdio.h>
#include <stdlib.h>
#include <pwd.h>
#include <iostream>
#include <queue>
#include "smalloc.h"
using namespace option;
using namespace std;
t_option_file::t_option_file(string _in_file,string _out_file) : in_file(_in_file),out_file(_out_file) {
	read_options();
};

void t_option_file::get_options(t_option_module &opt_module) {
	string opt_module_id=opt_module.get_id();
	if (modules.find(opt_module_id)!=modules.end()) {
		cout << "found opt_module_id:" << opt_module_id << endl;
		opt_module=*(modules[opt_module_id]);
		cout << "in map " << dynamic_cast<my_test_options*>(modules[opt_module_id])->number << " out map " << dynamic_cast<my_test_options&>(opt_module).number << endl;
	} else {
		cout << "register new option block:" << opt_module_id << endl;
		t_option_module* popt_module =  opt_module.new_block();
		modules[opt_module_id]=popt_module;
		popt_module->register_options(blocks[opt_module_id]);
		opt_module=*popt_module;
	};
};

t_option_file::~t_option_file() {
	int i;
	write_options();
	/*for (i=0; (i<ninp); i++) {
		sfree(inp[i].name);
		sfree(inp[i].value);
	}
	sfree(inp); */
};

void t_option_file::read_options() {
	//* spaeter liste von bloecken lesen */
	// jetzt so tun, als ob erstes block-label schon gelesen wurde:

	//inp=read_inpfile(const_cast<char*>(in_file.c_str()),&ninp);
	FILE *infile=NULL;
	infile=fopen(in_file.c_str(),"r");
	//read blocks
	bool bEOF=false;
	string block_name="test_options";
	bEOF=t_option_block("").read_lines(infile,block_name);//get first block name
	while(!bEOF) {
		if (blocks.find(block_name)==blocks.end()) {
			//create new block
			blocks[block_name]=t_option_block(block_name);
			bEOF=blocks[block_name].read_lines(infile,block_name);
		} else {
			cerr << "block : "<< block_name << " has already been read in option file" << endl;
			exit(1);
		}
	}
	cerr << "ending reading of inputs before block " << block_name << endl;
};

void nice_header (FILE *out,const char *fn)
{
	char   *unk = "onbekend";
	time_t clock;
	char   *user=NULL;
	int    gh;
	uid_t  uid;
	char   buf[256];

	struct passwd *pw;
	/* Print a nice header above the file */
	clock = time (0);
	fprintf (out,"%c\n",COMMENTSIGN);
	fprintf (out,"%c\tFile '%s' was generated\n",COMMENTSIGN,fn ? fn : unk);

	uid = getuid();
	pw  = getpwuid(uid);
	gh  = gethostname(buf,255);
	user= pw->pw_name;

	fprintf (out,"%c\tBy user: %s (%d)\n",COMMENTSIGN,
		 user ? user : unk,(int) uid);
	fprintf(out,"%c\tOn host: %s\n",COMMENTSIGN,(gh == 0) ? buf : unk);

	fprintf (out,"%c\tAt date: %s",COMMENTSIGN,ctime(&clock));
	fprintf (out,"%c\n",COMMENTSIGN);
}
bool comp_lines_by_entry_nr(const t_option_line &lin1,const t_option_line &lin2) {
	return lin1.entry_number<lin2.entry_number;
};
void t_option_block::write_lines(FILE *out) {
	bool bHaltOnUnknown = false;
	priority_queue<t_option_line,vector<t_option_line>,greater<t_option_line> > sorted_lines;

	int mm=-1;
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		mm=max(mm,it->second.entry_number);
		//    cerr << it->second.name <<" (" << it->second.entry_number << ")" <<endl;
	};
	for (t_map_lines::iterator it=lines.begin();it!=lines.end();++it) {
		t_option_line &line=it->second;
		if (line.entry_number<0)
			line.entry_number=mm++;
		sorted_lines.push(line);
	};
	while(!sorted_lines.empty()) {
		const t_option_line &line=sorted_lines.top();
		if (line.bKnown) {
			if(line.name[0]==';' || (line.name.length()>2 && line.name[1]==';'))
	fprintf(out,"%-24s\n",line.name.c_str());
			else
	fprintf(out,"%-24s = %s\n",line.name.c_str(),line.value.c_str());
		} else {
			fprintf(stderr,"WARNING: Unknown left-hand '%s' in parameter file\n",
				line.name.c_str());
			if (bHaltOnUnknown) {
	exit(1);
			};
		};
		sorted_lines.pop();
	};
};
void t_option_file::write_options() {
	//  write_inpfile(const_cast<char*>(out_file.c_str()),ninp,inp,FALSE);
	FILE *out;
	int  i;

	out=fopen(out_file.c_str(),"w");
	nice_header(out,out_file.c_str());


	for (t_map_blocks::iterator it=blocks.begin();it!=blocks.end();++it) {
	//  sort_inp(ninp,inp);
		t_option_block &ablock=it->second;
		fprintf(out,"[ %s ]\n;%s\n",ablock.name().c_str(),ablock.descr().c_str());

		ablock.write_lines(out);
	};
	fclose(out);
};







/* ################## block and line implementation ##################*/
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

bool t_option_block::read_lines(FILE *in,string &next_block) {
	char      buf[STRLEN],lbuf[STRLEN],rbuf[STRLEN];
	char      *ptr,*cptr;
	t_inpfile *inp=NULL;
	int       nin,lc,i,j,k;

	/* leftovers from the GROMACS code */
	char *fn="option_file"; //used in fprintf(debug,...
	FILE *debug = NULL;

	nin = lc  = 0;
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
			lines[name]=t_option_line(name,string(rbuf));

		}
	}
			}
		}
	} while (ptr);

	next_block="ENDE";
	return true; //EOF
};


/*void my_test_basic::options::construct_option_table() {
	//  add_line(new SOPT("title",title,"my_title","this is important"));
	};*/

void t_option_module::register_options(t_option_block& block) {
	if (block.name()!=id) {
		cerr << "try to register file_block " << block.name().c_str() << " with "
	 << "option module "<< id << ". That should not have happend. Bailing out" << endl;
		exit(1);
	};
	block.set_descr(module_description);
};
void t_option_block::reg_opt(const base_opt& opt) {
	t_option_line &line=lines[opt.name()];
	line.bKnown=true;
	line.name=opt.name();
	line.comment=opt.comment();
	line.entry_number=number_of_known_entries++;
	opt.set_value(line.value);
};


void my_test_basic::options::register_options(t_option_block& block) {
	t_option_module::register_options(block);
	block.reg_opt(COPT("VARIOUS TEST OPTIONS"));
	block.reg_opt(SOPT("title",title,"","this is the title"));
	block.reg_opt(IOPT("bla",number,3,"this defines the number of runs"));

	/* char      *tmp;   //used by STYPE macro
	string     dummy; //used by CTYPE/ CCTYPE macros
	 CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
	STYPE ("title",	title,	NULL);
	CTYPE ("Preprocessor - specify a full path if necessary.");
	STYPE ("cpp",		cpp,	"cpp");
	STYPE ("include",	include,	NULL);
	STYPE ("define",	define,	NULL);
	ITYPE("bla",number,3);*/
};

using namespace option;
int main(int argc, char **argv) {
		t_option_file option_file("in.mdp","out.mdp");
		my_test_options op1;
		my_test_options op2;
		option_file.get_options(op1);
		option_file.get_options(op2);
		cout << "op1" << op1.cpp << " op2 " << op2.cpp << endl;


};
