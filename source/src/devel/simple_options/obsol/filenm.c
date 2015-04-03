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

#include "filenm.h"


/* Use bitflag ... */
#define IS_SET(fn) ((fn.flag & ffSET) != 0)
#define IS_OPT(fn) ((fn.flag & ffOPT) != 0)
#define IS_MULT(fn) ((fn.flag & ffMULT) != 0)
#define UN_SET(fn) (fn.flag = (fn.flag & ~ffSET))
#define DO_SET(fn) (fn.flag = (fn.flag |  ffSET))

enum { eftASC, eftBIN, eftXDR, eftGEN, eftNR };

#define NTPSS asize(tpss)

typedef struct {
  int  ftype;
  char *ext;
  char *defnm;
  char *defopt;
  const char *descr;
  int  ntps;
  const int  *tps;
} t_deffile;


void pr_def(FILE *fp,int ftp)
{
  t_deffile *df;
  char *s=NULL,*flst;
  const char *ext,*desc;

  df=&(deffile[ftp]);
  /* find default file extension and \tt-ify description */
  flst="";
  if (df->ntps) {
    ext = deffile[df->tps[0]].ext;
    desc= strdup(df->descr);
    s = strstr(desc,": ")+1;
    if (s) {
      s[0] = '\0';
      s++;
      snew(flst,strlen(s)+6);
      strcpy(flst, " \\tt ");
      strcat(flst, s);
    }
  } else {
    ext = df->ext;
    desc= df->descr;
  }
  /* now skip dot */
  if (ext[0])
    ext++;
  else
    ext="";
  /* set file contents type */
  switch (df->ftype) {
  case eftASC: s="Asc";
    break;
  case eftBIN: s="Bin";
    break;
  case eftXDR: s="xdr";
    break;
  case eftGEN: s="";
    break;
  default:
    gmx_fatal(FARGS,"Unimplemented filetype %d %d",ftp,df->ftype);
  }
  fprintf(fp,"\\tt %8s & \\tt %3s & %3s & \\tt %2s & %s%s \\\\[-0.1ex]\n",
	  df->defnm, ext, s, df->defopt ? df->defopt : "",
	  check_tex(desc),check_tex(flst));
}

void pr_fns(FILE *fp,int nf,t_filenm tfn[])
{
  int  i,j,f;
  char buf[256],*wbuf,opt_buf[32];
#define OPTLEN 4
#define NAMELEN 14
  fprintf(fp,"%6s %12s  %-12s %s\n",
	  "Option","Filename","Type","Description");
  fprintf(fp,"------------------------------------------------------------\n");
  for(i=0; (i<nf); i++) {
    for(f=0; (f<tfn[i].nfiles); f++) {
      sprintf(buf, "%4s %14s  %-12s ",
	      (f==0) ? tfn[i].opt : "",tfn[i].fns[f],
	      (f==0) ? fileopt(tfn[i].flag,opt_buf,32) : "");
      if ( f < tfn[i].nfiles-1 )
	fprintf(fp, "%s\n", buf);
    }
    if (tfn[i].nfiles > 0) {
      strcat(buf, deffile[tfn[i].ftp].descr);
      if ( (strlen(tfn[i].opt)>OPTLEN) &&
	   (strlen(tfn[i].opt)<=
	    ((OPTLEN+NAMELEN)-strlen(tfn[i].fns[tfn[i].nfiles-1]))) ) {
	for(j=strlen(tfn[i].opt);
	    j<strlen(buf)-(strlen(tfn[i].opt)-OPTLEN)+1; j++)
	  buf[j]=buf[j+strlen(tfn[i].opt)-OPTLEN];
      }
      wbuf=wrap_lines(buf,78,35,FALSE);
      fprintf(fp,"%s\n",wbuf);
      sfree(wbuf);
    }
  }
  fprintf(fp,"\n");
  fflush(fp);
}

void pr_fopts(FILE *fp,int nf,t_filenm tfn[], int shell)
{
  int i,j;

  switch (shell) {
  case eshellCSH:
    for(i=0; (i<nf); i++) {
      fprintf(fp," \"n/%s/f:*.",tfn[i].opt);
      if (deffile[tfn[i].ftp].ntps) {
	fprintf(fp,"{");
	for(j=0; j<deffile[tfn[i].ftp].ntps; j++) {
	  if (j>0)
	    fprintf(fp,",");
	  fprintf(fp,"%s",deffile[deffile[tfn[i].ftp].tps[j]].ext+1);
	}
	fprintf(fp,"}");
      } else
	fprintf(fp,"%s",deffile[tfn[i].ftp].ext+1);
      fprintf(fp,"{");
      for(j=0; j<NZEXT; j++)
	fprintf(fp,",%s",z_ext[j]);
      fprintf(fp,"}/\"");
    }
    break;
  case eshellBASH:
    for(i=0; (i<nf); i++) {
	fprintf(fp,"%s) COMPREPLY=( $(compgen -X '!*.",tfn[i].opt);
      if (deffile[tfn[i].ftp].ntps) {
	fprintf(fp,"+(");
	for(j=0; j<deffile[tfn[i].ftp].ntps; j++) {
	  if (j>0)
	    fprintf(fp,"|");
	  fprintf(fp,"%s",deffile[deffile[tfn[i].ftp].tps[j]].ext+1);
	}
	fprintf(fp,")");
      } else
	fprintf(fp,"%s",deffile[tfn[i].ftp].ext+1);
      fprintf(fp,"*(");
      for(j=0; j<NZEXT; j++) {
	if (j>0)
	  fprintf(fp,"|");
	fprintf(fp,"%s",z_ext[j]);
      }
      fprintf(fp,")' -f $c ; compgen -S '/' -X '.*' -d $c ));;\n");
    }
    break;
  case eshellZSH:
    for(i=0; (i<nf); i++) {
      fprintf(fp,"- 'c[-1,%s]' -g '*.",tfn[i].opt);
      if (deffile[tfn[i].ftp].ntps) {
	fprintf(fp,"(");
	for(j=0; j<deffile[tfn[i].ftp].ntps; j++) {
	  if (j>0)
	    fprintf(fp,"|");
	  fprintf(fp,"%s",deffile[deffile[tfn[i].ftp].tps[j]].ext+1);
	}
	fprintf(fp,")");
      } else
	fprintf(fp,"%s",deffile[tfn[i].ftp].ext+1);
      fprintf(fp,"(");
      for(j=0; j<NZEXT; j++)
	fprintf(fp,"|%s",z_ext[j]);
      fprintf(fp,") *(/)' ");
    }
    break;
  }
}

static void check_opts(int nf,t_filenm fnm[])
{
  int       i;
  t_deffile *df;

  for(i=0; (i<nf); i++) {
    df=&(deffile[fnm[i].ftp]);
    if (fnm[i].opt == NULL) {
      if (df->defopt == NULL)
	gmx_fatal(FARGS,"No default cmd-line option for %s (type %d)\n",
		    deffile[fnm[i].ftp].ext,fnm[i].ftp);
      else
	fnm[i].opt=df->defopt;
    }
  }
}

int fn2ftp(char *fn)
{
  int  i,len;
  char *feptr,*eptr;

  if (!fn)
    return efNR;

  len=strlen(fn);
  if ((len >= 4) && (fn[len-4] == '.'))
    feptr=&(fn[len-4]);
  else
    return efNR;

  for(i=0; (i<efNR); i++)
    if ((eptr=deffile[i].ext) != NULL)
      if (strcasecmp(feptr,eptr)==0)
	break;

  return i;
}

static void set_extension(char *buf,int ftp)
{
  int len,extlen;
  t_deffile *df;

  /* check if extension is already at end of filename */
  df=&(deffile[ftp]);
  len=strlen(buf);
  extlen = strlen(df->ext);
  if ((len <= extlen) || (strcasecmp(&(buf[len-extlen]),df->ext) != 0))
    strcat(buf,df->ext);
}

static void add_filenm(t_filenm *fnm, char *filenm)
{
  srenew(fnm->fns, fnm->nfiles+1);
  fnm->fns[fnm->nfiles] = strdup(filenm);
  fnm->nfiles++;
}

static void set_grpfnm(t_filenm *fnm,char *name,bool bCanNotOverride)
{
  char buf[256],buf2[256];
  int  i,type;
  bool bValidExt;
  int  nopts;
  const int  *ftps;

  nopts = deffile[fnm->ftp].ntps;
  ftps  = deffile[fnm->ftp].tps;
  if ((nopts == 0) || (ftps == NULL))
    gmx_fatal(FARGS,"nopts == 0 || ftps == NULL");

  bValidExt = FALSE;
  if (name && (bCanNotOverride || (default_file_name == NULL))) {
    strcpy(buf,name);
    /* First check whether we have a valid filename already */
    type = fn2ftp(name);
    for(i=0; (i<nopts) && !bValidExt; i++)
      if (type == ftps[i])
	bValidExt = TRUE;
  } else
    /* No name given, set the default name */
    strcpy(buf,ftp2defnm(fnm->ftp));

  if (!bValidExt && (fnm->flag & ffREAD)) {
    /* for input-files only: search for filenames in the directory */
    for(i=0; (i<nopts) && !bValidExt; i++) {
      type = ftps[i];
      strcpy(buf2,buf);
      set_extension(buf2,type);
      if (fexist(buf2)) {
	bValidExt = TRUE;
	strcpy(buf,buf2);
      }
    }
  }

  if (!bValidExt)
    /* Use the first extension type */
    set_extension(buf,ftps[0]);

  add_filenm(fnm, buf);
}

static void set_filenm(t_filenm *fnm,char *name,bool bCanNotOverride)
{
  /* Set the default filename, extension and option for those fields that
   * are not already set. An extension is added if not present, if fn = NULL
   * or empty, the default filename is given.
   */
  char buf[256];
  int  i,len,extlen;

  if ((fnm->ftp < 0) || (fnm->ftp >= efNR))
    gmx_fatal(FARGS,"file type out of range (%d)",fnm->ftp);

  if ((fnm->flag & ffREAD) && name && fexist(name)) {
    /* check if filename ends in .gz or .Z, if so remove that: */
    len    = strlen(name);
    for (i=0; i<NZEXT; i++) {
      extlen = strlen(z_ext[i]);
      if (len > extlen)
	if (strcasecmp(name+len-extlen,z_ext[i]) == 0) {
	  name[len-extlen]='\0';
	  break;
	}
    }
  }

  if (deffile[fnm->ftp].ntps)
    set_grpfnm(fnm,name,bCanNotOverride);
  else {
    if ((name != NULL) && (bCanNotOverride || (default_file_name == NULL)))
      strcpy(buf,name);
    else
      strcpy(buf,deffile[fnm->ftp].defnm);
    set_extension(buf,fnm->ftp);

    add_filenm(fnm, buf);
  }
}

static void set_filenms(int nf,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nf); i++)
    if (!IS_SET(fnm[i]))
      set_filenm(&(fnm[i]),fnm[i].fn,FALSE);
}

void print_tty_formatted(FILE *out, int nldesc, char **desc,int indent)
{
  char *buf;
  char *temp;
  int i,j;

  /* Just to be sure */
  j=0;
  for(i=0; (i<nldesc); i++)
    j+=strlen(desc[i])+10;
  snew(buf,j);
  for(i=0; (i<nldesc); i++) {
    if ((strlen(buf)>0) &&
	(buf[strlen(buf)-1] !=' ') && (buf[strlen(buf)-1] !='\n'))
      strcat(buf," ");
    temp=check_tty(desc[i]);
    strcat(buf,temp);
    sfree(temp);
  }
  /* Make lines of at most 79 characters */
  temp = wrap_lines(buf,78,indent,FALSE);
  fprintf(out,"%s\n",temp);
  sfree(temp);
  sfree(buf);
}

void write_man(FILE *out,char *mantp,
	       char *program,
	       int nldesc,char **desc,
	       int nfile,t_filenm *fnm,
	       int npargs,t_pargs *pa,
	       int nbug,char **bugs,bool bHidden)
{
  char    *pr;
  int     i,npar;
  t_pargs *par;

  /* Don't write hidden options to completions, it just
   * makes the options more complicated for normal users
   */

  if (bHidden) {
    npar=npargs;
    par=pa;
  }
  else {
    snew(par,npargs);
    npar=0;
    for(i=0;i<npargs;i++)
      if (!is_hidden(&pa[i])) {
	par[npar]=pa[i];
	npar++;
      }
  }

  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else
    pr+=1;
  if (strcmp(mantp,"tex")==0)
    write_texman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"nroff")==0)
    write_nroffman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"ascii")==0)
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,TRUE);
  if (strcmp(mantp,"help")==0)
    write_ttyman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs,FALSE);
  if (strcmp(mantp,"html")==0)
    write_htmlman(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"py")==0)
    write_py(out,pr,nldesc,desc,nfile,fnm,npar,par,nbug,bugs);
  if (strcmp(mantp,"xml")==0)
    write_xmlman(out,pr,nldesc,desc,nfile,fnm,npargs,pa,nbug,bugs);
  if (strcmp(mantp,"completion-zsh")==0)
    write_zshcompl(out,nfile,fnm,npar,par);
  if (strcmp(mantp,"completion-bash")==0)
    write_bashcompl(out,nfile,fnm,npar,par);
  if (strcmp(mantp,"completion-csh")==0)
    write_cshcompl(out,nfile,fnm,npar,par);

  if (!bHidden)
    sfree(par);
}

void write_ttyman(FILE *out,
			 char *program,
			 int nldesc,char **desc,
			 int nfile,t_filenm *fnm,
			 int npargs,t_pargs *pa,
			 int nbug,char **bugs,bool bHeader)
{
  int i;
  char buf[256];
  char *tmp;

  if (bHeader) {
    fprintf(out,"%s\n\n",check_tty(program));
    fprintf(out,"%s\n%s\n",GromacsVersion(),mydate(buf,255));
  }
  if (nldesc > 0) {
    fprintf(out,"DESCRIPTION\n-----------\n");
    print_tty_formatted(out,nldesc,desc,0);
  }
  if (nbug > 0) {
    fprintf(out,"\n");
    fprintf(out,"KNOWN BUGS\n----------\n");
    for(i=0; i<nbug; i++) {
      snew(tmp,strlen(bugs[i])+3);
      strcpy(tmp,"* ");
      strcpy(tmp+2,check_tty(bugs[i]));
      fprintf(out,"%s\n",wrap_lines(tmp,78,2,FALSE));
      sfree(tmp);
    }
  }
  if (nfile > 0) {
    fprintf(out,"\n");
    pr_fns(out,nfile,fnm);
  }
  if (npargs > 0) {
    print_pargs(out,npargs,pa);
  }
}


void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
		     bool bKeep)
{
  int  i,j;
  bool *bRemove;

  check_opts(nf,fnm);

  for(i=0; (i<nf); i++)
    UN_SET(fnm[i]);

  if (*argc > 1) {
    snew(bRemove,(*argc)+1);
    i=1;
    do {
      for(j=0; (j<nf); j++) {
	if (strcmp(argv[i],fnm[j].opt) == 0) {
	  DO_SET(fnm[j]);
	  bRemove[i]=TRUE;
	  i++;
	  /* check if we are out of arguments for this option */
	  if ( (i >= *argc) || (argv[i][0] == '-') )
	    set_filenm(&fnm[j],fnm[j].fn,FALSE);
	  /* sweep up all file arguments for this option */
	  while ((i < *argc) && (argv[i][0] != '-')) {
	    set_filenm(&fnm[j],argv[i],TRUE);
	    bRemove[i]=TRUE;
	    i++;
	    /* only repeat for 'multiple' file options: */
	    if ( ! IS_MULT(fnm[j]) )
	      break;
	  }

	  break; /* jump out of 'j' loop */
	}
      }
      /* No file found corresponding to option argv[i] */
      if (j == nf)
	i++;
    } while (i < *argc);

    if (!bKeep) {
      /* Remove used entries */
      for(i=j=0; (i<=*argc); i++) {
	if (!bRemove[i])
	  argv[j++]=argv[i];
      }
      (*argc)=j-1;
    }
    sfree(bRemove);
  }

  set_filenms(nf,fnm);
}

char *opt2fn(char *opt,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0) {
      return fnm[i].fns[0];
    }

  fprintf(stderr,"No option %s\n",opt);

  return NULL;
}

int opt2fns(char **fns[], char *opt,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0) {
      *fns = fnm[i].fns;
      return fnm[i].nfiles;
    }

  fprintf(stderr,"No option %s\n",opt);
  return 0;
}

char *ftp2fn(int ftp,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return fnm[i].fns[0];

  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}

int ftp2fns(char **fns[], int ftp,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp) {
      *fns = fnm[i].fns;
      return fnm[i].nfiles;
    }

  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return 0;
}

bool ftp2bSet(int ftp,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return (bool) IS_SET(fnm[i]);

  fprintf(stderr,"ftp2bSet: No filetype %s\n",deffile[ftp].ext);

  return FALSE;
}

bool opt2bSet(char *opt,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      return (bool) IS_SET(fnm[i]);

  fprintf(stderr,"No option %s\n",opt);

  return FALSE;
}

char *opt2fn_null(char *opt,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0) {
      if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
	return NULL;
      else
	return fnm[i].fns[0];
    }
  fprintf(stderr,"No option %s\n",opt);
  return NULL;
}

char *ftp2fn_null(int ftp,int nfile,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp) {
      if (IS_OPT(fnm[i]) && !IS_SET(fnm[i]))
	return NULL;
      else
	return fnm[i].fns[0];
    }
  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}

static void add_filters(char *filter,int *n,int nf,const int ftp[])
{
  char buf[8];
  int  i;

  sprintf(filter,"*.{");
  for(i=0; (i<nf); i++) {
    sprintf(buf,"%s",ftp2ext(ftp[i]));
    if (*n > 0)
      strcat(filter,",");
    strcat(filter,buf);
    (*n) ++;
  }
  strcat(filter,"}");
}

char *ftp2filter(int ftp)
{
  int    n;
  static char filter[128];

  filter[0] = '\0';
  n         = 0;
  switch (ftp) {
  case efENX:
    add_filters(filter,&n,NENXS,enxs);
    break;
  case efTRX:
    add_filters(filter,&n,NTRXS,trxs);
    break;
  case efTRN:
    add_filters(filter,&n,NTRNS,trns);
    break;
  case efSTO:
    add_filters(filter,&n,NSTOS,stos);
    break;
  case efSTX:
    add_filters(filter,&n,NSTXS,stxs);
    break;
  case efTPX:
    add_filters(filter,&n,NTPXS,tpxs);
    break;
  default:
    sprintf(filter,"*%s",ftp2ext(ftp));
    break;
  }
  return filter;
}

bool is_optional(t_filenm *fnm)
{
  return ((fnm->flag & ffOPT) == ffOPT);
}

bool is_output(t_filenm *fnm)
{
  return ((fnm->flag & ffWRITE) == ffWRITE);
}

bool is_set(t_filenm *fnm)
{
  return ((fnm->flag & ffSET) == ffSET);
}


