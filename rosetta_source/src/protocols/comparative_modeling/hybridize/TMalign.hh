// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief A reimplementation of TM-align algorithm
/// @author Yifan Song

#ifndef INCLUDED_protocols_comparative_modeling_hybridize_TMalign_hh
#define INCLUDED_protocols_comparative_modeling_hybridize_TMalign_hh

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a>b?b:a
//#define MAXLEN 10000                        //maximum length of filenames

//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
//#include <time.h>
#include <string>
//#include <malloc.h>

//#include <sstream>
//#include <iostream>
//#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace std;
using core::Size;
	//using namespace numeric;

class TMalign {
/*
===============================================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   Center for Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
                                                       
           
   Please report bugs and questions to yangji@umich.edu or or zhng@umich.edu
===============================================================================
*/
public:


char version[20];                          //version 
	
	
//global variables
double D0_MIN;                             //for d0
double Lnorm;                              //normalization length
double score_d8, d0, d0_search, dcu0;      //for TMscore search
std::vector < std::vector < double > > score;                  //Input score table for dynamic programming
std::vector < std::vector < bool > >   path;                             //for dynamic programming  
std::vector < std::vector < double > > val;                              //for dynamic programming  
int    xlen, ylen, minlen;                 //length of proteins
std::vector < numeric::xyzVector < core::Real > > xa, ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
//in general, ya is regarded as native structure --> superpose xa onto ya
std::vector < int >   xresno, yresno;                   //residue numbers, used in fragment gapless threading 
vector < numeric::xyzVector <core::Real> > xtm, ytm;                       //for TMscore search engine
vector < numeric::xyzVector <core::Real> > xt;                               //for saving the superposed version of r_1 or xtm
std::string   seqx, seqy;                       //for the protein sequence 
std::vector < int >   secx, secy;                       //for the secondary structure 
vector < numeric::xyzVector <core::Real> > r1, r2;                         //for Kabsch rotation 
numeric::xyzVector <core::Real> t;
numeric::xyzMatrix <core::Real> u;                      //Kabsch translation vector and rotation matrix

int n_ali8_;
std::vector < int > m1_, m2_;
double d0_out_;
	
//argument variables
//char out_reg[MAXLEN];
double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
bool o_opt, a_opt, u_opt, d_opt, v_opt;
double TM3, TM4, TM5;



void PrintErrorAndQuit(string sErrorString)
{
	cout << sErrorString << endl;
	exit(1);
}

template < class A > void ResizeArray( A & array, int Narray1, int Narray2)
{
	array.resize(Narray1);
	for(int i=0; i<Narray1; i++) array[i].resize(Narray2);
}

/*
template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
	*array=new A* [Narray1];
	for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};

template <class A> void DeleteArray(A *** array, int Narray)
{
	for(int i=0; i<Narray; i++)
		if(*(*array+i)) delete [] *(*array+i);
	if(Narray) delete [] (*array);
	(*array)=NULL;
};
*/

char AAmap(string AA)
{
	char A=' ';
	if(     AA.compare("BCK")==0)   A='X';
	else if(AA.compare("GLY")==0)   A='G';
	else if(AA.compare("ALA")==0)   A='A';
	else if(AA.compare("SER")==0)   A='S';
	else if(AA.compare("CYS")==0)   A='C';
	else if(AA.compare("VAL")==0)   A='V';     
	else if(AA.compare("THR")==0)   A='T';
	else if(AA.compare("ILE")==0)   A='I';
	else if(AA.compare("PRO")==0)   A='P';
	else if(AA.compare("MET")==0)   A='M';
	else if(AA.compare("ASP")==0)   A='D';
	else if(AA.compare("ASN")==0)   A='N';
	else if(AA.compare("LEU")==0)   A='L';
	else if(AA.compare("LYS")==0)   A='K';
	else if(AA.compare("GLU")==0)   A='E';
	else if(AA.compare("GLN")==0)   A='Q';
	else if(AA.compare("ARG")==0)   A='R';
	else if(AA.compare("HIS")==0)   A='H';
	else if(AA.compare("PHE")==0)   A='F';
	else if(AA.compare("TYR")==0)   A='Y';
	else if(AA.compare("TRP")==0)   A='W';    
	else if(AA.compare("CYX")==0)   A='C';
	else
		A='Z'; //ligand
	
	return A;
}

void AAmap3(char A, char AA[3])
{
	if     ( A=='X')   strcpy(AA, "BCK");
	else if( A=='G')   strcpy(AA, "GLY");
	else if( A=='A')   strcpy(AA, "ALA");
	else if( A=='S')   strcpy(AA, "SER");
	else if( A=='C')   strcpy(AA, "CYS");
	else if( A=='V')   strcpy(AA, "VAL");
	else if( A=='T')   strcpy(AA, "THR");
	else if( A=='I')   strcpy(AA, "ILE");
	else if( A=='P')   strcpy(AA, "PRO");
	else if( A=='M')   strcpy(AA, "MET");
	else if( A=='D')   strcpy(AA, "ASP");
	else if( A=='N')   strcpy(AA, "ASN");
	else if( A=='L')   strcpy(AA, "LEU");
	else if( A=='K')   strcpy(AA, "LYS");
	else if( A=='E')   strcpy(AA, "GLU");
	else if( A=='Q')   strcpy(AA, "GLN");
	else if( A=='R')   strcpy(AA, "ARG");
	else if( A=='H')   strcpy(AA, "HIS");
	else if( A=='F')   strcpy(AA, "PHE");
	else if( A=='Y')   strcpy(AA, "TYR");
	else if( A=='W')   strcpy(AA, "TRP");
	else if( A=='C')   strcpy(AA, "CYX");
	else
		strcpy(AA, "UNK");           
}

/*
void get_xyz(string line, double *x, double *y, double *z, char *resname, int *no)
{
	char cstr[50];    
	
	strcpy(cstr, (line.substr(30, 8)).c_str());
	sscanf(cstr, "%lf", x);
	
	strcpy(cstr, (line.substr(38, 8)).c_str());
	sscanf(cstr, "%lf", y);  
	
	strcpy(cstr, (line.substr(46, 8)).c_str());
	sscanf(cstr, "%lf", z);
	
	strcpy(cstr, (line.substr(17, 3)).c_str());
	*resname=AAmap(cstr);
	
	strcpy(cstr, (line.substr(22, 4)).c_str());
	sscanf(cstr, "%d", no);
}
*/
	
	/*
int get_PDB_len(char *filename)
{
	int i=0;   
	string line;
	string atom ("ATOM "); 
	
	
	ifstream fin (filename);
	if (fin.is_open())
	{
		while ( fin.good() )
		{
			getline(fin, line);
			if(line.compare(0, atom.length(), atom)==0)
			{
				if( line.compare(12, 4, "CA  ")==0 ||\
				   line.compare(12, 4, " CA ")==0 ||\
				   line.compare(12, 4, "  CA")==0 )
				{
					if( line.compare(16, 1, " ")==0 ||\
					   line.compare(16, 1, "A")==0 )
					{                                  
						i++;
					}
				}                  
			}            
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	}
	
	return i;    
}
int read_PDB(char *filename, double **a, char *seq, int *resno)
{
	int i=0;
	string line, str;    
	string atom ("ATOM "); 
	
	
	ifstream fin (filename);
	if (fin.is_open())
	{
		while ( fin.good() )
		{
			getline(fin, line);
			if(line.compare(0, atom.length(), atom)==0)
			{
				if( line.compare(12, 4, "CA  ")==0 ||\
				   line.compare(12, 4, " CA ")==0 ||\
				   line.compare(12, 4, "  CA")==0 )
				{
					if( line.compare(16, 1, " ")==0 ||\
					   line.compare(16, 1, "A")==0 )
					{  
						get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
						i++;
					}
				}                  
			}            
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	} 
	seq[i]='\0';   
	
	return i;
}
	 */

int read_pose(core::pose::Pose const & pose, std::list <core::Size> const & residue_list, std::vector < numeric::xyzVector < core::Real > > & a, string & seq, std::vector < int > & resno)
{
	int i=0;
	seq.resize(residue_list.size());
	for (std::list<core::Size>::const_iterator it = residue_list.begin();
		 it != residue_list.end();
		 ++it) {
		core::Size ires = *it;
		if (!pose.residue_type(ires).is_protein()) continue;
		
        a[i] = pose.residue(ires).xyz("CA");
		//a[i][0] = pose.residue(ires).xyz("CA").x();
		//a[i][1] = pose.residue(ires).xyz("CA").y();
		//a[i][2] = pose.residue(ires).xyz("CA").z();
		seq[i] = pose.residue_type(ires).name1();
		resno[i] = pose.pdb_info()->number(ires);
		i++;
	}
	//std::cout << seq << std::endl;
	return residue_list.size();
}

	/*
int get_ligand_len(char *filename)
{
	int i=0;
	char cstr[100];    
	string line;    
	string atom ("HETATM "); 
	string finish ("END"); 
	
	
	ifstream fin (filename);
	if (fin.is_open())
	{
		while ( fin.good() )
		{
			getline(fin, line);
			if(line.compare(0, atom.length(), atom)==0)
			{
				strcpy(cstr, (line.substr(12, 4)).c_str()); 
				
				if(!strstr(cstr, "H"))
				{
					if( line.compare(16, 1, " ")==0 ||\
					   line.compare(16, 1, "A")==0 )
					{
						i++;
					}
				}                  
			}
			else if(line.compare(0, finish.length(), finish)==0) 
			{
				break;
			}          
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	} 
	
	return i;
}


int read_ligand(char *filename, double **a, char *seq, int *resno)
{
	int i=0;
	char cstr[100];    
	string line, str;    
	string atom ("HETATM "); 
	string finish ("END"); 
	
	
	ifstream fin (filename);
	if (fin.is_open())
	{
		while ( fin.good() )
		{
			getline(fin, line);
			if(line.compare(0, atom.length(), atom)==0)
			{
				strcpy(cstr, (line.substr(12, 4)).c_str()); 
				
				if(!strstr(cstr, "H"))
				{
					if( line.compare(16, 1, " ")==0 ||\
					   line.compare(16, 1, "A")==0 )
					{                                  
						get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
						i++;
					}
				}                  
			}
			else if(line.compare(0, finish.length(), finish)==0) 
			{
				break;
			}          
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);     
	} 
	seq[i]='\0';   
	
	return i;
}
*/

double dist(numeric::xyzVector<core::Real> x, numeric::xyzVector<core::Real> y)
{
	//double d1=x[0]-y[0];
	//double d2=x[1]-y[1];
	//double d3=x[2]-y[2];	
	
	//return (d1*d1 + d2*d2 + d3*d3);
	return x.distance(y);
}

	/*
double dot(std::vector < double > a, std::vector < double > b)
{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}
	 */

void transform(numeric::xyzVector <core::Real> const t, numeric::xyzMatrix <core::Real> const u, numeric::xyzVector < core::Real > x, numeric::xyzVector < core::Real > & x1)
{
	x1 = t + u*x;
	//x1[0]=t[0]+dot(&u[0][0], x);
	//x1[1]=t[1]+dot(&u[1][0], x);
	//x1[2]=t[2]+dot(&u[2][0], x);
}

void do_rotation(std::vector < numeric::xyzVector<core::Real> > const x, std::vector < numeric::xyzVector<core::Real> > & x1, int len, numeric::xyzVector<core::Real> const t, numeric::xyzMatrix<core::Real> const u)
{
	for(int i=0; i<len; i++)
	{
		transform(t, u, x[i], x1[i]);
	}    
}

/*
void output_align1(int *invmap0, int len)
{
	for(int i=0; i<len; i++)
	{
		if(invmap0[i]>=0)
		{
			cout << invmap0[i]+1 << " ";
		}			
		else
			cout << invmap0[i] << " ";
		
	}	
	cout << endl << endl;	
}


int output_align(int *invmap0, int len)
{
	int n_ali=0;
	for(int i=0; i<len; i++)
	{		
		cout <<  invmap0[i] << " ";
		n_ali++;		
	}	
	cout << endl << endl;	
	
	return n_ali;
}
*/
/*    Please note this fucntion is not a correct implementation of 
 *     the N-W dynamic programming because the score tracks back only 
 *     one layer of the matrix. This code was exploited in TM-align 
 *     because it is about 1.5 times faster than a complete N-W code
 *     and does not influence much the final structure alignment result.
 */
void NWDP_TM(int const len1, int const len2, double const gap_open, vector < int > & j2i)
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: score[1:len1, 1:len2], and gap_open
	//Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	
	double h, v, d;
	
	//initialization
	val[0][0]=0;
	for(Size i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}
	
	for(Size j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2]
	}      
	
	
	//decide matrix and path
	for(Size i=1; i<=len1; i++)	
	{	
		for(Size j=1; j<=len2; j++)
		{
			d=val[i-1][j-1]+score[i][j]; //diagonal
			
			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				
			
			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;
			
			
			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j
	
	//trace back to extract the alignment
	int i=len1;
	int j=len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;
			
			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;
			
			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}

void NWDP_TM(std::vector < numeric::xyzVector<core::Real> > const x, std::vector < numeric::xyzVector<core::Real> > const y, int const len1, int const len2, numeric::xyzVector <core::Real> const t, numeric::xyzMatrix <core::Real> const u, double d02, double gap_open, vector < int > & j2i)
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
	//Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	
	int i, j;
	double h, v, d;
	
	//initialization
	val[0][0]=0;
	for(i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}
	
	for(j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2]
	}
	numeric::xyzVector <core::Real> xx;
	double dij;
	
	
	//decide matrix and path
	for(i=1; i<=len1; i++)	
	{	
		transform(t, u, x[i-1], xx);
		for(j=1; j<=len2; j++)
		{
			//d=val[i-1][j-1]+score[i][j]; //diagonal
			dij=dist(xx, y[j-1]);    					
			d=val[i-1][j-1] +  1.0/(1+dij/d02);
			
			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				
			
			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;
			
			
			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j
	
	//trace back to extract the alignment
	i=len1;
	j=len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;
			
			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;
			
			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}

//+ss
void NWDP_TM( std::vector < int > const secx, std::vector < int > const secy, int const len1, int const len2, double gap_open, vector < int > & j2i)
{
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: secondary structure secx, secy, and gap_open
	//Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	
	int i, j;
	double h, v, d;
	
	//initialization
	val[0][0]=0;
	for(i=0; i<=len1; i++)
	{
		val[i][0]=0;
		path[i][0]=false; //not from diagonal
	}
	
	for(j=0; j<=len2; j++)
	{
		val[0][j]=0;
		path[0][j]=false; //not from diagonal
		j2i[j]=-1;	//all are not aligned, only use j2i[1:len2]
	}      
	
	//decide matrix and path
	for(i=1; i<=len1; i++)	
	{	
		for(j=1; j<=len2; j++)
		{
			//d=val[i-1][j-1]+score[i][j]; //diagonal			
			if(secx[i-1]==secy[j-1])
			{
				d=val[i-1][j-1] + 1.0;
			}
			else
			{
				d=val[i-1][j-1];
			}
			
			//symbol insertion in horizontal (= a gap in vertical)
			h=val[i-1][j];
			if(path[i-1][j]) //aligned in last position
				h += gap_open;				
			
			//symbol insertion in vertical
			v=val[i][j-1];
			if(path[i][j-1]) //aligned in last position
				v += gap_open;
			
			
			if(d>=h && d>=v)
			{
				path[i][j]=true; //from diagonal
				val[i][j]=d;
			}
			else 
			{
				path[i][j]=false; //from horizontal
				if(v>=h)
					val[i][j]=v;
				else					
					val[i][j]=h;
			}
		} //for i
	} //for j
	
	//trace back to extract the alignment
	i=len1;
	j=len2;
	while(i>0 && j>0)
	{
		if(path[i][j]) //from diagonal
		{
			j2i[j-1]=i-1;		
			i--;
			j--;
		}
		else 			
		{
			h=val[i-1][j];
			if(path[i-1][j]) h +=gap_open;
			
			v=val[i][j-1];
			if(path[i][j-1]) v +=gap_open;
			
			if(v>=h)
				j--;
			else
				i--;
		}
	}	
}

void
convert_xyz_to_vector(numeric::xyzVector <core::Real> const x,
					  vector <core::Real> & xx)
{
	xx.resize(3);
	xx[0]= x.x();
	xx[1]= x.y();
	xx[2]= x.z();
}
void
convert_xyz_to_matrix(numeric::xyzMatrix <core::Real> const x,
					  vector <vector <core::Real> > & xx)
{
	xx.resize(3);
	for (Size i = 0; i<3; ++i)
	{
		xx[i].resize(3);
		for (Size j = 0; j<3; ++j)
		{
			xx[i][j]= x(i+1,j+1);
		}
	}
}
void
convert_vector_to_xyz(vector <core::Real> const x,
						  numeric::xyzVector <core::Real> & xx)
{
	xx.x() = x[0];
	xx.y() = x[1];
	xx.z() = x[2];
}
void
convert_matrix_to_xyz(
					  vector <vector <core::Real> > const x,
					  numeric::xyzMatrix <core::Real> & xx)
{
	for (Size i = 0; i<3; ++i)
	{
		for (Size j = 0; j<3; ++j)
		{
			xx(i+1,j+1) = x[i][j];
		}
	}
}

void
convert_xyz_v_to_vectors(vector < numeric::xyzVector <core::Real> > const x,
				   vector < vector < double > > & xx)
{
	xx.resize(x.size());
	for (Size i=0;i<xx.size();++i) {
		convert_xyz_to_vector(x[i],xx[i]);
	}
}

// temporary fix before I have time to convert the Kabsch function -Y Song
bool Kabsch(vector < numeric::xyzVector <core::Real> > const x, 
			vector < numeric::xyzVector <core::Real> > const y, 
			int const n, 
			int const mode,
			double *rms,
			numeric::xyzVector <core::Real> & t,
			numeric::xyzMatrix <core::Real> & u
			)
{
	vector < vector < double > > xx;
	convert_xyz_v_to_vectors(x,xx);
	
	vector < vector < double > > yy;
	convert_xyz_v_to_vectors(y,yy);

	vector < double > tt;
	convert_xyz_to_vector(t,tt);
	
	vector < vector < double > > uu;
	convert_xyz_to_matrix(u,uu);
	
	bool retval =
	Kabsch(xx,yy,n,mode,rms,tt,uu);
	convert_vector_to_xyz(tt,t);
	convert_matrix_to_xyz(uu,u);
	
	return retval;
}
	
/**************************************************************************
 Implemetation of Kabsch algoritm for finding the best rotation matrix
 ---------------------------------------------------------------------------
 x    - x(i,m) are coordinates of atom m in set x            (input)
 y    - y(i,m) are coordinates of atom m in set y            (input)
 n    - n is number of atom pairs                            (input)
 mode  - 0:calculate rms only                                (input)
 1:calculate rms,u,t                                 (takes longer)
 rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
 u    - u(i,j) is   rotation  matrix for best superposition  (output)
 t    - t(i)   is translation vector for best superposition  (output)
 **************************************************************************/
bool Kabsch(vector < vector < double > > const x, 
			vector < vector < double > > const y, 
			int const n, 
			int const mode,
			double *rms,
			vector < double > & t,
			vector < vector < double > > & u
			)
{	
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;  
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3=1.73205080756888, tol=0.01;
	int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
	int ip2312[]={1, 2, 0, 1};
	
	int a_failed=0, b_failed=0;
	double epsilon=0.00000001;
	
	//initializtation
	*rms=0;
	rms1=0;
	e0=0;
	for(i=0; i<3; i++)
	{
		xc[i]=0.0;
		yc[i]=0.0;
		t[i]=0.0;
		for(j=0; j<3; j++)
		{
			u[i][j]=0.0;
			r[i][j]=0.0;
			a[i][j]=0.0;
			if(i==j)
			{
				u[i][j]=1.0;
				a[i][j]=1.0;
			}
		}
	}
	
	if(n<1)
	{
		return false;
	} 
	
	//compute centers for vector sets x, y
	for(i=0; i<n; i++)
	{
		xc[0] += x[i][0];
		xc[1] += x[i][1];
		xc[2] += x[i][2];
		
		yc[0] += y[i][0];
		yc[1] += y[i][1];
		yc[2] += y[i][2];
	}
	for(i=0; i<3; i++)
	{
		xc[i] = xc[i]/n;
		yc[i] = yc[i]/n;        
	}
	
	//compute e0 and matrix r
	for(m=0; m<n; m++)
	{
		for (i=0; i<3; i++)
		{
			e0 += (x[m][i]-xc[i])*(x[m][i]-xc[i])+\
			(y[m][i]-yc[i])*(y[m][i]-yc[i]);
			d = y[m][i] - yc[i];
			for(j=0; j<3; j++)
			{
				r[i][j] += d*(x[m][j] - xc[j]);
			}
		}        
	}
	//compute determinat of matrix r
	det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )\
	- r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )\
	+ r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] ); 
	sigma=det;
	
	//compute tras(r)*r
	m = 0;
	for(j=0; j<3; j++)
	{
		for (i=0; i<=j; i++)
		{            
			rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j];
			m++;
		}
	}
	
	double spur=(rr[0]+rr[2]+rr[5]) / 3.0;
	double cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5])\
					- rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
	det = det*det; 
	
	for (i=0; i<3; i++)
	{
		e[i]=spur;
	}
	
	if(spur>0) 
	{
		d = spur*spur;
		h = d - cof;
		g = (spur*cof - det)/2.0 - spur*h;
		
		if(h>0)
		{
			sqrth = sqrt(h);
			d = h*h*h - g*g;
			if(d<0.0) d=0.0;
			d = atan2( sqrt(d), -g ) / 3.0;			
			cth = sqrth * cos(d);
			sth = sqrth*sqrt3*sin(d);
			e[0]= (spur + cth) + cth;
			e[1]= (spur - cth) + sth;            
			e[2]= (spur - cth) - sth;
			
			if(mode!=0)
			{//compute a                
				for(l=0; l<3; l=l+2)
				{
					d = e[l];  
					ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
					ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
					ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
					ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
					ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];                
					ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1]; 
					
					if(fabs(ss[0])<=epsilon) ss[0]=0.0;
					if(fabs(ss[1])<=epsilon) ss[1]=0.0;
					if(fabs(ss[2])<=epsilon) ss[2]=0.0;
					if(fabs(ss[3])<=epsilon) ss[3]=0.0;
					if(fabs(ss[4])<=epsilon) ss[4]=0.0;
					if(fabs(ss[5])<=epsilon) ss[5]=0.0;
					
					if( fabs(ss[0]) >= fabs(ss[2]) )
					{
						j=0;                    
						if( fabs(ss[0]) < fabs(ss[5])) 
						{
							j = 2;
						}
					}
					else if ( fabs(ss[2]) >= fabs(ss[5]) )
					{
						j = 1;
					}
					else
					{
						j = 2;
					}  
					
					d = 0.0;
					j = 3 * j;
					for(i=0; i<3; i++)
					{
						k=ip[i+j];
						a[i][l] = ss[k];
						d = d + ss[k]*ss[k];						
					} 
					
					
					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if( d > epsilon ) d = 1.0 / sqrt(d);
					else d=0.0;
					for(i=0; i<3; i++)
					{
						a[i][l] = a[i][l] * d;
					}               
				}//for l
				
				d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2]))
				{
					m1=2;
					m=0;
				}
				else
				{
					m1=0;
					m=2;                
				}
				p=0;
				for(i=0; i<3; i++)
				{
					a[i][m1] = a[i][m1] - d*a[i][m];
					p = p + a[i][m1]*a[i][m1];
				}
				if( p <= tol )
				{
					p = 1.0;
					for(i=0; i<3; i++)
					{
						if (p < fabs(a[i][m])) 
						{
							continue;
						}
						p = fabs( a[i][m] );
						j = i;                    
					}
					k = ip2312[j];
					l = ip2312[j+1];
					p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] ); 
					if( p > tol )  
					{
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m]/p;
						a[l][m1] =  a[k][m]/p;                                                       
					}   
					else
					{//goto 40
						a_failed=1;
					}     
				}//if p<=tol
				else
				{
					p = 1.0 / sqrt(p);
					for(i=0; i<3; i++)
					{
						a[i][m1] = a[i][m1]*p;
					}                                  
				}//else p<=tol  
				if(a_failed!=1)
				{
					a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
					a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
					a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];       
				}                                   
			}//if(mode!=0)       
		}//h>0
		
		//compute b anyway
		if(mode!=0 && a_failed!=1)//a is computed correctly
		{
			//compute b
			for(l=0; l<2; l++)
			{
				d=0.0;
				for(i=0; i<3; i++)
				{
					b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
					d = d + b[i][l]*b[i][l];
				}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if( d > epsilon ) d = 1.0 / sqrt(d);
				else d=0.0;
				for(i=0; i<3; i++)
				{
					b[i][l] = b[i][l]*d;
				}                
			}            
			d =b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
			p=0.0;
			
			for(i=0; i<3; i++)
			{
				b[i][1] = b[i][1] - d*b[i][0];
				p += b[i][1]*b[i][1];
			}
			
			if( p <= tol )
			{
				p = 1.0;
				for(i=0; i<3; i++)
				{
					if(p<fabs(b[i][0]))
					{
						continue;
					}
					p = fabs( b[i][0] );
					j=i;
				}
				k = ip2312[j];
				l = ip2312[j+1];
				p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] ); 
				if( p > tol )  
				{
					b[j][1] = 0.0;
					b[k][1] = -b[l][0]/p;
					b[l][1] =  b[k][0]/p;        
				}
				else
				{
					//goto 40
					b_failed=1;
				}                
			}//if( p <= tol )
			else
			{
				p = 1.0 / sqrt(p);
				for(i=0; i<3; i++)
				{
					b[i][1]=b[i][1]*p;
				}
			}            
			if(b_failed!=1)
			{
				b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
				b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
				b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0]; 
				//compute u
				for(i=0; i<3; i++)
				{
					for(j=0; j<3; j++)
					{
						u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1]\
						+ b[i][2]*a[j][2];
					}
				}
			}
			
			//compute t
			for(i=0; i<3; i++)
			{
				t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1])\
				- u[i][2]*xc[2];
			}            
		}//if(mode!=0 && a_failed!=1)
	}//spur>0
	else //just compute t and errors
	{
		//compute t
		for(i=0; i<3; i++)
		{
			t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
		}
	}//else spur>0 
	
	//compute rms
	for(i=0; i<3; i++)
	{
		if( e[i] < 0 ) e[i] = 0;
		e[i] = sqrt( e[i] );           
	}            
	d = e[2];
	if( sigma < 0.0 )
	{
		d = - d;
	}
	d = (d + e[1]) + e[0];
	rms1 = (e0 - d) - d; 
	if( rms1 < 0.0 ) rms1 = 0.0;  
	
	*rms=rms1;
	
	return true;
	
}

/*
void load_PDB_allocate_memory(char *xname, char *yname)
{    
    //------get length first------>
	xlen=get_PDB_len(xname);
	ylen=get_PDB_len(yname);
	minlen=min(xlen, ylen);

    //------allocate memory for x and y------>
	NewArray(&xa, xlen, 3);
    seqx   = new char[xlen+1];
	secx   = new int[xlen];
	xresno = new int[xlen];

	NewArray(&ya, ylen, 3);
	seqy    = new char[ylen+1];
	yresno  = new int[ylen];
	secy    = new int[ylen];

    
    
    //------load data------>  
    read_PDB(xname, xa, seqx, xresno);
    read_PDB(yname, ya, seqy, yresno);
    
    
    //------allocate memory for other temporary varialbes------>
 	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);
	NewArray(&xtm, minlen, 3);
	NewArray(&ytm, minlen, 3);
	NewArray(&xt, xlen, 3);

	NewArray(&score, xlen+1, ylen+1);
	NewArray(&path, xlen+1, ylen+1);
	NewArray(&val, xlen+1, ylen+1);  	
}
*/
void load_pose_allocate_memory(core::pose::Pose pose1, core::pose::Pose pose2, std::list <core::Size> residue_list1, std::list <core::Size> residue_list2)
{    
	//------get length first------>
	xlen=residue_list1.size();
	ylen=residue_list2.size();
	minlen=min(xlen, ylen);
	
	//------allocate memory for x and y------>
	xa.resize(xlen);
	seqx.resize(xlen);
	secx.resize(xlen);
	xresno.resize(xlen);
	
	ya.resize(ylen);
	seqy.resize(ylen);
	yresno.resize(ylen);
	secy.resize(ylen);
	
	
	
	//------load data------>  
	read_pose(pose1, residue_list1, xa, seqx, xresno);
	read_pose(pose2, residue_list2, ya, seqy, yresno);
	
	
	//------allocate memory for other temporary varialbes------>
	r1.resize(minlen);
	r2.resize(minlen);
	xtm.resize(minlen);
	ytm.resize(minlen);
	xt.resize(xlen);
	//NewArray(&r1, minlen, 3);
	//NewArray(&r2, minlen, 3);
	//NewArray(&xtm, minlen, 3);
	//NewArray(&ytm, minlen, 3);
	//NewArray(&xt, xlen, 3);
	
	ResizeArray(score, xlen+1, ylen+1);
	ResizeArray(path, xlen+1, ylen+1);
	ResizeArray(val, xlen+1, ylen+1);  	
}


	/*
void free_memory()
{
	//DeleteArray(&path, xlen+1);
	//DeleteArray(&val, xlen+1);
	//DeleteArray(&score, xlen+1);
	DeleteArray(&xa, xlen);
	DeleteArray(&xt, xlen);
	DeleteArray(&ya, ylen);
	DeleteArray(&r1, minlen);
	DeleteArray(&r2, minlen);
	DeleteArray(&xtm, minlen);
	DeleteArray(&ytm, minlen);
    
    
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    delete [] xresno;
    delete [] yresno;
    
}
*/

//     1, collect those residues with dis<d;
//     2, calculate TMscore
int score_fun8( vector < numeric::xyzVector <core::Real> > const xa, 
                vector < numeric::xyzVector <core::Real> > const ya, 
                int const n_ali,
                double const d,
                vector <int> & i_ali, 
                double *score1,
                int score_sum_method
              )
{
    double score_sum=0, di;
    double d_tmp=d*d;
	double d02=d0*d0;
	double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

    while(1)
    {
		n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {				
                if(di<=score_d8_cut)
                {					
                    score_sum += 1/(1+di/d02);
                }                
            }
            else
            {
				score_sum += 1/(1+di/d02);
            }
        }
        //there are not enough feasible pairs, reliefe the threshold 		
        if(n_cut<3 && n_ali>3)
        {
			inc++;
			double dinc=(d+inc*0.5);
			d_tmp = dinc * dinc;
        }
        else
        {
            break;
        }

    }  

    *score1=score_sum/Lnorm;

    return n_cut;
}

// TMscore search engine
// input:   two aligned vector sets: x, y
//          scale parameter d0
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
//                                  
//          
// output:  the best rotaion matrix t0, u0 that results in highest TMscore
double TMscore8_search( vector < numeric::xyzVector <core::Real> > const  xtm, 
                        vector < numeric::xyzVector <core::Real> > const  ytm,
                        int Lali, 
                        numeric::xyzVector <core::Real> & t0,
                        numeric::xyzMatrix <core::Real> & u0,
                        int const simplify_step,
                        int const score_sum_method,
                        double *Rcomm
                       )
{   
    int i, m;
    double score_max, score, rmsd;    
    const int kmax=Lali;    
    vector <int> k_ali(kmax);
	int ka, k;
    numeric::xyzVector<core::Real> t;
    numeric::xyzMatrix<core::Real> u;
	double d;
	

	//iterative parameters
	int n_it=20;            //maximum number of iterations
    const int n_init_max=6; //maximum number of different fragment length 
    vector <int> L_ini(n_init_max);  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
    int L_ini_min=4;
    if(Lali<4) L_ini_min=Lali;   
    int n_init=0, i_init;      
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=(int) (Lali/pow(2.0, (double) i));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }
    
    score_max=-1;
    //find the maximum score starting from local structures superposition
    vector <int> i_ali(kmax);
	int n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment
    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;
      
        i=0;   
        while(1)
        {
            //extract the fragment starting from position i 
            ka=0;
            for(k=0; k<L_frag; k++)
            {
				int kk=k+i;
                r1[k]=xtm[kk]; 
                
                r2[k]=ytm[kk];  
                
                k_ali[ka]=kk;
                ka++;
            }
            
            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
            if(i_init==0)
            {
                *Rcomm=sqrt(rmsd/Lali);
            }
            do_rotation(xtm, xt, Lali, t, u);
            
            //get subsegment of this fragment
            d=d0_search-1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
            if(score>score_max)
            {
                score_max=score;
                
                //save the rotation matrix
				t0 = t;
				u0 = u;
				/*
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
				 */
            }
            
            //try to extend the alignment iteratively            
            d=d0_search+1;
            for(int it=0; it<n_it; it++)            
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k]=xtm[m];  
                    
                    r2[k]=ytm[m];  
                    
                    k_ali[ka]=m;
                    ka++;
                } 
                //extract rotation matrix based on the fragment                
                Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
					t0 = t;
					u0 = u;
					/*
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }
					 */
                }
                
                //check if it converges                 
			
                if(n_cut==ka)
                {				
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k])
						{
							break;
						}
                    }
                    if(k==n_cut)
                    {						
                        break; //stop iteration
                    }
                }                                                               
            } //for iteration            

			if(i<iL_max)
			{
				i=i+simplify_step; //shift the fragment		
				if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
			}
			else if(i>=iL_max)
			{
				break;
			}
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8          
// output:  the best rotaion matrix t, u that results in highest TMscore
double detailed_search( vector < numeric::xyzVector <core::Real> > const x,
                        vector < numeric::xyzVector <core::Real> > const y, 
                        int const, 
                        int const y_len, 
                        vector < int > const invmap0,
                        numeric::xyzVector <core::Real> & t,
                        numeric::xyzMatrix <core::Real> & u,
                        int simplify_step,
                        int score_sum_method                        
                       )
{
    //x is model, y is template, try to superpose onto y
    int i, j, k;     
    double tmscore;
    double rmsd;



    k=0;
    for(i=0; i<y_len; i++) 
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
			xtm[k] = x[j];
			ytm[k] = y[i];
            k++;
        }
    }
    
    //detailed search 40-->1
    tmscore=TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd);  
   
    return tmscore;
}



//compute the score quickly in three iterations
double get_score_fast(vector < numeric::xyzVector <core::Real> > const x, vector < numeric::xyzVector <core::Real> > const y, int const , int const y_len, vector < int > const invmap)
{
    double rms, tmscore, tmscore1, tmscore2;
    int i, j, k;

    k=0;
    for(j=0; j<y_len; j++)
    {
        i=invmap[j];
        if(i>=0)
        {
            r1[k]=x[i];
            r2[k]=y[j];
            xtm[k]=x[i];
            ytm[k]=y[j];
            
            k++;
        }
        else if(i!=-1)
        {
            PrintErrorAndQuit("Wrong map!\n");
        }       
    }
    Kabsch(r1, r2, k, 1, &rms, t, u);
    
    //evaluate score   
    double di;
	const int len=k;
    vector <double> dis(len);    
	double d00=d0_search;
	double d002=d00*d00;
	double d02=d0*d0;
	
    int n_ali=k;
	numeric::xyzVector <core::Real> xrot;
	tmscore=0;
	for(k=0; k<n_ali; k++)
	{
        transform(t, u, xtm[k], xrot);        
        di=dist(xrot, ytm[k]);
        dis[k]=di;
        tmscore += 1/(1+di/d02);
    }
	
   
   
   //second iteration 
    double d002t=d002;
    while(1)
    {
		j=0;
        for(k=0; k<n_ali; k++)
        {            
            if(dis[k]<=d002t)
            {
                r1[j]=xtm[k];
                r2[j]=ytm[k];
                
                j++;
            }
        }
        //there are not enough feasible pairs, relieve the threshold 
        if(j<3 && n_ali>3)
        {
            d002t += 0.5;
        }
        else
        {
            break;
        }
    }
    
    if(n_ali!=j)
    {
        Kabsch(r1, r2, j, 1, &rms, t, u);
    	tmscore1=0;
    	for(k=0; k<n_ali; k++)
    	{
            transform(t, u, xtm[k], xrot);        
            di=dist(xrot, ytm[k]);
            dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }
        
        //third iteration
        d002t=d002+1;
       
        while(1)
        {
			j=0;
            for(k=0; k<n_ali; k++)
            {            
                if(dis[k]<=d002t)
                {
                    r1[j]=xtm[k];
                    r2[j]=ytm[k];
                                        
                    j++;
                }
            }
            //there are not enough feasible pairs, relieve the threshold 
            if(j<3 && n_ali>3)
            {
                d002t += 0.5;
            }
            else
            {
                break;
            }
        }

        //evaluate the score
        Kabsch(r1, r2, j, 1, &rms, t, u);
        tmscore2=0;
        for(k=0; k<n_ali; k++)
        {
            transform(t, u, xtm[k], xrot);
            di=dist(xrot, ytm[k]);
            tmscore2 += 1/(1+di/d02);
        }    
    }
    else
    {
        tmscore1=tmscore;
        tmscore2=tmscore;
    }
    
      
    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;

    
    return tmscore; // no need to normalize this score because it will not be used for latter scoring
}


//perform gapless threading to find the best initial alignment
//input: x, y, x_len, y_len
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial( vector < numeric::xyzVector <core::Real> > const x, 
                    vector < numeric::xyzVector <core::Real> > const y, 
                    int const x_len,
                    int const y_len, 
                    vector < int > & y2x
                   )
{
    int min_len=getmin(x_len, y_len);
    if(min_len<=5) PrintErrorAndQuit("Sequence is too short <=5!\n");
    
    int min_ali= min_len/2;              //minimum size of considered fragment 
    if(min_ali<=5)  min_ali=5;    
    int n1, n2;
    n1 = -y_len+min_ali; 
    n2 = x_len-min_ali;

    int i, j, k, k_best;
    double tmscore, tmscore_max=-1;

    k_best=n1;
    for(k=n1; k<=n2; k++)
    {
        //get the map
        for(j=0; j<y_len; j++)
        {
            i=j+k;
            if(i>=0 && i<x_len)
            {
                y2x[j]=i;
            }
            else
            {
                y2x[j]=-1;
            }
        }
        
        //evaluate the map quickly in three iterations
		//this is not real tmscore, it is used to evaluate the goodness of the initial alignment
        tmscore=get_score_fast(x, y, x_len, y_len, y2x); 
        if(tmscore>=tmscore_max)
        {
            tmscore_max=tmscore;
            k_best=k;
        }
    }
    
    //extract the best map
    k=k_best;
    for(j=0; j<y_len; j++)
    {
        i=j+k;
        if(i>=0 && i<x_len)
        {
            y2x[j]=i;
        }
        else
        {
            y2x[j]=-1;
        }
    }    

    return tmscore_max;
}

void smooth(vector < int > & sec, int const len)
{
	int i, j;
	//smooth single  --x-- => -----
	for(i=2; i<len-2; i++)
	{
		if(sec[i]==2 || sec[i]==4)
		{
			j=sec[i];
			if(sec[i-2] != j)
			{
				if(sec[i-1] != j)
				{
					if(sec[i+1] != j)
					{
						if(sec[i+2] != j) 
						{
							sec[i]=1;
						}
					}
				}
			}
		}
	}

	//   smooth double 
	//   --xx-- => ------

	for(i=0; i<len-5; i++)
	{
		//helix
		if(sec[i] != 2)
		{
			if(sec[i+1] != 2)
			{
				if(sec[i+2] == 2)
				{
					if(sec[i+3] == 2)
					{
						if(sec[i+4] != 2)
						{
							if(sec[i+5] != 2)
							{
								sec[i+2]=1;
								sec[i+3]=1;
							}
						}
					}
				}
			}
		}

		//beta
		if(sec[i] != 4)
		{
			if(sec[i+1] != 4)
			{
				if(sec[i+2] ==4)
				{
					if(sec[i+3] == 4)
					{
						if(sec[i+4] != 4)
						{
							if(sec[i+5] != 4)
							{
								sec[i+2]=1;
								sec[i+3]=1;
							}
						}
					}
				}
			}
		}
	}

	//smooth connect
	for(i=0; i<len-2; i++)
	{		
		if(sec[i] == 2)
		{
			if(sec[i+1] != 2)
			{
				if(sec[i+2] == 2)
				{
					sec[i+1]=2;
				}
			}
		}
		else if(sec[i] == 4)
		{
			if(sec[i+1] != 4)
			{
				if(sec[i+2] == 4)
				{
					sec[i+1]=4;
				}
			}
		}
	}

}

int sec_str(double dis13, double dis14, double dis15, double dis24, double dis25, double dis35)
{
	int s=1;
	
	double delta=2.1;
	if(fabs(dis15-6.37)<delta)
	{
		if(fabs(dis14-5.18)<delta)
		{
			if(fabs(dis25-5.18)<delta)
			{
				if(fabs(dis13-5.45)<delta)
				{
					if(fabs(dis24-5.45)<delta)
					{
						if(fabs(dis35-5.45)<delta)
						{
							s=2; //helix						
							return s;
						}
					}
				}
			}
		}
	}

	delta=1.42;
	if(fabs(dis15-13)<delta)
	{
		if(fabs(dis14-10.4)<delta)
		{
			if(fabs(dis25-10.4)<delta)
			{
				if(fabs(dis13-6.1)<delta)
				{
					if(fabs(dis24-6.1)<delta)
					{
						if(fabs(dis35-6.1)<delta)
						{
							s=4; //strand
							return s;
						}
					}
				}
			}
		}
	}

	if(dis15 < 8)
	{
		s=3; //turn
	}	  


	return s;
}


//1->coil, 2->helix, 3->turn, 4->strand
void make_sec(std::vector < numeric::xyzVector <core::Real> > const x, int const len, vector < int > & sec)
{
    int j1, j2, j3, j4, j5;
	double d13, d14, d15, d24, d25, d35;
    for(int i=0; i<len; i++)
    { 	
		sec[i]=1;
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;		
        
        if(j1>=0 && j5<len)
		{
			d13=sqrt(dist(x[j1], x[j3]));
			d14=sqrt(dist(x[j1], x[j4]));
			d15=sqrt(dist(x[j1], x[j5]));
			d24=sqrt(dist(x[j2], x[j4]));
			d25=sqrt(dist(x[j2], x[j5]));
			d35=sqrt(dist(x[j3], x[j5]));
			sec[i]=sec_str(d13, d14, d15, d24, d25, d35);			
		}    
    } 
	smooth(sec, len);
}




//get initial alignment from secondary structure alignment
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ss(  std::vector < numeric::xyzVector <core::Real> > const x, 
					  std::vector < numeric::xyzVector <core::Real> > const y, 
					  int const x_len,
					  int const y_len, 
						std::vector < int > & y2x
                      )
{
    //assign secondary structures
	make_sec(x, x_len, secx);
	make_sec(y, y_len, secy);

	double gap_open=-1.0;
	NWDP_TM(secx, secy, x_len, y_len, gap_open, y2x);    
}


// get_initial5 in TMalign
//get initial alignment of local structure superposition
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
bool get_initial_local(  std::vector < numeric::xyzVector <core::Real> > const x, 
						 std::vector < numeric::xyzVector <core::Real> > const y, 
						 int const x_len,
						 int const y_len, 
					   std::vector < int > & y2x
						 )
{
    double GL, rmsd;    
    numeric::xyzVector <core::Real> t;
    numeric::xyzMatrix <core::Real> u;

	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;
	
	double GLmax=0;
	int n_frag=20; //length of fragment for superposition
	int ns=20; //tail length to discard
	vector < int > invmap(y_len+1);
      
	int aL=getmin(x_len, y_len);
	if(aL>250)
	{
		n_frag=50;
	}
	else if(aL>200)
	{
		n_frag=40;
	}
	else if(aL>150)
	{
		n_frag=30;
	}
	else
	{
		n_frag=20;
	}
	
	int smallest=aL/3; // I change here from aL/2 to aL/3

	if(n_frag>smallest) n_frag=smallest; 
	if(ns>smallest) ns=smallest;
    
	int m1=x_len-n_frag-ns;
	int m2=y_len-n_frag-ns;

	bool flag=false;

	for(int ii=0; ii<y_len; ii++) 
	{                
		y2x[ii]=-1;
	}

	int count=0;
	for(int i=ns-1; i<m1; i=i+n_frag) //index starts from 0, different from FORTRAN
	{
		for(int j=ns-1; j<m2; j=j+n_frag)
		{
			for(int k=0; k<n_frag; k++) //fragment in y
			{ 
				r1[k]=x[k+i];  
                r2[k]=y[k+j];  
			}


			Kabsch(r1, r2, n_frag, 1, &rmsd, t, u);			
			count++;

			double gap_open=0.0;			
			NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open, invmap);
			GL=get_score_fast(x, y, x_len, y_len, invmap);
			if(GL>GLmax)
			{
				GLmax=GL;
				for(int ii=0; ii<y_len; ii++) 
				{                
					y2x[ii]=invmap[ii];
				}
				flag=true;
			}
		}
	}	


	//delete [] invmap;
	return flag;

}



//with invmap(i) calculate score(i,j) using RMSD rotation
void score_matrix_rmsd(  std::vector < numeric::xyzVector <core::Real> > const x, 
						 std::vector < numeric::xyzVector <core::Real> > const y, 
						 int const x_len,
						 int const y_len,
					   std::vector < int > const y2x
						 )
{
	numeric::xyzVector <core::Real> t;
	numeric::xyzMatrix <core::Real> u;
	double rmsd, dij;
	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;

	numeric::xyzVector <core::Real> xx;
	int i, k=0;
	for(int j=0; j<y_len; j++)
	{
		i=y2x[j];
		if(i>=0)
		{
			r1[k]=x[i];  
			r2[k]=y[j];  
			
			k++;
		}
	}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);
	//do_rotation(x, xt, x_len, t, u);
	
	
	for(int ii=0; ii<x_len; ii++)
	{		
		transform(t, u, x[ii], xx);
		for(int jj=0; jj<y_len; jj++)
		{
			//dij=dist(&xt[ii][0], &y[jj][0]);   
			dij=dist(xx, y[jj]); 
			score[ii+1][jj+1] = 1.0/(1+dij/d02);
			//	cout << ii+1 << " " << jj+1 << " " << score[ii+1][jj+1]<< endl;
		}
	}		
}


void score_matrix_rmsd_sec(  vector < numeric::xyzVector <core::Real> > const x, 
							 vector < numeric::xyzVector <core::Real> > const  y, 
							 int const x_len,
							 int const y_len,
							 vector < int > const y2x
							 )
{
	numeric::xyzVector <core::Real> t;
	numeric::xyzMatrix <core::Real> u;
	double rmsd, dij;
	double d01=d0+1.5;
	if(d01 < D0_MIN) d01=D0_MIN;
	double d02=d01*d01;

	numeric::xyzVector <core::Real> xx;
	int i, k=0;
	for(int j=0; j<y_len; j++)
	{
		i=y2x[j];
		if(i>=0)
		{
			r1[k]=x[i];  
			r2[k]=y[j];  
			
			k++;
		}
	}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);

	
	for(int ii=0; ii<x_len; ii++)
	{		
		transform(t, u, x[ii], xx);
		for(int jj=0; jj<y_len; jj++)
		{
			dij=dist(xx, y[jj]); 
			if(secx[ii]==secy[jj])
			{
				score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
			}
			else
			{
				score[ii+1][jj+1] = 1.0/(1+dij/d02);
			}		
		}
	}		
}


//get initial alignment from secondary structure and previous alignments
//input: x, y, x_len, y_len
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
void get_initial_ssplus( std::vector < numeric::xyzVector <core::Real> > const x, 
						 std::vector < numeric::xyzVector <core::Real> > const y, 
						 int const x_len,
						 int const y_len,
						 std::vector < int > & y2x0,
						std::vector < int > & y2x
						 )
{

	//create score matrix for DP
	score_matrix_rmsd_sec(x, y, x_len, y_len, y2x0);
	
	double gap_open=-1.0;
	NWDP_TM(x_len, y_len, gap_open, y2x);
}


void find_max_frag(std::vector < numeric::xyzVector <core::Real> > const x, std::vector < int > const resno, int const len, int *start_max, int *end_max)
{
	int r_min, fra_min=4;           //minimum fragment for search
	double d;
	int start;
	int Lfr_max=0, flag;

	r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
	if(r_min > fra_min) r_min=fra_min;
	
	int inc=0;
	double dcu0_cut=dcu0*dcu0;;
	double dcu_cut=dcu0_cut;

	while(Lfr_max < r_min)
	{		
		Lfr_max=0;			
		int j=1;    //number of residues at nf-fragment
		start=0;
		for(int i=1; i<len; i++)
		{			
			d = dist(x[i-1], x[i]);
			flag=0;
			if(dcu_cut>dcu0_cut)
			{
				if(d<dcu_cut)
				{
					flag=1;
				}
			}
			else if(resno[i] == (resno[i-1]+1)) //necessary??
			{
				if(d<dcu_cut)
				{
					flag=1;
				}
			}

			if(flag==1)
			{
				j++;

				if(i==(len-1))
				{
					if(j > Lfr_max) 
					{
						Lfr_max=j;
						*start_max=start;
						*end_max=i;						
					}
					j=1;
				}
			}
			else
			{
				if(j>Lfr_max) 
				{
					Lfr_max=j;
					*start_max=start;
					*end_max=i-1;										
				}

				j=1;
				start=i;
			}
		}// for i;
		
		if(Lfr_max < r_min)
		{
			inc++;
			double dinc=pow(1.1, (double) inc) * dcu0;
			dcu_cut= dinc*dinc;
		}
	}//while <;	
}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, x_len, y_len
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double get_initial_fgt( std::vector < numeric::xyzVector < core::Real > > const x, 
						std::vector < numeric::xyzVector < core::Real > > const y, 
						int const x_len,
						int const y_len, 
						vector < int > const xresno,
						vector < int > const yresno,
						vector < int > & y2x
						)
{
	int fra_min=4;           //minimum fragment for search
	int fra_min1=fra_min-1;  //cutoff for shift, save time

	int xstart=0, ystart=0, xend=0, yend=0;

	find_max_frag(x, xresno, x_len,  &xstart, &xend);
	find_max_frag(y, yresno, y_len, &ystart, &yend);


	int Lx = xend-xstart+1;
	int Ly = yend-ystart+1;
	vector < int > ifr, y2x_;
	int L_fr=getmin(Lx, Ly);
	ifr.resize(L_fr);
	y2x_.resize(y_len+1);

	//select what piece will be used (this may araise ansysmetry, but
	//only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
	//if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

	if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	{		
		for(int i=0; i<L_fr; i++)
		{
			ifr[i]=xstart+i;
		}
	}
	else if(Lx>Ly || (Lx==Ly && x_len>y_len))
	{		
		for(int i=0; i<L_fr; i++)
		{
			ifr[i]=ystart+i;
		}	
	}

	
	int L0=getmin(x_len, y_len); //non-redundant to get_initial1
	if(L_fr==L0)
	{
		int n1= (int)(L0*0.1); //my index starts from 0
		int n2= (int)(L0*0.89);

		int j=0;
		for(int i=n1; i<= n2; i++)
		{
			ifr[j]=ifr[i];
			j++;
		}
		L_fr=j;
	}


	//gapless threading for the extracted fragment
	double tmscore, tmscore_max=-1;

	if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	{
		int L1=L_fr;
	    int min_len=getmin(L1, y_len);    
		int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		if(min_ali<=fra_min1)  min_ali=fra_min1;    
		int n1, n2;
		n1 = -y_len+min_ali; 
		n2 = L1-min_ali;

		int i, j, k;
		for(k=n1; k<=n2; k++)
		{
			//get the map
			for(j=0; j<y_len; j++)
			{
				i=j+k;
				if(i>=0 && i<L1)
				{				
					y2x_[j]=ifr[i];
				}
				else
				{
					y2x_[j]=-1;
				}
			}

			//evaluate the map quickly in three iterations
			tmscore=get_score_fast(x, y, x_len, y_len, y2x_);

			if(tmscore>=tmscore_max)
			{
				tmscore_max=tmscore;
				for(j=0; j<y_len; j++)
				{
					y2x[j]=y2x_[j];
				}
			}
		}
	}
	else
	{
		int L2=L_fr;
	    int min_len=getmin(x_len, L2);    
		int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		if(min_ali<=fra_min1)  min_ali=fra_min1;    
		int n1, n2;
		n1 = -L2+min_ali; 
		n2 = x_len-min_ali;

		int i, j, k;	

		for(k=n1; k<=n2; k++)
		{
			//get the map
			for(j=0; j<y_len; j++)
			{
				y2x_[j]=-1;
			}

			for(j=0; j<L2; j++)
			{
				i=j+k;
				if(i>=0 && i<x_len)
				{
					y2x_[ifr[j]]=i;
				}
			}
        
			//evaluate the map quickly in three iterations
			tmscore=get_score_fast(x, y, x_len, y_len, y2x_);
			if(tmscore>=tmscore_max)
			{
				tmscore_max=tmscore;
				for(j=0; j<y_len; j++)
				{
					y2x[j]=y2x_[j];
				}
			}
		}
	}    

	//delete [] ifr;
	//delete [] y2x_;
    return tmscore_max;
}





//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double DP_iter( vector < numeric::xyzVector < core::Real > > const x,
                vector < numeric::xyzVector < core::Real > > const y, 
                int const x_len, 
                int const y_len, 
                numeric::xyzVector < core::Real > t,
                numeric::xyzMatrix < core::Real > u,
				vector < int > & invmap0,
				int const g1,
				int const g2,
				int const iteration_max                                   
				)
{
    double gap_open[2]={-0.6, 0};
    double rmsd; 
    vector < int > invmap(y_len+1);
    
    int iteration, i, j, k;
    double tmscore, tmscore_max, tmscore_old=0;    
    int score_sum_method=8, simplify_step=40;
    tmscore_max=-1;

	//double d01=d0+1.5;
    double d02=d0*d0;
    for(int g=g1; g<g2; g++)
    {
        for(iteration=0; iteration<iteration_max; iteration++)
        {           
			NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open[g], invmap);
            
            k=0;
            for(j=0; j<y_len; j++) 
            {
                i=invmap[j];

                if(i>=0) //aligned
                {
                    xtm[k]=x[i];
                    ytm[k]=y[j];

                    k++;
                }
            }

            tmscore=TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd);

           
            if(tmscore>tmscore_max)
            {
                tmscore_max=tmscore;
                for(i=0; i<y_len; i++) 
                {                
                    invmap0[i]=invmap[i];                                      
                }				
            }
    
            if(iteration>0)
            {
                if(fabs(tmscore_old-tmscore)<0.000001)
                {     
                    break;       
                }
            }
            tmscore_old=tmscore;
        }// for iteration           
        
    }//for gapopen
    
    
    //delete []invmap;
    return tmscore_max;
}

/*
void output_superpose(char *xname,
					  char *yname,
					  int x_len,
					  int y_len,
					  numeric::xyzVector < core::Real > t,
					  numeric::xyzMatrix < core::Real > u,
					  double rmsd,
					  double d0_out,
					  int m1[], 
					  int m2[],
					  int n_ali8,
					  double seq_id,
					  double TM_0,
					  double Lnorm_0,
					  double d0_0
					 )
{
	int i, j, j1;
	double dis2;

	int max=2000;


	//aligned region
	FILE *fp=fopen(out_reg, "w");
	fprintf(fp, "load inline\n");
	fprintf(fp, "select atomno<%d\n", max);
	fprintf(fp, "wireframe .45\n");
	fprintf(fp, "select none\n");
	fprintf(fp, "select atomno>%d\n", max);
	fprintf(fp, "wireframe .20\n");
	fprintf(fp, "color white\n");


	do_rotation(xa, xt, x_len, t, u);
	for(i=0; i<n_ali8; i++)
	{
		j=m1[i];
		j1=m2[i];
		dis2=sqrt(dist(xt[j], ya[j1])) ;
		if(dis2<=d0_out)
		{
			fprintf(fp, "select atomno=%d\n", j+1);
			fprintf(fp, "color red\n");

			fprintf(fp, "select atomno=%d\n", max+j1+1);
			fprintf(fp, "color red\n");
		}
	}

	fprintf(fp, "select all\n");
	fprintf(fp, "exit\n");
	fprintf(fp, "REMARK TM-align Version %s\n", version);
	fprintf(fp, "REMARK Structure A:%s   Size=%d\n", yname, y_len);
	fprintf(fp, "REMARK Structure B:%s   Size=%d\n", xname, x_len);
	fprintf(fp, "REMARK TM-score is normalized by %.2f, d0=%.2f\n",  Lnorm_0, d0_0);
	fprintf(fp, "REMARK Aligned length=%d, RMSD=%.2f, TM-score=%.5f, ID=%.3f\n", n_ali8, rmsd, TM_0, seq_id);


	char AA[3];
	//superposed structure B
	for(i=0; i<n_ali8; i++)
	{
		j=m1[i];
		AAmap3(seqx[j], AA);
		fprintf(fp, "ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f\n", j+1, AA, xresno[j], xt[j][0], xt[j][1], xt[j][2]);
	}
	fprintf(fp, "TER\n");

	for(i=1; i<n_ali8; i++)
	{
		j=m1[i-1]+1;
		j1=m1[i]+1;
		fprintf(fp, "CONECT%5d%5d\n", j, j1);
	}
	//structure A
	for(i=0; i<n_ali8; i++)
	{
		j=m2[i];
		AAmap3(seqy[j], AA);
		fprintf(fp, "ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f\n", j+max+1, AA, yresno[j], ya[j][0], ya[j][1], ya[j][2]);
	}
	fprintf(fp, "TER\n");
	for(i=1; i<n_ali8; i++)
	{
		j=max+m2[i-1]+1;
		j1=max+m2[i]+1;
		fprintf(fp, "CONECT%5d%5d\n", j, j1);
	}

	fclose(fp);








	//all regions
	char str[3000];
	sprintf(str, "%s_all", out_reg);
	fp=fopen(str, "w");
	fprintf(fp, "load inline\n");
	fprintf(fp, "select atomno<%d\n", max);
	fprintf(fp, "wireframe .45\n");
	fprintf(fp, "select none\n");
	fprintf(fp, "select atomno>%d\n", max);
	fprintf(fp, "wireframe .20\n");
	fprintf(fp, "color white\n");

	for(i=0; i<n_ali8; i++)
	{
		j=m1[i];
		j1=m2[i];
		dis2=sqrt(dist(xt[j], ya[j1])) ;
		if(dis2<=d0_out)
		{
			fprintf(fp, "select atomno=%d\n", j+1);
			fprintf(fp, "color red\n");

			fprintf(fp, "select atomno=%d\n", max+j1+1);
			fprintf(fp, "color red\n");
		}
	}


	fprintf(fp, "select all\n");
	fprintf(fp, "exit\n");
	fprintf(fp, "REMARK TM-align Version %s\n", version);
	fprintf(fp, "REMARK Structure A:%s   Size=%d\n", yname, y_len);
	fprintf(fp, "REMARK Structure B:%s   Size=%d\n", xname, x_len);
	fprintf(fp, "REMARK TM-score is normalized by %.2f, d0=%.2f\n",  Lnorm_0, d0_0);
	fprintf(fp, "REMARK Aligned length=%d, RMSD=%.2f, TM-score=%.5f, ID=%.3f\n", n_ali8, rmsd, TM_0, seq_id);



	//superposed structure B
	for(i=0; i<x_len; i++)
	{
		j=i;
		AAmap3(seqx[j], AA);
		fprintf(fp, "ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f\n", j+1, AA, xresno[j], xt[j][0], xt[j][1], xt[j][2]);
	}
	fprintf(fp, "TER\n");

	for(i=1; i<x_len; i++)
	{
		fprintf(fp, "CONECT%5d%5d\n", i, i+1);
	}
	//structure A
	for(i=0; i<y_len; i++)
	{
		j=i;
		AAmap3(seqy[j], AA);
		fprintf(fp, "ATOM  %5d  CA  %3s %5d    %8.3f%8.3f%8.3f\n", j+max+1, AA, yresno[j], ya[j][0], ya[j][1], ya[j][2]);
	}
	fprintf(fp, "TER\n");
	for(i=1; i<y_len; i++)
	{
		fprintf(fp, "CONECT%5d%5d\n", max+i, max+i+1);
	}

	fclose(fp);
}
*/
void alignment2AtomMap(
					   core::pose::Pose const & pose,
					   core::pose::Pose const & ref_pose,
					   std::list <core::Size> const & residue_list,
					   std::list <core::Size> const & ref_residue_list,
					   core::Size & n_mapped_residues,
					   core::id::AtomID_Map< core::id::AtomID > & atom_map)
{
	double seq_id;          
    int k;
    double d;
	
	do_rotation(xa, xt, xlen, t, u);
	
	seq_id=0;

	n_mapped_residues = 0;
	for(k=0; k<n_ali8_; k++)
	{
		d=sqrt(dist(xt[m1_[k]], ya[m2_[k]]));
		//std::cout << m1_[k] << " " << m2_[k] << " " << d << " " << d0_out_ << std::endl;
		if(d<d0_out_)
		{
			std::list<core::Size>::const_iterator it1 = residue_list.begin(); advance(it1, m1_[k]);
			core::Size ires = *it1;
			core::id::AtomID const id1( pose.residue_type(ires).atom_index("CA"), ires );

			std::list<core::Size>::const_iterator it2 = ref_residue_list.begin(); advance(it2, m2_[k]);
			core::Size jres = *it2;
			core::id::AtomID const id2( ref_pose.residue_type(jres).atom_index("CA"), jres );
			
			//std::cout << pose.residue_type(ires).name3() << " " << ires << " " << ref_pose.residue_type(jres).name3() << " " << jres << " " << std::endl;
			atom_map[ id1 ] = id2;
			++n_mapped_residues;
		}
	}
}
	
void alignment2strings(string & seqxA,
					   string & seqyA,
					   string & seqM
)
{
	double seq_id;          
    int i, j, k;
    double d;

	int ali_len=xlen+ylen; //maximum length of alignment
	seqM.resize(ali_len);
	seqxA.resize(ali_len);
	seqyA.resize(ali_len);
	
	do_rotation(xa, xt, xlen, t, u);
	
	seq_id=0;
	int kk=0, i_old=0, j_old=0;

	for(k=0; k<n_ali8_; k++)
	{
		for(i=i_old; i<m1_[k]; i++)
		{
			//align x to gap
			seqxA[kk]=seqx[i];
			seqyA[kk]='-';
			seqM[kk]=' ';					
			kk++;
		}
		
		for(j=j_old; j<m2_[k]; j++)
		{
			//align y to gap
			seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
		}
		
		seqxA[kk]=seqx[m1_[k]];
		seqyA[kk]=seqy[m2_[k]];
		if(seqxA[kk]==seqyA[kk])
		{
			seq_id++;
		}
		d=sqrt(dist(xt[m1_[k]], ya[m2_[k]]));
		if(d<d0_out_)
		{
			seqM[kk]=':';
		}
		else
		{
			seqM[kk]='.';
		} 
		kk++;  
		i_old=m1_[k]+1;
		j_old=m2_[k]+1;
	}
	
	//tail
	for(i=i_old; i<xlen; i++)
	{
		//align x to gap
		seqxA[kk]=seqx[i];
		seqyA[kk]='-';
		seqM[kk]=' ';					
		kk++;
	}    
	for(j=j_old; j<ylen; j++)
	{
		//align y to gap
		seqxA[kk]='-';
		seqyA[kk]=seqy[j];
		seqM[kk]=' ';
		kk++;
	}
	
    seqxA.resize(kk);
    seqyA.resize(kk);
    seqM.resize(kk);
	
	seq_id=seq_id/( n_ali8_+0.00000001); //what did by TMalign, but not reasonable, it should be n_ali8    
	
	//std::cout << seqxA << std::endl;
	//std::cout << seqM << std::endl;
	//std::cout << seqyA << std::endl;
}

//output the final results
	/*
void output_results(char *xname,
					 char *yname,
					 int x_len,
					 int y_len,
					 numeric::xyzVector <core::Real> const t,
					 numeric::xyzMatrix <core::Real> const u,
					 double TM1,
					 double TM2,
					 double rmsd,
					 double d0_out,
					 int m1[], 
					 int m2[],
					 int n_ali8,
					 int n_ali,
					 double TM_0,
					 double Lnorm_0,
					 double d0_0
					 )
{
    double seq_id;          
    int i, j, k;
    double d;
    int ali_len=x_len+y_len; //maximum length of alignment
	char *seqM, *seqxA, *seqyA;
	seqM=new char[ali_len];
	seqxA=new char[ali_len];
	seqyA=new char[ali_len];
	

	do_rotation(xa, xt, x_len, t, u);

	seq_id=0;
	int kk=0, i_old=0, j_old=0;
	for(k=0; k<n_ali8; k++)
	{
		for(i=i_old; i<m1[k]; i++)
		{
			//align x to gap
			seqxA[kk]=seqx[i];
			seqyA[kk]='-';
			seqM[kk]=' ';					
			kk++;
		}

		for(j=j_old; j<m2[k]; j++)
		{
			//align y to gap
			seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
		}

		seqxA[kk]=seqx[m1[k]];
		seqyA[kk]=seqy[m2[k]];
		if(seqxA[kk]==seqyA[kk])
		{
			seq_id++;
		}
		d=sqrt(dist(&xt[m1[k]][0], &ya[m2[k]][0]));
		if(d<d0_out)
		{
			seqM[kk]=':';
		}
		else
		{
			seqM[kk]='.';
		} 
		kk++;  
		i_old=m1[k]+1;
		j_old=m2[k]+1;
	}

	//tail
	for(i=i_old; i<x_len; i++)
	{
		//align x to gap
		seqxA[kk]=seqx[i];
		seqyA[kk]='-';
		seqM[kk]=' ';					
		kk++;
	}    
	for(j=j_old; j<y_len; j++)
	{
		//align y to gap
		seqxA[kk]='-';
		seqyA[kk]=seqy[j];
		seqM[kk]=' ';
		kk++;
	}
 
    seqxA[kk]='\0';
    seqyA[kk]='\0';
    seqM[kk]='\0';
	

	seq_id=seq_id/( n_ali8+0.00000001); //what did by TMalign, but not reasonable, it should be n_ali8    




 



	
	cout <<endl;	
	cout << " *****************************************************************************" << endl
		 << " * TM-align (Version "<< version <<"): A protein structural alignment algorithm     *" << endl
		 << " * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)       *" << endl
		 << " * Please email your comments and suggestions to Yang Zhang (zhng@umich.edu) *" << endl
		 << " *****************************************************************************" << endl;	


	
	printf("\nName of structure A: %s\n", yname); 
	printf("Name of structure B: %s (to be superimposed onto structure A)\n", xname); 
	printf("Length of structure A: %d residues\n", y_len);
	printf("Length of structure B: %d residues\n\n", x_len);

	printf("Aligned length=%d, RMSD=%6.2f, Seq_ID=n_identical/n_aligned=%4.3f\n\n", n_ali8, rmsd, seq_id); 
	
	printf("TM-score=%6.5f (if normalized by length of structure B, i.e., LN=%d, d0=%.2f)\n", TM2, x_len, d0B);
	printf("TM-score=%6.5f (if normalized by length of structure A, i.e., LN=%d, d0=%.2f)\n", TM1, y_len, d0A);
	if(a_opt)
	{
		double L_ave=(x_len+y_len)*0.5;
		printf("TM-score=%6.5f (if normalized by average length of two structures, i.e., LN=%.2f, d0=%.2f)\n", TM3, L_ave, d0a);
	}
	if(u_opt)
	{		
		printf("TM-score=%6.5f (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
	}
	if(d_opt)
	{		
		printf("TM-score=%6.5f (if scaled by user-specified d0=%.2f, and LN=%.2f)\n", TM5, d0_scale, Lnorm_0);
	}
	printf("(You should use TM-score normalized by length of the reference protein)\n");

   

    printf("\n----- The rotation matrix to rotate Structure B to Structure A -----\n");
    
    printf("i\t%18s %15s %15s %15s\n", "t[i]", "u[i][0]", "u[i][1]", "u[i][2]");
    for(k=0; k<3; k++)
    {
        printf("%d\t%18.10f %15.10f %15.10f %15.10f\n",\
                k, t[k], u[k][0], u[k][1] ,u[k][2]);
    }
    printf("\nCode for rotating Structure B from (x,y,z) to (X,Y,Z):\n");
    printf("for(k=0; k<L; k++)\n");
    printf("{\n");
    printf("   X[k] = t[0] + u[0][0]*x[k] + u[0][1]*y[k] + u[0][2]*z[k]\n");
    printf("   Y[k] = t[1] + u[1][0]*x[k] + u[1][1]*y[k] + u[1][2]*z[k]\n");
    printf("   Z[k] = t[2] + u[2][0]*x[k] + u[2][1]*y[k] + u[2][2]*z[k]\n");    
    printf("}\n"); 
    
    
    //output structure alignment
    printf("\n(\":\" denotes residue pairs of d < %4.1f Angstrom, ", d0_out);
    printf("\".\" denotes other aligned residues)\n");
    printf("%s\n", seqxA);
    printf("%s\n", seqM);
    printf("%s\n", seqyA);

	cout << endl;




	if(o_opt)
	{
		output_superpose(xname, yname, x_len, y_len, t, u, rmsd, d0_out, m1, m2, n_ali8, seq_id, TM_0, Lnorm_0, d0_0);
	}


	//delete [] seqM;
	//delete [] seqxA;
	//delete [] seqyA;
    
}
*/
void parameter_set4search(int xlen, int ylen)
{
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN=0.5; 
	dcu0=4.25;                       //update 3.85-->4.25
	
	Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
	if(Lnorm<=19)                    //update 15-->19
	{
		d0=0.168;                   //update 0.5-->0.168
	}
	else
	{
		d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
	}
	D0_MIN=d0+0.8;              //this should be moved to above
	d0=D0_MIN;                  //update: best for search    
	
	
	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;
	
	
	score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final(double len)
{
	D0_MIN=0.5; 
	
	Lnorm=len;            //normaliz TMscore by this in searching
	if(Lnorm<=21)         
	{
		d0=0.5;          
	}
	else
	{
		d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
	}
	if(d0<D0_MIN) d0=D0_MIN;   
	
	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  
	
}


void parameter_set4scale(int len, double d_s)
{
	
	d0=d_s;          
	Lnorm=len;            //normaliz TMscore by this in searching
	
	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  
	
}

int apply(core::pose::Pose const & pose1, core::pose::Pose const & pose2)
{
	std::list <core::Size> residue_list1;
	std::list <core::Size> residue_list2;
	for ( Size ires=1; ires<= pose1.total_residue(); ++ires ) {
		if ( !pose1.residue(ires).is_protein() ) continue;
		residue_list1.push_back(ires);
	}
	for ( Size ires=1; ires<= pose2.total_residue(); ++ires ) {
		if ( !pose2.residue(ires).is_protein() ) continue;
		residue_list2.push_back(ires);
	}

	return apply(pose1, pose2, residue_list1, residue_list2);
}
	
core::Real TMscore(Size length) {
	d0_out_=5.0;  
	Size simplify_step=1;
	Size score_sum_method=0;
	
	//normalized by user assigned length		
	parameter_set4final(length);		
	d0A=d0;
	//d0_0=d0A;
	//Lnorm_0=Lnorm_ass;
	core::Real rmsd;
	return TMscore8_search(xtm, ytm, n_ali8_, t, u, simplify_step, score_sum_method, &rmsd);	
}
	
int apply(core::pose::Pose const & pose1, core::pose::Pose const & pose2, std::list <core::Size> residue_list1, std::list <core::Size> residue_list2)
{
	/*********************************************************************************/
	/*                                load data                                      */ 
	/*********************************************************************************/
	load_pose_allocate_memory(pose1, pose2, residue_list1, residue_list2);
	
	
	/*********************************************************************************/
	/*                                parameter set                                  */ 
	/*********************************************************************************/
	parameter_set4search(xlen, ylen);          //please set parameters in the function
	int simplify_step     = 40;               //for similified search engine
	int score_sum_method  = 8;                //for scoring method, whether only sum over pairs with dis<score_d8
	
	int i;
	vector < int > invmap0(ylen+1); 
	vector < int > invmap(ylen+1); 
	double TM, TMmax=-1;
	for(i=0; i<ylen; i++)
	{
		invmap0[i]=-1;
	}	
	
	
	double ddcc=0.4;
	if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
	
	
	/*********************************************************************************/
	/*         get initial alignment with gapless threading                          */ 
	/*********************************************************************************/
	get_initial(xa, ya, xlen, ylen, invmap0);
	//find the max TMscore for this initial alignment with the simplified search_engin
	TM=detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method);
	if(TM>TMmax)
	{
		TMmax=TM;
	}           
	//run dynamic programing iteratively to find the best alignment
	TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
	if(TM>TMmax)
	{        
		TMmax=TM;
		for(int i=0; i<ylen; i++)
		{
			invmap0[i]=invmap[i];
		}
	}
	
	/*********************************************************************************/
	/*         get initial alignment based on secondary structure                    */ 
	/*********************************************************************************/	
	get_initial_ss(xa, ya, xlen, ylen, invmap);
	TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
	if(TM>TMmax)
	{
		TMmax=TM;
		for(int i=0; i<ylen; i++)
		{
			invmap0[i]=invmap[i];
		}
	} 
	if(TM > TMmax*0.2)
	{
		TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
		if(TM>TMmax)
		{
			TMmax=TM;
			for(int i=0; i<ylen; i++)
			{
				invmap0[i]=invmap[i];
			}
		}   
	}
	//output_align(invmap0, ylen);
	
	
	
	/*********************************************************************************/
	/*         get initial alignment based on local superposition                    */ 
	/*********************************************************************************/	
	//=initial5 in original TM-align
	if(get_initial_local(xa, ya, xlen, ylen, invmap))
	{		
		TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
		if(TM>TMmax)
		{
			TMmax=TM;
			for(int i=0; i<ylen; i++)
			{
				invmap0[i]=invmap[i];
			}
		} 
		if(TM > TMmax*ddcc)
		{
			TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 2);
			if(TM>TMmax)
			{
				TMmax=TM;
				for(int i=0; i<ylen; i++)
				{
					invmap0[i]=invmap[i];
				}
			}   
		}
		//output_align(invmap0, ylen);
	}    
	else
	{
		cout << endl << endl << "Warning: initial alignment from local superposition fail!" << endl << endl <<endl;		
	}
	
	
	/*********************************************************************************/
	/*    get initial alignment based on previous alignment+secondary structure      */ 
	/*********************************************************************************/	
	//=initial3 in original TM-align
	get_initial_ssplus(xa, ya, xlen, ylen, invmap0, invmap);
	TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);
	if(TM>TMmax)
	{
		TMmax=TM;
		for(i=0; i<ylen; i++)
		{
			invmap0[i]=invmap[i];
		}
	} 
	if(TM > TMmax*ddcc)
	{
		TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30);
		if(TM>TMmax)
		{
			TMmax=TM;
			for(i=0; i<ylen; i++)
			{
				invmap0[i]=invmap[i];
			}
		}   
	}
	//output_align(invmap0, ylen);   
	
	/*********************************************************************************/
	/*        get initial alignment based on fragment gapless threading              */ 
	/*********************************************************************************/	    
	//=initial4 in original TM-align
	get_initial_fgt(xa, ya, xlen, ylen, xresno, yresno, invmap);
	TM=detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method);	
	if(TM>TMmax)
	{
		TMmax=TM;
		for(i=0; i<ylen; i++)
		{
			invmap0[i]=invmap[i];
		}
	} 
	if(TM > TMmax*ddcc)
	{
		TM=DP_iter(xa, ya, xlen, ylen, t, u, invmap, 1, 2, 2);
		if(TM>TMmax)
		{         
			TMmax=TM;
			for(i=0; i<ylen; i++)
			{
				invmap0[i]=invmap[i];
			}
		}   
	} 
	//output_align(invmap0, ylen); 
	
	//*********************************************************************************//
	//     The alignment will not be changed any more in the following                 //
	//*********************************************************************************//
	//check if the initial alignment is generated approately	
	bool flag=false;
	for(i=0; i<ylen; i++)
	{
		if(invmap0[i]>=0)
		{
			flag=true;
			break;			
		}			
	}		
	if(!flag) 
	{
		cout << "There is no alignment between the two proteins!" << endl;
		cout << "Program stop with no result!" << endl;
		return 1;
	}
	//cout << "final alignment" << endl;
	//output_align(invmap0, ylen);
	
	
	
	//*********************************************************************************//
	//       Detailed TMscore search engine  --> prepare for final TMscore             //
	//*********************************************************************************//       
	//run detailed TMscore search engine for the best alignment, and 
	//extract the best rotation matrix (t, u) for the best alginment
	simplify_step=1;
	score_sum_method=8;
	TM=detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method);
	
	//select pairs with dis<d8 for final TMscore computation and output alignment
	n_ali8_=0;
	int k=0;
	int n_ali=0;
	double d;
	m1_.resize(xlen); //alignd index in x
	m2_.resize(ylen); //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k=0;
	for(int j=0; j<ylen; j++)
	{
		i=invmap0[j];
		if(i>=0)//aligned
		{
			n_ali++;        
			d=sqrt(dist(xt[i], ya[j]));
			if(d <= score_d8)
			{
				m1_[k]=i;
				m2_[k]=j;
				
				xtm[k]=xa[i];
				ytm[k]=ya[j];
				
				k++;
			}
		}
	}
	n_ali8_=k;
	
	
	
	//*********************************************************************************//
	//                               Final TMscore                                     //
	//                     Please set parameters for output                            //
	//*********************************************************************************//
	double rmsd, TM1, TM2;
	d0_out_=5.0;  
	simplify_step=1;
	score_sum_method=0;
	
	numeric::xyzVector <core::Real> t0;
	numeric::xyzMatrix <core::Real> u0;
	//double t0[3], u0[3][3];
	double d0_0, TM_0;
	double Lnorm_0=ylen;
	
	
	//normalized by length of structure A
	parameter_set4final(Lnorm_0);
	d0A=d0;
	d0_0=d0A;
	TM1=TMscore8_search(xtm, ytm, n_ali8_, t0, u0, simplify_step, score_sum_method, &rmsd);
	TM_0=TM1;
	
	//normalized by length of structure B
	parameter_set4final(xlen+0.0);
	d0B=d0;
	TM2=TMscore8_search(xtm, ytm, n_ali8_, t, u, simplify_step, score_sum_method, &rmsd);
	
	/*
	std::cout << d0_out_ << std::endl;
	if(a_opt)
	{
		//normalized by average length of structures A, B
		Lnorm_0=(xlen+ylen)*0.5;
		parameter_set4final(Lnorm_0);
		d0a=d0;
		d0_0=d0a;
		TM3=TMscore8_search(xtm, ytm, n_ali8_, t, u, simplify_step, score_sum_method, &rmsd);
		TM_0=TM3;
	}
	if(u_opt)
	{	
		//normalized by user assigned length		
		parameter_set4final(Lnorm_ass);		
		d0u=d0;		
		d0_0=d0u;
		Lnorm_0=Lnorm_ass;
		TM4=TMscore8_search(xtm, ytm, n_ali8_, t, u, simplify_step, score_sum_method, &rmsd);	
		TM_0=TM4;
	}
	if(d_opt)
	{	
		//scaled by user assigned d0
		parameter_set4scale(ylen, d0_scale);
		d0_out_=d0_scale;
		d0_0=d0_scale;
		//Lnorm_0=ylen;
		Lnorm_d0=Lnorm_0;
		TM5=TMscore8_search(xtm, ytm, n_ali8_, t, u, simplify_step, score_sum_method, &rmsd);	
		TM_0=TM5;
	}
	
	std::cout << d0_out_ << std::endl;
	*/
	//output_results(xname, yname, xlen, ylen, t0, u0, TM1, TM2, rmsd, d0_out, m1, m2, n_ali8, n_ali, TM_0, Lnorm_0, d0_0);
	
	
	
	//*********************************************************************************//
	//                            Done! Free memory                                    //
	//*********************************************************************************//           
	//free_memory();
	//delete [] invmap0;
	//delete [] invmap;
	//delete [] m1_;
	//delete [] m2_;
	return 0;
}
	
};
}
}
}
#endif
