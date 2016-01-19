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

#ifndef INCLUDED_protocols_hybridization_TMalign_hh
#define INCLUDED_protocols_hybridization_TMalign_hh

// Class headers
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <core/id/AtomID_Map.hh>


namespace protocols {
namespace hybridization {

using core::Size;

class TMalign {

public:

	TMalign();

	void PrintErrorAndQuit(std::string sErrorString);

	template < class A > void
	ResizeArray( A & array, int Narray1, int Narray2);

	int read_pose(
		core::pose::Pose const & pose,
		std::list <core::Size> const & residue_list,
		std::vector < numeric::xyzVector < core::Real > > & a,
		std::string & seq,
		std::vector < int > & resno);

	double dist(numeric::xyzVector<core::Real> x, numeric::xyzVector<core::Real> y);

	void transform(
		numeric::xyzVector <core::Real> const & t,
		numeric::xyzMatrix <core::Real> const & u,
		numeric::xyzVector < core::Real > x,
		numeric::xyzVector < core::Real > & x1);

	void do_rotation(
		std::vector < numeric::xyzVector<core::Real> > const & x,
		std::vector < numeric::xyzVector<core::Real> > & x1,
		int len,
		numeric::xyzVector<core::Real> const & t,
		numeric::xyzMatrix<core::Real> const & u);

	//    Please note this function is not a correct implementation of
	//     the N-W dynamic programming because the score tracks back only
	//     one layer of the matrix. This code was exploited in TM-align
	//     because it is about 1.5 times faster than a complete N-W code
	//     and does not influence much the final structure alignment result.
	void NWDP_TM(Size const len1, Size const len2, double const gap_open, std::vector < int > & j2i);

	void NWDP_TM(
		std::vector < numeric::xyzVector<core::Real> > const & x,
		std::vector < numeric::xyzVector<core::Real> > const & y,
		int const len1, int const len2,
		numeric::xyzVector <core::Real> const & t, numeric::xyzMatrix <core::Real> const & u,
		double d02, double gap_open, std::vector < int > & j2i);
		//NW dynamic programming for alignment
		//not a standard implementation of NW algorithm
		//Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
		//Output: j2i[1:len2] \in {1:len1} U {-1}
		//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical



	void NWDP_TM(
		std::vector < int > const & secx,
		std::vector < int > const & secy,
		int const len1, int const len2,
		double gap_open, std::vector < int > & j2i);
		//NW dynamic programming for alignment
		//not a standard implementation of NW algorithm
		//Input: secondary structure secx, secy, and gap_open
		//Output: j2i[1:len2] \in {1:len1} U {-1}
		//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical

	void
	convert_xyz_to_vector(
		numeric::xyzVector <core::Real> const & x,
		std::vector <core::Real> & xx);

	void
	convert_xyz_to_matrix(
		numeric::xyzMatrix <core::Real> const & x,
		std::vector <std::vector <core::Real> > & xx);

	void
	convert_vector_to_xyz(
		std::vector <core::Real> const & x,
		numeric::xyzVector <core::Real> & xx);

	void
	convert_matrix_to_xyz(
		std::vector <std::vector <core::Real> > const & x,
		numeric::xyzMatrix <core::Real> & xx);

	void
	convert_xyz_v_to_vectors(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < std::vector < double > > & xx);


	// wrapper is a temp fix -yfsong
	bool Kabsch(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const n,
		int const mode,
		double *rms,
		numeric::xyzVector <core::Real> & t,
		numeric::xyzMatrix <core::Real> & u ); 


	// Implemetation of Kabsch algoritm for finding the best rotation matrix
	// ---------------------------------------------------------------------------
	// x    - x(i,m) are coordinates of atom m in set x            (input)
	// y    - y(i,m) are coordinates of atom m in set y            (input)
	// n    - n is number of atom pairs                            (input)
	// mode  - 0:calculate rms only                                (input)
	// 1:calculate rms,u,t                                 (takes longer)
	// rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
	// u    - u(i,j) is   rotation  matrix for best superposition  (output)
	// t    - t(i)   is translation vector for best superposition  (output)
	bool Kabsch(
		std::vector < std::vector < double > > const & x,
		std::vector < std::vector < double > > const & y,
		int const n,
		int const mode,
		double *rms,
		std::vector < double > & t,
		std::vector < std::vector < double > > & u );


	void load_pose_allocate_memory(core::pose::Pose const & pose1, core::pose::Pose const & pose2, std::list <core::Size>  & residue_list1, std::list <core::Size>  & residue_list2);

	//     1, collect those residues with dis<d;
	//     2, calculate TMscore
	int score_fun8(
		std::vector < numeric::xyzVector <core::Real> > const & xa,
		std::vector < numeric::xyzVector <core::Real> > const & ya,
		int const n_ali,
		double const d,
		std::vector <int> & i_ali,
		double *score1,
		int score_sum_method );


	// TMscore search engine
	// input:   two aligned vector sets: x, y
	//          scale parameter d0
	//          simplify_step: 1 or 40 or other integers
	//          score_sum_method: 0 for score over all pairs
	//                            8 for socre over the pairs with dist<score_d8
	// output:  the best rotaion matrix t0, u0 that results in highest TMscore
	double TMscore8_search(
		std::vector < numeric::xyzVector <core::Real> > const  & xtm,
		std::vector < numeric::xyzVector <core::Real> > const  & ytm,
		int Lali,
		numeric::xyzVector <core::Real> & t0,
		numeric::xyzMatrix <core::Real> & u0,
		int const simplify_step,
		int const score_sum_method,
		double *Rcomm );

	//Comprehensive TMscore search engine
	// input:   two vector sets: x, y
	//          an alignment invmap0[] between x and y
	//          simplify_step: 1 or 40 or other integers
	//          score_sum_method: 0 for score over all pairs
	//                            8 for socre over the pairs with dist<score_d8
	// output:  the best rotaion matrix t, u that results in highest TMscore
	double detailed_search(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const,
		int const y_len,
		std::vector < int > const & invmap0,
		numeric::xyzVector <core::Real> & t,
		numeric::xyzMatrix <core::Real> & u,
		int simplify_step,
		int score_sum_method);

	//compute the score quickly in three iterations
	double get_score_fast(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const , int const y_len, std::vector < int > const & invmap);

	//perform gapless threading to find the best initial alignment
	//input: x, y, x_len, y_len
	//output: y2x0 stores the best alignment: e.g.,
	//y2x0[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0
	//the jth element in y is aligned to a gap in x if i==-1
	double get_initial(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > & y2x );

	void smooth(std::vector < int > & sec, int const len);


	int sec_str(double dis13, double dis14, double dis15, double dis24, double dis25, double dis35);


	//1->coil, 2->helix, 3->turn, 4->strand
	void make_sec(std::vector < numeric::xyzVector <core::Real> > const & x, int const len, std::vector < int > & sec);


	//get initial alignment from secondary structure alignment
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g.,
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0
	//the jth element in y is aligned to a gap in x if i==-1
	void get_initial_ss(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len, int const y_len,
		std::vector < int > & y2x );


	// get_initial5 in TMalign
	//get initial alignment of local structure superposition
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g.,
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0
	//the jth element in y is aligned to a gap in x if i==-1
	bool get_initial_local(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > & y2x );


	//with invmap(i) calculate score(i,j) using RMSD rotation
	void score_matrix_rmsd(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > const & y2x );

	void score_matrix_rmsd_sec(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > const & y2x );


	//get initial alignment from secondary structure and previous alignments
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g.,
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0
	//the jth element in y is aligned to a gap in x if i==-1
	void get_initial_ssplus(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < numeric::xyzVector <core::Real> > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > & y2x0,
		std::vector < int > & y2x );


	void find_max_frag(
		std::vector < numeric::xyzVector <core::Real> > const & x,
		std::vector < int > const & resno,
		int const len, int *start_max, int *end_max);


	//perform fragment gapless threading to find the best initial alignment
	//input: x, y, x_len, y_len
	//output: y2x0 stores the best alignment: e.g.,
	//y2x0[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0
	//the jth element in y is aligned to a gap in x if i==-1
	double get_initial_fgt(
		std::vector < numeric::xyzVector < core::Real > > const & x,
		std::vector < numeric::xyzVector < core::Real > > const & y,
		int const x_len,
		int const y_len,
		std::vector < int > const & xresno,
		std::vector < int > const & yresno,
		std::vector < int > & y2x );

	//heuristic run of dynamic programing iteratively to find the best alignment
	//input: initial rotation matrix t, u
	//       vectors x and y, d0
	//output: best alignment that maximizes the TMscore, will be stored in invmap
	double DP_iter(
		std::vector < numeric::xyzVector < core::Real > > const & x,
		std::vector < numeric::xyzVector < core::Real > > const & y,
		int const x_len, int const y_len,
		numeric::xyzVector < core::Real > t,
		numeric::xyzMatrix < core::Real > u,
		std::vector < int > & invmap0,
		int const g1, int const g2, int const iteration_max );


// Adding signature that will work with Pyrosetta ( ie. can't pass in n_mapped_resiudes by refs from python)
	void alignment2AtomMap(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::id::AtomID_Map< core::id::AtomID > & atom_map);


	void alignment2AtomMap(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size & n_mapped_residues,
		core::id::AtomID_Map< core::id::AtomID > & atom_map);

	void alignment2AtomMap(
		core::pose::Pose const & pose,
		core::pose::Pose const & ref_pose,
		std::list <core::Size> const & residue_list,
		std::list <core::Size> const & ref_residue_list,
		core::Size & n_mapped_residues,
		core::id::AtomID_Map< core::id::AtomID > & atom_map);

	void alignment2strings(
		std::string & seqxA,
		std::string & seqyA,
		std::string & seqM );

	void parameter_set4search(int xlen, int ylen);

	void parameter_set4final(double len);

	void parameter_set4scale(int len, double d_s);

	core::Real TMscore(Size length);

	int apply(core::pose::Pose const & pose1, core::pose::Pose const & pose2);

	int apply(core::pose::Pose const & pose1, core::pose::Pose const & pose2, std::list <core::Size> residue_list1, std::list <core::Size> residue_list2);


private:
	double D0_MIN;
	double Lnorm;                                           //normalization length
	double score_d8, d0, d0_search, dcu0;                   //for TMscore search
	std::vector < std::vector < double > > score;           //Input score table for dynamic programming
	std::vector < std::vector < bool > >   path;            //for dynamic programming
	std::vector < std::vector < double > > val;             //for dynamic programming
	int xlen, ylen, minlen;                                 //length of proteins
	std::vector < numeric::xyzVector<core::Real> > xa, ya;  //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
	std::vector < int >   xresno, yresno;                   //residue numbers, used in fragment gapless threading
	std::vector < numeric::xyzVector <core::Real> > xtm, ytm;    //for TMscore search engine
	std::vector < numeric::xyzVector <core::Real> > xt;          //for saving the superposed version of r_1 or xtm
	std::string   seqx, seqy;                               //for the protein sequence
	std::vector < int >   secx, secy;                       //for the secondary structure
	std::vector < numeric::xyzVector <core::Real> > r1, r2;      //for Kabsch rotation
	numeric::xyzVector <core::Real> t;
	numeric::xyzMatrix <core::Real> u;                      //Kabsch translation vector and rotation matrix

	int n_ali8_;
	std::vector < int > m1_, m2_;
	double d0_out_;

	double d0A, d0B;

}; // class TMalign

}  // namespace hybridization

}  // namespace protocols

#endif
