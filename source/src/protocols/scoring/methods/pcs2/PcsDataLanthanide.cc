// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsDataLanthanide.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsDataLanthanide.hh>

// Package headers
#include <protocols/scoring/methods/pcs2/PcsTensor.hh>
#include <protocols/scoring/methods/pcs2/PcsInputLine.hh>

// Project headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/PCS.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


// Numeric headers

// Objexx headers

// C++ headers

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR_PcsDataLanthanide( "protocols.scoring.methods.pcs.PcsDataLanthanide" );

PcsDataLanthanide::~PcsDataLanthanide(){
}

PcsDataLanthanide::PcsDataLanthanide(PcsDataLanthanide const &other):
	filename_(other.filename_),  weight_(other.weight_), svd_s_(other.svd_s_)
{

	// TR_PcsDataLanthanide << " () called" << std::endl;

	n_pcs_ = other.n_pcs_;
	A_index_ = other.A_index_;
	// fstyle_A_ = other.fstyle_A_;
	cstyle_A_ = other.cstyle_A_;
	// fstyle_b_ = other.fstyle_b_;
	cstyle_b_ = other.cstyle_b_;
	cstyle_b_individual_scale_ = other.cstyle_b_individual_scale_;
	individual_scale_ = other.individual_scale_;
	normalization_1_ = other.normalization_1_;
	normalization_2_ = other.normalization_2_;
	normalization_3_ = other.normalization_3_;
	normalization_factor_ = other.normalization_factor_;
	normalization_factor_inversed_ = other.normalization_factor_inversed_;

}

core::Real
PcsDataLanthanide::get_weight() const{
	return weight_;
}


core::Real
PcsDataLanthanide::get_individual_scale() const{
	return individual_scale_;
}


core::Real
PcsDataLanthanide::get_normalization_factor() const{
	return normalization_factor_;
}


core::Real
PcsDataLanthanide::get_normalization_factor_inversed() const{
	return normalization_factor_inversed_;
}


PcsDataLanthanide &
PcsDataLanthanide::operator=( PcsDataLanthanide const & other )
{

	// TR_PcsDataLanthanide << " =  called" << std::endl;

	if ( this != &other ) {
		n_pcs_ = other.n_pcs_;
		A_index_ = other.A_index_;
		//  fstyle_A_ = other.fstyle_A_;
		cstyle_A_ = other.cstyle_A_;
		//  fstyle_b_ = other.fstyle_b_;
		cstyle_b_ = other.cstyle_b_;
		cstyle_b_individual_scale_ = other.cstyle_b_individual_scale_;
		individual_scale_ = other.individual_scale_;
		svd_s_ = other.svd_s_;
		normalization_1_ = other.normalization_1_;
		normalization_2_ = other.normalization_2_;
		normalization_3_ = other.normalization_3_;
		normalization_factor_ = other.normalization_factor_;
		normalization_factor_inversed_ = other.normalization_factor_inversed_;
	}
	return *this;
}

void
PcsDataLanthanide::update_my_A_matrix_for_svd(utility::vector1< utility::vector1<core::Real> > & A_all){
	core::Size i, j;
	for ( i = 1; i <= A_index_.size(); ++i ) {
		for ( j = 1; j <= 5; ++j ) {
			cstyle_A_[i][j] = A_all[A_index_[i]][j];
		}
	}
	svd_s_.set_matrix_A(cstyle_A_);
	//svd_s_.set_matrix_A_point(&cstyle_A_);
}

void
PcsDataLanthanide::update_my_A_matrix_for_cstyle(utility::vector1< utility::vector1<core::Real> > & A_all){
	core::Size i, j;
	for ( i = 1; i <= A_index_.size(); ++i ) {
		for ( j = 1; j <= 5; ++j ) {
			cstyle_A_[i][j] = A_all[A_index_[i]][j];
		}
	}
}


core::Real
PcsDataLanthanide::calculate_cost_only_with_svd(){

	core::Real score;

	if ( weight_ == 0 ) { //This is the individual weight for each lanthanide
		return 0;
	}

	svd_s_.run_decomp_svd();
	score = svd_s_.run_score_svd_without_solving();
	return(score*weight_/normalization_factor_);
}

void
PcsDataLanthanide::retrieve_tensor_from_svd(PcsTensor &pcs_t){

	// FArray1D< core::Real > const & X(svd_s_.get_svd_solution());
	utility::vector1<core::Real> const & f (svd_s_.get_svd_solution());
	pcs_t.reset_tensor(f[1], f[2], f[3], f[4], f[5]);
	// pcs_t.reset_tensor(X(1), X(2), X(3), X(4), X(5));

}

core::Real
PcsDataLanthanide::calculate_tensor_and_cost_with_svd(PcsTensor &pcs_t){

	core::Real score;

	//#define FAST_SCORING
	//FAST SCORING is actually not faster (nor slower)
	//FAST SCORING gives the PCS score without actually solving Ax = b
	//(but still, the matrix has to be decomposed.
	//So FAST_SCORING doesn't give the tensor parameters...
	//FAST SCORING might become handy when working on unassigned data, that's why I keep it.

	if ( weight_ == 0 ) { //This is the individual weight for each lanthanide
		return 0;
	}

	svd_s_.run_decomp_svd();

#ifdef FAST_SCORING
	score = svd_s_.run_score_svd_without_solving();
	pcs_t.reset_tensor(0, 0, 0, 0, 0); //Junk tensor to avoid warning at compilation
#else
	svd_s_.run_solve_svd();
	score = svd_s_.run_score_svd_on_matrix(cstyle_A_);
	// The 2 following calls are not that necessary for just scoring
	// FArray1D<core::Real> const & f (svd_s_.get_svd_solution());

	utility::vector1<core::Real> const & f (svd_s_.get_svd_solution());
	// pcs_t.reset_tensor(f(1), f(2), f(3), f(4), f(5));
	pcs_t.reset_tensor(f[1], f[2], f[3], f[4], f[5]);
#endif

	//Weight is the individual weight for each lanthanide

	return(score*weight_/normalization_factor_);
}

void
PcsDataLanthanide::calculate_tensor_only_with_svd(PcsTensor &pcs_t){

	// core::Real score;
	svd_s_.run_decomp_svd();
	svd_s_.run_solve_svd();
	// FArray1D<core::Real> f (svd_s_.get_svd_solution());
	// pcs_t.reset_tensor(f(1), f(2), f(3), f(4), f(5));
	utility::vector1<core::Real> const & f (svd_s_.get_svd_solution());
	pcs_t.reset_tensor(f[1], f[2], f[3], f[4], f[5]);
}


core::Size
PcsDataLanthanide::get_n_pcs() const{
	return (n_pcs_);
}

const utility::vector1<core::Size> &
PcsDataLanthanide::get_A_index() const{
	return(A_index_);
}

const utility::vector1<core::Real> &
PcsDataLanthanide::get_cstyle_b() const{
	return(cstyle_b_);
}


const utility::vector1<core::Real> &
PcsDataLanthanide::get_cstyle_b_individual_scale() const{
	return(cstyle_b_individual_scale_);
}

const utility::vector1< utility::vector1<core::Real> > &
PcsDataLanthanide::get_cstyle_A() const{
	return(cstyle_A_);
}

/*
const FArray1D< core::Real > &
PcsDataLanthanide::get_fstyle_b() const{
return(fstyle_b_);
}

const FArray2D< core::Real > &
PcsDataLanthanide::get_fstyle_A() const{
return(fstyle_A_);
}
*/

bool
do_I_skip(PcsInputLine & pcs_i_l, core::Size start, core::Size end){
	core::Size residue_number(pcs_i_l.get_residue_num());

	if ( start == 0 ) {
		return false;
	}

	if ( (residue_number >= start)&&(residue_number <= end) ) {
		return false;
	}
	return true;
}

core::Size
reduced_size(utility::vector1< PcsInputLine > const & pcs_i_l_v, core::Size start, core::Size end){
	core::Size i;
	core::Size n_pcs_reduced(0);

	for ( i = 1; i <= pcs_i_l_v.size(); ++i ) {
		PcsInputLine pcs_i_l (pcs_i_l_v[i]);
		if ( do_I_skip(pcs_i_l, start, end) ) {
			continue;
		}
		n_pcs_reduced ++;
	}
	// std::cerr<< n_pcs_reduced << " instead of " << pcs_i_l_v.size() << std::endl;
	return(n_pcs_reduced);
}


PcsDataLanthanide::PcsDataLanthanide(std::string filename,
	core::Real const weight,
	utility::vector1< PcsInputLine > & pcs_i_l_v,
	core::Size start,
	core::Size end,
	core::Real individual_scale):
	filename_(filename), weight_(weight), svd_s_(basic::svd::SVD_Solver(reduced_size(pcs_i_l_v, start, end), 5))
{

	// TR_PcsDataLanthanide << " constructor called" << std::endl;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size i;
	core::Size M, N;
	// utility::vector1<PcsInputLine>::iterator it;
	individual_scale_ = individual_scale;

	n_pcs_ =  reduced_size(pcs_i_l_v, start, end);
	A_index_.resize(n_pcs_);

	M = n_pcs_;
	N = 5;

	svd_s_ = basic::svd::SVD_Solver(M, N);
	utility::vector1<core::Real> vec_temp;
	// fstyle_A_.dimension(M, N);
	// fstyle_b_.dimension(M);

	cstyle_b_.resize(M);
	cstyle_b_individual_scale_.resize(M);

	cstyle_A_.resize(M);
	for ( i = 1; i <= M; i++ ) {
		(cstyle_A_[i]).resize(N);
	}

	normalization_1_ = 0;
	core::Real sum (0);

	core::Size reduced_i(1);
	for ( i = 1; i <= pcs_i_l_v.size(); ++i ) {
		PcsInputLine & pcs_i_l = pcs_i_l_v[i];
		if ( do_I_skip(pcs_i_l, start, end) ) {
			continue;
		}
		vec_temp.push_back(pcs_i_l.get_PCS_experimental());
		//  fstyle_b_(i) = pcs_i_l.get_PCS_experimental();
		cstyle_b_[reduced_i] = pcs_i_l.get_PCS_experimental();

		if ( individual_scale >= 0.0 ) {
			if ( i == 1 ) {
				std::cerr << "WARNING: INDIVIDUAL_NOARMALIZATION ON (message2) with " << individual_scale << std::endl;
			}
			cstyle_b_individual_scale_[reduced_i] = 1.0 / (pcs_i_l.get_PCS_experimental() * pcs_i_l.get_PCS_experimental() + individual_scale)  ;
			normalization_1_ += pcs_i_l.get_PCS_experimental() * pcs_i_l.get_PCS_experimental() * cstyle_b_individual_scale_[reduced_i];

		} else {
			normalization_1_ += pcs_i_l.get_PCS_experimental() * pcs_i_l.get_PCS_experimental();
		}

		sum += pcs_i_l.get_PCS_experimental();
		reduced_i++;
	}
	svd_s_.set_vector_b(vec_temp);
	normalization_1_ = sqrt(normalization_1_);

	normalization_3_ = sqrt(normalization_1_ / M);

	core::Real average(sum/M);

	normalization_2_ = 0;
	reduced_i = 1;
	for ( i = 1; i <= pcs_i_l_v.size(); ++i ) {
		if ( do_I_skip(pcs_i_l_v[i], start, end) ) {
			continue;
		}
		normalization_2_ += (cstyle_b_[reduced_i] - average) *  (cstyle_b_[reduced_i] - average);
		reduced_i++;
	}
	normalization_2_ = sqrt( normalization_2_ / M );

	if ( option[ basic::options::OptionKeys::PCS::normalization_id ].user() ) {
		core::Size norma_id (option[ basic::options::OptionKeys::PCS::normalization_id ]());

		switch (norma_id){
		case 1 : {
			TR_PcsDataLanthanide << "Using Normalization '1': " << normalization_1_ << ". Normalize each data set by SQRT( SUM( PCSexp(i)^2 ) ) for " << filename_ << std::endl;
			normalization_factor_ = normalization_1_;
			break;
		}
		case 2 : {
			TR_PcsDataLanthanide << "Using Normalization '2': " << normalization_2_ << ". Normalize each data set by Standard Deviation for " << filename_ << std::endl;
			normalization_factor_ = normalization_2_;
			break;
		}
		case 3 : {
			TR_PcsDataLanthanide << "Using Normalization '3': " << normalization_3_ << ". Normalize each data set by SQRT( SUM( PCSexp(i)^2 ) / N) for " << filename_ << std::endl;
			normalization_factor_ = normalization_3_;
			break;
		}
		default : {
			TR_PcsDataLanthanide << "Normalization '"<< norma_id << "' id not recognized " << std::endl;
			utility_exit_with_message("You should use a valid normalization id ('1' or '2' or '3')");
			break;
		}
		}
	} else {
		normalization_factor_ = 1.0;
	}

	normalization_factor_inversed_ = 1.0/normalization_factor_;

	// std::cerr << "OK" << std::endl;
}

std::string
PcsDataLanthanide::get_filename() const {
	return (filename_);
}

void
PcsDataLanthanide::set_A_index(core::Size j, core::Size index){
	A_index_[j] = index;
}

std::ostream &
operator<<(std::ostream& out, const PcsDataLanthanide &me){
	core::Size i; //, j;
	out << std::setprecision(4) << std::endl;
	out << "  Filename : " << me.get_filename() << std::endl;
	out << "  Number of pcs : " << me.get_n_pcs() << std::endl;
	out << "  b vector (pcs values)  : " << std::endl;

	for ( i = 1; i <= me.get_n_pcs(); ++i ) {
		out << "  " << std::setw(4) << i << ":"<< std::setw(10) <<me.cstyle_b_[i]  << std::endl;
	}
	out << std::endl;

	out << "  A_index  : " << std::endl;
	for ( i = 1; i <= me.get_n_pcs(); ++i ) {
		out << "  " << std::setw(4) << i << ":"  << std::setw(4) <<me.A_index_[i];
		out << std::endl;
	}
	/*
	out << "  Matrix A:"<< std::endl;
	for( i = 1; i <= me.n_pcs_ ; ++i ){
	out << "  " << std::setw(4) << i << ":";
	for (j = 1; j <= 5; ++j){
	out << " " << std::setw(12) << me.fstyle_A_( i, j);
	}
	out << std::endl;
	}
	*/
	return out;
}

}//namespcacs PCS
}//namespace methods
}//namespace scoring
}//namespace protocols
