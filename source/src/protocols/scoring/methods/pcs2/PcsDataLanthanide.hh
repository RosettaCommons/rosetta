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
/// @file protocols/scoring/methods/pcs2/PcsDataLanthanide.hh
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

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsDataLanthanide_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsDataLanthanide_hh


// Package headers
#include <protocols/scoring/methods/pcs2/PcsTensor.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsInputLine.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/svd/SVD_Solver.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// c++ headers

//#define LOGPCS

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsDataLanthanide{
private:
	std::string const filename_;
	core::Real const weight_;

	core::Size n_pcs_;
	utility::vector1<core::Size> A_index_; //index on the giant matrix A to build all the smalls matrix A_
	//// ObjexxFCL::FArray2D< core::Real > fstyle_A_; //We are going to SVD Ax = b
	//// ObjexxFCL::FArray1D< core::Real > fstyle_b_; //I should make this one const
	utility::vector1< utility::vector1<core::Real> > cstyle_A_;
	utility::vector1<core::Real> cstyle_b_;
	utility::vector1<core::Real> cstyle_b_individual_scale_;
	basic::svd::SVD_Solver svd_s_;
	core::Real normalization_1_; // SQRT(SUMi( PCS(calc,i)^2 ) )
	core::Real normalization_2_; // Standard deviation
	core::Real normalization_3_; // SQRT(SUMi( PCS(calc,i)^2 )/N )
	core::Real normalization_factor_;
	core::Real normalization_factor_inversed_;
	core::Real individual_scale_;

public:


	void
	update_my_A_matrix_for_cstyle(utility::vector1< utility::vector1<core::Real> > & A_all);

	PcsDataLanthanide(); //construct

	~PcsDataLanthanide(); //destruct

	PcsDataLanthanide(PcsDataLanthanide const &other); //copy

	PcsDataLanthanide &
	operator=( PcsDataLanthanide const & other ); //=

	PcsDataLanthanide(std::string, core::Real const weight, utility::vector1< PcsInputLine > & pcs_i_l, core::Size start, core::Size end, core::Real individual_scale);

	/// @brief Set a value of the A_index_ vector
	void
	set_A_index(core::Size j, core::Size n_pcs_spin_);

	/// @brief update the A matrix given A_all matrix.
	/// Dimensions of A_all >= dimension of A.
	/// A_all is common to all the lanthanide sharing the same center
	void
	update_my_A_matrix_for_svd(utility::vector1< utility::vector1<core::Real> > & A_all);

	/// @brief give me the weight associated with this PCS data
	core::Real
	get_weight() const;

	core::Real
	get_individual_scale() const;


	/// @brief give me the normalization factor associated with this PCS data
	core::Real
	get_normalization_factor() const;

	/// @brief give me the normalization factor associated with this PCS data
	core::Real
	get_normalization_factor_inversed() const;


	/// @brief give me the filename associated with this PCS data
	std::string
	get_filename() const;

	/// @brief give me the number of PCS data
	core::Size
	get_n_pcs() const;

	/// @brief Give me the A_index_ vector
	utility::vector1<core::Size> const &
	get_A_index() const;

	/// @brief return the b vector in FArray1D format
	/*
	ObjexxFCL::FArray1D< core::Real > const &
	get_fstyle_b() const;

	ObjexxFCL::FArray2D< core::Real > const &
	get_fstyle_A() const;
	*/

	const utility::vector1< utility::vector1<core::Real> > &
	get_cstyle_A() const;

	const utility::vector1<core::Real> &
	get_cstyle_b() const;

	const utility::vector1<core::Real> &
	get_cstyle_b_individual_scale() const;

	/// @Print me
	friend
	std::ostream &
	operator << ( std::ostream& out, const PcsDataLanthanide &me );

	/// @This return the score and populate the PcsTensor with svd.
	core::Real
	calculate_tensor_and_cost_with_svd(PcsTensor &pcs_t);

	core::Real
	calculate_cost_only_with_svd();


	/// @This populate the PcsTensor with svd.
	void
	calculate_tensor_only_with_svd(PcsTensor &pcs_t);

	void
	retrieve_tensor_from_svd(PcsTensor &pcs_t);


	/*
	core::Real
	calculate_tensor_and_cost_with_svd_precalc(PcsTensor &pcs_t);
	*/
};

bool
do_I_skip(PcsInputLine & pcs_i_l, core::Size start, core::Size end);


}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
