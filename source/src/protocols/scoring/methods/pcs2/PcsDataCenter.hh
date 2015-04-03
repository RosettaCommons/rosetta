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
 /// @file protocols/scoring/methods/pcs2/PcsDataCenter.hh
 ///
 /// @brief Hold the PCS data on which the SVD will be applyed
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

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenter_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenter_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsInputCenter.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsDataLanthanide.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/constants.hh>

// ObjexxFCL headers

// c++ headers

#ifdef WIN32
	#include <protocols/scoring/methods/pcs2/PcsInputLine.hh>
#endif


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

class PcsDataCenter : public utility::pointer::ReferenceCount{
private:
	core::Size n_lanthanides_;
	core::Size n_pcs_spin_;
	utility::vector1<PcsInputLine> PCS_data_line_all_spin_;
	utility::vector1<PcsDataLanthanide> PCS_data_per_lanthanides_all_;
	utility::vector1< utility::vector1<core::Real> > A_all_;
	utility::vector1<core::Real> X_all_;
	utility::vector1<core::Real> Y_all_;
	utility::vector1<core::Real> Z_all_;

public:
	PcsDataCenter(); //construct

	virtual ~PcsDataCenter(); //destruct

	PcsDataCenter(PcsDataCenter const &other); // copy

	PcsDataCenter &
	operator=( PcsDataCenter const & src ); // =

	PcsDataCenter(PcsInputCenter & pcs_i_c, core::Size start, core::Size end, core::Real individual_scale);

	/// @brief Give me the number of lanthanides for this center
	core::Size
	get_n_lanthanides() const;

	/// @brief Give me the matrix A_all_
	utility::vector1< utility::vector1<core::Real> >  const &
	get_A_all() const;
	/*
	/// @brief Give me the vector r5_all_
	utility::vector1<core::Real> const &
	get_r5_all() const;
	*/
	/// @brief Give me the vector X_all_
	utility::vector1<core::Real> const &
	get_X_all() const;

	/// @brief Give me the vector Y_all_
	utility::vector1<core::Real> const &
	get_Y_all() const;

	/// @brief Give me the vector Z_all_
	utility::vector1<core::Real> const &
	get_Z_all() const;

	/// @brief Give the index number of the PcsInputLine given
	core::Size
	where_is_line(PcsInputLine & pcs_i_l);

	/// @brief This is called each time the pose is changed
	void
	update_X_Y_Z_all(core::pose::Pose const & pose);

	/// @brief Call update_my_A_matrix for all lanthanide data.
	/// X Y Z are the new coordinate of the center
	void
	update_matrix_A_all(core::Real const X,
											core::Real const Y,
											core::Real const Z);

	/// @brief Call update_my_A_matrix for all lanthanide data.
	/// X Y Z are the new coordinate of the center
	/// It also update individual smaller matrice for svd
	void
	update_matrix_A_all_for_svd(core::Real const X,
															core::Real const Y,
															core::Real const Z);

	void
	update_matrix_A_all_for_cstyle(core::Real const X,
																 core::Real const Y,
																 core::Real const Z);


	/// @brief Give me the vector PCS_data_per_lanthanides_all_
	utility::vector1<PcsDataLanthanide>&
	get_pcs_data_per_lanthanides_all();

	/// @brief Give me the vector PCS_data_line_all_spin_
	const utility::vector1<PcsInputLine> &
	get_PCS_data_line_all_spin() const;

	/// @brief Give me the vector PCS_data_per_lanthanides_all_ (const version)
	const utility::vector1<PcsDataLanthanide>&
	get_pcs_data_per_lanthanides_all() const;

	/// @brief Print me
	friend std::ostream &
	operator<<(std::ostream& out, const PcsDataCenter & me);


private:
	void
	update_matrix_A();

	void
	update_matrix_A_cstyle();

};


core::Real
fill_A_line_fast(utility::vector1<core::Real> & A_line,
								 core::Real const xM,
								 core::Real const yM,
								 core::Real const zM,
								 core::Real const x,
								 core::Real const y,
								 core::Real const z
);

void
fill_A_line_slow(utility::vector1<core::Real> & A_line,
								 core::Real const xM,
								 core::Real const yM,
								 core::Real const zM,
								 core::Real const x,
								 core::Real const y,
								 core::Real const z
);


}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

static const core::Real OTHER_FACT_USI_PRECALC_FOR_A_3( (10000.0/12.0/ core::Real( numeric::constants::d::pi ) ) * 3.0 );

#endif
