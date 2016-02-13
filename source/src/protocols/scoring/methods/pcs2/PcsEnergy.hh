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
/// @file protocols/scoring/methods/pcs2/PcsEnergy.hh
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


#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsEnergy_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsEnergy_hh

// Package headers
//#include <protocols/scoring/methods/pcs2/PcsEnergy.fwd.hh>

#include <protocols/scoring/methods/pcs2/PcsDataCenter.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenterManager.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsDataCenterManagerSingleton.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsTensor.fwd.hh>
#include <protocols/scoring/methods/pcs2/PcsEnergy.fwd.hh>
#include <protocols/scoring/methods/pcs2/GridSearchIteratorCA.fwd.hh>

// Project headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/optimization/Minimizer.fwd.hh>

// Utility headers

// Numeric headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>


// Objexx headers

// C++ headers


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsEnergy : public core::scoring::methods::WholeStructureEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent;

public:

	PcsEnergy(); //Construct

	~PcsEnergy(); //Destruct

	PcsEnergy & // =
	operator=(PcsEnergy const & other);

	PcsEnergy(PcsEnergy const & other); //copy

	virtual core::scoring::methods::EnergyMethodOP
	clone() const; // Clone

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const;

	/// @brief This is called to start the PCS machinerie and get the score (set in totals)
	void
	finalize_total_energy(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const;

	/// @brief Return the PCS score given the pose, the given PcsDataCenter, and the lanthanide number
	core::Real
	calculate_pcs_score_on_PCS_data_center_CA( core::pose::Pose & pose,
		bool print_to_tracer,
		PcsDataCenter & pcs_d,
		core::Size i_multi_data,
		GridSearchIteratorCA & grid_it ) const;


	/// @brief Return the PCS score given the pose, the given PcsDataCenter, and the lanthanide number,
	/// return the vec of best score, vec of best tensor, and the vec of best x-y-z coordinate
	core::Real
	CA_search_scores_and_tensors(
		utility::vector1<core::Real> & vec_best_score,
		utility::vector1<PcsTensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		core::pose::Pose const & pdb,
		PcsDataCenter & pcs_d,
		core::Size i_multi_data,
		GridSearchIteratorCA & grid_it
	) const;

	core::Real
	CA_search_scores_and_tensors_with_svd(utility::vector1<core::Real> & vec_best_score,
		utility::vector1<PcsTensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		core::pose::Pose const & /* pdb*/,
		PcsDataCenter & pcs_d_c,
		core::Size,
		GridSearchIteratorCA & grid_it) const;


	core::Real
	minimize_tensors_from_PCS_data(
		utility::vector1<PcsTensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		PcsDataCenter const & pcs_d,
		core::optimization::Minimizer & minimizer,
		utility::vector1<core::Real> & vect_to_opti
	) const;

	core::Real
	minimize_tensors_from_PCS_data_with_svd( utility::vector1<PcsTensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		PcsDataCenter const & /*pcs_d_c*/,
		core::optimization::Minimizer & minimizer,
		utility::vector1<core::Real> & vect_to_opti
	) const;

	core::Real
	minimize_tensors_fix_from_PCS_data(
		utility::vector1<PcsTensor> & vec_best_tensor,
		PcsDataCenter const & pcs_d/*,
		core::Real xM,
		core::Real yM,
		core::Real zM*/
	) const;


	/* core::Real
	calculate_single_score_and_tensor_from_PCS_data_per_lanthanides(
	PcsTensor & PCS_t,
	PcsDataLanthanide & pcs_d_p_l
	) const; */


	PcsDataCenterManager &
	PCS_multi_data_from_pose(core::pose::Pose & pose) const;

	PcsDataCenterManagerSingleton &
	PCS_multi_data_from_noone() const;

	void
	dump_PCS_info(
		utility::vector1<PcsTensor> const & vec_tensor,
		numeric::xyzVector< core::Real > const & best_coo,
		PcsDataCenter const &pcs_d
	) const;
	virtual
	core::Size version() const;
};

} //PCS
} //methods
} //scoring
} //core

#endif
