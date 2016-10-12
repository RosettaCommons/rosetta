// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_ncbb_SecStructMinimizeMultiFunc_hh
#define INCLUDED_protocols_ncbb_SecStructMinimizeMultiFunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/id/TorsionID.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

//C++ headers
#include <string>

namespace protocols {
namespace ncbb {

/// @brief class for fitting a length of secondary structure keeping the same dihedrals
class SecStructMinimizeMultiFunc : public core::optimization::Multifunc
{
public: // Creation

	/// @brief Destructor
	///
	~SecStructMinimizeMultiFunc() override;

	/// @brief Constructor
	///
	SecStructMinimizeMultiFunc(
		core::pose::Pose &pose,
		core::scoring::ScoreFunction &scorefxn_in,
		core::optimization::MinimizerMap & min_map,
		std::string const alpha_beta_pattern,
		std::string const dihedral_pattern
	);

public: // Methods

	/// @brief Calculate function value (rms squared)
	///

	core::Real
	operator ()( core::optimization::Multivec const & vars ) const override;

	/// @brief Calculate function value derivatives
	///

	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const override;

	/// @brief give short set of torsions from full dofs
	core::optimization::Multivec
	dofs_to_vars( core::optimization::Multivec const & dofs ) const;

	/// @brief give full dofs from short set of torsions
	core::optimization::Multivec
	vars_to_dofs( core::optimization::Multivec const & vars ) const;


	/// @brief Error state reached -- derivative does not match gradient
	///

	void
	dump( core::optimization::Multivec const & vars, core::optimization::Multivec const & vars2 ) const override;

private:

	/*************************************************************
	PRIVATE DATA
	*************************************************************/

	void get_dofs_for_pose0( );
	void get_dofs_map( );
	core::Real dofs_for_pose0( core::Size const i_dof ) const { return dofs_for_pose0_[i_dof]; }
	core::optimization::Multivec
	dEddofs_to_dEdvars( core::optimization::Multivec const & dEdtors ) const;

	// Map between BB torsionID <-> min_map DOF_ID
	// Ah! and I need a map to a vector of core::Size for bb to dof...
	std::map< core::Size, utility::vector1< core::Size > > map_BB_to_DOF_;
	std::map< core::Size, core::Size > map_DOF_to_BB_;

	void
	setup_minimization_graph( core::pose::Pose & pose, core::scoring::ScoreFunction const & sfxn, core::optimization::MinimizerMap const & min_map ) const;

	bool uniq_refers_to_beta( char const uniq ) const;

	core::pose::Pose & pose_;

	core::scoring::ScoreFunction & scorefxn_;

	core::optimization::MinimizerMap & min_map_;

	// Reference pose during minimization( set as initial structure )
	core::pose::Pose & pose0_;
	core::optimization::Multivec dofs_for_pose0_;

	std::string alpha_beta_pattern_;

	/// @brief The pattern at which dihedrals are applied
	///
	std::string dihedral_pattern_;
	std::map< core::Size, utility::vector1< core::id::TorsionID > > vars_index_to_torsion_id_;
	core::Size nvar_;


}; // SecStructMinimizeMultiFunc


} // namespace ncbb
} // namespace protocols


#endif // INCLUDED_protocols_ncbb_SecStructMinimizeMultiFunc_HH
