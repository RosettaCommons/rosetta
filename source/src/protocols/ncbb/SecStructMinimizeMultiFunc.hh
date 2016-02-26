// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

using namespace core;
using namespace core::optimization;

/// @brief class for fitting a length of secondary structure keeping the same dihedrals
class SecStructMinimizeMultiFunc : public Multifunc
{
public: // Creation

	/// @brief Destructor
	///
	virtual ~SecStructMinimizeMultiFunc();

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
	virtual
	Real
	operator ()( Multivec const & vars ) const;

	/// @brief Calculate function value derivatives
	///
	virtual
	void
	dfunc( Multivec const & vars, Multivec & dE_dvars ) const;

	/// @brief give short set of torsions from full dofs
	Multivec
	dofs_to_vars( Multivec const & dofs ) const;

	/// @brief give full dofs from short set of torsions
	Multivec
	vars_to_dofs( Multivec const & vars ) const;


	/// @brief Error state reached -- derivative does not match gradient
	///
	virtual
	void
	dump( Multivec const & vars, Multivec const & vars2 ) const;

private:

	/*************************************************************
	PRIVATE DATA
	*************************************************************/

	void get_dofs_for_pose0( );
	void get_dofs_map( );
	Real dofs_for_pose0( Size const i_dof ) const { return dofs_for_pose0_[i_dof]; }
	Multivec
	dEddofs_to_dEdvars( Multivec const & dEdtors ) const;

	// Map between BB torsionID <-> min_map DOF_ID
	// Ah! and I need a map to a vector of Size for bb to dof...
	std::map< Size, utility::vector1< Size > > map_BB_to_DOF_;
	std::map< Size, Size > map_DOF_to_BB_;

	void
	setup_minimization_graph( pose::Pose & pose, scoring::ScoreFunction const & sfxn, MinimizerMap const & min_map ) const;

	bool uniq_refers_to_beta( char const uniq ) const;

	core::pose::Pose & pose_;

	core::scoring::ScoreFunction & scorefxn_;

	core::optimization::MinimizerMap & min_map_;

	// Reference pose during minimization( set as initial structure )
	pose::Pose & pose0_;
	Multivec dofs_for_pose0_;

	std::string alpha_beta_pattern_;

	/// @brief The pattern at which dihedrals are applied
	///
	std::string dihedral_pattern_;
	std::map< Size, utility::vector1< core::id::TorsionID > > vars_index_to_torsion_id_;
	Size nvar_;


}; // SecStructMinimizeMultiFunc


} // namespace ncbb
} // namespace protocols


#endif // INCLUDED_protocols_ncbb_SecStructMinimizeMultiFunc_HH
