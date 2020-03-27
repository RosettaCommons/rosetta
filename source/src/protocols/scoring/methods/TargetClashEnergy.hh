// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/TargetClashEnergy.hh
/// @brief  for de novo folding in the context of the target
/// @author Longxing (longxing@uw.edu)


#ifndef INCLUDED_protocols_scoring_methods_TargetClashEnergy_hh
#define INCLUDED_protocols_scoring_methods_TargetClashEnergy_hh


// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/scoring/methods/TargetClashEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// context hashing
#include <numeric/geometry/hashing/MinimalClashHash.hh>

#include <utility/vector1.hh>


// Project headers

// Utility headers


namespace protocols {
namespace scoring {
namespace methods {


class TargetClashEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:


	TargetClashEnergy() = delete;
	TargetClashEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	void initiate_voxel( );

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;
	void
	setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;
	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:
	// const-ref to scoring database
	std::string target_clash_pdb_;
	bool voxel_initialized_;
	core::Real clash_atom_scale_;
	std::string clash_check_resn_;
	numeric::geometry::hashing::MinimalClashHashOP context_clash_hash_;
	core::Size version() const override;
};


} // methods
} // scoring
} // protocols

#endif
