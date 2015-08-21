// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SurfaceEnergy.hh
/// @brief  Header file for SurfaceEnergy; Not really a CD1B energy; this class is needed for optE
/// @author Ron Jacak


#ifndef INCLUDED_core_pack_interaction_graph_SurfaceEnergy_hh
#define INCLUDED_core_pack_interaction_graph_SurfaceEnergy_hh

// Unit Headers
#include <core/pack/interaction_graph/SurfaceEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/scoring/InterchainPotential.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace pack {
namespace interaction_graph {


//----------------------------------------------------------------------------//
//------------------------ SurfaceEnergy Class ----------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines a (pseudo) context-dependent one-body surface energy.  Really, this class
/// is only being used as a hack for the optE protocol so that the non-PD surface
/// energy can be optimized together with the other PD-terms.
///
class SurfaceEnergy : public scoring::methods::ContextDependentOneBodyEnergy  {
public:
	typedef scoring::methods::ContextDependentOneBodyEnergy  parent;

public:

	SurfaceEnergy();

	virtual
	scoring::methods::EnergyMethodOP
	clone() const;

	virtual
	void
	setup_for_scoring( pose::Pose & pose, scoring::ScoreFunction const & ) const;

	virtual
	void
	residue_energy( conformation::Residue const & rsd, pose::Pose const &, scoring::EnergyMap & ) const;

	virtual
	void
	finalize_total_energy( pose::Pose & pose, scoring::ScoreFunction const &, scoring::EnergyMap & totals ) const;

	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {};


private:

	// const-ref to scoring database
	//InterchainPotential const & potential_;
	virtual
	core::Size version() const;

};


} // namespace interaction_graph
} // namespace pack
} // namespace core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
