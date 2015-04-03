// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HPatchEnergy.hh
/// @brief  HPatch energy method, a score based on exact SASA calculations and a new patch identification algorithm, header file
/// @author Ron Jacak (ronj@email.unc.edu)

#ifndef INCLUDED_core_pack_interaction_graph_HPatchEnergy_hh
#define INCLUDED_core_pack_interaction_graph_HPatchEnergy_hh

// Unit headers
#include <core/pack/interaction_graph/HPatchEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace interaction_graph {


//----------------------------------------------------------------------------//
//---------------------- HPatchEnergy Class -----------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines a (pseudo) context-dependent one-body surface energy.  Really, this class
/// is only being used as a hack for the optE protocol so that the non-PD surface
/// energy can be optimized together with the other PD-terms.
/// The difference from this energy method from the plain SurfaceEnergy method is that
/// it calculates the patch area using methods in sasa.cc instead of using average
/// values. This new method also uses a new approach for finding which residues to include
/// in a patch, not just all residues within 10A.
///
class HPatchEnergy : public scoring::methods::ContextDependentOneBodyEnergy  {
public:
	typedef scoring::methods::ContextDependentOneBodyEnergy  parent;

public:
	HPatchEnergy();

	virtual
	scoring::methods::EnergyMethodOP clone() const;

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
virtual
core::Size version() const;


private:

};

} // namespace interaction_graph
} // namespace pack
} // namesapce core


#endif // INCLUDED_core_pack_interaction_graph_HPatchEnergy_HH
