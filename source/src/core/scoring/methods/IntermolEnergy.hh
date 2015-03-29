// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/IntermolEnergy.hh
/// @brief  Radius of gyration score for METHODS, to match Rosetta++
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_methods_IntermolEnergy_hh
#define INCLUDED_core_scoring_methods_IntermolEnergy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class IntermolEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:

	/// @brief Defines the cost of instantiating a new chain.
	IntermolEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

/////////////////////////////////
	void
	eval_atom_derivative(
	  id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2 ) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:

Size
get_num_chains_frozen( pose::Pose const & pose ) const;

private:

	virtual
	core::Size version() const;

	Real const penalty_at_1M_;
	Real const log_conc_;
};


}
}
}

#endif // INCLUDED_core_scoring_methods_IntermolEnergy_HH
