// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RG_LocalEnergy.hh
/// @brief  only calculates RG the length of the repeat. Uses code from RG_Energy_fast
/// @author TJ Brunette


#ifndef INCLUDED_core_scoring_methods_RG_LocalEnergy_hh
#define INCLUDED_core_scoring_methods_RG_LocalEnergy_hh


// Package headers
#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <basic/datacache/CacheableData.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class RG_Local_MinData:  public RG_MinData {
public:
	RG_Local_MinData():RG_MinData(){}
};

typedef utility::pointer::shared_ptr< RG_Local_MinData > RG_Local_MinDataOP;

class RG_LocalEnergy: public RG_Energy_Fast  {

public:
	/// @brief Defines a center of mass based RG calculation that is O(n) rather
	/// than O(n^2).
	RG_LocalEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;


	core::Real
	calculate_rg_score( pose::Pose const & pose ) const;

	core::Real
	calculate_rg_score(
		pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues) const;

	// derivatives
	virtual void setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	virtual void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:
	RG_Local_MinData const & mindata_from_pose( pose::Pose const & ) const;
	RG_Local_MinData & nonconst_mindata_from_pose( pose::Pose & ) const;
	core::Size lastRes_;
	core::Size firstRes_;
};


}
}
}

#endif // INCLUDED_core_scoring_methods_RG_Energy_Fast_HH
