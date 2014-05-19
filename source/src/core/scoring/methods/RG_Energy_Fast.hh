// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RG_Energy_Fast.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author James Thompson


#ifndef INCLUDED_core_scoring_methods_RG_Energy_Fast_hh
#define INCLUDED_core_scoring_methods_RG_Energy_Fast_hh


// Package headers
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

class RG_MinData:  public basic::datacache::CacheableData {
public:
	RG_MinData(): com(0,0,0), rg(0), nres_scored(0) {}

	basic::datacache::CacheableDataOP clone() const {
		return new RG_MinData(*this);
	}

	numeric::xyzVector< core::Real > com;
	core::Real rg;
	core::Size nres_scored;
};

typedef utility::pointer::owning_ptr< RG_MinData > RG_MinDataOP;

class RG_Energy_Fast : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;
public:

	/// @brief Defines a center of mass based RG calculation that is O(n) rather
	/// than O(n^2).
	RG_Energy_Fast();

	/// clone
	virtual
	EnergyMethodOP
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


	virtual
	core::Size version() const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}
    
    ///@brief function to get membrane protein TMH centers of mass
    utility::vector1<core::Size>
    get_tmh_CoMs(pose::Pose const& pose) const;
    
private:
	RG_MinData const & mindata_from_pose( pose::Pose const & ) const;
	RG_MinData & nonconst_mindata_from_pose( pose::Pose & ) const;
};


}
}
}

#endif // INCLUDED_core_scoring_methods_RG_Energy_Fast_HH
