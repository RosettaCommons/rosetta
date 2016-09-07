// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/WaterAdductHBondPotential.hh
/// @brief
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_WaterAdductHBondPotential_hh
#define INCLUDED_core_scoring_WaterAdductHBondPotential_hh

#include <core/scoring/WaterAdductHBondPotential.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class WaterAdductHBondPotential  : public utility::pointer::ReferenceCount{
public:
	typedef conformation::Residue Residue;

public:
	/// ctor
	WaterAdductHBondPotential();
	~WaterAdductHBondPotential() override;
	///
	Real
	water_adduct_hbond_score(
		Residue const & rsd1,
		Residue const & rsd2
	) const;


	Real
	h2o_hbond_score_1way(
		Residue const & h2o_rsd,
		Residue const & other_rsd
	) const;

	/// @details Fills the hbond set with hydrogen bonds to water; clears the original hbonds in the set
	/// and also resets the hbond_options in the input hbond set.
	void
	fill_h2o_hbond_set(
		pose::Pose const & pose,
		hbonds::HBondSet & hbond_set
	) const;

	void
	get_residue_residue_h2o_hbonds_1way(
		// input
		Residue const & don_rsd,
		Residue const & acc_rsd,
		int const & don_nb,
		int const & acc_nb,
		// output
		hbonds::HBondSet & hbond_set
	) const;

private:
	scoring::hbonds::HBondOptionsOP   hbondoptions_;
	scoring::hbonds::HBondDatabaseCOP hb_database_;
};


} // scoring
} // core

#endif
