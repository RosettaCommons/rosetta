// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/methods/SplitUnfoldedTwoBodyEnergyCreator.cc
/// @briefEnergy creator for the split unfolded two body energy method
/// @author Riley Simmons-Edler (rse231@nyu.edu)


#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergyCreator.hh>
#include <core/scoring/methods/SplitUnfoldedTwoBodyEnergy.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>

#include <iostream>

namespace core {
namespace scoring {
namespace methods {


methods::EnergyMethodOP SplitUnfoldedTwoBodyEnergyCreator::create_energy_method(const methods::EnergyMethodOptions & options) const
{
	if ( options.has_method_weights( split_unfolded_two_body ) ) {
				utility::vector1<Real> const & v = options.method_weights( split_unfolded_two_body );
				assert( v.size() == scoring::n_score_types );
				EnergyMap e;
				for ( Size ii = 1; ii < scoring::n_score_types; ++ii ) {
					e[(ScoreType)ii]=v[ii];
				}
			//using the same type option as unfolded state energy since those two need to match when using the split unfolded energy(since unfolded state holds the one body component).
				return SplitUnfoldedTwoBodyEnergyOP( new SplitUnfoldedTwoBodyEnergy( options.split_unfolded_label_type(), options.split_unfolded_value_type(), options.unfolded_energies_type(), e ) );
	}
	return SplitUnfoldedTwoBodyEnergyOP( new SplitUnfoldedTwoBodyEnergy( options.split_unfolded_label_type(), options.split_unfolded_value_type(), options.unfolded_energies_type() ) );
}

ScoreTypes SplitUnfoldedTwoBodyEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( split_unfolded_two_body );
	sts.push_back( fa_atr_ref );
	sts.push_back( fa_rep_ref );
	sts.push_back( fa_sol_ref );
	sts.push_back( fa_elec_ref );
	sts.push_back( hbond_ref );
	sts.push_back( dslf_fa13_ref );
	return sts;
}


}
}
}
