// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FreeMoietyEnergy.cc
/// @brief  FreeMoiety energy method implementation
/// @author Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <core/scoring/methods/FreeMoietyEnergy.hh>
#include <core/scoring/methods/FreeMoietyEnergyCreator.hh>

// Package Headers
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

using namespace core::pose::full_model_info;

static basic::Tracer TR( "core.scoring.methods.FreeMoietyEnergy", basic::t_info );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the FreeMoietyEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FreeMoietyEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new FreeMoietyEnergy;
}

ScoreTypes
FreeMoietyEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( free_suite );
	sts.push_back( free_2HOprime );
	sts.push_back( free_side_chain );
	return sts;
}

/// ctor
FreeMoietyEnergy::FreeMoietyEnergy() :
	parent( new FreeMoietyEnergyCreator ),
	free_suite_bonus_( -1.0 ), // this is ad hoc for now.
	free_2HOprime_bonus_( -1.0 ), // this is ad hoc for now.
	free_sugar_bonus_( basic::options::option[ basic::options::OptionKeys::score::free_sugar_bonus ] ), // this is -1.0 by default (also ad hoc)
	pack_phosphate_penalty_( basic::options::option[ basic::options::OptionKeys::score::pack_phosphate_penalty ] ),
	free_side_chain_bonus_( basic::options::option[ basic::options::OptionKeys::score::free_side_chain_bonus ] )
{}

FreeMoietyEnergy::~FreeMoietyEnergy() {}

/// clone
core::scoring::methods::EnergyMethodOP
FreeMoietyEnergy::clone() const
{
	return new FreeMoietyEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

void
FreeMoietyEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const{
	nonconst_full_model_info( pose ); // does the setup.
}

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
FreeMoietyEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{

	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();
	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();

	if ( rsd.has_variant_type( "VIRTUAL_PHOSPHATE" ) ){
		Size const seqpos_in_full_model = res_list[ rsd.seqpos() ];
		if ( seqpos_in_full_model > 1 && !cutpoint_open_in_full_model.has_value( seqpos_in_full_model - 1 ) ){
			emap[ free_suite ] += free_suite_bonus_;
		}
	}
	if ( rsd.has_variant_type( "FIVE_PRIME_PACKABLE_PHOSPHATE" ) ){ // this always comes with a virtual phosphate
		emap[ free_suite ] += pack_phosphate_penalty_;
	}
	if ( rsd.has_variant_type( "THREE_PRIME_PACKABLE_PHOSPHATE" ) ){
		emap[ free_suite ] += pack_phosphate_penalty_;
	}
	if ( rsd.has_variant_type( "FIVE_PRIME_PHOSPHATE" ) ){
		// weird, I know. these variants should go on as intermediates in stepwise assembly, replacing
		// virtual phosphates. Ideally the bonus & the penalty should cancel, but somehow that is too stringent
		// and would end up disallowing phosphates from ever being instantiated.
		emap[ free_suite ] += free_suite_bonus_ + pack_phosphate_penalty_;
	}
	if ( rsd.has_variant_type( "THREE_PRIME_PHOSPHATE" ) ){
		emap[ free_suite ] += pack_phosphate_penalty_;
	}
	if ( rsd.has_variant_type( "VIRTUAL_RIBOSE" ) )	emap[ free_suite ] += free_sugar_bonus_;

	if ( rsd.has_variant_type( "VIRTUAL_O2PRIME_HYDROGEN" ) )	emap[ free_2HOprime ] += free_2HOprime_bonus_;

	if ( rsd.has_variant_type( "VIRTUAL_SIDE_CHAIN" ) ) 	emap[ free_side_chain ] += free_side_chain_bonus_ * rsd.nchi() ;

}


/// @brief FreeMoietyEnergy is context independent; indicates that no context graphs are required
void
FreeMoietyEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
FreeMoietyEnergy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core

