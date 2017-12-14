// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/SpecialRotamerRotSetOps.cc
/// @brief  combine two rotamer sets
/// @author Summer Thyme, sthyme@gmail.com, Jan 2010

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/SpecialRotamerRotSetOps.hh>

//Project headers

#include <core/pose/variant_util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

static basic::Tracer tr( "protocols.toolbox.RotamerSetOperations.SpecialRotamerRotSetOps" );

SpecialRotamerRSO::SpecialRotamerRSO( core::Size seqpos )
: parent(),
	seqpos_(seqpos)
{
	new_rots_.clear();
}

SpecialRotamerRSO::SpecialRotamerRSO( SpecialRotamerRSO const & /*src*/ ) = default;


SpecialRotamerRSO::~SpecialRotamerRSO()= default;

core::pack::rotamer_set::RotamerSetOperationOP
SpecialRotamerRSO::clone() const{
	return core::pack::rotamer_set::RotamerSetOperationOP( new SpecialRotamerRSO( *this ) );
}


/// @brief all this alter_rotamer_set function does
/// is to append a rotamer_set that is a member
/// of this RSO to the input rotamer_set
/// Making sure that the variant type is not
/// on the rotamers of the original rotamer set
void
SpecialRotamerRSO::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::pack::task::PackerTask const &,
	utility::graph::GraphCOP,
	core::pack::rotamer_set::RotamerSet & rotamer_set
)
{
	using namespace core::pack::rotamer_set;
	if ( new_rots_.empty() ) {
		std::cerr << "Warning!! The rotamer_set (new_rots_) you are trying to append is empty" << std::endl;
		return;
	}
	if ( rotamer_set.resid() == seqpos_ ) {

		// Make two passes over the rotamer set. In the first pass, look for all of the existing
		// SPECIAL_ROT rotamers (i.e. "variant pink") and for each one, create a new rotamer
		// that does not have the variant that will be a replacement for the SPECIAL_ROT rotamer.
		// Record which rotamers should be replaced. Then add all of the new rotamers to the
		// Rotamer set using the add_rotamer_into_existing_group function.
		// Finally, iterate across the rotamer set a second time, and mark all the original
		// SPECIAL_ROT rotamers for deletion in a rotamers2dropb vector. Give this vector
		// to the rotamer_set so that it can delete them en masse.

		core::pack::rotamer_set::Rotamers variant_rotamers;
		std::set< core::conformation::ResidueCOP > rotamers_2_drop;
		for ( core::Size r(1); r<=rotamer_set.num_rotamers(); ++r ) {
			if ( ! rotamer_set.rotamer(r)->has_variant_type( core::chemical::SPECIAL_ROT ) ) continue;
			rotamers_2_drop.insert( rotamer_set.rotamer(r) );
			core::conformation::ResidueOP variant_rot( core::pose::remove_variant_type_from_residue( *rotamer_set.rotamer(r), core::chemical::SPECIAL_ROT, pose ) );
			variant_rotamers.push_back( variant_rot );
		}
		for ( auto const & variant_rot : variant_rotamers ) {
			rotamer_set.add_rotamer_into_existing_group( *variant_rot );
		}

		utility::vector1< bool > rotamers2dropb( rotamer_set.num_rotamers(), false );
		for ( core::Size r(1); r <= rotamer_set.num_rotamers(); ++r ) {
			if ( rotamers_2_drop.count( rotamer_set.rotamer(r) ) ) {
				//tr << "Temp: dropping rotamer " << r << "; has SPECIAL_ROT?: " << rotamer_set.rotamer(r)->has_variant_type( core::chemical::SPECIAL_ROT ) << std::endl;
				rotamers2dropb[ r ] = true;
			}
		}
		rotamer_set.drop_rotamers( rotamers2dropb );

		// add the new rotamers
		for ( core::Size i(1); i <= new_rots_.size(); ++i ) {
			rotamer_set.add_rotamer_into_existing_group( *new_rots_[i] );
		}
	}
}


// this function takes the input rotamers from the motif screen
void
SpecialRotamerRSO::set_new_rots(
	core::pack::rotamer_set::Rotamers & new_rots
)
{
	for ( core::Size i(1); i <= new_rots.size(); ++i ) {
		new_rots_.push_back( new_rots[i] );
	}
}

} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols
