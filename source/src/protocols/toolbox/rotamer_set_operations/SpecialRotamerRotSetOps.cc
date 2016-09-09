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

#include <core/pose/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.RotamerSetOperations.SpecialRotamerRotSetOps" );

SpecialRotamerRSO::SpecialRotamerRSO( core::Size seqpos )
: parent(),
	seqpos_(seqpos)
{
	new_rots_.clear();
}

SpecialRotamerRSO::SpecialRotamerRSO( SpecialRotamerRSO const & src )
: parent( src ),
	seqpos_(src.seqpos_),
	new_rots_(src.new_rots_)
{}


SpecialRotamerRSO::~SpecialRotamerRSO(){}

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
	//mjo commenting out 'sfxn' because it is unused and causes a warning
	core::scoring::ScoreFunction const & /*sfxn*/,
	//mjo commenting out 'ptask' because it is unused and causes a warning
	core::pack::task::PackerTask const & /*ptask*/,
	//mjo commenting out 'packer_neighbor_graph' because it is unused and causes a warning
	utility::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
)
{
	using namespace core::pack::rotamer_set;
	if ( new_rots_.empty() ) {
		std::cerr << "Warning!! The rotamer_set (new_rots_) you are trying to append is empty" << std::endl;
		return;
	}
	if ( rotamer_set.resid() == seqpos_ ) {
		// make sure to remove any rotamers from original set that were already set as variant pink
		core::pack::rotamer_set::Rotamers variant_rotamers;
		utility::vector1< bool > rotamers2dropb( rotamer_set.num_rotamers(), false );
		for ( core::Size r(1); r<=rotamer_set.num_rotamers(); ++r ) {
			if ( ! rotamer_set.rotamer(r)->has_variant_type( core::chemical::SPECIAL_ROT ) ) continue;
			rotamers2dropb[r] = true;
			core::conformation::ResidueOP variant_rot( core::pose::remove_variant_type_from_residue( *rotamer_set.rotamer(r), core::chemical::SPECIAL_ROT, pose ) );
			variant_rotamers.push_back( variant_rot );
			rotamer_set.add_rotamer( *variant_rot );
			rotamers2dropb.push_back( false );
		}
		rotamer_set.drop_rotamers( rotamers2dropb );

		// add the new rotamers
		for ( core::Size i(1); i <= new_rots_.size(); ++i ) {
			rotamer_set.add_rotamer( *new_rots_[i] );

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
