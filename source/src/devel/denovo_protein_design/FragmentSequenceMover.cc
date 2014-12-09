// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/DenovoProteinDesign/FragmentSequenceMover.cc
/// @brief FragmentSequenceMover methods implemented
/// @author


// Unit Headers
#include <devel/denovo_protein_design/FragmentSequenceMover.hh>


// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>


// Utility Headers
#include <utility/assert.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "devel.DenovoProteinDesign.FragmentSequenceMover" );

namespace devel {
namespace denovo_protein_design {

///@details
void FragmentSequenceMover::apply( core::pose::Pose & pose ){

	// how do we respect the movemap?

	core::Size FramePosition = numeric::random::random_range( fragset_->min_pos() , fragset_->max_pos());

	core::fragment::FrameList FrameList;

	core::Size NumberofFrames = fragset_->frames( FramePosition, FrameList);
	core::Size FrameNumber = numeric::random::random_range( 1, NumberofFrames );
	core::Size NumberofFragmentsinthisFrame = FrameList[ FrameNumber ]->nr_frags();
	core::Size FragNumber = numeric::random::random_range( 1, NumberofFragmentsinthisFrame );

	FrameList[ FrameNumber ]->apply( FragNumber, pose ); // we can make this frame apply movemap aware
	// but we may try to insert in a non insertable region

	TR << " SEQUENCE OF FRAG INSERTION " << FrameList[ FrameNumber ]->fragment( FragNumber ).sequence() << std::endl;
	std::string new_seq = FrameList[ FrameNumber ]->fragment( FragNumber ).sequence();
	TR << "SIZE STRING " << new_seq.size() << std::endl;
	core::chemical::ResidueTypeSetCAP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	for ( Size i=1; i<= new_seq.size(); ++i ) {
		core::conformation::Residue const & old_rsd( pose.residue( FramePosition + i - 1 ) );
		for ( Size j=1; j<=old_rsd.type().variant_types().size(); j++ ) {
			std::cout << old_rsd.type().variant_types()[j] << std::endl;
		}
		// get all residue types with same AA
		char aa = new_seq[i-1];
		core::chemical::AA this_aa = core::chemical::aa_from_oneletter_code( aa );
		core::chemical::ResidueTypeCOPs const & rsd_types( residue_set->aa_map( this_aa ) );
		core::conformation::ResidueOP new_rsd( 0 );
		// now look for a rsdtype with same variants
		for ( Size j=1; j<= rsd_types.size(); ++j ) {
			core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
			if ( old_rsd.type().variants_match( new_rsd_type ) ) {
				new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() );
				break;
			}
		}

		pose.replace_residue( FramePosition + i - 1, *new_rsd, false );
	}

}//apply

///@brief
FragmentSequenceMover::FragmentSequenceMover(
) : Mover()
{
	Mover::type( "FragmentSequenceMover" );
}

///@brief
FragmentSequenceMover::FragmentSequenceMover( core::fragment::FragSetCOP fragset, core::kinematics::MoveMapCOP movemap
																							) : Mover(), fragset_(fragset), movemap_(movemap)
{
	Mover::type( "FragmentSequenceMover" );
}





FragmentSequenceMover::~FragmentSequenceMover(){}

}//DenovoProteinDesign
}//devel

