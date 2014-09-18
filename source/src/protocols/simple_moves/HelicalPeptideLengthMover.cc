// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/HelicalPeptideLengthMover.cc
/// @brief Appends or prepends a single residue to a helical peptide that's attached to the rest of the pose by a jump
/// @authors Jacob Corn

// Unit Headers
#include <protocols/simple_moves/HelicalPeptideLengthMover.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <core/id/AtomID.hh>

#include <core/kinematics/MoveMap.hh>

//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/after_opts.hh>
//#include <basic/options/util.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <numeric/constants.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.HelicalPeptideLengthMover" );

// C++ Headers
#include <sstream>
#include <fstream>

// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

HelicalPeptideLengthMover::HelicalPeptideLengthMover( core::Size const seqpos, std::string const resname, residue_addition_type const addition_type ) :
	protocols::moves::Mover(),
	seqpos_(seqpos),
	resname_(resname),
	addition_type_( addition_type )
{
	protocols::moves::Mover::type( "HelicalPeptideLengthMover" );
}

/// @brief Places and minimizes a PeptideStaple (i+4 or i+7, depending on staple_gap) at seqpos in pose
void HelicalPeptideLengthMover::apply( core::pose::Pose & pose )
{
	if( pose.num_jump() < 1 ) {
		TR << "Pose does not contain a jump! Aborting!" << std::endl;
		return;
	}
	else if( pose.num_jump() > 1 ) {
		TR << "Pose contains more than one jump! Aborting!" << std::endl;
		return;
	}

/*
	if( pose.conformation().chain_begin(seqpos) > seqpos_ ) {
		TR << seqpos_ << " is before the beginning of chain 2! Aborting!" << std::endl;
		return;
	}
	else if( pose.conformation().chain_end(1) < seqpos _ ) {
		TR << seqpos_ << " is after the end of chain2! Aborting!" << std::endl;
		return;
	}
*/
	core::conformation::ResidueCOP residue = create_residue_from_resname_( resname_ );
	add_residue_( pose, seqpos_, residue, addition_type_ );
}

core::conformation::ResidueCOP HelicalPeptideLengthMover::create_residue_from_resname_( std::string const resname )
{
	core::conformation::ResidueCOP residue;

	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueType const default_restype( residue_set->name_map( resname ) );
	core::chemical::ResidueType const & min_restype1 = residue_set->get_residue_type_with_variant_removed( default_restype, core::chemical::UPPER_TERMINUS );
	core::chemical::ResidueType const & min_restype2 = residue_set->get_residue_type_with_variant_removed( min_restype1, core::chemical::LOWER_TERMINUS ) ;

	residue = core::conformation::ResidueFactory::create_residue( min_restype2 );

	return residue;
}


void HelicalPeptideLengthMover::add_residue_( core::pose::Pose & pose, core::Size const seqpos, core::conformation::ResidueCOP residue, residue_addition_type const append_prepend )
{
	utility::vector1< core::Real > helical_torsions;
	helical_torsions.push_back( -57.8 ); // phi
	helical_torsions.push_back( -47 ); // psi
	helical_torsions.push_back( 180 ); // omega
	if( append_prepend == append ) {
		pose.conformation().safely_append_polymer_residue_after_seqpos( *residue, seqpos-1, true /*build_idea_geometry*/ );
		pose.set_phi( seqpos, -57.8 );
		pose.set_psi( seqpos, -47.0 );
		pose.set_omega( seqpos, 180 );
	}
	else if ( append_prepend == prepend ) {
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *residue, seqpos, true /*build_ideal_geometry*/ );
		pose.set_phi( seqpos, -57.8 );
		pose.set_psi( seqpos, -47.0 );
		pose.set_omega( seqpos, 180 );

	}
	else if ( append_prepend == insert ) {
		TR << "Residue insertion is not supported!" << std::endl;
		return;
	}
	else {
		TR << "Unrecognized residue insertion type! Aborting residue addition!" << std::endl;
		return;
	}
}


} // moves
} // protocols
