// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/membrane/AddMPLigandMover.cc
///
/// @brief  Add "single" ligand to to membrane pose
/// @details  Accommodate membrane protein ligand in the membrane framework by
///    reorganizing the current foldtree. Resulting foldtree will
///    keep the membrane attached to the COM and ligand to the closest
///    binding pocket residue, provided in the constructor.
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// #RosettaMPMover

// Unit Headers
#include <protocols/membrane/AddMPLigandMover.hh>
#include <protocols/membrane/AddMPLigandMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.membrane.AddMPLigandMover" );

namespace protocols {
namespace membrane {

////////////////////
/// Constructors ///
////////////////////

/// @brief Add membrane protein ligand mover
/// @details Attach ligand downstream in the foldtree
/// for refinement at the last residue as a default.
/// DO NOT USE
AddMPLigandMover::AddMPLigandMover() :
	protocols::moves::Mover(),
	closest_rsd_( 0 ),
	ligand_seqpos_( 0 )
{}

/// @brief Add Membrane protein ligand mover (custom)
/// @details Attach ligand downstream in the foldtree of the
/// closest residue to the binding pocket
AddMPLigandMover::AddMPLigandMover( core::Size closest_rsd, core::Size ligand_seqpos ) :
	protocols::moves::Mover(),
	closest_rsd_( closest_rsd ),
	ligand_seqpos_( ligand_seqpos )
{}

/// @brief Copy Constructor
/// @details Mkae a deep copy of this mover
AddMPLigandMover::AddMPLigandMover( AddMPLigandMover const & src ) :
	protocols::moves::Mover( src ),
	closest_rsd_( src.closest_rsd_ ),
	ligand_seqpos_( src.ligand_seqpos_ )
{}

/// @brief Destructor
AddMPLigandMover::~AddMPLigandMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
AddMPLigandMover::clone() const {
	return ( protocols::moves::MoverOP( new AddMPLigandMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
AddMPLigandMover::fresh_instance() const {
	return ( protocols::moves::MoverOP( new AddMPLigandMover() ) );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
AddMPLigandMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// Read in closest residue option
	if ( tag->hasOption( "closest_rsd" ) ) {
		closest_rsd_ = tag->getOption< Size >( "closest_rsd" );
	}

	// Read in sequence position of the ligand
	if ( tag->hasOption( "ligand_seqpos" ) ) {
		ligand_seqpos_ = tag->getOption< Size >( "ligand_seqpos" );
	}

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
AddMPLigandMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMPLigandMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
AddMPLigandMoverCreator::keyname() const {
	return AddMPLigandMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
AddMPLigandMoverCreator::mover_name() {
	return "AddMPLigandMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Mover Apply Method
void
AddMPLigandMover::apply( core::pose::Pose & pose ) {

	using namespace core::kinematics; 

	// Check the pose is a membrane framework pose
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot perform add ligand to membrane pose operation because this pose is not a membrane pose!" );
	}

	TR << "Closest Residue: " << closest_rsd_ << " Ligand Seqpos: " << ligand_seqpos_ << std::endl;

	// Check the closest rsd and ligand seqpos parameters are valid
	// Note - default constructor parameters will make at least one of these statements
	// fail so please don't use it!
	if ( closest_rsd_ <= 0 || closest_rsd_ > pose.total_residue() ) {
		utility_exit_with_message( "User specified closest residue to ligand is out of bounds.");
	}

	if ( ligand_seqpos_ <= 0 || ligand_seqpos_ > pose.total_residue() ) {
		utility_exit_with_message( "User specified ligand position in sequence is out of bounds.");
	}

	if ( closest_rsd_ == ligand_seqpos_ ) {
		utility_exit_with_message( "Cannot sepcify closest residue as the ligand. Self-attachemnt is not valid");
	}

	// Get the following parameters from the foldtree: Residue COM, MP rsd position
	core::Size mp_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );
	core::Size rsd_com( residue_center_of_mass( pose, 1, pose.total_residue() ) ); // COM calc includes the ligand!

	// Create a new simple foldtree (assumes ligand at the end!)
	FoldTree ft;
	ft.simple_tree( pose.total_residue() ); // Exclude Ligand and MEM
	ft.new_jump( rsd_com, mp_rsd, rsd_com );
	ft.new_jump( closest_rsd_, ligand_seqpos_, closest_rsd_ );
	ft.reorder( rsd_com );
	pose.fold_tree( ft );
	pose.fold_tree().show( std::cout );

	// Update jump number in membrane info (will always be 1 here)
	pose.conformation().membrane_info()->set_membrane_jump( 1 );

}

/// @brief Show the name of this mvoer
std::string
AddMPLigandMover::get_name() const {
	return "AddMPLigandMover";
}

} // membrane
} // protocols


