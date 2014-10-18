// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/symmetry/SetupForMembraneSymmetry.hh
///
/// @brief		Setup a Symmetric Membrane Protein Using the Membrane Framework
/// @details	The setup for membrane symmetry class first adds a membrane residue
///				to the asymmetric unit of a protein, creates a symmetric complex, and
///				then adds the remainder of the membrane framework to capture the entire
///				symmetric system. This should work for both relax and docking protocols
///				currently in Rosetta.
///
///				Last Modified: 10/3/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_cc
#define INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_cc

// Unit Headers
#include <protocols/membrane/symmetry/SetupForMembraneSymmetry.hh>
#include <protocols/membrane/symmetry/SetupForMembraneSymmetryCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <core/conformation/symmetry/SymmData.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/geometry/util.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utiility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

static thread_local basic::Tracer TR( "protocols.membrane.symmetry.SetupForMembraneSymmetry" );

namespace protocols {
namespace membrane {
namespace symmetry {

using namespace protocols::moves;
using namespace core::pose;
	
/// @brief Default Constructor for Setup for MP Symm Mover
SetupForMembraneSymmetry::SetupForMembraneSymmetry() :
	Mover(),
	symmdef_file_( "" )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor for Setup for MP Symm Mover
SetupForMembraneSymmetry::SetupForMembraneSymmetry( std::string symmdef_file ) :
	Mover(),
	symmdef_file_( symmdef_file )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a Deep Copy of this setup for membrane symmetry mover
SetupForMembraneSymmetry::SetupForMembraneSymmetry( SetupForMembraneSymmetry const & src ) :
	Mover( src ),
	symmdef_file_( src.symmdef_file_ )
{}

/// @brief Default Destructor
SetupForMembraneSymmetry::~SetupForMembraneSymmetry() {}

/// @brief
MoverOP
SetupForMembraneSymmetry::clone() const {
	return ( MoverOP( new SetupForMembraneSymmetry( *this ) ) );
}

/// @brief
MoverOP
SetupForMembraneSymmetry::fresh_instance() const {
	return MoverOP( new SetupForMembraneSymmetry() );
}

/// @brief
void
SetupForMembraneSymmetry::parse_my_tag(
	 utility::tag::TagCOP tag,
	 basic::datacache::DataMap &,
	 protocols::filters::Filters_map const &,
	 protocols::moves::Movers_map const &,
	 core::pose::Pose const &
) {

	// Read in symmetry definition file
	if ( tag->hasOption( "symmdef_file" ) ) {
		symmdef_file_ = tag->getOption< std::string >( "symmdef_file" );
	}
}
	
/// @brief Create a new copy of this mover
MoverOP
SetupForMembraneSymmetryCreator::create_mover() const {
	return MoverOP( new SetupForMembraneSymmetry );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SetupForMembraneSymmetryCreator::keyname() const {
	return SetupForMembraneSymmetryCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SetupForMembraneSymmetryCreator::mover_name() {
	return "SetupForMembraneSymmetry";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply setup for membrane symmetry move
void
SetupForMembraneSymmetry::apply( Pose & pose ) {
	
	using namespace core::conformation::symmetry;
	using namespace protocols::simple_moves::symmetry;
	using namespace protocols::membrane;
	
	// If user has yet to provide a symmetry definition file, fail with no return
	if ( symmdef_file_ == "" ) {
		utility_exit_with_message("Cannot setup a membrane protein for symmetry without starting symmetry definition file. Did you forget to pass -symmetry::symmetry_definition?");
	}
	
	// Add Membrane Residue to Pose (Part of the Asymmetric Unit)
	core::SSize membrane_pos = setup_asymm_membrane( pose );
	
	// Modify Symmdata to point to the membrane residue in the master
	// asymmetric unit
	SymmDataOP symmdef( new SymmData() );
	symmdef->read_symmetry_data_from_file( symmdef_file_ );
	symmdef->set_anchor_residue( utility::to_string( membrane_pos ) );
	
	// Setup for Symmetry (Create symmetric unit)
	SetupForSymmetryMoverOP symm_mover( new SetupForSymmetryMover( symmdef ) );
	symm_mover->apply( pose );
	
	// Position the membrane based on the full symmetric subunit??
	// needs to answer this question - leaving it out for now
	
	// Add the membrnae component to the full symmetric unit
	AddMembraneMoverOP add_memb( new AddMembraneMover( pose.total_residue() ) );
	add_memb->apply( pose );
	
}

/// @brief Get the name of this mover
std::string
SetupForMembraneSymmetry::get_name() const {
	return "SetupForMembraneSymmetry";
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Helper Methods - add membrane residue (pre-symmdef)
core::SSize
SetupForMembraneSymmetry::setup_asymm_membrane( Pose & pose ) {
				
	TR << "Adding a membrane residue representing the master membrane" << std::endl;
	
	using namespace protocols::membrane::geometry;
	using namespace core::kinematics;
	using namespace core::conformation;
	using namespace core::chemical;
	
	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set(
	  ChemicalManager::get_instance()->residue_type_set( pose.is_fullatom() ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
	  );
	
	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
	ResidueType const & membrane( *rsd_type_list[1] );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );
	
	// Compute residue COM of the subunit
	core::SSize rsd_com = residue_center_of_mass( pose, 1, pose.total_residue() );
	
	// Append residue by jump, don't create a new chain
	pose.append_residue_by_jump( *rsd, rsd_com, "", "", false );
	FoldTreeOP ft( new FoldTree( pose.fold_tree() ) );
	ft->reorder( rsd_com );
	pose.fold_tree( *ft );
	
	pose.fold_tree().show( std::cout );
	
	// Updating Chain Record in PDB Info
	char curr_chain = pose.pdb_info()->chain( pose.total_residue()-1 );
	pose.pdb_info()->chain( pose.total_residue(), curr_chain );
	pose.pdb_info()->obsolete(false);
	
	return pose.total_residue();
}

/// @brief Register options from the command line
void
SetupForMembraneSymmetry::register_options() {
	
	using namespace basic::options;
	option.add_relevant( OptionKeys::symmetry::symmetry_definition );

}

/// @brief Initialize option from the command line
void
SetupForMembraneSymmetry::init_from_cmd() {

	using namespace basic::options;
	
	if ( option[ OptionKeys::symmetry::symmetry_definition].user() ) {
		symmdef_file_ = option[ OptionKeys::symmetry::symmetry_definition ]();
	}
}

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetry_cc
