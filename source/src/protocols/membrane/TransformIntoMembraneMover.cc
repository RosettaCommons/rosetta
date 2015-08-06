// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///	@file		protocols/membrane/TransformIntoMembraneMover.cc
/// @brief		Transform a pose into a membrane coordinate frame
/// @author		Rebecca Faye Alford (rfalford12@gmail.com)
/// @author     JKLeman (julia.koehler1982@gmail.com)
///				CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///				MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
/// Last Modified: 6/11/15
/// #RosettaMPMover

#ifndef INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc
#define INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc

// Unit Headers
#include <protocols/membrane/TransformIntoMembraneMover.hh> 
#include <protocols/membrane/TransformIntoMembraneMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/util.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh> 
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.TransformIntoMembraneMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
		
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Transform the protein into a default membrane
/// @details Transform the protein into default membrane, protein
/// embedding computed from structure and spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover() :
	jump_( 0 ), 
	new_mem_cntr_( 0, 0, 0 ), 
	new_mem_norm_( 0, 0, 1 ), 
	current_embedding_( new EmbeddingDef( new_mem_norm_, new_mem_cntr_ ) ), 
	use_default_membrane_( false ), 
	user_defined_membrane_( false )
{}

// TODO: use and test this constructor
// Use custom jump to transform protein into membrane
// Using user-specified jump, transform the downstream partner
// into default membrane, partner embedding computed from structure & spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover( core::Size jump ) :
	jump_( jump ), 
	new_mem_cntr_( 0, 0, 0 ), 
	new_mem_norm_( 0, 0, 1 ), 
	current_embedding_( new EmbeddingDef( new_mem_norm_, new_mem_cntr_ ) ),  
	use_default_membrane_( false ), 
	user_defined_membrane_( false )
{}

/// @brief Transform the protein with a user-specified protein embedding into
/// a default membrane
/// @details Transform the protein with a user-defined embedding (might have
/// been optimized before) into the default membrane
TransformIntoMembraneMover::TransformIntoMembraneMover( EmbeddingDefOP current_embedding ) : 
	jump_( 0 ), 
	new_mem_cntr_( 0, 0, 0 ), 
	new_mem_norm_( 0, 0, 1 ), 
	current_embedding_( current_embedding ), 
	use_default_membrane_( false ), 
	user_defined_membrane_( false )
{}

/// @brief Transform the protein into user-specified membrane coordinates
/// @details Transform the protein into a user-defined membrane, protein
/// embedding is computed from structure and spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover( Vector new_mem_cntr, Vector new_mem_norm ) :
	jump_( 0 ), 
	new_mem_cntr_( new_mem_cntr ), 
	new_mem_norm_( new_mem_norm ), 
	current_embedding_( 0 ),
	use_default_membrane_( false ), 
	user_defined_membrane_( true )
{
	core::Vector zero_vec( 0, 0, 0 ); 
	current_embedding_ = EmbeddingDefOP( new EmbeddingDef( zero_vec, zero_vec ) ); 
}

/// @brief Transform the protein into user-specified membrane coordinates
/// @details Transform the protein into a user-defined membrane, protein
/// embedding is computed from structure and spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover( EmbeddingDefOP current_embedding, Vector new_mem_cntr, Vector new_mem_norm ) :
	jump_( 0 ), 
	new_mem_cntr_( new_mem_cntr ), 
	new_mem_norm_( new_mem_norm ), 
	current_embedding_( current_embedding ), 
	use_default_membrane_( false ), 
	user_defined_membrane_( true )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TransformIntoMembraneMover::TransformIntoMembraneMover( TransformIntoMembraneMover const & src ) :
	protocols::moves::Mover( src ),
	jump_( src.jump_ ),
	new_mem_cntr_( src.new_mem_cntr_ ),
	new_mem_norm_( src.new_mem_norm_ ),
	current_embedding_( src.current_embedding_ ),
	use_default_membrane_( src.use_default_membrane_ )
{}

/// @brief Destructor
TransformIntoMembraneMover::~TransformIntoMembraneMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
TransformIntoMembraneMover::clone() const {
	return ( protocols::moves::MoverOP( new TransformIntoMembraneMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
TransformIntoMembraneMover::fresh_instance() const {
	return protocols::moves::MoverOP( new TransformIntoMembraneMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
TransformIntoMembraneMover::parse_my_tag(
	 utility::tag::TagCOP tag,
	 basic::datacache::DataMap &,
	 protocols::filters::Filters_map const &,
	 protocols::moves::Movers_map const &,
	 core::pose::Pose const &
	 ) {

	// Read in jump option
	if ( tag->hasOption( "jump") ) {
		jump_ = tag->getOption< core::Size >( "jump" ); 
	}

	// Read in option use default membrane
	if ( tag->hasOption( "use_default_membrane" ) ) {
		use_default_membrane_ = tag->getOption< bool >( "use_default_membrane" ); 
	}

	// User defined membrane
	if ( tag->hasOption( "user_defined_membrane" ) ) {
		user_defined_membrane_ = tag->getOption< bool >( "user_defined_membrane" ); 
	}
	
	// Reading in center/normal pair 
	read_center_normal_from_tag( new_mem_cntr_, new_mem_norm_, tag ); 
	
} // parse my tag

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
TransformIntoMembraneMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TransformIntoMembraneMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
TransformIntoMembraneMoverCreator::keyname() const {
	return TransformIntoMembraneMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
TransformIntoMembraneMoverCreator::mover_name() {
	return "TransformIntoMembraneMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Use the default membrane (cntr 0,0,0 and normal 0,0,1) instead
///			of the membrane from the MEM coordinates stored in MembraneInfo
void TransformIntoMembraneMover::use_default_membrane( bool truefalse ) {
	use_default_membrane_ = truefalse;
}

/// @brief Get the name of this Mover (TransformIntoMembraneMover)
std::string
TransformIntoMembraneMover::get_name() const {
	return "TransformIntoMembraneMover";
}

/// @brief Move the pose into membrane coordinate frame
void
TransformIntoMembraneMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace core::kinematics;
	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
    
    TR << "Transforming the pose into new membrane coordinates" << std::endl;
	
	// Initialize options from JD2 and this mover via commandline
	register_options();
	init_from_cmd();
	
	// Initial checks
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane protein transformation to a non membrane pose! Initialize the membrane framework using AddMembraneMover first." );
	}
	
	// starting point is default membrane, this is overwritten here with the MEM
	// info from the pose
	if ( use_default_membrane_ == false && user_defined_membrane_ == false ) {

		new_mem_cntr_ = pose.conformation().membrane_info()->membrane_center();
		new_mem_norm_ = pose.conformation().membrane_info()->membrane_normal();
	}

	// for user-defined membrane, set membrane position
	if ( user_defined_membrane_ == true ){

		SetMembranePositionMoverOP setmem( new SetMembranePositionMover( new_mem_cntr_, new_mem_norm_ ) );
		setmem->apply( pose );
	}
	
	// Initialize jump for transformation
	// initial jump should be zero so we now it's bogus; the user can use any
	// other jump
     if ( jump_ == 0 ) {
         jump_ = pose.conformation().membrane_info()->membrane_jump();
     }
         
	// Making a copy of the foldtree to reset it to after the mover
	FoldTreeOP orig_ft = FoldTreeOP( new FoldTree( pose.fold_tree() ) );
	
	// Check that the membrane is at the root of the foldtree
	if ( !is_membrane_fixed( pose ) ) {
		if ( is_membrane_moveable_by_itself( pose ) ) {
			
			TR << "Membrane is currently not fixed, but is an independently moveable branch of the foldtree. Temporarily adjusting " << std::endl;
			TR << "such that the membrane is the root then will set the original foldtree back afterward" << std::endl;
			
			FoldTreeOP curr_ft = FoldTreeOP( new FoldTree( pose.fold_tree() ) );
			Size membrane_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );
			curr_ft->reorder( membrane_rsd );
			pose.fold_tree( *curr_ft );
			
			TR << "Foldtree Reorder complete" << std::endl;
			pose.fold_tree().show( std::cout );
			
		} else {
			utility_exit_with_message( "Membrane residue is not fixed and also not an independent branch of the foldtree, therefore a direct reorder is unsafe to be used by this mover. Please make your foldtree smarter - see RosettaMP framework documentation" );
		}
	}
	
	// recompute embedding from pose and topology if it hasn't been set
	if ( current_embedding_->normal().length() == 0 ){
		current_embedding_ = compute_structure_based_embedding( pose );
	}
	
	TR.Debug << "...transforming..." << std::endl;
		  
	// translate and rotate pose into membrane
	TranslationRotationMoverOP rt( new TranslationRotationMover( current_embedding_->center(), current_embedding_->normal(), new_mem_cntr_, new_mem_norm_, jump_ ) );
	rt->apply( pose );
	
	// reset foldtree
	TR << "Restoring the original foldtree" << std::endl;
	pose.fold_tree( *orig_ft );
	
} // apply

/////////////////////
/// Setup Methods ///
/////////////////////


/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp, seutp options: center, normal, spanfile and
void TransformIntoMembraneMover::register_options() {
	
	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );
	
}

/// @brief Initialize Mover options from the comandline
/// @details Initialize mover settings from the commandline
/// mainly in the mp, setup group: center, normal,
/// spanfile
void TransformIntoMembraneMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Read center & normal options from the commnadline
	read_center_normal_from_cmd( new_mem_cntr_, new_mem_norm_ ); 

}// init from cmd


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc
