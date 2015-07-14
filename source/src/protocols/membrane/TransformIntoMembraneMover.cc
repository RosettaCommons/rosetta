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
#include <core/conformation/membrane/types.hh>

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

/// @brief Transform the protein into a defalt membrane
/// @details Transform the protein into default membrane, current protein
/// embedding computed from structure and spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover() :
    protocols::moves::Mover(),
    jump_( 1 ),
    mem_center_( mem_center ),
    mem_normal_( mem_normal ),
    embedding_( 0 ),
    keep_current_protein_embedding_( false )
{}

/// @brief Use custom jump to transform protein into membrane
/// @details Using user-specified jump, transform the protein into default
/// membrane, current protein computed from structure & spanfile
TransformIntoMembraneMover::TransformIntoMembraneMover( core::Size jump ) :
    protocols::moves::Mover(),
    jump_( jump ),
    mem_center_( mem_center ),
    mem_normal_( mem_normal ),
    embedding_( 0 ),
    keep_current_protein_embedding_( false )
{}
    
/// @brief Transform the protein into a defalt membrane and user specified embedding
/// @details Transform the protein with a user-defined embedding (might have
/// been optimized before) into the default membrane
TransformIntoMembraneMover::TransformIntoMembraneMover( EmbeddingDefOP new_protein_embedding ) :
	protocols::moves::Mover(),
    jump_( 1 ),
    mem_center_( mem_center ),
    mem_normal_( mem_normal ),
    embedding_( new_protein_embedding ),
    keep_current_protein_embedding_( false )
{}

/// @brief Transform the protein into user-specified mmebrane coordinates
/// @details Transform the protein with a user-defined embedding (might have
/// been optimized before and is saved as MEM) into a user-defined membrane
/// which is ARGV for the constructor
TransformIntoMembraneMover::TransformIntoMembraneMover(
		Vector center,
		Vector normal,
		bool keep_current_protein_embedding
		) :
    protocols::moves::Mover(),
    jump_( 1 ),
    mem_center_( center ),
    mem_normal_( normal ),
    embedding_( 0 ),
    keep_current_protein_embedding_( keep_current_protein_embedding )
{}


/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TransformIntoMembraneMover::TransformIntoMembraneMover( TransformIntoMembraneMover const & src ) :
	protocols::moves::Mover( src ),
    jump_( src.jump_ ),
    mem_center_( src.mem_center_ ),
    mem_normal_( src.mem_normal_ ),
    embedding_( src.embedding_ ),
    keep_current_protein_embedding_( src.keep_current_protein_embedding_ )
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
	
	// Read in membrane center & normal
	if ( tag->hasOption( "center" ) ) {
		std::string center = tag->getOption< std::string >( "center" );
		utility::vector1< std::string > str_cen = utility::string_split_multi_delim( center, ":,'`~+*&|;." );
		
		if ( str_cen.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			mem_center_.x() = std::atof( str_cen[1].c_str() );
			mem_center_.y() = std::atof( str_cen[2].c_str() );
			mem_center_.z() = std::atof( str_cen[3].c_str() );
		}
	}

	if ( tag->hasOption( "normal" ) ) {
		std::string normal = tag->getOption< std::string >( "normal" );
		utility::vector1< std::string > str_norm = utility::string_split_multi_delim( normal, ":,'`~+*&|;." );
		
		if ( str_norm.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			mem_normal_.x() = std::atof( str_norm[1].c_str() );
			mem_normal_.y() = std::atof( str_norm[2].c_str() );
			mem_normal_.z() = std::atof( str_norm[3].c_str() );
		}
	}
	
	// If option specified, use the current membrane coordinates
	// instead of the current protein embedding
	if ( tag->hasOption( "keep_current_protein_embedding" ) ) {
		keep_current_protein_embedding_ = tag->getOption< bool >( "keep_current_protein_embedding" );
	}
}

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

/// @brief Get the name of this Mover (TransformIntoMembraneMover)
std::string
TransformIntoMembraneMover::get_name() const {
	return "TransformIntoMembraneMover";
}

/// @brief Move the pose into membrane coordinate frame
void
TransformIntoMembraneMover::apply( Pose & pose ) {
	
	using namespace core::kinematics;
	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
    
    TR << "Transforming the pose into a new set of membranr coordinates" << std::endl; 
	
	// Initialize options from JD2 and this mover via commandline
	register_options();
	init_from_cmd();
	
	// Initial checks
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane protein transformation to a non membrane pose! Initialize the membrane framework using AddMembraneMover first" );
	}
	
	// Initialize jump for transformation
     if ( jump_ == 1 ) {
         jump_ = pose.conformation().membrane_info()->membrane_jump();
     }
         
	// Making a copy of the foldtree just in case
	FoldTreeOP ft_copy = FoldTreeOP( new FoldTree( pose.fold_tree() ) );
	
	// Check the setup of the foldtree is as expected
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
			utility_exit_with_message( "Membrane residue is not fixed and also not an independent branch of the foldtree, therefore a direct reorder is unsafe for reorder and to be used by this mover. Please make your foldtree smarter - see RosettaMP framework documentation" );
		}
	}
	
	// Declaring data for the previous old center / normal pair
	Vector old_center, old_normal;
	
	if ( !keep_current_protein_embedding_ ) {
	
		TR << "Using the current structure based protein embedding for transformation" << std::endl;
		// If user did not provide embedding, calculate directly from the pose
		if ( embedding_ == 0 ) {
			embedding_ = compute_structure_based_embedding( pose );
		}
	
		// Grab the current position of the protein (normal/center) prior to
		// transformation
		old_center = embedding_->center();
		old_normal = embedding_->normal();
	
	} else {
		
		TR << "Using the current membrane position for transformation" << std::endl;
		// Otherwise, read in the preivous center/normal pair from the
		// membrane residue
		old_center = pose.conformation().membrane_info()->membrane_center();
		old_normal = pose.conformation().membrane_info()->membrane_normal();
	
	}
	
	TR << "Transforming pose into a membrane coordinate frame" << std::endl;
		  
	// translate and rotate pose into membrane
	TranslationRotationMoverOP rt( new TranslationRotationMover( old_center, old_normal, mem_center_, mem_normal_, jump_ ) );
	rt->apply( pose );

	TR << "Restoring the original foldtree" << std::endl;
	pose.fold_tree( *ft_copy );
	
} // apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp, seutp options: center, normal, spanfile and
void
TransformIntoMembraneMover::register_options() {
	
	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );
	option.add_relevant( OptionKeys::mp::transform::keep_current_protein_embedding );
	
}

/// @brief Initialize Mover options from the comandline
/// @details Initialize mover settings from the commandline
/// mainly in the mp, setup group: center, normal,
/// spanfile
void
TransformIntoMembraneMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Read in Center Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::mp::setup::center ].user() ) {
		mem_center_.x() = option[ OptionKeys::mp::setup::center ]()[1];
		mem_center_.y() = option[ OptionKeys::mp::setup::center ]()[2];
		mem_center_.z() = option[ OptionKeys::mp::setup::center ]()[3];
	}
	
	// Read in Normal Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::mp::setup::normal ].user() ) {
		mem_normal_.x() = option[ OptionKeys::mp::setup::normal ]()[1];
		mem_normal_.y() = option[ OptionKeys::mp::setup::normal ]()[2];
		mem_normal_.z() = option[ OptionKeys::mp::setup::normal ]()[3];
	}
	
	// Should I use the current membrane coordinates as the initial criteria for transformation?
	if ( option[ OptionKeys::mp::transform::keep_current_protein_embedding ].user() ) {
		keep_current_protein_embedding_ = option[ OptionKeys::mp::transform::keep_current_protein_embedding ]();
	}
	
}// init from cmd


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc
