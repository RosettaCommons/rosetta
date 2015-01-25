// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/TransformInfoMembraneMoverCreator.hh
/// @brief      Transform pose into membrane coordinates (Rosetta Scripts Hook)
/// @details	Requires a MembraneInfo object with all of its associated
///				information. This can be done by calling AddMembraneMover
///				beforehand. pose.conformation().is_membrane() should always
///				return true.
///				CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///				MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc
#define INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc

// Unit Headers
#include <protocols/membrane/TransformIntoMembraneMover.hh> 
#include <protocols/membrane/TransformIntoMembraneMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/TranslationRotationMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh> 
#include <core/conformation/membrane/SpanningTopology.hh> 
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/util.hh>

// Package Headers
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh> 
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/conformation/membrane/types.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>
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

/// @brief Default Constructor
/// @details uses default membrane coords: center(0,0,0), normal(0,0,15)
TransformIntoMembraneMover::TransformIntoMembraneMover() :
	protocols::moves::Mover(),
	fullatom_( true ), 
	mem_center_( mem_center ),
	mem_normal_( mem_normal )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor - mainly for PyRosetta
/// @details user can specify desired membrane coords
///		requires center and normal of the membrane into which the protein should be
///		moved to (NOT the coordinates of the protein relative to the membrane!!!)
TransformIntoMembraneMover::TransformIntoMembraneMover(
	 Vector center,
	 Vector normal,
	 std::string spanfile
	 ) :
	 protocols::moves::Mover(),
	 fullatom_( true ), 
	 mem_center_( center ),
	 mem_normal_( normal ),
	 spanfile_( spanfile )
{
	register_options();
	init_from_cmd();
}
				
/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TransformIntoMembraneMover::TransformIntoMembraneMover( TransformIntoMembraneMover const & src ) :
	protocols::moves::Mover( src ), 
	fullatom_( src.fullatom_ ), 
	mem_center_( src.mem_center_ ),
	mem_normal_( src.mem_normal_ ),
	spanfile_( src.spanfile_ )
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
	
	// Read in boolean fulltom
	if ( tag->hasOption( "fullatom" ) ) {
		fullatom_ = tag->getOption< bool >( "fullatom" );
	}
	
	// Read in spanfile information
	if ( tag->hasOption( "spanfile" ) ) {
		spanfile_ = tag->getOption< std::string >( "spanfile" );
	}
	
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
	
	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
	
	TR << "Transforming pose into membrane coordinates" << std::endl;

	// reorder foldtree
	core::kinematics::FoldTree foldtree = pose.fold_tree();
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( foldtree );
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);
	
	// get current protein embedding from pose
	// uses however many spans are in the topology object held in MembraneInfo
	EmbeddingDefOP embedding = compute_structure_based_embedding( pose );
	Vector old_center = embedding->center();
	Vector old_normal = embedding->normal();
	
	// get membrane jump number: membrane is upstream, pose downstream
	SSize jumpnum = pose.conformation().membrane_info()->membrane_jump();
	
	// translate and rotate pose into membrane
	TranslationRotationMoverOP rt( new TranslationRotationMover( old_center, old_normal, mem_center_, mem_normal_, jumpnum ) );
	rt->apply( pose );

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// membrane_new, seutp options: center, normal, spanfile and
void
TransformIntoMembraneMover::register_options() {
	
	using namespace basic::options;

	option.add_relevant( OptionKeys::membrane_new::setup::center );
	option.add_relevant( OptionKeys::membrane_new::setup::normal );
	option.add_relevant( OptionKeys::membrane_new::setup::spanfiles );
	
}

/// @brief Initialize Mover options from the comandline
/// @details Initialize mover settings from the commandline
/// mainly in the membrane_new, setup group: center, normal,
/// spanfile
void
TransformIntoMembraneMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Set user-defined fullatom param
	if ( option[ OptionKeys::in::file::fullatom ].user() ) {
		fullatom_ = option[ OptionKeys::in::file::fullatom ]();
	}
	
	// Read in User-Provided spanfile
	if ( option[ OptionKeys::membrane_new::setup::spanfiles ].user() ) {
		spanfile_ = option[ OptionKeys::membrane_new::setup::spanfiles ]()[1];
	}

	// Read in Center Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::membrane_new::setup::center ].user() ) {
		mem_center_.x() = option[ OptionKeys::membrane_new::setup::center ]()[1];
		mem_center_.y() = option[ OptionKeys::membrane_new::setup::center ]()[2];
		mem_center_.z() = option[ OptionKeys::membrane_new::setup::center ]()[3];
	}
	
	// Read in Normal Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::membrane_new::setup::normal ].user() ) {
		mem_normal_.x() = option[ OptionKeys::membrane_new::setup::normal ]()[1];
		mem_normal_.y() = option[ OptionKeys::membrane_new::setup::normal ]()[2];
		mem_normal_.z() = option[ OptionKeys::membrane_new::setup::normal ]()[3];
	}
	
}// init from cmd


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_cc
