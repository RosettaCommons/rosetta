// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	    protocols/membrane/VisualizeMembraneMover.hh
///
/// @brief      Visualize Membrane Planes with Virtual Residues
/// @details    Add a set of virtual residues as a third chain to the
///				membrane pose. This tool is strictly for visualization of
///				the implicit membrane and should not be present in modeling.
///				Last Modified: 6/19/14
///
/// @author		Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_cc
#define INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_cc

// Unit Headers
#include <protocols/membrane/visualize/VisualizeMembraneMover.hh>
#include <protocols/membrane/visualize/VisualizeMembraneMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

static thread_local basic::Tracer TR( "protocols.membrane.visualize.VisualizeMembraneMover" );

namespace protocols {
namespace membrane {
namespace visualize {

using namespace core;
using namespace core::pose;

////////////////////
/// Constructors ///
////////////////////

/// @brief	  Defualt Constructor
/// @details  Construct membrane residues with spacing = 5,
///	          width = 100
VisualizeMembraneMover::VisualizeMembraneMover() :
	protocols::moves::Mover(),
	spacing_( 5 ),
	width_( 100 ),
	thickness_( 12.5 )
{
	register_options();
	init_from_cmd();
}

/// @brief    Construct with User specified spacing & width
/// @details  Construct membranes with a given spacing and
///			  width in angstroms
VisualizeMembraneMover::VisualizeMembraneMover( Real spacing, Real width, Real thickness ) :
	protocols::moves::Mover(),
	spacing_( spacing ),
	width_( width ),
	thickness_( thickness )
{}

/// @brief Copy Constructor
/// @details Creates a deep copy of the visualize membrane mover class
VisualizeMembraneMover::VisualizeMembraneMover( VisualizeMembraneMover const & src ) :
	protocols::moves::Mover( src ),
	spacing_( src.spacing_ ),
	width_( src.width_ ),
	thickness_( src.thickness_ )
{}

/// @brief Assignment Operator
/// @details Overloads "=" assignemnt for deep copying
VisualizeMembraneMover &
VisualizeMembraneMover::operator=( VisualizeMembraneMover const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new VisualizeMembraneMover( *this ) );
	
}

/// @brief Destructor
VisualizeMembraneMover::~VisualizeMembraneMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
VisualizeMembraneMover::clone() const {
	return ( new VisualizeMembraneMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
VisualizeMembraneMover::fresh_instance() const {
	return new VisualizeMembraneMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
VisualizeMembraneMover::parse_my_tag(
  utility::tag::TagCOP tag,
  basic::datacache::DataMap &,
  protocols::filters::Filters_map const &,
  protocols::moves::Movers_map const &,
  core::pose::Pose const &
  ) {
  
	// Read in spacing option
	if ( tag->hasOption( "spacing" ) ) {
		spacing_ = tag->getOption< Real >( "spacing" );
	}
	
	// Read in Width Option
	if ( tag->hasOption( "width" ) ) {
		width_ = tag->getOption< Real >( "width" );
	}
  
	// Read in thickness option
	if ( tag->hasOption( "thickness" ) ) {
		thickness_ = tag->getOption< Real >( "thickness" );
	}
	
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
VisualizeMembraneMoverCreator::create_mover() const {
	return new VisualizeMembraneMover;
}
	
/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
VisualizeMembraneMoverCreator::keyname() const {
	return VisualizeMembraneMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
VisualizeMembraneMoverCreator::mover_name() {
	return "VisualizeMembraneMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

void
VisualizeMembraneMover::apply( Pose & pose ) {

	using namespace core::kinematics;
	
	TR << "Adding membrane planes represented as virtual residues to pose" << std::endl;

	// Check that I am a membrane pose
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Cannot visualize a non-membrane pose!");
	}
	
	// Determine if the pose is fullatom (needed for residue typesets)
	bool fullatom = pose.is_fullatom();
	
	utility::vector1< ResidueOP > membrane_residues;
	
	// IF first residue, create new chain
	bool is_first( true );
	
	// Compute residues in the top plane
	for ( Real i = -width_; i <= width_; i += spacing_ ) {
		for ( Real j = -width_; j <= width_; j += spacing_ ) {

			// Create initial position, apply rotation
			Vector pt_upper( i, j, thickness_ );
			Vector pt_lower( i, j, -thickness_ );
	
			// Create a new residue, append to list
			ResidueOP membrane_upper = create_membrane_virtual( pt_upper, fullatom );
			ResidueOP membrane_lower = create_membrane_virtual( pt_lower, fullatom );
			
			// Append residue to list
			membrane_residues.push_back( membrane_upper );
			membrane_residues.push_back( membrane_lower );
			
		}
	}
	
	// Append Residues to the pose
	for ( Size i = 1; i <= membrane_residues.size(); ++i ) {
		if ( is_first ) {
			pose.append_residue_by_jump( *membrane_residues[i], pose.total_residue(), "", "", true );
			is_first = false;
		} else {
			pose.append_residue_by_jump( *membrane_residues[i], pose.total_residue(), "", "", false );
		}
	}
}

std::string
VisualizeMembraneMover::get_name() const {
	return "VisualizeMembraneMover";
}

//////////////////////
/// Helper Methods ///
//////////////////////

/// @brief Register Options with JD2
void
VisualizeMembraneMover::register_options() {
	
	using namespace basic::options;
	
	option.add_relevant( OptionKeys::membrane_new::visualize::spacing );
	option.add_relevant( OptionKeys::membrane_new::visualize::width );
	
}
	
/// @brief Initialize Options from the Command Line
/// @details Options allowed are vrt spacing and plane width
void
VisualizeMembraneMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Specify spacing option
	if ( option[ OptionKeys::membrane_new::visualize::spacing ].user() ) {
		spacing_ = option[ OptionKeys::membrane_new::visualize::spacing ]();
	}
	
	// Specify width option
	if ( option[ OptionKeys::membrane_new::visualize::width ].user() ) {
		width_ = option[ OptionKeys::membrane_new::visualize::width ]();
	}
	
	// Specify Thickness Option
	if ( option[ OptionKeys::membrane_new::visualize::thickness ].user() ) {
		thickness_ = option[ OptionKeys::membrane_new::visualize::thickness ]();
	}
	
}

ResidueOP
VisualizeMembraneMover::create_membrane_virtual( Vector pos, bool fullatom ) {
	
	using namespace core::conformation;
	using namespace core::chemical;
	
	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCAP const & residue_set(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
		);
	
	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
	ResidueType const & membrane( *rsd_type_list[1] );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );
	
	// setup membrane thicnkess
	Vector tk(0, 1.0, 30.0);
	
	// Fill the Residue with normal/center info
	rsd->set_xyz(1, pos);
	rsd->set_xyz(3, tk);

	return rsd;
}


} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeMembraneMover_cc
