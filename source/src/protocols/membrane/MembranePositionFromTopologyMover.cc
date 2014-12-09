// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	    protocols/membrane/MembranePositionFromTopologyMover.cc
///
/// @brief      Computes and sets the initial position of the membrane
/// @details	Computes and sets the initial position of the membrane from
///				sequence or structure (can be specified by the user at construction
///				or as a setup cmd flag).
///				CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///
///				NOTE: Requires a membrane pose!
///				NOTE: sequence not yet implemented
///				Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc

// Unit Headers
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMoverCreator.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Residue.hh>
#include <protocols/membrane/geometry/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.membrane.MembranePositionMoverFromTopologyMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;

/// @brief Defualt Constructor
/// @details Compute the embedding of the pose based on xyz coordinates
/// and spanning topology provided in MembraneInfo
MembranePositionFromTopologyMover::MembranePositionFromTopologyMover() :
	protocols::moves::Mover(),
	structure_based_( true )
{}

/// @brief Custom Constructor - for Pyrosetta
/// @details Compute the embedding of the pose - if structure_based is
/// true do this based on xyz coordinates. If structure_based is false,
/// compute based on sequence.
MembranePositionFromTopologyMover::MembranePositionFromTopologyMover( bool structure_based ) :
	protocols::moves::Mover(),
	structure_based_( structure_based )
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover
MembranePositionFromTopologyMover::MembranePositionFromTopologyMover( MembranePositionFromTopologyMover const & src ) :
	protocols::moves::Mover( src ),
	structure_based_( src.structure_based_ )
{}

/// @brief Destructor
MembranePositionFromTopologyMover::~MembranePositionFromTopologyMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MembranePositionFromTopologyMover::clone() const {
	return ( protocols::moves::MoverOP( new MembranePositionFromTopologyMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MembranePositionFromTopologyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MembranePositionFromTopologyMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MembranePositionFromTopologyMover::parse_my_tag(
   utility::tag::TagCOP tag,
   basic::datacache::DataMap &,
   protocols::filters::Filters_map const &,
   protocols::moves::Movers_map const &,
   core::pose::Pose const &
   ) {

	if ( tag->hasOption( "structure_based" ) ) {
		structure_based_ = tag->getOption< bool >("structure_based");
	}

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MembranePositionFromTopologyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MembranePositionFromTopologyMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MembranePositionFromTopologyMoverCreator::keyname() const {
	return MembranePositionFromTopologyMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MembranePositionFromTopologyMoverCreator::mover_name() {
	return "MembranePositionFromTopologyMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Update Membrane position in pose
/// @details Compute membrane posiiton based on sequence or structure
/// and then call pose.update_membrane_position() to update the membrane position
void
MembranePositionFromTopologyMover::apply( Pose & pose ) {

	// Check pose is a membrane pose
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Warning! Pose is not a membrane pose. Cannot perform mpframework operation on a non membrane pose~" );
	}

	TR << "Computing initial membrane position from structure..." << std::endl;

	using namespace core;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	// Initialize starting vectors
	Vector normal( 0, 0, 0 );
	Vector center( 0, 0, 0 );

	// Compute position from pose
	compute_structure_based_membrane_position( pose, center, normal );

	// Update membrane position - shift normal along center
	pose.conformation().update_membrane_position( center, center + normal );

	TR << "Done" << std::endl;

}

/// @brief Get the name of this mover
std::string
MembranePositionFromTopologyMover::get_name() const {
	return "MembranePositionFromTopologyMover";
}

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_cc
