// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/conformation/membrane/MembraneInfo.cc
///
/// @brief		MembraneInfo - Membrane Pose Definition Object
/// @details	The membrane info object is responsible for storing non-coordinate
///				derived information, extending the traditional definiton of a pose
///				to describe a membrane protein in Rosetta. This information includes:
///				 - the resiue number of the membrane virtual residue containing
///					the posiiton of the membrane
///				 - membrane spanning topology
///				 - membrane lipophilicity info (user-specified)
///
///				This object belongs to the conformation and should be accessed
///				via pose.conformation().membrane().
///
///	     		Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <core/conformation/membrane/MembraneInfo.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

#include <core/conformation/membrane/MembranePlanes.hh>

// Package Headers
#include <core/conformation/Conformation.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/id/AtomID.hh>

#include <core/conformation/membrane/Exceptions.hh>

#include <core/conformation/util.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

static basic::Tracer TR("core.membrane.MembraneInfo");

namespace core {
namespace conformation {
namespace membrane {

using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace core::kinematics;
	
////////////////////
/// Constructors ///
////////////////////

/// @brief Default Constructor
/// @details Creates a default copy of the membrane info object
/// with the membrane virtual at nres = 1, thicnkess = 15A, steepenss
/// of fullatom membrane transition = 10A and no spanning topology or
/// lips info defined
MembraneInfo::MembraneInfo() :
	utility::pointer::ReferenceCount(),
	thickness_( 15.0 ),
	steepness_( 10.0 ),
	membrane_rsd_num_( 1 ),
	membrane_jump_( 2 ),
	view_in_pymol_( false ),
	plane_info_( 0 )
{}

/// @brief Custom Constructor - Membrane pos & topology
/// @details Creates a default copy of the membrane info object
/// with the membrane virtual at nres = membrane_pos, thicnkess = 15A, steepenss
/// of fullatom membrane transition = 10A and spanning topology defined
/// by API provided vector1 of spanning topology obejcts (by chain)
MembraneInfo::MembraneInfo(
	core::Size membrane_pos,
	SpanningTopologyOP topology,
	core::SSize membrane_jump,
	bool view_in_pymol
	) :
	utility::pointer::ReferenceCount(),
	thickness_( 15.0 ),
	steepness_( 10.0 ),
	membrane_rsd_num_( membrane_pos ),
	membrane_jump_( membrane_jump ),
	spanning_topology_( topology ),
	view_in_pymol_( view_in_pymol ),
	plane_info_( 0 )
{}

/// @brief Custom Constructor - Membrane pos, topology & lips
/// @details Creates a default copy of the membrane info object
/// with the membrane virtual at nres = membrane_pos, thicnkess = 15A, steepenss
/// of fullatom membrane transition = 10A, spanning topology defined
/// by API provided vector1 of spanning topology obejcts (by chain), and
/// lipid accessibility defined by a vector1 of lipid acc objects.
MembraneInfo::MembraneInfo(
	core::Size membrane_pos,
	SpanningTopologyOP topology,
	LipidAccInfoOP lips,
	core::SSize membrane_jump,
	bool view_in_pymol
	) :
	utility::pointer::ReferenceCount(),
	thickness_( 15.0 ),
	steepness_( 10.0 ),
	membrane_rsd_num_( membrane_pos ),
	membrane_jump_( membrane_jump ),
	lipid_acc_data_( lips ),
	spanning_topology_( topology ),
	view_in_pymol_( view_in_pymol ),
	plane_info_( 0 )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this object
MembraneInfo::MembraneInfo( MembraneInfo const & src ) :
	utility::pointer::ReferenceCount(),
	thickness_( src.thickness_ ),
	steepness_( src.steepness_ ),
	membrane_rsd_num_( src.membrane_rsd_num_ ),
	membrane_jump_( src.membrane_jump_ ),
	lipid_acc_data_( src.lipid_acc_data_ ),
	spanning_topology_( src.spanning_topology_ ),
	view_in_pymol_( src.view_in_pymol_ ),
	plane_info_( src.plane_info_ )
{}

/// @brief Assignment Operator
/// @details Create a deep copy of this object while overloading the assignemnt
/// operator "="
MembraneInfo &
MembraneInfo::operator=( MembraneInfo const & src ) {
	
	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new MembraneInfo( *this ) );
}

/// @brief Destructor
MembraneInfo::~MembraneInfo() {}

/// @brief  Generate string representation of MembraneInfo for debugging purposes.
void
MembraneInfo::show(std::ostream & output ) const {
	
	// Show generic membrane info
	output << "MembraneInfo: Information about this Membrane Protein" << std::endl;
	output << "Membrane Residue Num: " << membrane_rsd_num_ << std::endl;
	output << "Membrane Fold Tree Jump: " << membrane_jump_ << std::endl;
	output << "Membrane Thicnkess: " << thickness_ << std::endl;
	output << "Membrane Steepness: " << steepness_ << std::endl;
	output << "Membrane Spanning Topology " << std::endl;

	// SHow spanning topology object
	spanning_topology_->show();
	
	// Skipping lips for now, will go back to it
	
}


////////////////////////////
/// Membrane Data Access ///
////////////////////////////

/// @brief Return position of membrane residue
/// @details Return the residue number of the membrane virtual
/// residue in the pose
core::Size
MembraneInfo::membrane_rsd_num() const { return membrane_rsd_num_; }
	
/// @brief Return membrane thickness
/// @details Return the membrane thicnkess used by the
/// fullatom energy method (centroid is hard coded for now
core::Real
MembraneInfo::membrane_thickness() const {
	return thickness_;
}

/// @brief Return the membrane transition steepness
/// @details Return the steepness betwen isotropic and ansitropic
/// layers of the membrane used by the fullatom energy methods
core::Real
MembraneInfo::membrane_steepness() const {
	return steepness_;
}

/// @brief Return a list of membrane spanning topology objects
/// @details Return a vector1 of spanning topology objects defining
/// the starting and ending position of membrane spans per chain.
SpanningTopologyOP
MembraneInfo::spanning_topology(){
	return spanning_topology_;
}

/// @brief Return a list of lipid accessibility objects
/// @details Return a vector1 of lipid accessibility info objects
/// describing lipid exposre of individual residues in the pose
LipidAccInfoOP
MembraneInfo::lipid_acc_data() const {
	return lipid_acc_data_;
}

////////////////////////////////////
/// Membrane Base Fold Tree Info ///
////////////////////////////////////

/// @brief Get the number of the membrane jump
/// @details Get a core::SSize (int) denoting the number of the fold tree jump
/// relating the membrane residue to the rest of the pose
core::SSize
MembraneInfo::membrane_jump() const { return membrane_jump_; }

/// @brief Check membrane fold tree
/// @details Check that the membrane jump num is a jump point and located at the root
/// of the fold tree in addition to maintaining a reasonable fold tree
bool
MembraneInfo::check_membrane_fold_tree( FoldTree const & ft_in ) const {

	// Check regular fold tree
	if (! ft_in.check_fold_tree() ) {
		return false;
	}
	
	// Check memrbane residue is a jump point
	if (! ft_in.is_jump_point( membrane_rsd_num_ ) ) {
		return false;
	}
	
	// Otherwise, looks reasonable!
	return true;
	
}

/////////////////////
/// Visualizaiton ///
/////////////////////

/// @brief Check membrane planes are initialized for visualization
bool
MembraneInfo::view_in_pymol() const {
	return plane_info_ != 0;
}

/// @brief Setup Planes Info Object
void
MembraneInfo::setup_plane_visualization(
	utility::vector1< Size > top_points,
	utility::vector1< Size > bottom_points
	) {
	
	using namespace core::conformation::membrane;
	plane_info_ = new MembranePlanes( top_points, bottom_points );
}

/// @brief Membrane Planes Points
/// @details Return object containing membrane planes info. Initialized at setup
MembranePlanesOP
MembraneInfo::plane_info() {
	return plane_info_;
}
	
} // membrane
} // conformation
} // core
