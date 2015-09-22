// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/conformation/membrane/MembraneInfo.cc
///
/// @brief  Information about the membrane bilayer and its relationship with the protein(s)
/// @details MembraneInfo is responsible for describing attributes of the
///    membrane bilayer *and* the position & orientation of the bilayer
///    in 3D space. All of the data members work together to accomplish
///    this representation. And the players are:
///     = A pointer to the Pose Conformation, which includes an MEM residue
///     = Topology of transmembrane spans
///     = Per-residue lipophilicity
///     = Membrane thickness & steepness - derived from chemical profiles
///
///    This object is a member of Conformation and should be accessed by
///    pose.conofrmation().membrane_info(). DO NOT access the MEM residue
///    outside of the framework!!!
///
///    Last Updated: 7/23/15
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

// Unit Headers
#include <core/conformation/membrane/MembraneInfo.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <core/conformation/membrane/MembraneParams.hh>

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

static THREAD_LOCAL basic::Tracer TR( "core.conformation.membrane.MembraneInfo" );

namespace core {
namespace conformation {
namespace membrane {

/// @brief Create a default version of MembraneInfo (DONT USE)
/// @details Initializes all data members to dummy or empty values
/// Use the fully specified constructors instead. MembraneInfo is a
/// data container but is NOT responsible for initialization.
MembraneInfo::MembraneInfo() :
	conformation_( *(new core::conformation::Conformation()) ),
	thickness_( 0 ),
	steepness_( 0 ),
	membrane_rsd_num_( 0 ),
	membrane_jump_( 0 )
{}

/// @brief Create MembraneInfo from initialized data
/// @details Creates a MembraneInfo object by linking the conformation
/// to the pose, specify the  membrane residue number, membrane jump number,
/// spanning topology object and optional lipophilicity data. Thickness and
/// steepness are currently constants
MembraneInfo::MembraneInfo(
	core::conformation::Conformation & conformation,
	core::Size membrane_pos,
	core::SSize membrane_jump,
	SpanningTopologyOP topology
) :
	conformation_( conformation ),
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_rsd_num_( membrane_pos ),
	membrane_jump_( membrane_jump ),
	lipid_acc_data_( 0 ),
	spanning_topology_( topology )
{}

/// @brief Create MembraneInfo from initialized data with lipophilicity
/// @details Creates a MembraneInfo object by linking the conformation
/// to the pose, specify the  membrane residue number, membrane jump number,
/// spanning topology object and optional lipophilicity data. Thickness and
/// steepness are currently constants
MembraneInfo::MembraneInfo(
	core::conformation::Conformation & conformation,
	core::Size membrane_pos,
	core::SSize membrane_jump,
	LipidAccInfoOP lips,
	SpanningTopologyOP topology
) :
	conformation_( conformation ),
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_rsd_num_( membrane_pos ),
	membrane_jump_( membrane_jump ),
	lipid_acc_data_( lips ),
	spanning_topology_( topology )
{}

/// @brief Create a deep copy of all data in this object.
MembraneInfo::MembraneInfo( MembraneInfo const & src ) :
	utility::pointer::ReferenceCount(),
	conformation_( src.conformation_ ),
	thickness_( src.thickness_ ),
	steepness_( src.steepness_ ),
	membrane_rsd_num_( src.membrane_rsd_num_ ),
	membrane_jump_( src.membrane_jump_ ),
	lipid_acc_data_( src.lipid_acc_data_ ),
	spanning_topology_( src.spanning_topology_ )
{}

/// @brief create a deep copy of all data in thsi object upon assignment
MembraneInfo &
MembraneInfo::operator=( MembraneInfo const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Make a deep copy of everything
	this->conformation_ = src.conformation_;
	this->thickness_ = src.thickness_;
	this->steepness_ = src.steepness_;
	this->membrane_rsd_num_ = src.membrane_rsd_num_;
	this->membrane_jump_ = src.membrane_jump_;
	this->lipid_acc_data_ = src.lipid_acc_data_;
	this->spanning_topology_ = src.spanning_topology_;

	return *this;
}

/// @brief Destructor
MembraneInfo::~MembraneInfo() {}

/// @brief Generate a string representation of information represented by ths MembraneInfo
void
MembraneInfo::show(std::ostream & output ) const {

	// Show generic membrane info
	output << "MembraneInfo: Information about this Membrane Protein" << std::endl;
	output << "Membrane Residue Num: " << membrane_rsd_num_ << std::endl;
	output << "Membrane Fold Tree Jump: " << membrane_jump_ << std::endl;
	output << "Membrane Thickness: " << thickness_ << std::endl;
	output << "Membrane Steepness: " << steepness_ << std::endl;
	output << "Membrane Spanning Topology " << std::endl;

	// Grab membrane center/normal
	core::Vector center( membrane_center() );
	core::Vector normal( membrane_normal() );

	// Show Current Membrane Position
	output << "Membrane Center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;
	output << "Membrane Normal: " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;

	// SHow spanning topology object
	spanning_topology_->show();

	// TODO
	// Skipping lips for now, will go back to it

}

// Chemical Information about this Membrane

/// @brief Effective thickness of the membrane
/// @details For IMM = default is 15. Otherwise, defined as the distance between the g3p
/// linkers in the phospholipid
core::Real
MembraneInfo::membrane_thickness() const {
	return thickness_;
}

/// @brief Steepness of hydrophobic -> hydrophillic transition
/// @details For IMM - default is 10. Otherwise, caluclated as the slope in the polarity
/// gradient from low to high charge density from the chemical profile.
core::Real
MembraneInfo::membrane_steepness() const {
	return steepness_;
}

// membrane position & orientation

/// @brief Membrane center
/// @details Returns the xyzVector describing the center of the membrane
/// This is the same as the MPct atom of the membrane (MEM) residue.
core::Vector
MembraneInfo::membrane_center() const  {

	return conformation_.residue( membrane_rsd_num() ).xyz( membrane::center );
}

/// @brief Membrane normal
/// @details Returns the membrane normal, which describes the membrane
/// orientation. This is the same as the xyzVector in the MPnm atom
/// in the membrane residue.
core::Vector
MembraneInfo::membrane_normal() const {

	core::Vector normal_tracked = conformation_.residue( membrane_rsd_num() ).xyz( membrane::normal );
	core::Vector normal = normal_tracked - membrane_center();

	return normal.normalize();
}

/// @brief Is residue in the membrane? Takes CA coordinate
/// @details Uses the thickness stored in MembraneInfon and the residue_z_position
bool
MembraneInfo::in_membrane( core::Size resnum ) const {

	bool in_mem( false );

	if ( residue_z_position( resnum ) >= -membrane_thickness() &&
			residue_z_position( resnum ) <= membrane_thickness() ) {
		in_mem = true;
	}

	return in_mem;
} // in membrane?

/// @brief Compute residue position relative to membrane normal
/// @details Calculate the z coordinate of the residue, projected onto
/// the membrane normal axis. Objective is to maintain correct coordinates
/// in relative coordinate frame.
core::Real
MembraneInfo::residue_z_position( core::Size resnum ) const {

	// Compute z_position
	core::Vector const & xyz( conformation_.residue( resnum ).atom( "CA" ).xyz() );

	// membrane normal is normalized to 15, that's why dividing it here by 15
	core::Vector normalized_to_1( membrane_normal() );
	normalized_to_1.normalize();
	core::Real result = dot( xyz - membrane_center(), normalized_to_1 );
	return result;
}

/// @brief Compute atom position relative to membrane normal
/// @details Calculate the z coordinate of the atom, projected onto
/// the membrane normal axis. Objective is to maintain correct coordinates
/// in relative coordinate frame.
core::Real
MembraneInfo::atom_z_position( core::Size resnum, core::Size atomnum ) const {

	// Compute z_position
	core::Vector const & xyz( conformation_.residue( resnum ).atom( atomnum ).xyz() );

	// membrane normal is normalized to 15, that's why dividing it here by 15
	core::Vector normalized_to_1( membrane_normal() );
	normalized_to_1.normalize();
	core::Real result = dot( xyz - membrane_center(), normalized_to_1 );
	return result;
}

/// @brief Sequence position of the membrane residue
/// @details Return the residue number of MEM (rsd.seqpos()) in the pose
core::Size
MembraneInfo::membrane_rsd_num() const {
	return membrane_rsd_num_;
}

/// @brief Indeitifier for the membrane jump
/// @details Returns an integer (core::Size) denoting the jump number in the foldtree
/// representing the jump relating the membrane residue to the rest of the molecule
core::SSize
MembraneInfo::membrane_jump() const { return membrane_jump_; }

/// @brief Allow a protocol to set a new jump number for the membrane jump
/// @details Set the membrane jump number (core::SSize)
void
MembraneInfo::set_membrane_jump( core::SSize jumpnum ) {

	TR << "Setting a new membrane jump number in MembraneInfo to " << jumpnum << "." << std::endl;
	TR << "Use with caution!" << std::endl;
	membrane_jump_ = jumpnum;
}

/// @brief Somewhat weak check that a membrane foldtree is valid. Use checks in
/// protocols/membrane/util.hh instead!
bool
MembraneInfo::check_membrane_fold_tree( core::kinematics::FoldTree const & ft_in ) const {

	// Check regular fold tree
	if ( ! ft_in.check_fold_tree() ) {
		return false;
	}

	// Check memrbane residue is a jump point
	if ( ! ft_in.is_jump_point( membrane_rsd_num_ ) ) {
		return false;
	}

	// Otherwise, looks reasonable!
	return true;
}

// topology of TM spans and lipophilicity

/// @brief Transmembrane spaning topology
/// @details Return a SpanningTopology object, which includes a
/// list of Span objects, describing the start and end sequence
/// positions of each transmembrane span
SpanningTopologyOP
MembraneInfo::spanning_topology() const {
	return spanning_topology_;
}

/// @brief Does this MembraneInfo includes lipophilicity information?
bool
MembraneInfo::include_lips() const {

	// If lips initialized, return true
	if ( lipid_acc_data_ != 0 ) {
		return true;
	}

	return false;
}

/// @brief Per-residue lipophilicity (probability of exposure to lipid)
/// @details Returns a LipidAccInfo describing per residue probability
/// of exposure to lipid. Data calcualted via the run_lips.pl script
/// and provided by the user on the commandline if applicable
LipidAccInfoOP
MembraneInfo::lipid_acc_data() const {
	return lipid_acc_data_;
}


/// @brief Show MembraneInfo method for pyrosetta
std::ostream & operator << ( std::ostream & os, MembraneInfo const & mem_info )
{

	// Grab membrane position from the pose
	core::Vector center( mem_info.membrane_center() );
	core::Vector normal( mem_info.membrane_normal() );

	// Grab a const version of spanning topology
	os << "Membrane residue located at position " << mem_info.membrane_rsd_num();
	os << "Membrane Position: " << "center = " << center.to_string() << "; normal = " << normal.to_string() << std::endl;
	os << "Number of transmembrane spans: " << mem_info.spanning_topology()->nspans() << std::endl;

	return os;
}

} // membrane
} // conformation
} // core
