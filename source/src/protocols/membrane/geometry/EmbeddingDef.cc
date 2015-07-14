// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/geometry/EmbeddingDef.cc
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author		Rebecca Faye Alford (rfalford12@gmail.com)
///	Last Modified: 6/11/15
/// #RosettaMPData

// Unit headers
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/conformation/membrane/types.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <iostream>

static thread_local basic::Tracer TR( "protocols.membrane.geometry.EmbeddingDef" );

namespace protocols {
namespace membrane {
namespace geometry {

using namespace core;
using namespace core::conformation::membrane;
using namespace protocols::membrane; 

/// @brief Default Constructor
EmbeddingDef::EmbeddingDef() :
	utility::pointer::ReferenceCount()
{
	center_ = mem_center;
	normal_ = mem_normal;
}

/// @brief Standard Constructor
EmbeddingDef::EmbeddingDef( core::Vector const center, core::Vector const normal ) :
    utility::pointer::ReferenceCount(),
    center_( center ),
	normal_( normal )
{}

/// @brief Constructor from pose, two residue numbers, and bool if in positive z-direction
EmbeddingDef::EmbeddingDef( core::pose::Pose const & pose, core::Size start, core::Size end, bool pos_z ) :
utility::pointer::ReferenceCount()
{
	center_ = mem_center;
	normal_ = mem_normal;
	
	if ( pos_z == false ) {
		from_span( pose, start, end );
	}
	else {
		from_span_positive_z(pose, start, end);
	}
}

/// @brief Copy Constructor
EmbeddingDef::EmbeddingDef( EmbeddingDef const & src ) :
    utility::pointer::ReferenceCount(),
    center_( src.center_ ),
	normal_( src.normal_ )
{}

/// @brief Assignment Operator
EmbeddingDef & EmbeddingDef::operator = ( EmbeddingDef const & src ) {

	// Abort self-assignment.
	if ( this == &src ) { return *this; }

	// Deep Copy of the data
	this->center_ = src.center_;
	this->normal_ = src.normal_;

	return *this;
}

/// @brief Destructor
EmbeddingDef::~EmbeddingDef() {}

/// @brief Standard Rosetta Show Method for Debugging
void
EmbeddingDef::show( std::ostream & ) const {
	TR << "Embedding: center: " << center_.to_string() << ", normal: " << normal_.to_string() << std::endl;
}

/// @brief Check reasonable range of vectors in embedding object
void EmbeddingDef::check_range() const{
    check_vector( normal_ );
    check_vector( center_ );
}

/// @brief Translate by center and normal
void EmbeddingDef::translate_by( EmbeddingDef const & translation ){
    center_ += translation.center();
    normal_ += translation.normal();
    check_range();
}

/// @brief Set New Normal
void EmbeddingDef::set_normal( core::Vector normal ) {
    normal_ = normal;
    check_range();
}

/// @brief Set New Center
void EmbeddingDef::set_center( core::Vector center ) {
    center_ = center;
    check_range();
}

/// @brief Invert normal
void EmbeddingDef::invert() {
	if ( normal_.length() > 0 ) {
		normal_.negate();
	}
}

/// @brief Access Normal Param
core::Vector EmbeddingDef::normal() const {
	return normal_;
}

/// @brief Access center param
core::Vector EmbeddingDef::center() const {
	return center_;
}

/// @brief Check Object Equality
bool EmbeddingDef::equals( EmbeddingDef const & other ) const {
    
    if ( normal_ != other.normal() ) return false;
    if ( center_ != other.center() ) return false;
    
    return true;
}

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this.
void EmbeddingDef::from_span( core::pose::Pose const & pose, core::Size const start, core::Size const end ) {

    TR << "Computing membrane embedding from TMspan " << start << " to " << end << std::endl;
    
	// get CA atom positions of anchor residues
	core::Vector pos1 = pose.residue( start ).atom( 2 ).xyz();
	core::Vector pos2 = pose.residue( end ).atom( 2 ).xyz();
	
	// check TMspan
	// the reason this is not an exception is that we need this functionality
	// when transforming a protein into the membrane
	if ( pos1.z() < 0 && pos2.z() < 0 ){
		TR << "WARNING: If your starting PDB is already translated into the " << std::endl;
		TR << "membrane, then your TMspan does not span the membrane!" << std::endl;
	}
	else if ( pos1.z() > 0 && pos2.z() > 0 ){
		TR << "WARNING: If your starting PDB is already translated into the " << std::endl;
		TR << "membrane, then your TMspan does not span the membrane!" << std::endl;
	}
	
	// compute center
	core::Vector center = 0.5 * ( pos1 + pos2 );

	// compute normal
	core::Vector normal = pos2 - pos1;
	normal.normalize( mem_thickness );
	
	center_.assign( center.x(), center.y(), center.z() );
	normal_.assign( normal.x(), normal.y(), normal.z() );
	
}// from span

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this. Normal always shows in positive z-direction!
void EmbeddingDef::from_span_positive_z( core::pose::Pose const & pose, core::Size start, core::Size end ) {

	TR << "Computing membrane embedding from TMspan " << start << " to " << end << std::endl;

	// get CA atom positions of anchor residues
	core::Vector pos1 = pose.residue( start ).atom( 2 ).xyz();
	core::Vector pos2 = pose.residue( end ).atom( 2 ).xyz();
	
	// check TMspan
	// the reason this is not an exception is that we need this functionality
	// when transforming a protein into the membrane
	if ( pos1.z() < 0 && pos2.z() < 0 ){
	TR << "WARNING: If your starting PDB is already translated into the " << std::endl;
	TR << "membrane, then your TMspan does not span the membrane!" << std::endl;
	}
	else if ( pos1.z() > 0 && pos2.z() > 0 ){
	TR << "WARNING: If your starting PDB is already translated into the " << std::endl;
	TR << "membrane, then your TMspan does not span the membrane!" << std::endl;
	}

	// compute center
	core::Vector center = 0.5 * ( pos1 + pos2 );

	// compute normal in direction of increasing z-axis
	core::Vector normal;
	if ( pos1.z() > pos2.z() ){
		normal = pos1 - pos2;
	}
	else {
		normal = pos2 - pos1;
	}

	normal.normalize( mem_thickness );

	center_.assign( center.x(), center.y(), center.z() );
	normal_.assign( normal.x(), normal.y(), normal.z() );

}// from span, positive z direction

} // geometry
} // membrane
} // protocols

