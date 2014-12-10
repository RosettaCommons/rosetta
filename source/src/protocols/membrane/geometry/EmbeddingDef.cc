// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/geometry/EmbeddingDef.cc
///
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
///
/// @note       Last Modified: 6/17/14
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)

// Unit headers
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.membrane.geometry.EmbeddingDef" );

namespace protocols {
namespace membrane {
namespace geometry {

using namespace core; 

/// @brief Default Constructor
EmbeddingDef::EmbeddingDef() :
	utility::pointer::ReferenceCount()
{
	center_.assign( 0, 0, 0 );
	normal_.assign( 0, 0, 1 );
}

/// @brief Standard Constructor
EmbeddingDef::EmbeddingDef( core::Vector const center, core::Vector const normal ) :
    utility::pointer::ReferenceCount(),
	normal_( normal ),
    center_( center )
{}

/// @brief Constructor from pose and two residue numbers
EmbeddingDef::EmbeddingDef( core::pose::PoseOP pose, core::Size start, core::Size end ) :
    utility::pointer::ReferenceCount()
{
	EmbeddingDef( *pose, start, end );
}

/// @brief Constructor from pose and two residue numbers
EmbeddingDef::EmbeddingDef( core::pose::Pose & pose, core::Size start, core::Size end ) :
utility::pointer::ReferenceCount()
{
	normal_.assign( 0, 0, 0 );
	center_.assign( 0, 0, 0 );
	
	from_span( pose, start, end );
}

/// @brief Copy Constructor
EmbeddingDef::EmbeddingDef( EmbeddingDef const & config ) :
    utility::pointer::ReferenceCount(),
    normal_( config.normal_ ),
    center_( config.center_ )
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
EmbeddingDef::show( std::ostream & output ) const {
	output << "Embedding: center: " << center_.to_string() << ", normal: " << normal_.to_string() << std::endl;
}

/// @brief Check reasonable range of vectors in embedding object
void EmbeddingDef::check_range() const{
    check_vector( normal_ );
    check_vector( center_ );
}

/// @brief Translate by center and normal
void EmbeddingDef::translate_by( EmbeddingDefOP translation ){
    center_ += translation->center();
    normal_ += translation->normal();
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

/// @brief Access Normal Param
core::Vector EmbeddingDef::normal() const {
	return normal_;
}

/// @brief Access center param
core::Vector EmbeddingDef::center() const {
	return center_;
}

/// @brief Check Object Equality
bool EmbeddingDef::equals( EmbeddingDef & other ) {
    
    if ( normal_ != other.normal() ) return false;
    if ( center_ != other.center() ) return false;
    
    return true;
}

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this.
void EmbeddingDef::from_span( core::pose::Pose & pose, core::Size start, core::Size end ) {

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
	normal.normalize(15);

	center_.assign( center.x(), center.y(), center.z() );
	normal_.assign( normal.x(), normal.y(), normal.z() );
	
}// from span

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this.
void EmbeddingDef::from_span( core::pose::PoseOP pose, core::Size start, core::Size end ) {
	from_span( *pose, start, end );
}

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this. Normal always shows in positive z-direction!
void EmbeddingDef::from_span_positive_z( core::pose::Pose & pose, core::Size start, core::Size end ) {

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
	normal.normalize(15);

	center_.assign( center.x(), center.y(), center.z() );
	normal_.assign( normal.x(), normal.y(), normal.z() );

}// from span, positive z direction

/// @brief Embedding object from span
/// @details Takes the CA coords of two residues and calculates center and normal
///				from this. Normal always shows in positive z-direction!
void EmbeddingDef::from_span_positive_z( core::pose::PoseOP pose, core::Size start, core::Size end ) {
	from_span_positive_z( *pose, start, end);
}

} // geometry
} // membrane
} // protocols

