// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingDef.cc
///
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
///
/// @note       Last Modified: 6/17/14
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)

// Unit headers
#include <core/membrane/geometry/EmbeddingDef.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "core.membrane.geometry.EmbeddingDef" );

namespace core {
namespace membrane {
namespace geometry {

/// @brief Default Constructor
EmbeddingDef::EmbeddingDef() :
	utility::pointer::ReferenceCount()
{
	Vector normal(0, 0, 1);
	Vector center(0, 0, 0);
	EmbeddingDef( normal, center );
}

/// @brief Standard Constructor
EmbeddingDef::EmbeddingDef( core::Vector const normal, core::Vector const center ) :
    utility::pointer::ReferenceCount(),
    normal_( normal ),
    center_( center )
{}

/// @brief Copy Constructor
EmbeddingDef::EmbeddingDef( EmbeddingDef const & config ) :
    utility::pointer::ReferenceCount(),
    normal_( config.normal_ ),
    center_( config.center_ )
{}

/// @brief Destructor
EmbeddingDef::~EmbeddingDef() {}

/// @brief Standard Rosetta Show Method for Debugging
void
EmbeddingDef::show( std::ostream & output ) const {
    
    output << "Embedding Setup " << std::endl;
    output << "normal=(" << normal_.x() << "," << normal_.y() << "," << normal_.z() << ")" << std::endl;
    output << "center=(" << center_.x() << "," << center_.y() << "," << center_.z() << ")" << std::endl;
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
core::Vector EmbeddingDef::normal() const { return normal_; }

/// @brief Access center param
core::Vector EmbeddingDef::center() const { return center_; }

/// @brief Check Object Equality
bool EmbeddingDef::equals( EmbeddingDef & other ) {
    
    if ( normal_ != other.normal() ) return false;
    if ( center_ != other.center() ) return false;
    
    return true;
}

/// @brief Check reasonable range of vector
void EmbeddingDef::check_vector( core::Vector vector ) const {
    TR << "Checking vector " << std::endl;
    
    // warn if vector is origin
    if ( vector.to_string() == "(0, 0, 0)"){
        TR << "WARNING: your vector is (0, 0, 0)!" << std::endl;
    }
    
    // Fail if vector is out of range
    if ( vector.x() < -1000 || vector.x() > 1000 ||
        vector.y() < -1000 || vector.y() > 1000 ||
        vector.z() < -1000 || vector.z() > 1000 ){
        
        throw new conformation::membrane::EXCN_Illegal_Arguments("Unreasonable range for center or normal! Check your input vectors!");
    }
}// check_vector


} // geometry
} // membrane
} // core

