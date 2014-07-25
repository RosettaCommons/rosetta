// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/geometry/EmbeddingDef.fwd.hh
///
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
///
/// @note       Last Modified: 6/17/14
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_EmbeddingDef_hh
#define INCLUDED_core_membrane_geometry_EmbeddingDef_hh

// Unit headers
#include <core/membrane/geometry/EmbeddingDef.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace geometry {
    
    /// @brief Embedding Setup Helper Class
class EmbeddingDef : public utility::pointer::ReferenceCount {
	
public:
	
	/// @brief Default Constructor
	EmbeddingDef();
	
	/// @brief Standard Constructor
	EmbeddingDef( core::Vector const normal, core::Vector const center );
	
	/// @brief Copy Constructor
	EmbeddingDef( EmbeddingDef const & EmbeddingDef );
	
	/// @brief Destructor
	~EmbeddingDef();
	
	/// @brief Standard Rosetta Show Method for Debugging
	virtual void show( std::ostream & output=std::cout ) const;
	
	/// @brief Check reasonable range of vectors in embedding object
	void check_range() const;
	
	/// @brief Access Normal Param
	core::Vector normal() const;
	
	/// @brief Access center param
	core::Vector center() const;
	
	/// @brief Translate by center and normal
	void translate_by( EmbeddingDefOP translation );
	
	/// @brief Set Center Param
	void set_center( core::Vector center );
	
	/// @brief Set Normal Param
	void set_normal( core::Vector normal );
	
	/// @brief Equals method
	bool equals( EmbeddingDef & other );
	
private: // data
	
	core::Vector normal_;
	core::Vector center_;
	
private: // function
	
	void check_vector( core::Vector vector ) const;
	
};
    
} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_EmbeddingDef_hh

