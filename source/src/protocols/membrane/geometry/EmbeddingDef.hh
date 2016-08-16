// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/membrane/geometry/EmbeddingDef.hh
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// Last Modified: 6/11/15
/// #RosettaMPData

#ifndef INCLUDED_protocols_membrane_geometry_EmbeddingDef_hh
#define INCLUDED_protocols_membrane_geometry_EmbeddingDef_hh

// Unit headers
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {
namespace geometry {

/// @brief Embedding Setup Helper Class
class EmbeddingDef : public utility::pointer::ReferenceCount {

public:

	/// @brief Default Constructor
	EmbeddingDef();

	/// @brief Standard Constructor
	EmbeddingDef( core::Vector const center, core::Vector const normal );

	/// @brief Constructor from pose, two residue numbers, and bool if in positive z-direction
	EmbeddingDef( core::pose::Pose const & pose, core::Size const start, core::Size const end, bool pos_z=false );

	/// @brief Copy Constructor
	EmbeddingDef( EmbeddingDef const & src );

	/// @brief Assignment Operator
	EmbeddingDef & operator = ( EmbeddingDef const & src );

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

	/// @brief Invert normal
	void invert();

	/// @brief Translate by center and normal
	void translate_by( EmbeddingDef const & translation );

	/// @brief Set Center Param
	void set_center( core::Vector center );

	/// @brief Set Normal Param
	void set_normal( core::Vector normal );

	/// @brief Equals method
	bool equals( EmbeddingDef const & other ) const;

	/// @brief From span
	void from_span( core::pose::Pose const & pose, core::Size const start, core::Size const end );

	/// @brief From span
	void from_span_positive_z( core::pose::Pose const & pose, core::Size const start, core::Size const end );


private: // data

	core::Vector center_;
	core::Vector normal_;

};

} // geometry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_geometry_EmbeddingDef_hh

