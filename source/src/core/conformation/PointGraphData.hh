// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/PointGraphData.hh
/// @brief  classes to work with UpperEdgeGraph for fast neighbor detection
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_conformation_PointGraphData_hh
#define INCLUDED_core_conformation_PointGraphData_hh

// Project Headers
#include <core/types.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

namespace core {
namespace conformation {

class PointGraphVertexData
{
public:
	PointGraphVertexData() : xyz_( 0.0, 0.0, 0.0 ) {}
	PointGraphVertexData( numeric::xyzVector< core::Real > const & coors ) : xyz_( coors ) {}

	/// @brief Get a non-const reference to xyz data in order to set the data by reference.
	numeric::xyzVector< core::Real > & xyz() { return xyz_; }
	numeric::xyzVector< core::Real > const & xyz() const { return xyz_; }

	static int const NUM_EDGES_TO_RESERVE = 50;

private:
	numeric::xyzVector< core::Real > xyz_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class PointGraphEdgeData
{
public:
	PointGraphEdgeData() : dsq_( 0.0 ) {}

	/// @brief inputs and outputs are distances squared
	PointGraphEdgeData( platform::Real d2 ) : dsq_( d2 ) {}

	core::Real & dsq() { return dsq_; }
	core::Real dsq() const { return dsq_; }

private:
	core::Real dsq_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}

#endif
