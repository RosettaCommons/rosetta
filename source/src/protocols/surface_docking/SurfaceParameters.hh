// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//      rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//          under license.
// (c) The Rosetta software is developed by the contributing members of the
//          Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions
//          about this can be
// (c) addressed to University of Washington UW TechTransfer,
//                                            email: license@u.washington.edu.

/// @file   SurfaceParameters.hh
/// @brief
/// @author Robin A Thottungal (rathottungal@gmail.com)
#ifndef INCLUDED_protocols_surface_docking_SurfaceParameters_hh
#define INCLUDED_protocols_surface_docking_SurfaceParameters_hh

// Project header
#include <core/pose/datacache/CacheableDataType.hh> //dono if I need this
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.hh>
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/exit.hh>


// C++ headers
#include <iostream>
#include <map>
#include <string>


// for storing plane equation, ax+by+cz+d=0
typedef utility::vector1< core::Real > Plane;
using namespace core;

namespace protocols {
namespace surface_docking {

//class SurfaceParameters: public utility::pointer::ReferenceCount  {

class SurfaceParameters: public basic::datacache::CacheableData   {

public:
	///@brief default constructor to initialize all values
	SurfaceParameters(); //:

	SurfaceParameters(std::string strSURFA0,
			std::string strSURFA1, std::string strSURFA2);

	// Require for CacheableData Copy ctor
	// @details Copy constructors must copy all data, not just some...
	SurfaceParameters( SurfaceParameters const & src );

	// Require for CacheableData, Leaving this makes this class abstract class
	//~SurfaceParameters();


	basic::datacache::CacheableDataOP
		clone() const
		{
			return new  SurfaceParameters( *this );
		}

	void GenerateSurfaceParameters( Vector SurfCG );

	Vector CalcNormalVector ( Vector Apoint , Vector Bpoint , Vector Cpoint );

	Vector CalcAxisIntersect ( Vector point1 , Vector point2 , Vector Bvector ,
															Vector Cvector );

	Plane GeneratePlane( Vector Apoint ,Vector Bpoint , Vector Cpoint );


	Vector PlanePointIntersection
		( Plane plane_abcd , Vector point_outofplane , Vector normal_plane );

	//overloading to use in surfaceOrientMover
	Vector PlanePointIntersection( Vector Point );

	Vector SplitSurfaceVectorString( std::string surfVectString );
	// Need to write a print operator for the class
// Not a good idea keeping the members private, temporary fix
public:

	// Holds the xyz co-ordinates that form the AB and AC vector
	Vector SURFA0;
	Vector SURFA1;
	Vector SURFA2;

	//string to hold each line that is used to write in the output pdb
	std::string strSURFA0;
	std::string strSURFA1;
	std::string strSURFA2;

	// A,B,C points of the unit cell, translated to the centroid of the surface
	Vector A,B,C; // its xyz coordinates not a vector
	// AB & AC Vector
	Vector vecAB,vecAC;

	Vector SurfaceCG; //Centroid of the surface

	// Parameters that relates to the plane
	Vector surfaceNormalVec;
	Vector unitsurfaceNormalVec;
	Vector surfaceAntiNormalVec;
	Plane surfacePlane; // ax+by+cz+d=0

	// Vector along which the slide into surface need to occur
	Vector slideaxis;
};

} // namespace surface_docking
} // namespace protocols


#endif
