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
/// @author Robin A Thottungal (raugust1@jhu.edu)
/// @author Michael Pacella (mpacella88@gmail.com)

// Unit Headers
#include <protocols/surface_docking/SurfaceParameters.hh>

// Numeric Headers
#include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/exit.hh>
// Basic Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.surfaceDocking.SurfaceParameters");

namespace protocols {
namespace surface_docking {

using namespace core;
using namespace numeric;
	
	SurfaceParameters::SurfaceParameters(
		xyzVector< core::Real > SURFA0,
		xyzVector< core::Real > SURFA1,
		xyzVector< core::Real > SURFA2 )
	{
	// error check for correct length 
	SURFA0_ = SURFA0;
	SURFA1_ = SURFA1;
	SURFA2_ = SURFA2;
	TR<<"Surface Vectors value inside constructor:"<<std::endl;
	TR<<"SURFA0:"<<SURFA0<<std::endl;
	TR<<"SURFA1:"<<SURFA1<<std::endl;
	TR<<"SURFA2:"<<SURFA2<<std::endl;
	
	// Generating the 2 basis vectors of the surface crystal system
	vecAB_ = SURFA1_-SURFA0_;
	vecAC_ = SURFA2_-SURFA0_;
	
	// Calculate the normal vector to the surface (slide_axis), the surface docking protocol
	// will determine which side of the surface the protein is on and assign the direction
	//accordingly
	slide_axis_ = vecAB_.cross(vecAC_);
	}

SurfaceParameters::SurfaceParameters( SurfaceParameters const & src ):
	utility::pointer::ReferenceCount()
	{
	SURFA0_=src.SURFA0_;
	SURFA1_=src.SURFA1_;
	SURFA2_=src.SURFA2_;
	vecAB_=src.vecAB_;
	vecAC_=src.vecAB_;
	slide_axis_=src.slide_axis_;
	}
	
SurfaceParametersOP SurfaceParameters::clone() const
	{
		return new SurfaceParameters( *this );
	}

SurfaceParameters::~SurfaceParameters(){}



}// namespace surface_docking
}// namespace protocols
