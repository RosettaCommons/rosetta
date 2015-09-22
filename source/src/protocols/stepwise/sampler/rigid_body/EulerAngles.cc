// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rigid_body/EulerAngles.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/rigid_body/EulerAngles.hh>
#include <numeric/conversions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyzMatrix.hh>

#include <cmath>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.swa.rna.EulerAngles" );

/////////////////////////////////////////////////////////////////
// Yet another rotation matrix / euler angle object.
//
// Includes nice functions set up by Parin Sripakdeevong for
//  "floating base" modeler in RNA stepwise assembly.
//
// Probably will re-use for aptamer modeling.
//
/////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

//Constructor
EulerAngles::EulerAngles():
	alpha_( 0.0 ),
	beta_( 0.0 ),
	gamma_( 0.0 ),
	z_( 1.0 )
{}

//Constructor
EulerAngles::EulerAngles( numeric::xyzMatrix< core::Real > const & rotation_matrix )
{
	initialize_from_rotation_matrix( rotation_matrix );
}

//Destructor
EulerAngles::~EulerAngles()
{}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
EulerAngles::set_beta( Real const setting ){
	beta_ = setting;
	z_ = cos( beta_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
EulerAngles::set_z( Real const setting ){
	z_ = setting;
	runtime_assert( (z_ >= -1.0 ) && ( z_ <= 1.0 ) );
	beta_ = acos( z_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
EulerAngles::initialize_from_rotation_matrix( numeric::xyzMatrix< core::Real > const & rotation_matrix ){

	Real determinant = rotation_matrix.det();
	runtime_assert(  std::abs( determinant - 1.0 ) < 0.000001 );

	numeric::xyzMatrix< core::Real > const & M = rotation_matrix;
	const Real DEGS_PER_RAD = 180. / numeric::NumericTraits < Real > ::pi();

	alpha_ = atan2( M.xz(),  - M.yz() ) * DEGS_PER_RAD;
	z_     = M.zz();
	gamma_ = atan2( M.zx(), M.zy() ) * DEGS_PER_RAD ;  //tan2(y,x)=gamma

	beta_ = acos( z_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
EulerAngles::convert_to_rotation_matrix( numeric::xyzMatrix< core::Real > & rotation_matrix ){

	//Probably could save time if determine x and y and take cross product to determine z
	//Should determine using both ways and check for consistency

	rotation_matrix.xx( cos( alpha_ )*cos( gamma_ ) - sin( alpha_ )*cos( beta_ )*sin( gamma_ ) );
	rotation_matrix.xy(  - cos( alpha_ )*sin( gamma_ ) - sin( alpha_ )*cos( beta_ )*cos( gamma_ ) );
	rotation_matrix.xz( sin( alpha_ )*sin( beta_ ) );
	rotation_matrix.yx( sin( alpha_ )*cos( gamma_ ) + cos( alpha_ )*cos( beta_ )*sin( gamma_ ) );
	rotation_matrix.yy(  - sin( alpha_ )*sin( gamma_ ) + cos( alpha_ )*cos( beta_ )*cos( gamma_ ) ); //Found bug on Feb 13, 2010...previously had cos(gamma_) instead of sin(gamma_)
	rotation_matrix.yz(  - cos( alpha_ )*sin( beta_ ) );
	rotation_matrix.zx( sin( beta_ ) *sin( gamma_ ) );
	rotation_matrix.zy( sin( beta_ ) *cos( gamma_ ) );
	rotation_matrix.zz( cos( beta_ ) );

	Real determinant = rotation_matrix.det();
	runtime_assert(  std::abs( determinant - 1.0 ) < 0.000001 );
}


} //rigid_body
} //sampler
} //stepwise
} //protocols
