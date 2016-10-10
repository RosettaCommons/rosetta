// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/bb_sampler/BBDihedralSampler.cc
/// @brief This class functions to hold, access, and set independent and dependent dihedral data.
///   It can act as a base class for particular types of data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)  and Jason W. Labonte (JWLabonte@jhu.edu)

#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace simple_moves {
namespace bb_sampler {

using core::pose::Pose;
using core::Size;
using core::Real;

////////////////////////////////////////
BBDihedralSamplerBase::BBDihedralSamplerBase():
	torsion_type_(core::id::phi_dihedral),
	sampling_type_(probability)

{

}

BBDihedralSamplerBase::BBDihedralSamplerBase( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type ):
	torsion_type_(torsion_type),
	sampling_type_(sampling_type)
{

}

BBDihedralSamplerBase::~BBDihedralSamplerBase(){}

BBDihedralSamplerBase::BBDihedralSamplerBase( BBDihedralSamplerBase const & src ):
	torsion_type_(src.torsion_type_),
	sampling_type_(src.sampling_type_)
{

}

BBDihedralSamplerBaseOP
BBDihedralSamplerBase::clone() const {
	utility_exit_with_message( "clone has been called on a BBDihedralSampler which has not overridden the base class implementation.");
	return BBDihedralSamplerBaseOP(NULL);
}


//////////////////// 1D ////////////////
BBDihedralSampler::BBDihedralSampler():
	BBDihedralSamplerBase()
{

}

BBDihedralSampler::BBDihedralSampler( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type ):
	BBDihedralSamplerBase( torsion_type, sampling_type )
{

}

BBDihedralSampler::~BBDihedralSampler(){}

BBDihedralSampler::BBDihedralSampler( BBDihedralSampler const & src ):
	BBDihedralSamplerBase( src )
{

}

BBDihedralSamplerOP
BBDihedralSampler::clone() const {
	utility_exit_with_message( "clone has been called on a BBDihedralSampler which has not overridden the base class implementation.");
	return BBDihedralSamplerOP(NULL);
}


//////////////////// 2D ////////////////////

BBDihedralSampler2D::BBDihedralSampler2D():
	BBDihedralSamplerBase()
{

}

BBDihedralSampler2D::BBDihedralSampler2D( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type ):
	BBDihedralSamplerBase( torsion_type, sampling_type )
{

}

BBDihedralSampler2D::~BBDihedralSampler2D(){}

BBDihedralSampler2D::BBDihedralSampler2D( BBDihedralSampler2D const & src ):
	BBDihedralSamplerBase( src )
{

}

BBDihedralSampler2DOP
BBDihedralSampler2D::clone() const {
	utility_exit_with_message( "clone has been called on a BBDihedralSampler which has not overridden the base class implementation.");
	return BBDihedralSampler2DOP(NULL);
}


//////////////////// 3D ////////////////////

BBDihedralSampler3D::BBDihedralSampler3D():
	BBDihedralSamplerBase()
{

}

BBDihedralSampler3D::BBDihedralSampler3D( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type ):
	BBDihedralSamplerBase( torsion_type, sampling_type )
{

}

BBDihedralSampler3D::~BBDihedralSampler3D(){}

BBDihedralSampler3D::BBDihedralSampler3D( BBDihedralSampler3D const & src ):
	BBDihedralSamplerBase( src )
{

}

BBDihedralSampler3DOP
BBDihedralSampler3D::clone() const {
	utility_exit_with_message( "clone has been called on a BBDihedralSampler which has not overridden the base class implementation.");
	return BBDihedralSampler3DOP(NULL);
}


/////////////////// ND /////////////////

BBDihedralSamplerND::BBDihedralSamplerND():
	BBDihedralSamplerBase()
{

}

BBDihedralSamplerND::BBDihedralSamplerND( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type ):
	BBDihedralSamplerBase( torsion_type, sampling_type )
{

}

BBDihedralSamplerND::~BBDihedralSamplerND(){}

BBDihedralSamplerND::BBDihedralSamplerND( BBDihedralSamplerND const & src ):
	BBDihedralSamplerBase( src )
{

}

BBDihedralSamplerNDOP
BBDihedralSamplerND::clone() const {
	utility_exit_with_message( "clone has been called on a BBDihedralSampler which has not overridden the base class implementation.");
	return BBDihedralSamplerNDOP(NULL);
}

} //protocols
} //pose
} //core
