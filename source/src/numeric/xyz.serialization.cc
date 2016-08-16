// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyz.serialization.cc
/// @brief  functions for serializing xyzVector and xyzMatrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifdef SERIALIZATION

// Package Headers
#include <numeric/xyz.serialization.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

namespace numeric {

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzVector< float > const & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzVector< double > const & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzVector< float > & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzVector< double > & xyz
)
{
	arch( xyz.x(), xyz.y(), xyz.z() );
}

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzMatrix< float > const & xyz
)
{
	arch( xyz.xx(), xyz.xy(), xyz.xz() );
	arch( xyz.yx(), xyz.yy(), xyz.yz() );
	arch( xyz.zx(), xyz.zy(), xyz.zz() );
}


template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzMatrix< double > const & xyz
)
{
	arch( xyz.xx(), xyz.xy(), xyz.xz() );
	arch( xyz.yx(), xyz.yy(), xyz.yz() );
	arch( xyz.zx(), xyz.zy(), xyz.zz() );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzMatrix< float > & xyz
)
{
	arch( xyz.xx(), xyz.xy(), xyz.xz() );
	arch( xyz.yx(), xyz.yy(), xyz.yz() );
	arch( xyz.zx(), xyz.zy(), xyz.zz() );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzMatrix< double > & xyz
)
{
	arch( xyz.xx(), xyz.xy(), xyz.xz() );
	arch( xyz.yx(), xyz.yy(), xyz.yz() );
	arch( xyz.zx(), xyz.zy(), xyz.zz() );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::xyzVector< float  > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::xyzVector< double > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::xyzMatrix< float  > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::xyzMatrix< double > );

}

#endif
