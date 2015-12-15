// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/HomogeneousTransform.srlz.cc
/// @brief  functions for serializing HomogeneousTransforms
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifdef SERIALIZATION

// Unit headers
#include <numeric/HomogeneousTransform.srlz.hh>

// Package Headers
#include <numeric/HomogeneousTransform.hh>

// Utility headers
#include <utility/serialization/serialization.hh>

namespace numeric {

template < class Archive, class T >
void save_to_archive(
	Archive & arch,
	numeric::HomogeneousTransform< T > const & ht
)
{
	arch( ht.xx(), ht.xy(), ht.xz() );
	arch( ht.yx(), ht.yy(), ht.yz() );
	arch( ht.zx(), ht.zy(), ht.zz() );
	arch( ht.px(), ht.py(), ht.pz() );
}

template < class Archive, class T >
void load_from_archive(
	Archive & arch,
	numeric::HomogeneousTransform< T > & ht
)
{
	arch( ht.xx(), ht.xy(), ht.xz() );
	arch( ht.yx(), ht.yy(), ht.yz() );
	arch( ht.zx(), ht.zy(), ht.zz() );
	arch( ht.px(), ht.py(), ht.pz() );
}

template < class Archive >
void
save(
	Archive & arch,
	numeric::HomogeneousTransform< float > const & ht
)
{
	save_to_archive( arch, ht );
}

template < class Archive >
void
save(
	Archive & arch,
	numeric::HomogeneousTransform< double > const & ht
)
{
	save_to_archive( arch, ht );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::HomogeneousTransform< float > & ht
)
{
	load_from_archive( arch, ht );
}

template < class Archive >
void
load(
	Archive & arch,
	numeric::HomogeneousTransform< double > & ht
)
{
	load_from_archive( arch, ht );
}


EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::HomogeneousTransform< float  > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::HomogeneousTransform< double > );

}

#endif
