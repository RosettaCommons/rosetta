// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyz.serialization.hh
/// @brief  functions for serializing xyzVector and xyzMatrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_numeric_xyz_serialization_HH
#define INCLUDED_numeric_xyz_serialization_HH

#ifdef SERIALIZATION

// Package headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// This file declares save/load functions (at global scope) for serializing and
// deserializing numerix::xyzVector and numeric::xyzMatrix objects.

namespace numeric {

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzVector< float > const & xyz
);

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzVector< double > const & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzVector< float > & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzVector< double > & xyz
);

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzMatrix< float > const & xyz
);

template < class Archive >
void
save(
	Archive & arch,
	numeric::xyzMatrix< double > const & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzMatrix< float > & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	numeric::xyzMatrix< double > & xyz
);

}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_xyz_serialization_HH
