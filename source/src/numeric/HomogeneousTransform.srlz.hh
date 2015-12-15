// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/HomogeneousTransform.srlz.hh
/// @brief  functions for serializing xyzVector and xyzMatrix
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_numeric_HomogeneousTransform_SRLZ_HH
#define INCLUDED_numeric_HomogeneousTransform_SRLZ_HH

#ifdef SERIALIZATION

// Package headers
#include <numeric/HomogeneousTransform.fwd.hh>

// This file declares save/load functions for serializing and  deserializing
// HomogeneousTransform objects.

namespace numeric {

template < class Archive >
void
save(
	Archive & arch,
	HomogeneousTransform< float > const & xyz
);

template < class Archive >
void
save(
	Archive & arch,
	HomogeneousTransform< double > const & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	HomogeneousTransform< float > & xyz
);

template < class Archive >
void
load(
	Archive & arch,
	HomogeneousTransform< double > & xyz
);

}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_xyz_serialization_HH
