// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/membrane/types.hh
/// @brief  Default values for membrane code
/// @author JKLeman (julia.koehler1982@gmail.com)
/// @details Should supersede files in Rosetta/main/database/membrane

#ifndef INCLUDED_core_conformation_membrane_types_hh
#define INCLUDED_core_conformation_membrane_types_hh

// Package headers
#include <numeric/xyzVector.hh>
#include <core/types.hh>

#include <string>

namespace core {
namespace conformation {
namespace membrane {

// These are the default values that should be used in the membrane code!

// membrane center and normal
static numeric::xyzVector< Real > const mem_center( 0, 0, 0 );
static numeric::xyzVector< Real > const mem_normal( 0, 0, 15 );

// half of the membrane thickness: should match with normal.z
static Real const mem_thickness( 15 );

// anchor point of membrane residue - first residue of the protein
static Size const mem_jump( 2 );
static Size const mem_anchor( 1 );

// values for the energy function
static Real const mem_steepness( 10 );


} // namespace membrane
} // namespace conformation
} // namespace core


#endif // INCLUDED_core_id_types_HH
