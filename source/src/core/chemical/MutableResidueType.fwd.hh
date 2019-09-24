// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief declaration of MutableResidueType class
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_MutableResidueType_fwd_hh
#define INCLUDED_core_chemical_MutableResidueType_fwd_hh


// Unit headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <utility/vector1.fwd.hh>

namespace core {
namespace chemical {


class MutableResidueType;

typedef  utility::pointer::shared_ptr< MutableResidueType >  MutableResidueTypeOP;
typedef  utility::pointer::shared_ptr< MutableResidueType const >  MutableResidueTypeCOP;

typedef  utility::vector1< MutableResidueTypeOP >  MutableResidueTypeOPs;
typedef  utility::vector1< MutableResidueTypeCOP >  MutableResidueTypeCOPs;

} // chemical
} // core

#endif // INCLUDED_core_chemical_Residues_HH
