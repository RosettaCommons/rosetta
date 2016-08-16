// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief declaration of implementation class for abstract class Residue
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueType_fwd_hh
#define INCLUDED_core_chemical_ResidueType_fwd_hh


// Unit headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


// C++ headers

namespace core {
namespace chemical {


class ResidueType;

typedef  utility::pointer::weak_ptr< ResidueType >  ResidueTypeAP;
typedef  utility::pointer::weak_ptr< ResidueType const >  ResidueTypeCAP;
typedef  utility::pointer::shared_ptr< ResidueType >  ResidueTypeOP;
typedef  utility::pointer::shared_ptr< ResidueType const >  ResidueTypeCOP;
typedef  utility::vector1< ResidueTypeOP >  ResidueTypeOPs;
typedef  utility::vector1< ResidueTypeCAP >  ResidueTypeCAPs;
typedef  utility::vector1< ResidueTypeCOP >  ResidueTypeCOPs;

typedef utility::vector1< Size > AtomIndices;
typedef utility::vector1<Size> OrbitalIndices;

} // chemical
} // core

#endif // INCLUDED_core_chemical_Residues_HH
