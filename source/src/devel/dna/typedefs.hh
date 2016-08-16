// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file fwd.hh
/// @brief non-class and non-pointer typedefs for this namespace
/// @author ashworth

#ifndef INCLUDED_devel_dna_typedefs_hh
#define INCLUDED_devel_dna_typedefs_hh

#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <string>
#include <map>

namespace devel {
namespace dna {

typedef std::map< core::Size, core::chemical::ResidueTypeCOP > ResTypeSequence;
typedef utility::vector1< ResTypeSequence > ResTypeSequences;
typedef std::map< core::Size, core::chemical::ResidueTypeCOPs > TypesMap;
typedef utility::vector1< core::Size > Positions;
typedef utility::vector1< std::string > Strings;

}
}

#endif
