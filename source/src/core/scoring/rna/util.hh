// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/util.hh
/// @brief  Some useful nonmember functions

#ifndef INCLUDED_core_scoring_rna_util_HH
#define INCLUDED_core_scoring_rna_util_HH

#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace rna {


Size
rna_residue_name_to_num( char const c );

Size
protein_atom_name_to_num( std::string const & name );

Size
protein_atom_name_to_num( std::string const & name, std::string const & resname );

}
}
}

#endif

