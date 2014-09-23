// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/ddG/Mutation.fwd.hh
///
/// @brief      Informaiton about a Point Mutation
/// @details    Object for scoring information about a point mutation - includes
///				the position to mutate to, and the specified amino acid.
///				Last Modified: 7/18/14
///
/// @remarks	This is pretty general and can probably live somewhere in core pack?
///				Not sure if anyone needs anything like this
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_ddG_Mutation_fwd_hh
#define INCLUDED_protocols_membrane_ddG_Mutation_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace ddG {

class Mutation;
typedef utility::pointer::shared_ptr< Mutation > MutationOP;
typedef utility::pointer::shared_ptr< Mutation const > MutationCOP;

} // membrane
} // conformation
} // core

#endif // INCLUDED_protocols_membrane_ddG_Mutation_fwd_hh

