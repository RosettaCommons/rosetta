// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

#ifndef INCLUDED_devel_constrained_sequence_Symmetrizer_fwd_HH
#define INCLUDED_devel_constrained_sequence_Symmetrizer_fwd_HH
#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace matdes {

class Symmetrizer;

typedef utility::pointer::shared_ptr< Symmetrizer > SymmetrizerOP;
typedef utility::pointer::shared_ptr< Symmetrizer const > SymmetrizerCOP;


} // matdes
} // devel

#endif
