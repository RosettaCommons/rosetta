// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file devel/constrianed_sequence_design/ConstaintManager,fwd.hh
/// @brief
/// @author Javier Castellanos	(javiercv@uw.edu)
///

#ifndef INCLUDED_devel_constrained_sequence_SymmetricMultimerDesign_fwd_HH
#define INCLUDED_devel_constrained_sequence_SymmetricMultimerDesign_fwd_HH
#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace constrained_sequence_design {

class SymmetricMultimerDesign;

typedef utility::pointer::owning_ptr< SymmetricMultimerDesign > SymmetricMultimerDesignOP;
typedef utility::pointer::owning_ptr< SymmetricMultimerDesign const > SymmetricMultimerDesignCOP;


} // constrained_sequence_design
} // devel

#endif
