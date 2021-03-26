// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaMultipleSequenceAlignment.fwd.hh
/// @brief A class to store multiple sequence alignment data for input into trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_fwd_hh
#define INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_fwd_hh

#ifdef USE_TENSORFLOW

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace trRosetta {

class trRosettaMultipleSequenceAlignment;

using trRosettaMultipleSequenceAlignmentOP = utility::pointer::shared_ptr< trRosettaMultipleSequenceAlignment >;
using trRosettaMultipleSequenceAlignmentCOP = utility::pointer::shared_ptr< trRosettaMultipleSequenceAlignment const >;

} //trRosetta
} //protocols

#endif //USE_TENSORFLOW

#endif //INCLUDED_protocols_trRosetta_trRosettaMultipleSequenceAlignment_fwd_hh
