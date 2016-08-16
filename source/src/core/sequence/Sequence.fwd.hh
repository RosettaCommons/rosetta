// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Sequence.fwd.hh
/// @brief forward declaration for Sequence class.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_Sequence_fwd_hh
#define INCLUDED_core_sequence_Sequence_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace sequence {

class Sequence;
typedef utility::pointer::shared_ptr< Sequence > SequenceOP;
typedef utility::pointer::shared_ptr< Sequence const > SequenceCOP;

} // sequence
} // core

#endif
