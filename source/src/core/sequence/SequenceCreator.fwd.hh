// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/sequence/SequenceCreator.hh
/// @brief  Base class for SequenceCreators for the Sequence load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_sequence_SequenceCreator_fwd_hh
#define INCLUDED_core_sequence_SequenceCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace sequence {

/// @brief Abstract base class for a Sequence factory; the Creator class is responsible for
/// creating a particular Sequence class.
class SequenceCreator;

typedef utility::pointer::shared_ptr< SequenceCreator > SequenceCreatorOP;
typedef utility::pointer::shared_ptr< SequenceCreator const > SequenceCreatorCOP;

} //namespace sequence
} //namespace core

#endif
