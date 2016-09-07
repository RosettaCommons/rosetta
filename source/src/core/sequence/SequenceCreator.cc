// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/sequence/SequenceCreator.hh
/// @brief  Base class for SequenceCreators for the Sequence load-time factory registration scheme
/// @author James Thompson

// Unit Headers
#include <core/sequence/SequenceCreator.hh>

namespace core {
namespace sequence {

SequenceCreator::SequenceCreator() {}
SequenceCreator::~SequenceCreator() = default;

} // namespace sequence
} // namespace core
