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

#ifndef INCLUDED_core_sequence_SequenceCreator_hh
#define INCLUDED_core_sequence_SequenceCreator_hh

// Unit Headers
#include <core/sequence/SequenceCreator.fwd.hh>

// Package Headers
#include <core/sequence/Sequence.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace core {
namespace sequence {

/// @brief Abstract base class for a Mover factory; the Creator class is responsible for
/// creating a particular mover class.
class SequenceCreator : public utility::pointer::ReferenceCount
{
public:
	SequenceCreator();
	~SequenceCreator() override;

	virtual SequenceOP create_sequence() const = 0;
	virtual std::string keyname() const = 0;
};

} //namespace
} //namespace core

#endif
