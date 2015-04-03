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

#include <core/sequence/BasicSequenceCreators.hh>


#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceCoupling.hh>
#include <core/sequence/CompositeSequence.hh>
#include <core/sequence/ChemicalShiftSequence.hh>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

// class def for SimpleSequence
SimpleSequenceCreator::SimpleSequenceCreator() {}
SimpleSequenceCreator::~SimpleSequenceCreator() {}
SequenceOP SimpleSequenceCreator::create_sequence() const {
	return SequenceOP( new Sequence );
}

std::string SimpleSequenceCreator::keyname() const {
	return "sequence";
}

SequenceCouplingCreator::SequenceCouplingCreator() {}
SequenceCouplingCreator::~SequenceCouplingCreator() {}
SequenceOP SequenceCouplingCreator::create_sequence() const {
	return SequenceOP( new SequenceCoupling );
}

std::string SequenceCouplingCreator::keyname() const {
	return "sequence_coupling";
}
SequenceProfileCreator::SequenceProfileCreator() {}
SequenceProfileCreator::~SequenceProfileCreator() {}
SequenceOP SequenceProfileCreator::create_sequence() const {
	return SequenceOP( new SequenceProfile );
}

std::string SequenceProfileCreator::keyname() const {
	return "sequence_profile";
}

CompositeSequenceCreator::CompositeSequenceCreator() {}
CompositeSequenceCreator::~CompositeSequenceCreator() {}
SequenceOP CompositeSequenceCreator::create_sequence() const {
	return SequenceOP( new CompositeSequence );
}

std::string CompositeSequenceCreator::keyname() const {
	return "composite_sequence";
}

ChemicalShiftSequenceCreator::ChemicalShiftSequenceCreator() {}
ChemicalShiftSequenceCreator::~ChemicalShiftSequenceCreator() {}
SequenceOP ChemicalShiftSequenceCreator::create_sequence() const {
	return SequenceOP( new ChemicalShiftSequence );
}

std::string ChemicalShiftSequenceCreator::keyname() const {
	return "composite_sequence";
}

} //namespace sequence
} //namespace core
