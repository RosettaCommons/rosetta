// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableSequenceProfile.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_core_sequence_CacheableSequenceProfile_hh
#define INCLUDED_core_sequence_CacheableSequenceProfile_hh

// unit headers
#include <core/sequence/CacheableSequenceProfile.fwd.hh>
// package headers
#include <basic/datacache/CacheableData.hh>

// project headers
#include <core/sequence/SequenceProfile.hh>


namespace core {
namespace sequence {


// simple class for holding onto a SequenceProfile
class CacheableSequenceProfile : public basic::datacache::CacheableData {
	CacheableSequenceProfile( core::sequence::SequenceProfile prof ) {
		profile_ = prof;
	}

	virtual ~CacheableSequenceProfile() {};
	virtual basic::datacache::CacheableDataOP clone() const {
		return new CacheableSequenceProfile(*this);
	}

	core::sequence::SequenceProfile & profile() {
		return profile_;
	}

	core::sequence::SequenceProfile const & profile() const {
		return profile_;
	}
private:
	core::sequence::SequenceProfile profile_;
};


} // namespace sequence
} // namespace core

#endif /* INCLUDED_core_sequence_CacheableSequenceProfile_HH */
