// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author Javier Castellanos ( javiercv@uw.edu )

#ifndef INCLUDED_devel_constrained_sequence_design_SequenceConstraintSet_hh
#define INCLUDED_devel_constrained_sequence_design_SequenceConstraintSet_hh

// Package Header
#include <devel/constrained_sequence_design/SequenceConstraint.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.fwd.hh>

// Project Header
#include <utility/pointer/ReferenceCount.hh>

// STL Header
#include <set>

namespace devel {
namespace constrained_sequence_design {

class SequenceConstraintSet : public utility::pointer::ReferenceCount {
public:
	typedef std::set< SequenceConstraintOP >::const_iterator const_iterator;
	typedef std::set< SequenceConstraintOP >::iterator iterator;
public:
	SequenceConstraintSet();
	virtual ~SequenceConstraintSet();
	SequenceConstraintSet(const SequenceConstraintSet& cs);

	inline Size size() { return constraints_.size(); }
	inline iterator begin() { return constraints_.begin(); }
	inline iterator end() { return constraints_.end(); }

	inline void insert( const SequenceConstraintOP& c) { constraints_.insert(c) ;}
private:
	std::set < SequenceConstraintOP > constraints_;
};

} // namespace devel
} // namespace constrained_sequence_design

#endif
