// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/MultistateAggregateFunction.hh
/// @brief  declaration for class MultistateAggregateFunction to work with the PackDeamon classes
///         (not to be confused with J. Ashworth's MultiStateAggregateFunction class)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_MultistateAggregateFunction_hh
#define INCLUDED_protocols_pack_daemon_MultistateAggregateFunction_hh

// Unit headers
#include <protocols/pack_daemon/MultistateAggregateFunction.fwd.hh>

// Project headers
#include <core/types.hh>
#include <protocols/genetic_algorithm/Entity.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace pack_daemon {

class MultistateAggregateFunction : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;
	typedef utility::vector1< core::Real > StateEnergies;
	typedef utility::vector1< core::Size > StateIndices;

	typedef protocols::genetic_algorithm::Entity Entity;

public:
	MultistateAggregateFunction() : parent() {}
	virtual ~MultistateAggregateFunction();

	virtual core::Real   evaluate( StateEnergies const &, StateEnergies const &, Entity const &  ) = 0;
	virtual StateIndices select_relevant_states( StateEnergies const &, StateEnergies const &, Entity const & ) = 0;
};

}
}

#endif
