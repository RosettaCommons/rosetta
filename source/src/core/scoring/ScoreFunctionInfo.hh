// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_ScoreFunctionInfo_hh
#define INCLUDED_core_scoring_ScoreFunctionInfo_hh


// Unit Headers
#include <core/scoring/ScoreFunctionInfo.fwd.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {


/// Info on the scorefunction settings

/**
This object is handed to the pose's Energies object and stored
along with the cached energies. It is used in the next scoring
evaluation to decide whether it's safe to reuse cached energies.

It must describe the kinds of context that the Energies object provides
to the scoring function.  If any scoring function requires a tenA
neighbor graph, then that information is stored here.  If a scoring
function requires a different kind of neighborness graph (think
centroid mode), then this class will indicate the kind of storage.

Finally, it must describe the maximum atom-to-atom distance cutoff that characterizes
the energy function.

So all we need is a constructor and an operator==

**/

class ScoreFunctionInfo : public utility::pointer::ReferenceCount
{

public:


	/// default constructor -- fill this in
	ScoreFunctionInfo();
	ScoreFunctionInfo( ScoreFunctionInfo const & src );
	virtual ~ScoreFunctionInfo();


	/// copy constructor -- fill this in
	ScoreFunctionInfo( ScoreFunction const & scorefxn );

	void initialize_from( ScoreFunction const & scorefxn );

	/// comparison -- fill this in
	friend
	bool
	operator==( ScoreFunctionInfo const & a, ScoreFunctionInfo const & b ) /* PHIL */;

	/// comparison -- fill this in
	inline
	friend
	bool
	operator!=( ScoreFunctionInfo const & a, ScoreFunctionInfo const & b ) {
		return ! ( operator == ( a, b ));
	}

	EnergyMap const &
	scores_present() const {
		return scores_present_;
	}

	Distance
	max_atomic_interaction_distance() const {
		return max_atomic_interaction_distance_;
	}

	Distance
	max_context_neighbor_cutoff() const {
		return max_context_neighbor_cutoff_;
	}

	bool
	requires_context_graph( ContextGraphType cgt ) const;

private:
	ScoreFunctionInfo const & operator = ( ScoreFunctionInfo const & rhs ); // private, unimplemented; make sure this isn't called, or implement it

private:

	Distance max_atomic_interaction_distance_;
	Distance max_context_neighbor_cutoff_;

	utility::vector1< bool > context_graphs_required_;
	EnergyMap scores_present_;

	methods::EnergyMethodOptionsOP energy_method_options_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_ScoreFunctionInfo )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_ScoreFunction_HH
