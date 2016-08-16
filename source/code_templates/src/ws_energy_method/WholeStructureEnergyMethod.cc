// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

// Unit headers
#include <--path--/--class--.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/EnergyMap.hh> 

// Basic/Utility headers
#include <utility/vector1.hh> 
#include <core/types.hh> 
#include <basic/Tracer.hh> 

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::-class--();

// copy constructor (not needed unless you need deep copies)
//--class--::--class--( --class-- const & src );

// destructor (important for properly forward-declaring smart-pointer members)
--class--::~--class--();

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
ore::scoring::methods::EnergyMethodOP 
--class--::clone() const {

}

/// @brief Indicate required setup steps for scoring
void 
--class--::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {

}

/// @brief Is the score context dependent or context independent? 
void 
--class--::indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const {

}

/// @brief Indicates the current version of this score term
core::Size 
--class--::version() const {

}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void 
--class--::finalize_total_energy( 
	core::pose::Pose & pose, 
	ScoreFunction const &, 
	EnergyMap & totals 
	) const {

}

Distance 
--class--::atomic_interaction_cutoff() const {

} 

--end_namespace--

