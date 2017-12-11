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
#include <core/scoring/methods/OneBodyEnergy.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethod.hh> 
#include <core/scoring/EnergyMap.fwd.hh> 
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MinimizerMapBase.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <utility/vector1.hh>

static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

using namespace core::scoring;

--class--::--class--():
	methods::OneBodyEnergy( --class--::class_name() )
{

}

--class--::~--class--(){}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
methods::EnergyMethodOP
--class--::clone() const;

/// @brief Indicate required setup steps for scoring
void 
--class--::setup_for_scoring( core::pose::Pose & pose, ScoreFunction const & ) const;

/// @brief Is the score context dependent or context independent? 
void 
--class--::indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const;

/// @brief Indicates the current version of this score term
core::Size 
--class--::version() const;

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void 
--class--::finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;

Distance 
--class--::atomic_interaction_cutoff() const; 

--end_namespace--

