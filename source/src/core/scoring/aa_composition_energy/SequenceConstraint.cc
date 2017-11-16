// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/aa_composition_energy/SequenceConstraint.cc
/// @brief A base class for constraining sequences, analogous to a geometric constraint.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace aa_composition_energy {

static basic::Tracer TR( "core.scoring.aa_composition_energy.SequenceConstraint" );

/// @brief Constructor
///
SequenceConstraint::SequenceConstraint( core::scoring::ScoreType const & t ):
	core::scoring::constraints::Constraint( t ),
	dummy_atomid_(0,0)
	//TODO -- initialize variables here.
{}

/// @brief Copy constructor
///
SequenceConstraint::SequenceConstraint( SequenceConstraint const &src ):
	core::scoring::constraints::Constraint( src.score_type() ),
	dummy_atomid_(0,0)
	//TODO -- copy variables here.
{
}

/// @brief Destructor
///
SequenceConstraint::~SequenceConstraint() {}


} // aa_composition_energy
} // scoring
} // core
