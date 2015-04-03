// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MissingEnergy.cc
/// @brief  Cost of failing to close loops
/// @author Arvind Kannan


// Unit headers
#include <core/scoring/methods/MissingEnergy.hh>
#include <core/scoring/methods/MissingEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

// Utility headers
#include <utility/vector1.hh>

// C++
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.methods.MissingEnergy" );

/////////////////////////////////////////////////////////////////////////////////////
//
// Created as an attempted approach to raising the fraction of closed loops
// in SWM calculations.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MissingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MissingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MissingEnergy );
}

ScoreTypes
MissingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( missing_res );
	return sts;
}


/// c-tor
MissingEnergy::MissingEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MissingEnergyCreator ) )
{}

/// clone
methods::EnergyMethodOP
MissingEnergy::clone() const
{
	return methods::EnergyMethodOP( new MissingEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
MissingEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {

	Size const num_missing_residue_connections = pose::full_model_info::get_number_missing_residues_and_connections( pose );
	totals[ missing_res ] = num_missing_residue_connections;

} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
void
MissingEnergy::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
 	) const
{
	// no op.
} // eval atom derivative

core::Size
MissingEnergy::version() const
{
	return 2; // Initial versioning
}


} // methods
} // scoring
} // core
