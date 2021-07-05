// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/MissingEnergy.cc
/// @brief  Cost of failing to close loops
/// @author Arvind Kannan


// Unit headers
#include <core/energy_methods/MissingEnergy.hh>
#include <core/energy_methods/MissingEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/full_model_info/util.hh>

// Utility headers

// C++
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.energy_methods.MissingEnergy" );

/////////////////////////////////////////////////////////////////////////////////////
//
// Created as an attempted approach to raising the fraction of closed loops
// in SWM calculations.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the MissingEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
MissingEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< MissingEnergy >();
}

core::scoring::ScoreTypes
MissingEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( missing_res );
	return sts;
}


/// c-tor
MissingEnergy::MissingEnergy() :
	parent( utility::pointer::make_shared< MissingEnergyCreator >() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
MissingEnergy::clone() const
{
	return utility::pointer::make_shared< MissingEnergy >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
MissingEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {

	Size const num_missing_residue_connections = pose::full_model_info::get_number_missing_residues_and_connections( pose );
	totals[ core::scoring::missing_res ] = num_missing_residue_connections;

} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
void
MissingEnergy::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &,
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


} // scoring
} // core
