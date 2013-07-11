// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/IntermolEnergy.cc
/// @brief  Cost of bringing two chains together.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/methods/IntermolEnergy.hh>
#include <core/scoring/methods/IntermolEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>

// C++
#include <basic/Tracer.hh>

static basic::Tracer tr("core.scoring.methods.InterMolEnergy");

/////////////////////////////////////////////////////////////////////////////////////
//
// Created in attempt to fit 'Turner rules' for METHODS.
//  assumes 1 M standard state -- later will allow
//  specification of strand concentration(s) from
//  command line.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the IntermolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
IntermolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new IntermolEnergy;
}

ScoreTypes
IntermolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( intermol );
	return sts;
}


/// c-tor
IntermolEnergy::IntermolEnergy() :
	parent( new IntermolEnergyCreator )
{}

/// clone
methods::EnergyMethodOP
IntermolEnergy::clone() const
{
	return new IntermolEnergy;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {


	utility::vector1< Size > const chains = pose::full_model_info::figure_out_chains_from_full_model_info( pose );

	// now look to see which chains are represented in this pose...
	utility::vector1< Size > chains_in_pose;
	for ( Size k = 1; k <= pose.total_residue(); k++ ){
		if ( !chains_in_pose.has( chains[k] ) ) chains_in_pose.push_back( chains[k] );
	}

	Size const num_chains = chains_in_pose.size();

	tr.Debug << "Number of chains " <<   num_chains << std::endl;

	Size const num_chains_frozen = (num_chains > 0) ? (num_chains - 1) : 0;

	totals[ intermol ] = num_chains_frozen;

} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::eval_atom_derivative(
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
IntermolEnergy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core
