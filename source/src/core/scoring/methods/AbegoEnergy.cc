// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/AbegoEnergy.cc
/// @brief  ABEGO energy method class implementation
/// @author imv@uw.edu

// Unit Headers
#include <core/scoring/methods/AbegoEnergy.hh>
#include <core/scoring/methods/AbegoEnergyCreator.hh>

// Package Headers
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/EnergyMap.hh>
#include <basic/options/option.hh>
#include <core/scoring/P_AA_ABEGO3.hh>

// Utility headers
#include <numeric/conversions.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#include <basic/Tracer.hh>
static basic::Tracer TR( "core.scoring.AbegoEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the Abego class,
/// never an instance already in use
methods::EnergyMethodOP
AbegoEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new Abego );
}

ScoreTypes
AbegoEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( paa_abego3 );
	return sts;
}


/// ctor
Abego::Abego() :
	//parent( EnergyMethodCreatorOP( new AbegoEnergyCreator ) ),
	WholeStructureEnergy( EnergyMethodCreatorOP( new AbegoEnergyCreator ) ),
	paa_abego3_( ScoringManager::get_instance()->get_P_AA_ABEGO3() )
{}

/// clone
EnergyMethodOP
Abego::clone() const
{
	return EnergyMethodOP( new Abego );
}

void Abego::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const
{

	TR << "ABEGO::setup_for_scoring" << std::endl;

	energy_positive_sum_ = 0;
	energy_positive_sum_count_ = 0;
	energy_sum_ = 0;
	energy_sum_count_ = 0;

	utility::vector1< std::string > abego = abegoManager_.get_symbols(pose, 1);
	for ( core::Size chain(1), num_chains(pose.conformation().num_chains()); chain <= num_chains; ++chain ) {
		Size chain_begin = pose.conformation().chain_begin(chain);
		Size chain_end = pose.conformation().chain_end(chain);
		for ( Size triplet_center = chain_begin + 1; triplet_center < chain_end - 1; ++triplet_center ) {
			// ABEGO letters are assigned to each amino acid. The ABEGOManager will assign values to the first and last
			// residue in a chain, but these will be assigned codes 'G' and 'O' respectively, to match Gabe's method and
			// data file from which scores are loaded.
			// Note: For some reason ABEGOManageR::get_symbols(Pose& const, Size level) returns a vector of 1-character strings
			// instead of just a string... hence the [0] indexing into each item.
			char abego_previous = (triplet_center == chain_begin + 1)? 'G' : abego[triplet_center - 1][0];
			char abego_current = abego[triplet_center][0];
			char abego_next = (triplet_center == chain_end - 2)? 'O' : abego[triplet_center + 1][0];

			// All 3 residues must have a standard ABEGO type. '-' is assigned in all other cases.
			core::chemical::AA aa = pose.conformation().residue(triplet_center).aa();
			if ( abego_previous == '-' || abego_current == '-' || abego_next == '-' ) {
				continue;
			}

			// Skip noncanonical AA's
			if ( aa == core::chemical::aa_none || core::chemical::num_canonical_aas < aa ) {
				continue;
			}

			core::Real triplet_energy = paa_abego3_.P_AA_ABEGO3_energy(abego_previous, abego_current, abego_next, aa);

			if ( 0 < triplet_energy ) {
				energy_positive_sum_ += triplet_energy;
				++energy_positive_sum_count_;
			}

			energy_sum_ += triplet_energy;
			++energy_sum_count_;
		}
	}
}

void Abego::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & totals) const
{
	TR << "ABEGO::finalize_total_energy start" << std::endl;
	totals[ paa_abego3 ] = energy_sum_;
	TR << "ABEGO::finalize_total_energy end" << std::endl;
}

/// @brief Abego3 Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
Abego::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
Abego::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

