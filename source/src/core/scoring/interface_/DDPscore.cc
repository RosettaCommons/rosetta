// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/Interface_/DDPscore.cc
/// @brief Implementation of distance dependent interface score
/// @details The distance dependent score is a knowledge based potential specialized for
/// protein interfaces. This class loads and accesses the lookup tables.
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)

#include <utility/vector1.hh>
#include <core/scoring/interface_/DDPscore.hh>
#include <core/scoring/interface_/DDPscoreCreator.hh>


#include <core/scoring/ScoringManager.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>


#include <core/kinematics/Jump.hh>

namespace core {
namespace scoring {
namespace interface_ {


methods::EnergyMethodOP
DDPscoreCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DDPscore );
}

ScoreTypes
DDPscoreCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( interface_dd_pair );
	return sts;
}


DDPscore::DDPscore() :
	parent( methods::EnergyMethodCreatorOP( new DDPscoreCreator ) ),
	lookup_table_( core::scoring::ScoringManager::get_instance()->get_DDPLookupTable() )
{ }

methods::EnergyMethodOP DDPscore::clone() const
{
	return methods::EnergyMethodOP( new DDPscore(*this) );
}

void DDPscore::setup_for_scoring(pose::Pose&, const ScoreFunction& ) const
{
	//pose.update_residue_neighbors();
}

void DDPscore::setup_for_packing(
	pose::Pose&,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & ) const
{
	//pose.update_residue_neighbors();
}

void DDPscore::setup_for_derivatives(pose::Pose&, const ScoreFunction& ) const
{
	//pose.update_residue_neighbors();
}

bool DDPscore::defines_intrares_energy(core::scoring::EnergyMap const &) const
{
	return false;
}

void DDPscore::eval_intrares_energy(
	const core::conformation::Residue &,
	const core::pose::Pose &,
	const core::scoring::ScoreFunction &,
	core::scoring::EnergyMap &
) const {}


void DDPscore::residue_pair_energy(
	core::conformation::Residue const & rsd2,
	core::conformation::Residue const & rsd1,
	const core::pose::Pose &,
	const core::scoring::ScoreFunction &,
	core::scoring::EnergyMap & emap
) const
{
	debug_assert (rsd1.seqpos() != rsd2.seqpos()); //only call for distinct residues

	// Only score contacts across the interface
	if ( rsd1.chain() == rsd2.chain() )  return;

	core::Real distance = 1e3;
	for ( auto atom_it_1 = rsd1.atom_begin(), end1 = rsd1.heavyAtoms_end(); atom_it_1 != end1; ++atom_it_1 ) {
		for ( auto atom_it_2 = rsd2.atom_begin(), end2 = rsd2.heavyAtoms_end(); atom_it_2 != end2; ++atom_it_2 ) {
			if ( atom_it_1->xyz().distance(atom_it_2->xyz()) < distance ) {
				distance = atom_it_1->xyz().distance(atom_it_2->xyz());
			}
		}
	}

	if ( distance >= 10. || distance < 1.5 ||
			lookup_table_.get_potentials( rsd1.aa(), rsd2.aa(), distance ) > 0. ) {
		emap[ interface_dd_pair ] += 0.; // noop
	} else {
		emap[ interface_dd_pair ] += lookup_table_.get_potentials( rsd1.aa(), rsd2.aa(), distance );
	}
}


void DDPscore::indicate_required_context_graphs(utility::vector1< bool > & /*context_graphs_required*/ ) const
{
}

core::Distance DDPscore::atomic_interaction_cutoff() const{
	return 8.0;
}
core::Size
DDPscore::version() const
{
	return 1; // Initial versioning
}


} // InterfacePotentials
} // scoring
} // core
