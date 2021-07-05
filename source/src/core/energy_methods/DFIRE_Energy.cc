// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/DFIRE_Energy.cc
/// @author James Thompson

// Unit Headers
#include <core/energy_methods/DFIRE_Energy.hh>
#include <core/energy_methods/DFIRE_EnergyCreator.hh>

// Package Headers
#include <core/scoring/DenseEnergyContainer.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/methods/dfire/DFIRE_Potential.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/Methods.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "core.energy_methods.DFIRE_Energy" );

namespace core {
namespace energy_methods {

namespace dfire {

/// @details This must return a fresh instance of the DFIRE_Energy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
DFIRE_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< DFIRE_Energy >();
}

core::scoring::ScoreTypes
DFIRE_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( core::scoring::DFIRE );
	return sts;
}

/// ctor
DFIRE_Energy::DFIRE_Energy() :
	parent( utility::pointer::make_shared< DFIRE_EnergyCreator >() )
{
	potential_is_loaded_ = core::scoring::methods::dfire::get_DFIRE_potential().is_loaded();
}

core::scoring::methods::LongRangeEnergyType
DFIRE_Energy::long_range_type() const {
	return scoring::methods::DFIRE;
}

void
DFIRE_Energy::setup_for_scoring(
	pose::Pose & pose, core::scoring::ScoreFunction const &
) const {
	using namespace core::scoring::methods;
	//std::cout << "called setup_for_scoring" << std::endl;
	// Do we have a potential yet?
	if ( !potential_is_loaded_ ) {
		utility_exit_with_message("No potential loaded.");
	}

	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );
	// create a container
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;

	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		core::scoring::LREnergyContainerOP new_dec( new core::scoring::DenseEnergyContainer( pose.size(), core::scoring::gb_elec ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


/// clone
core::scoring::methods::EnergyMethodOP
DFIRE_Energy::clone() const
{
	return utility::pointer::make_shared< DFIRE_Energy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods
/////////////////////////////////////////////////////////////////////////////


void
DFIRE_Energy::residue_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::EnergyMap &
) const {
	//std::cout << "called residue_energy" << std::endl;
	return;
}

bool DFIRE_Energy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	return ( pose.residue_type(res1).is_protein() && pose.residue_type(res2).is_protein() );
}

void
DFIRE_Energy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const {
	//std::cout << "called eval_intrares_energy" << std::endl;
	return;
}

void
DFIRE_Energy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	//std::cout << "residue_pair_energy" << std::endl;
	emap[ core::scoring::DFIRE ] += core::scoring::methods::dfire::get_DFIRE_potential().eval_dfire_pair_energy( rsd1, rsd2 );
}

/// @brief Energy is context independent and thus indicates that no context
/// graphs need to be maintained.
void
DFIRE_Energy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const {}

core::Size
DFIRE_Energy::version() const {
	return 1; // Initial versioning
}


} // dfire
} // scoring
} // core
