// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DFIRE_Energy.cc
/// @author James Thompson

// Unit Headers
#include <core/scoring/methods/dfire/DFIRE_Energy.hh>
#include <core/scoring/methods/dfire/DFIRE_EnergyCreator.hh>

// Package Headers
#include <core/scoring/DenseEnergyContainer.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/dfire/DFIRE_Potential.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/Methods.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "core.scoring.methods.DFIRE_Energy" );

namespace core {
namespace scoring {
namespace methods {
namespace dfire {

/// @details This must return a fresh instance of the DFIRE_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
DFIRE_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DFIRE_Energy );
}

ScoreTypes
DFIRE_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( core::scoring::DFIRE );
	return sts;
}

/// ctor
DFIRE_Energy::DFIRE_Energy() :
	parent( methods::EnergyMethodCreatorOP( new DFIRE_EnergyCreator ) )
{
	potential_is_loaded_ = core::scoring::methods::dfire::get_DFIRE_potential().is_loaded();
}

methods::LongRangeEnergyType
DFIRE_Energy::long_range_type() const { 
	return scoring::methods::DFIRE; 
}

void
DFIRE_Energy::setup_for_scoring(
	pose::Pose & pose, ScoreFunction const &
) const {
	using namespace core::scoring::methods;
	//std::cout << "called setup_for_scoring" << std::endl;
	// Do we have a potential yet?
	if (!potential_is_loaded_) {
		utility_exit_with_message("No potential loaded.");
	}
	
	LongRangeEnergyType const & lr_type( long_range_type() );
	// create a container
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;

	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		LREnergyContainerOP new_dec( new DenseEnergyContainer( pose.total_residue(), gb_elec ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}


/// clone
EnergyMethodOP
DFIRE_Energy::clone() const
{
	return EnergyMethodOP( new DFIRE_Energy );
}

/////////////////////////////////////////////////////////////////////////////
// methods
/////////////////////////////////////////////////////////////////////////////

///
void
DFIRE_Energy::residue_energy(
	conformation::Residue const &,
	pose::Pose const &,
	EnergyMap &
) const {
	//std::cout << "called residue_energy" << std::endl;
	return;
}

bool DFIRE_Energy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	//std::cout << "testing " << res1 << "," << res2 << " ";
	//if (pose.residue_type(res1).is_protein() && pose.residue_type(res2).is_protein() ) {
	//	std::cout << "true!" << std::endl;
	//} else {
	//	std::cout << "false!" << std::endl;
	//}

	return ( pose.residue_type(res1).is_protein() && pose.residue_type(res2).is_protein() );
}

void
DFIRE_Energy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {
	//std::cout << "called eval_intrares_energy" << std::endl;
	return;
}

void
DFIRE_Energy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
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
} // methods
} // scoring
} // core
