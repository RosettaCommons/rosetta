// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/scoring/methods/EnsembleEnergy.cc
/// @brief  Function to subtract out individual reference energies during ensemble docking
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/methods/EnsembleEnergy.hh>
#include <protocols/scoring/methods/EnsembleEnergyCreator.hh>

// Package headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>

#include <core/pose/datacache/CacheableDataType.hh>

// Project headers
#include <core/pose/Pose.hh>


// Utility headers

#include <basic/Tracer.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.scoring.methods.EnsembleEnergy" );


// C++

namespace protocols {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the EnsembleEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
EnsembleEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new EnsembleEnergy );
}

core::scoring::ScoreTypes
EnsembleEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( dock_ens_conf );
	return sts;
}


/// c-tor
EnsembleEnergy::EnsembleEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new EnsembleEnergyCreator ) )
{}


/// clone
core::scoring::methods::EnergyMethodOP
EnsembleEnergy::clone() const
{
	return core::scoring::methods::EnergyMethodOP( new EnsembleEnergy() );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

// calculate the score of the residues to constrain to the interface
void
EnsembleEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core;
	using namespace core::scoring;
	Real unbound_score ( 0.0 );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		basic::datacache::CacheableStringFloatMapOP data
			= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringFloatMap >
			( pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		for ( std::map< std::string, float >::const_iterator iter = data->map().begin(),
				end = data->map().end(); iter != end; ++iter
				) {
			// skip score entry, as it gets confusing
			if ( iter->first == "dock_ens_conf1" || iter->first == "dock_ens_conf2" ) {
				//      TR << iter->first << ' ' << iter->second << std::endl;
				unbound_score += iter->second;
			}
		}

	}

	// TR << "unbound score: " << unbound_score << std::endl;

	emap[ dock_ens_conf ] = -unbound_score;
}

core::Size
EnsembleEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
