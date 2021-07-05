// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/EnvPairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/energy_methods/MembraneCbetaEnergy.hh>
#include <core/energy_methods/MembraneCbetaEnergyCreator.hh>

// Package headers
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/EnergyMap.hh>


// Utility headers

// C++

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the MembraneCbetaEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
MembraneCbetaEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< MembraneCbetaEnergy >();
}

core::scoring::ScoreTypes
MembraneCbetaEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( Mcbeta );
	return sts;
}


/// c-tor
MembraneCbetaEnergy::MembraneCbetaEnergy() :
	parent( utility::pointer::make_shared< MembraneCbetaEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_MembranePotential() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
MembraneCbetaEnergy::clone() const
{
	return utility::pointer::make_shared< MembraneCbetaEnergy >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
MembraneCbetaEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );

	//No need to calculate embedding...
}


///////////////////////////////////////
//
// ENV SCORE
void
MembraneCbetaEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const
{
	//Real env_score( 0.0 ),
	Real cb_score( 0.0 ); //, cb_score12( 0.0 ), cb_score( 0.0 );
	if ( rsd.aa() == core::chemical::aa_vrt ) {
		emap[ core::scoring::Mcbeta ] = cb_score;// * rsd_wt;
		return;
	}

	if ( MembraneTopology_from_pose( pose ).allow_scoring(rsd.seqpos()) ) {

		// potential_.evaluate_env( pose, rsd, env_score);
		potential_.evaluate_cbeta( pose, rsd,cb_score);

		//env_score *= 2.019;
		cb_score = 2.667 * ( cb_score ) * 0.3;
	}
	//bw is this something we like?
	//core::Real rsd_wt = core::scoring::methods::get_residue_weight_by_ss( pose.conformation().secstruct( rsd.seqpos() ) );

	///emap[ core::scoring::Menv   ] = env_score;// * rsd_wt;
	//std::cout << "CB " << MembraneTopology_from_pose( pose ).allow_scoring(rsd.seqpos()) << " "  << rsd.seqpos() << " " << cb_score << "\n";
	emap[ core::scoring::Mcbeta ] += cb_score;// * rsd_wt;
} // residue_energy

void
MembraneCbetaEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{
	potential_.finalize( pose );
}

/*
bool
MembraneCbetaEnergy::allow_scoring(Size const seqpos) const
{

return ((*( static_cast< core::scoring::MembraneTopology const * >( pose.data().get_const_ptr( basic::MEMBRANE_TOPOLOGY )() ))).allow_scoring(seqpos));
}
*/
core::scoring::MembraneTopology const &
MembraneCbetaEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
core::Size
MembraneCbetaEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
