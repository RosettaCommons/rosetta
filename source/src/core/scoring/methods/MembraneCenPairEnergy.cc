// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MembraneCenPairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/MembraneCenPairEnergy.hh>
#include <core/scoring/methods/MembraneCenPairEnergyCreator.hh>

// Package headers
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/Residue.hh>

//symmetry
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>


#include <core/scoring/EnergyMap.hh>
#include <core/scoring/MembraneTopology.hh>
#include <utility/vector1.hh>


// Utility headers


// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MembraneCenPairEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MembraneCenPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MembraneCenPairEnergy );
}

ScoreTypes
MembraneCenPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( Mpair );
	return sts;
}


/// c-tor
MembraneCenPairEnergy::MembraneCenPairEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MembraneCenPairEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_MembranePotential() )
{}


/// clone
EnergyMethodOP
MembraneCenPairEnergy::clone() const
{
	return EnergyMethodOP( new MembraneCenPairEnergy() );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
MembraneCenPairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
	potential_.compute_membrane_embedding( pose );


}


//////////////////////////////////////////////////////////////
//
//     MembraneCenTROID PAIR SCORE
//      and
//     "MembraneCenTROID PACK" SCORE (helps reproduce pairwise correlations
//                            between MembraneCentroids, observed in PDB)
//
void
MembraneCenPairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	Real pair_score( 0.0 ); //, cenpack_score( 0.0 );

	Size rsd1Seq(rsd1.seqpos()), rsd2Seq(rsd2.seqpos());
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;
		SymmetricConformation const & symm_conf (
			dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		if ( !symm_info->bb_is_independent(rsd1.seqpos()) ) {
			rsd1Seq = symm_info->bb_follows(rsd1.seqpos());
		}
		if ( symm_info->is_virtual(rsd1.seqpos()) ) {
			rsd1Seq = 0;
		}

		if ( !symm_info->bb_is_independent(rsd2.seqpos()) ) {
			rsd2Seq = symm_info->bb_follows(rsd2.seqpos());
		}
		if ( symm_info->is_virtual(rsd2.seqpos()) ) {
			rsd2Seq = 0;
		}

	}
	if ( rsd1Seq ==0 || rsd2Seq ==0 ) {
		return;
	}


	if ( MembraneTopology_from_pose( pose ).allow_scoring(rsd1Seq) && MembraneTopology_from_pose( pose ).allow_scoring(rsd2Seq) ) {
		/// assumes centroids are being used
		conformation::Atom const & cen1 ( rsd1.atom( rsd1.nbr_atom() ) ), cen2 (rsd2.atom( rsd2.nbr_atom() ) );
		Real const cendist = cen1.xyz().distance_squared( cen2.xyz() );

		/// accumulate total energies

		potential_.evaluate_pair( pose, rsd1, rsd2, cendist, pair_score);

		//if ( rsd1.aa() == chemical::aa_his && rsd2.aa() == chemical::aa_his && true /*replace with option[ no_his_his_pairE ]*/ ) {
		// pair_score = 0;
		//}

		//pair_score *= 2.019f;
		// cenpack_score *= 2.0f;

		//core::Real rsd_wt = 0.5 *
		// ( get_residue_weight_by_ss( pose.conformation().secstruct( rsd1.seqpos() ) ) +
		//   get_residue_weight_by_ss( pose.conformation().secstruct( rsd2.seqpos() ) )
		// );

		//Rosetta++ used the first residue's weight for both sides of the pair. I hate that. The above
		//comment is an example of an alternative we should probably test in the distant future.
		//bw is this something we like?
		// core::Real rsd_wt =  get_residue_weight_by_ss( pose.conformation().secstruct( rsd1.seqpos() )) ;
	}
	emap[ Mpair ]    += pair_score ; // * rsd_wt;
	//}
	// emap[ cenpack ] += cenpack_score;
}

void
MembraneCenPairEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{
	potential_.finalize( pose );
}
MembraneTopology const &
MembraneCenPairEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
///// @brief MembraneCenPairEnergy distance cutoff
Distance
MembraneCenPairEnergy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 6.0 from cutoffs in MembraneCentroid params files
	//  return 0.0; /// since all the cutoffs for MembraneCentroid mode are rolled into the MembraneCendist check
}
core::Size
MembraneCenPairEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
