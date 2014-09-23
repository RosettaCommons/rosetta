// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenPairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/util.hh>
#include <core/scoring/methods/CenPairEnergy.hh>
#include <core/scoring/methods/CenPairEnergyCreator.hh>

// Package headers
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



// Utility headers



// C++


namespace core {
namespace scoring {
namespace methods {

methods::EnergyMethodOP
CenPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new methods::CenPairEnergy;
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
ScoreTypes
CenPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( pair );
	sts.push_back( cenpack );
	return sts;
}


/// c-tor
CenPairEnergy::CenPairEnergy() :
	parent( methods::EnergyMethodCreatorOP( new CenPairEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_EnvPairPotential() )
{}


/// clone
EnergyMethodOP
CenPairEnergy::clone() const
{
	return new CenPairEnergy();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


///
void
CenPairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	potential_.compute_centroid_environment( pose );
	pose.update_residue_neighbors();
}


//////////////////////////////////////////////////////////////
//
//     CENTROID PAIR SCORE
//      and
//     "CENTROID PACK" SCORE (helps reproduce pairwise correlations
//                            between centroids, observed in PDB)
//
void
CenPairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ){
			return;
	}
	if(rsd1.aa()==core::chemical::aa_unk) return;
	if(rsd2.aa()==core::chemical::aa_unk) return;

	/// assumes centroids are being used
	conformation::Atom const & cen1 ( rsd1.atom( rsd1.nbr_atom() ) ), cen2 (rsd2.atom( rsd2.nbr_atom() ) );
	Real const cendist = cen1.xyz().distance_squared( cen2.xyz() );

	//fpd ignore cen-cen distances above 12.05A
	if (cendist > 12.05*12.05) return;

	/// accumulate total energies
	Real pair_score( 0.0 ), cenpack_score( 0.0 );
	potential_.evaluate_pair_and_cenpack_score( rsd1, rsd2, cendist,
		pair_score, cenpack_score );

	//if ( rsd1.aa() == chemical::aa_his && rsd2.aa() == chemical::aa_his && true /*replace with option[ no_his_his_pairE ]*/ ) {
	//	pair_score = 0;
	//}

	pair_score *= 2.019f;
	cenpack_score *= 2.0f;

	//core::Real rsd_wt = 0.5 *
	//	( get_residue_weight_by_ss( pose.conformation().secstruct( rsd1.seqpos() ) ) +
	//	  get_residue_weight_by_ss( pose.conformation().secstruct( rsd2.seqpos() ) )
	//	);

	//Rosetta++ used the first residue's weight for both sides of the pair. I hate that. The above
	//comment is an example of an alternative we should probably test in the distant future.
	core::Real rsd_wt =  get_residue_weight_by_ss( pose.conformation().secstruct( rsd1.seqpos() )) ;

	emap[ pair ]    += pair_score    * rsd_wt;
	emap[ cenpack ] += cenpack_score;
}

void
CenPairEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{
	potential_.finalize( pose );
}

/// @brief CenPairEnergy distance cutoff
Distance
CenPairEnergy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 6.0 from cutoffs in centroid params files
// 	return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}
core::Size
CenPairEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
