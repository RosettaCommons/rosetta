// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNP_LowResEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Kalli Kappel


// Unit headers
#include <core/scoring/rna/RNP_LowResEnergy.hh>
#include <core/scoring/rna/RNP_LowResPotential.hh>
#include <core/scoring/rna/RNP_LowResEnergyCreator.hh>

// Package headers
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//Auto Headers
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is a statistically derived low-resolution potential for RNA/protein interactions
// For RNA/protein modeling, this is meant to supplement the RNA low-res and protein low-res score
// functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "core.scoring.rna.RNP_LowResEnergy" );
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys::score;

namespace core {
namespace scoring {
namespace rna {

/// @details This must return a fresh instance of the RNP_LowResEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNP_LowResEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNP_LowResEnergy );
}

ScoreTypes
RNP_LowResEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rnp_base_pair );
	sts.push_back( rnp_stack );
	sts.push_back( rnp_pair );
	//sts.push_back( rnp_stack_xy );
	sts.push_back( rnp_aa_to_rna_backbone );
	return sts;
}


/// c-tor
RNP_LowResEnergy::RNP_LowResEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RNP_LowResEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_RNP_LowResPotential() )
{
	//std::cout << "Constructed the RNP energy" << std::endl;
}

//clone
methods::EnergyMethodOP
RNP_LowResEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNP_LowResEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
RNP_LowResEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & /*scfxn*/ ) const
{
	//////////////////////////////////////////////////////////////////
	// Need to make all of this smarter, faster, etc.
	// Just a test implementation to see if it's worth the effort to keep it!
	// This will fail if the protein residues are not in centroid mode!
	//////////////////////////////////////////////////////////////////

	//std::cout << "Setting up for scoring!" << std::endl;
	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );

	utility::vector1< bool > & is_interface_residues = rna_scoring_info.nonconst_is_interface();
	utility::vector1< bool > & is_buried_residues = rna_scoring_info.nonconst_is_buried();

	is_interface_residues.clear();
	is_buried_residues.clear();

	core::Real const INTERFACE_CUTOFF( 10.0 );
	core::Real const NEIGHBOR_CUTOFF( 10.0 );
	// Calculate the environment for each residue in the pose
	// this option defaults to false
	bool const use_actual_centroid( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() );
	
	// Figure out the interface residues
	// and buried or not for protein residues
	for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
	        core::Size interface_nbrs = 0;
		core::Size buried_nbrs = 0;

	        // Get the rsd_centroid
	        // use the rna base centroid for RNA
	        Vector rsd_centroid;
	        bool const is_rsd_protein( pose.residue(rsd).is_protein() );
		bool is_centroid = pose.residue(rsd).type().residue_type_set()->name() == core::chemical::CENTROID;
		if ( is_rsd_protein && !is_centroid && !use_actual_centroid ) {
			//TR << "Warning: rnp low res pair energy not computed b/c protein is not centroid" << std::endl;
			return;
		}
	        if (pose.residue(rsd).is_RNA()) {
	                rsd_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd), false /*verbos*/ );
	        } else {
	                if (!pose.residue(rsd).is_protein()) {
				// the residue is not RNA or protein
				// then say 0 interface and 0 buried
				is_interface_residues.push_back( false );
				is_buried_residues.push_back( false );
				// then continue to the next residue in the pose 
				continue;
			}
			if ( !use_actual_centroid ) {
	                	rsd_centroid = pose.residue(rsd).xyz("CEN");
			} else {
				rsd_centroid = pose.residue(rsd).actcoord();
			}
	        }

		// Loop through the rest of the residues in the pose
		// this is going to be really stupid if the pose is big though...
	        for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
	                if ( rsd == rsd2 ) continue;
	                Vector rsd2_centroid;
	                bool const is_rsd2_protein( pose.residue(rsd2).is_protein() );
	                // Only look at RNA/protein and protein/RNA distances here
	                //if ((is_rsd2_protein && is_rsd_protein) || (!is_rsd2_protein && !is_rsd_protein)) continue;
			if ( pose.residue(rsd).is_RNA() && pose.residue(rsd2).is_RNA() ) {
				// then this interaction doesn't count for "interface"
				// would just count for buried, but we dont care about "buried" for RNA
				continue;
			}
	
	                if (pose.residue(rsd2).is_RNA()) {
	                        rsd2_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd2), false /*verbos*/ );
	                } else {
	                        if (!pose.residue(rsd2).is_protein()) continue;
				if ( !use_actual_centroid ) {
	                        	rsd2_centroid = pose.residue(rsd2).xyz("CEN");
				} else {
					rsd2_centroid = pose.residue(rsd2).actcoord();
				}
	                        //rsd2_centroid = pose.residue(rsd2).actcoord();
	                }
	                core::Real distance = ( rsd_centroid - rsd2_centroid ).length();
	                if ((is_rsd2_protein && !is_rsd_protein) || (!is_rsd2_protein && is_rsd_protein)) {
	                	if ( distance < INTERFACE_CUTOFF ) {
	                        	++interface_nbrs;
	                	}
			}
			else if ( is_rsd2_protein && is_rsd_protein ) {
			//else if ((is_rsd2_protein && is_rsd_protein) || (!is_rsd2_protein && !is_rsd_protein)) {
				if ( distance < NEIGHBOR_CUTOFF ) {
					++buried_nbrs;
				}
			}
	        }
	        // Write out the number of neighbors for this residue
	        if ( interface_nbrs > 0 ) {
			is_interface_residues.push_back( true );
	        } else {
			is_interface_residues.push_back( false );
	        }
		// not buried if the residue is RNA
		if ( buried_nbrs > 10 ) {
			is_buried_residues.push_back( true );
		} else {
			is_buried_residues.push_back( false );
		}
	}
	//std::cout << "Done setting up for scoring" << std::endl;

}

//void
//RNP_LowResEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
//{
//}

//////////////////////////////////////////////////////////////////////////////////////////
void
RNP_LowResEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/, // need this back for rnp_pair
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// Only evaluate these score terms between RNA and protein residues
	if ( !(( rsd1.is_RNA() && rsd2.is_protein() ) || ( rsd1.is_protein() && rsd2.is_RNA() )) ) return;

	// this option is false by default
	bool const use_actual_centroid( basic::options::option[ basic::options::OptionKeys::score::FA_low_res_rnp_scoring ]() );

	// Just give a warning and return if the protein isn't coarse
	if ( rsd1.is_protein() ) {
		bool is_centroid = rsd1.type().residue_type_set()->name() == core::chemical::CENTROID;
		if ( !is_centroid && !use_actual_centroid ) {
			//TR << "Warning: rnp low res pair energy not computed b/c protein is not centroid" << std::endl;
			return;
		}
	} else { // rsd2 is protein
		bool is_centroid = rsd2.type().residue_type_set()->name() == core::chemical::CENTROID;
		if ( !is_centroid && !use_actual_centroid ) {
			//TR << "Warning: rnp low res pair energy not computed b/c protein is not centroid" << std::endl;
			return;
		}
	}

	Vector rna_centroid;
	Vector protein_centroid;
	numeric::xyzMatrix< core::Real > rna_base_coord_sys;
	
	if ( rsd1.is_RNA() ) {
		rna_centroid = chemical::rna::get_rna_base_centroid( rsd1 );
		rna_base_coord_sys = chemical::rna::get_rna_base_coordinate_system( rsd1, rna_centroid );
		if ( !use_actual_centroid ) {
			protein_centroid = rsd2.xyz( "CEN" );
		} else {
			protein_centroid = rsd2.actcoord();
		}
	} else {
		rna_centroid = chemical::rna::get_rna_base_centroid( rsd2 );
		rna_base_coord_sys = chemical::rna::get_rna_base_coordinate_system( rsd2, rna_centroid );
		if ( !use_actual_centroid ) {
			protein_centroid = rsd1.xyz( "CEN" );
		} else {
			protein_centroid = rsd1.actcoord();
		}
	}
	Vector x_rna = rna_base_coord_sys.col_x();
	Vector y_rna = rna_base_coord_sys.col_y();
	Vector z_rna = rna_base_coord_sys.col_z();

	Vector dist_rna_protein = protein_centroid - rna_centroid;
	Real const dist_x = dot_product( dist_rna_protein, x_rna );
	Real const dist_y = dot_product( dist_rna_protein, y_rna );
	Real const dist_z = dot_product( dist_rna_protein, z_rna );

	Real rnp_base_pair_score( 0.0 );
	Real const max_interaction_dist_x( 10.0 );
	Real const max_interaction_dist_y( 10.0 );
	if ( std::abs(dist_z) < 3.0 ) {
		if ( std::abs( dist_x ) < max_interaction_dist_x && std::abs( dist_y ) < max_interaction_dist_y ) {
			potential_.evaluate_rnp_base_pair_score( rsd1, rsd2, dist_x, dist_y, rnp_base_pair_score ); // just for now, obviously stupid
		}
	}
	emap[ rnp_base_pair ] += rnp_base_pair_score;

//	// Get the stack score
//	Real rnp_stack_xy_score( 0.0 );
//	if ( std::abs(dist_z) > 3.0 && std::abs(dist_z) < 6.5 ) {
//		potential_.evaluate_rnp_stack_xy_score( rsd1, rsd2, dist_x, dist_y, rnp_stack_xy_score );
//	}
//	emap[ rnp_stack_xy ] += rnp_stack_xy_score;

	// stack score (not xy)
	bool calculate_rnp_stack = false;
	// check protein residue type, stacking only occurs with certain restypes
	if ( rsd1.is_protein() ) {
		// check if this is TRP, TYR, PHE (perhaps expand this further...)
		if ( rsd1.name1() == 'W' || rsd1.name1() == 'Y' || rsd1.name1() == 'F' || rsd1.name1() == 'V' || rsd1.name1() == 'L') {
		//if ( rsd1.name1() == 'W' || rsd1.name1() == 'Y' || rsd1.name1() == 'F' ) {
			calculate_rnp_stack = true;
		}
	} else { // rsd2 is protein
		//if ( rsd2.name1() == 'W' || rsd2.name1() == 'Y' || rsd2.name1() == 'F' ) {
		if ( rsd2.name1() == 'W' || rsd2.name1() == 'Y' || rsd2.name1() == 'F' || rsd2.name1() == 'V' || rsd2.name1() == 'L') {
			calculate_rnp_stack = true;
		}
	}
	
	if ( calculate_rnp_stack ) {
		if ( std::abs( dist_z ) > 3.0 && std::abs( dist_z ) < 6.5 && (dist_x*dist_x + dist_y*dist_y) < 16.0 ) {
			emap[ rnp_stack ] += -1.0;
		}
	}

	// Get the rnp backbone score
	Vector rna_backbone_P_xyz;
	if ( rsd1.is_RNA() ) {
		// get the backbone phosphate position
		rna_backbone_P_xyz = rsd1.xyz( " P  " );
	} else {
		rna_backbone_P_xyz = rsd2.xyz( " P  " );
	}

	Real dist_to_backbone = (protein_centroid - rna_backbone_P_xyz).length();
	Real const min_backbone_interaction_dist = 3.0;
	Real const max_backbone_interaction_dist = 10.0;

	Real backbone_score = 0.0;

	if ( dist_to_backbone >= min_backbone_interaction_dist && dist_to_backbone <= max_backbone_interaction_dist ) {
		if ( rsd1.is_protein() ) {
			backbone_score = potential_.evaluate_rnp_aa_rna_backbone_score( rsd1, dist_to_backbone );
		} else {
			backbone_score = potential_.evaluate_rnp_aa_rna_backbone_score( rsd2, dist_to_backbone );
		}
	}
	
	emap[ rnp_aa_to_rna_backbone ] += backbone_score;

//	// Get the pair score
//	Real rnp_pair_score( 0.0 );
//	
//	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
//	utility::vector1< bool > const & is_interface = rna_scoring_info.is_interface();
//	utility::vector1< bool > const & is_buried = rna_scoring_info.is_buried();
//
//	bool const rsd1_is_interface = is_interface[ rsd1.seqpos() ];
//	bool const rsd1_is_buried = is_buried[ rsd1.seqpos() ];
//
//	bool const rsd2_is_interface = is_interface[ rsd2.seqpos() ];
//	bool const rsd2_is_buried = is_buried[ rsd2.seqpos() ];
//
//	//std::cout << "About to evaluate rnp pair score" << std::endl;
//	potential_.evaluate_rnp_pair_score( rsd1, rsd2, rsd1_is_interface, rsd1_is_buried, 
//			rsd2_is_interface, rsd2_is_buried, dist_rna_protein.length(), rnp_pair_score );
//	//std::cout << "Done evaluating rnp pair score" << std::endl;
//
//	emap[ rnp_pair ] += rnp_pair_score;
//
//	// TR << rsd1.name3()  << rsd1.seqpos() << "---" << rsd2.name3() << rsd2.seqpos() << ": " << (score1+score2) << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////////
//void
//RNP_LowResEnergy::finalize_total_energy(
//	pose::Pose & pose,
//	ScoreFunction const &,
//	EnergyMap &
//) const {
//
//}

//??
/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNP_LowResEnergy::atomic_interaction_cutoff() const
{
	// For testing 
	return 0.0; /// Uh, I don't know.
}

core::Size
RNP_LowResEnergy::version() const
{
	return 1; // Initial versioning
}

//etable::count_pair::CountPairFunctionCOP
//RNP_LowResEnergy::get_intrares_countpair(
//	conformation::Residue const &,
//	pose::Pose const &,
//	ScoreFunction const &
//) const
//{
//	utility_exit_with_message( "FA_ElecEnergy does not define intra - residue pair energies; do not call get_intrares_countpair()" );
//	return 0;
//}
//
//etable::count_pair::CountPairFunctionCOP
//RNP_LowResEnergy::get_count_pair_function(
//	Size const res1,
//	Size const res2,
//	pose::Pose const & pose,
//	ScoreFunction const &
//) const
//{
//}
//
//
//etable::count_pair::CountPairFunctionCOP
//RNP_LowResEnergy::get_count_pair_function(
//	conformation::Residue const & rsd1,
//	conformation::Residue const & rsd2
//) const
//{
//}


} //rna
} //scoring
} //core
