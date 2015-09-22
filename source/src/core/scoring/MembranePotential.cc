// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/MembranePotential.cc
///
/// @brief  Membrane Potential - Base Scoring Methods for LowRes Energy Function
/// @details Compute Low Res membrane energy terms: Menv, MPair, MCBeta and Membrane
///    penalties. Also contains pass-through methods for accessing and updating
///    mp framework supported data in a membrane conformation.
///    Last Modified: 3/11/14
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  Bjorn Wallner (Original)

#ifndef INCLUDE_core_scoring_MembranePotential_cc
#define INCLUDE_core_scoring_MembranePotential_cc

// Unit Headers
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/EnvPairPotential.hh>

// Project Headers
#include <core/scoring/MembraneTopology.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <basic/database/open.hh>

#include <core/pose/Pose.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Symmetry Package Headers
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Resource Manager Headers
#include <basic/resource_manager/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.MembranePotential" );

namespace core {
namespace scoring {

/// @brief Copy Constructor for a Membrane Embedding Object
MembraneEmbed::MembraneEmbed( MembraneEmbed const & src ) :
	CacheableData()
{
	depth_=src.depth_;
	center_=src.center_;
	penalty_=src.penalty_;
	normal_=src.normal_;
	spanning_=src.spanning_;
	calculated_ = src.calculated_;
}

/// @brief Initialize a Membrane Embedding Object
void
MembraneEmbed::initialize( pose::Pose const & pose )
{
	Size const nres( pose.total_residue() );
	depth_.resize(nres,0.0);
	std::fill(depth_.begin(),depth_.end(),30.0);
	center_.assign(0,0,0);
	normal_.assign(0,0,1);
}

////////////////////////// This is an Entire Constructor - Might want to Split for Aesthetics ///////////////////////////////////////

/// @brief Construct a Membrane Potential Object - Acts as the base layer for membrane protein scoring
MembranePotential::MembranePotential():
	//cems transition regions between environment bins
	//cems transition is from +/- sqrt(36+pad6) +/- sqrt(100+pad10) etc
	cen_dist5_pad( 0.5 ),
	// unused cen_dist6_pad( 0.6 ),
	cen_dist7_pad( 0.65 ),
	cen_dist10_pad( 1.0 ),
	cen_dist12_pad( 1.2 ),

	cen_dist5_pad_plus ( cen_dist5_pad  + 25.0 ),
	//cen_dist6_pad_plus( cen_dist6_pad + 36.0 ),
	cen_dist7_pad_plus ( cen_dist7_pad  + 56.25 ),
	cen_dist10_pad_plus( cen_dist10_pad + 100.0 ),
	// cen_dist12_pad_plus( cen_dist12_pad + 144.0 ),

	cen_dist5_pad_minus ( cen_dist5_pad  - 25.0 ),
	cen_dist7_pad_minus ( cen_dist7_pad  - 56.25 ),
	cen_dist10_pad_minus( cen_dist10_pad - 100.0 ),
	cen_dist12_pad_minus( cen_dist12_pad - 144.0 ),

	cen_dist5_pad_hinv ( 0.5 / cen_dist5_pad ),
	// cen_dist6_pad_hinv ( 0.5 / cen_dist6_pad ),
	cen_dist7_pad_hinv ( 0.5 / cen_dist7_pad ),
	cen_dist10_pad_hinv( 0.5 / cen_dist10_pad ),
	cen_dist12_pad_hinv( 0.5 / cen_dist12_pad )
	// calculated_(false)
{

	// Initialize Database Info
	Size const max_aa( 20 ); // just the standard aa's for now
	Size const env_log_table_cen6_bins( 15 );
	Size const env_log_table_cen10_bins( 40 );
	Size const pair_log_table_size( 5 );
	Size const cbeta_den_table_size( 45 );
	Size const max_mem_layers( 3 );
	Size const min_mem_layers( 2 );

	///////OPTIONS

	// This will be replaced by a scoring object that will be loaded via the resource manager
	membrane_normal_cycles_=basic::options::option[basic::options::OptionKeys::membrane::normal_cycles]();
	membrane_normal_magnitude_=numeric::conversions::radians(basic::options::option[ basic::options::OptionKeys::membrane::normal_mag ]());
	membrane_center_magnitude_=basic::options::option[ basic::options::OptionKeys::membrane::center_mag ]();
	smooth_move_frac_=basic::options::option[ basic::options::OptionKeys::membrane::smooth_move_frac ]();
	no_interpolate_Mpair_=basic::options::option[ basic::options::OptionKeys::membrane::no_interpolate_Mpair ]();
	Menv_penalties_=basic::options::option[ basic::options::OptionKeys::membrane::Menv_penalties ]();
	Membed_init_=basic::options::option[ basic::options::OptionKeys::membrane::Membed_init ]();
	memb_center_search_=basic::options::option[ basic::options::OptionKeys::membrane::center_search ]();
	memb_normal_search_=basic::options::option[ basic::options::OptionKeys::membrane::normal_search ]();
	membrane_center_max_delta_=basic::options::option[basic::options::OptionKeys::membrane::center_max_delta]();
	membrane_normal_start_angle_=basic::options::option[basic::options::OptionKeys::membrane::normal_start_angle]();
	membrane_normal_delta_angle_=basic::options::option[basic::options::OptionKeys::membrane::normal_delta_angle]();
	membrane_normal_max_angle_=basic::options::option[basic::options::OptionKeys::membrane::normal_max_angle]();

	/////

	// Initialize the Membrane Etable Data from within
	std::string tag,line;
	chemical::AA aa;

	{ // mem_env_cen6:
		mem_env_log6_.dimension( max_aa, max_mem_layers,env_log_table_cen6_bins );
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt" );
		for ( Size i=1; i<= max_aa; ++i ) {
			getline( stream, line );
			std::istringstream l(line);
			l >> tag >> aa;
			if ( l.fail() || tag != "MEM_ENV_LOG_CEN6:"  ) {
				utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt (cen6)");
			}
			for ( Size j=1; j<=max_mem_layers; ++j ) {
				for ( Size k=1; k<=env_log_table_cen6_bins; ++k ) {
					l >> mem_env_log6_(aa,j,k);
				}
			}
		}
	}
	{ // mem_env_cen10:
		mem_env_log10_.dimension( max_aa, max_mem_layers,env_log_table_cen10_bins );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt" );
		for ( Size i=1; i<= max_aa; ++i ) {
			getline( stream, line );
			std::istringstream l(line);
			l >> tag >> aa;
			if ( l.fail() || tag != "MEM_ENV_LOG_CEN10:"  ) {
				utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt (cen10)");
			}
			for ( Size j=1; j<=max_mem_layers; ++j ) {
				for ( Size k=1; k<=env_log_table_cen10_bins; ++k ) {
					l >> mem_env_log10_(aa,j,k);
				}
			}
		}
	}

	{ // cbeta_den_6/12
		mem_cbeta_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_den12_.dimension( cbeta_den_table_size );
		mem_cbeta_2TM_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_2TM_den12_.dimension( cbeta_den_table_size );
		mem_cbeta_4TM_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_4TM_den12_.dimension( cbeta_den_table_size );


		// Read in etables into private vars
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/memcbeta_den.txt" );

		{ // den6
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >>mem_cbeta_den6_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN6)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_2TM_den6_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_2TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN6)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_4TM_den6_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_4TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN6)");
			}
		}

		{
			{ // den12
				getline( stream, line );
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_den12_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN12)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_2TM_den12_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_2TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN12)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_4TM_den12_(i);
				}
				if ( l.fail() || tag != "MEMCBETA_4TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN12)");
			}

		}
	}

	{ // pair_log
		mem_pair_log_.dimension( min_mem_layers,pair_log_table_size,max_aa, max_aa );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/mem_pair_log.txt" );
		for ( Size i=1; i<= min_mem_layers; ++i ) {
			for ( Size j=1; j<= pair_log_table_size; ++j ) {
				for ( Size k=1; k<= max_aa; ++k ) {
					getline( stream, line );
					std::istringstream l(line);
					Size ii,jj;
					l >> tag >> ii >> jj >> aa;
					assert( Size(aa) == k );
					for ( Size n=1; n<=max_aa; ++n ) {
						l >> mem_pair_log_(i,j,aa,n);
					}
					if ( l.fail() || ii != i || jj != j || tag != "MEM_PAIR_LOG:"  ) utility_exit_with_message("scoring/score_functions/MembranePotential/mem_pair_log.txt");
				}

			}
		}
	}
}

/// @brief Fianlize Setup of MP Scoring
void
MembranePotential::finalize( pose::Pose & pose ) const
{
	CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;
	MembraneEmbed & membrane_embed( nonconst_MembraneEmbed_from_pose( pose ));
	membrane_embed.calculated() = false;
}

////////////////////////////////////////////////////////////////////////////////////

/// @brief Evaluate Membrane Environment Term (by residue)
void
MembranePotential::evaluate_env(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real & membrane_env_score
) const
{
	if ( MembraneEmbed_from_pose( pose ).spanning() ) {
		Real termini_pen(0);
		Real const MembraneDepth (MembraneEmbed_from_pose( pose ).depth(rsd.seqpos() ) );
		evaluate_env(pose,rsd,MembraneDepth,membrane_env_score);
		Vector const normal(MembraneEmbed_from_pose( pose ).normal());
		Vector const center(MembraneEmbed_from_pose( pose ).center());
		if ( Menv_penalties_ && ( rsd.seqpos()==1 || rsd.seqpos()==pose.total_residue() ) ) {
			Vector const & xyz( pose.residue(rsd.seqpos()).atom( 2 ).xyz() );
			Real depth=dot(xyz-center,normal)+30;
			if ( depth>18 &&
					depth<42 ) {
				termini_pen++;
			}
			membrane_env_score+=50*termini_pen;
		}
	} else {
		membrane_env_score=100;
	}
}

////////////////////////////////////////////////////////////////////////////////////

/// @brief Full Evaluate Membrane Env Method (given depth and const menv score)
void
MembranePotential::evaluate_env(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real const MembraneDepth,
	Real & membrane_env_score
) const {
	Real t2 = 2.0;
	Real t3 = 2.0;
	//int layer1, layer2, layer;
	Real f, z, zn, low;

	Real const env6_weight=1.0;
	Real const env10_weight=1.0;

	Real fcen6  ( cenlist_from_pose( pose ).fcen6( rsd.seqpos() ) );
	Real fcen10 ( cenlist_from_pose( pose ).fcen10( rsd.seqpos() ) );

	// in rare cases, the density is over 15 within 6
	if ( fcen6 > 15 ) fcen6 = 15 ;
	if ( fcen10 > 40 ) fcen10 = 40;

	if ( rsd.is_protein() ) {

		if ( ( MembraneDepth < 11.0 ) || ( MembraneDepth > 49.0 ) ) {
			//pure water layer
			int layer = 3;
			//B_layer = 1;  // set but never used ~Labonte
			Real score6 (env6_weight*mem_env_log6_( rsd.aa(), layer, static_cast< int >( fcen6 ) ));
			Real score10 (env10_weight* mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ) );
			membrane_env_score = score6 + score10;

		} else if ( ( MembraneDepth >= 11.0 && MembraneDepth <= 13.0 ) || ( MembraneDepth >= 47.0 && MembraneDepth <= 49.0 ) ) {
			//interpolate between water and interface phases
			int layer1 = 2; //interface layer
			int layer2 = 3; //water layer
			//B_layer = 1;  // set but never used ~Labonte

			if ( MembraneDepth <= 13.0 ) {
				low = 13.0;
			} else {
				low = 47.0;
			}
			z = 2*std::abs( (MembraneDepth - low) ) / t2;
			int s2 = 14;
			zn = std::pow( z, s2 );
			f = zn/(1 + zn);

			Real score6_layer2( env6_weight*mem_env_log6_( rsd.aa(), layer2, static_cast< int >( fcen6 ) ) );
			Real score10_layer2( env10_weight* mem_env_log10_( rsd.aa(), layer2, static_cast< int >( fcen10 ) ) );
			Real score6_layer1( env6_weight*mem_env_log6_( rsd.aa(), layer1, static_cast< int >( fcen6 ) ) );
			Real score10_layer1( env10_weight*mem_env_log10_( rsd.aa(), layer1, static_cast< int >( fcen10 ) ) );

			membrane_env_score = f * ( score6_layer2 + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );

			// amw what's the point of this? it immediately leaves scope and wasn't being used
			// in higher scope either.
			//int layer = ( MembraneDepth <= 12.0 || MembraneDepth >= 48.0 ) ? 2 : 3;

		} else if ( ( MembraneDepth > 13.0 && MembraneDepth < 17.0 ) || ( MembraneDepth > 43.0 && MembraneDepth < 47.0 ) ) {
			//pure interface phase
			int layer = 2; //interface layer
			//B_layer = 1;  // set but never used ~Labonte
			Real score6 ( env6_weight*mem_env_log6_( rsd.aa(), layer, static_cast< int >( fcen6 ) ) );
			Real score10 ( env10_weight*mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ) );
			membrane_env_score = score6 + score10;
		} else if ( ( MembraneDepth >= 17.0 && MembraneDepth <= 19.0 ) || ( MembraneDepth >= 41.0 && MembraneDepth <= 43.0 ) ) {
			//interpolate between interface and hydrophobic phases
			int layer1 = 1; //hydrophobic layer
			int layer2 = 2; //interface layer

			if ( MembraneDepth <= 19.0 ) {
				low = 19.0;
			} else {
				low = 41.0;
			}
			z = 2*std::abs( (MembraneDepth - low) ) / t3;
			int s3 = 14;
			zn = std::pow( z, s3 );
			f = zn/(1 + zn);

			Real score6_layer2(env6_weight*mem_env_log6_( rsd.aa(), layer2, static_cast< int >( fcen6 )));
			Real score10_layer2( env10_weight*mem_env_log10_( rsd.aa(), layer2, static_cast< int >( fcen10 ) ) );
			Real score6_layer1( env6_weight*mem_env_log6_( rsd.aa(), layer1, static_cast< int >( fcen6 ) ) );
			Real score10_layer1( env10_weight*mem_env_log10_( rsd.aa(), layer1, static_cast< int >( fcen10 ) ) );

			membrane_env_score = f * ( score6_layer2  + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );

			//int layer = ( MembraneDepth <= 18.0 || MembraneDepth >= 42.0 ) ? 2 : 1;

		} else {
			//pure hydrophobic phase
			int layer = 1;

			Real score6  (env6_weight *mem_env_log6_(  rsd.aa(), layer, static_cast< int >( fcen6  ) ));
			Real score10 (env10_weight*mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ));
			membrane_env_score = score6+score10;
		}

		membrane_env_score*=0.5; //bw membrane_embed_weight...

	} else { // amino acid check
		membrane_env_score = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
MembranePotential::evaluate_cbeta(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real & membrane_cb_score
) const
{
	membrane_cb_score=0;
	Real const fcen6 ( cenlist_from_pose( pose ).fcen6( rsd.seqpos() ) );
	Real const fcen12 ( cenlist_from_pose( pose ).fcen12( rsd.seqpos() ) );

	Real membrane_cb_score6,membrane_cb_score12;

	Size const TMHs (MembraneTopology_from_pose( pose ).tmh_inserted() );
	// interp1 rounds down to nearest (non-negative) integer.
	int const interp1 = static_cast< int >( fcen6 );
	int const interp3 = static_cast< int >( fcen12 );
	// note cen6 is always at least 1.0
	// lower bound
	Real const interp2 = fcen6-interp1;
	Real const interp4 = fcen12-interp3;
	if ( TMHs <= 2 ) {
		membrane_cb_score6 =
			( 1.0-interp2 ) * mem_cbeta_2TM_den6_( interp1 )+
			interp2         * mem_cbeta_2TM_den6_( interp1+1 );

	} else if ( TMHs <= 4 ) {
		membrane_cb_score6 =
			(1.0-interp2) * mem_cbeta_4TM_den6_( interp1 )+
			interp2       * mem_cbeta_4TM_den6_( interp1+1 );

	} else {
		membrane_cb_score6 =
			(1.0-interp2) * mem_cbeta_den6_( interp1 )+
			interp2       * mem_cbeta_den6_( interp1+1 );
	}
	membrane_cb_score12 =
		(1.0-interp4) * mem_cbeta_den12_( interp3 )+
		interp4       * mem_cbeta_den12_( interp3+1 );

	membrane_cb_score = (membrane_cb_score6+membrane_cb_score12);
}

///////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Evaluate Membrane Rsd Pair Term
void
MembranePotential::evaluate_pair(
	pose::Pose const & pose,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real & membrane_pair_score
) const
{

	membrane_pair_score = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	//CAR  no pair score if a disulfide
	if ( aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
			rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
			rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ <= 8 ) return;

	//$$$  we now try to find which bin the pair distance lies in
	//$$$  I note this could in principle be calculated and updatded
	//$$$  just like cen_dist is if there is a need for speed.
	//$$$  this function interpolates between bins.
	//$$$  An important(!) requirement on pair_log is that the
	//$$$  value should approach zero as the radius increases.
	//$$$  this fact permits us not to have to compute and score pairs are larger
	//$$$  than cen_dist > cutoff.

	int icon = 5;
	Real interp2( 0.0 );
	//Real interp1( 0.0);
	Real const MembraneDepth1 (MembraneEmbed_from_pose( pose ).depth(rsd1.seqpos() ) );
	Real const MembraneDepth2 (MembraneEmbed_from_pose( pose ).depth(rsd2.seqpos() ) );

	int hydro_layer=1;  //1 not_hydrophobic_core 2 hydrophobic core
	Real AverageDepth=(MembraneDepth1+MembraneDepth2)/2;
	if ( MembraneDepth1 > 18 &&
			MembraneDepth1 < 42 &&
			MembraneDepth2 >18 &&
			MembraneDepth2 <42 ) { //bw currently both residues have to be in the hydrophobic core
		hydro_layer=2;
	}

	if ( cendist > cen_dist10_pad_plus ) {
		icon = 4;
		interp2 = ( cendist + cen_dist12_pad_minus ) * cen_dist12_pad_hinv;
	} else {
		if ( cendist > cen_dist7_pad_plus ) {
			icon = 3;
			interp2 = ( cendist + cen_dist10_pad_minus ) * cen_dist10_pad_hinv;
		} else {
			if ( cendist > cen_dist5_pad_plus ) {
				icon = 2;
				interp2 = ( cendist + cen_dist7_pad_minus ) * cen_dist7_pad_hinv;
			} else {
				icon = 1;
				interp2 = ( cendist + cen_dist5_pad_minus ) * cen_dist5_pad_hinv;
			}
		}
	}
	if ( interp2 < 0.0 ) interp2 = 0.0;

	// note in theory this will never happen but in practice round off
	// error can cause problem
	if ( interp2 > 1.0 ) interp2 = 1.0;
	// handle last bin specially since icon+1 would be past array end
	Real f(0);
	if ( icon != 5 ) {
		if ( !no_interpolate_Mpair_ ) { //bw new mini specfic, true by default.

			if ( std::abs(AverageDepth - 18)<4 ) {
				f=1/(1+std::exp(1.5*(18-AverageDepth)));
				membrane_pair_score = ( ( 1.0f - interp2 ) * ((1-f)*mem_pair_log_( 1, icon  , aa1, aa2 ) + f*mem_pair_log_( 2, icon  , aa1, aa2 )) +
					(       interp2 ) *  ((1-f)*mem_pair_log_( 1, icon+1, aa1, aa2 ) + f*mem_pair_log_( 2, icon+1, aa1, aa2 )));
			} else if ( std::abs(AverageDepth - 42)<4 ) {

				f=1/(1+std::exp(1.5*(AverageDepth-42)));

				membrane_pair_score = ( ( 1.0f - interp2 ) * ((1-f)*mem_pair_log_( 1, icon  , aa1, aa2 ) + f*mem_pair_log_( 2, icon  , aa1, aa2 )) +
					(       interp2 ) *  ((1-f)*mem_pair_log_( 1, icon+1, aa1, aa2 ) + f*mem_pair_log_( 2, icon+1, aa1, aa2 )));

			}  else {

				membrane_pair_score = ( ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer, icon  , aa1, aa2 ) +
					(       interp2 ) *  mem_pair_log_( hydro_layer, icon+1, aa1, aa2 ));
			}
		} else {
			membrane_pair_score = ( ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer, icon  , aa1, aa2 ) +
				(       interp2 ) *  mem_pair_log_( hydro_layer, icon+1, aa1, aa2 ));
		}
	} else {
		membrane_pair_score =   ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer,icon  , aa1, aa2 );
	}
	membrane_pair_score*=2.019;//why is this multiplied by 2.019???
	return;
}

// duplicated in Embedding Factory and Embedding Residues
// part of it should be replaced
void
MembranePotential::compute_membrane_embedding(pose::Pose & pose) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	MembraneEmbed & membrane_embed(nonconst_MembraneEmbed_from_pose( pose ));

	//pba
	if ( !nonconst_cenlist_from_pose( pose ).calculated() ) {
		nonconst_cenlist_from_pose( pose ).initialize( pose );
	}

	//pba
	//read in spanfile
	MembraneTopology & topology( nonconst_MembraneTopology_from_pose(pose) );
	if ( !topology.initialized() ) {
		if ( option[in::file::spanfile].user() ) {
			std::string spanfile(option[OptionKeys::in::file::spanfile]());
			TR << "Reading spanfile " << spanfile << std::endl;
			topology.initialize(spanfile);
		} else {
			std::cerr << "spanfile missing ... " << std::endl;
		}
	}

	if ( !membrane_embed.calculated() ) {
		membrane_embed.initialize( pose );
		Vector init_normal,init_center;
		init_membrane_center_normal(pose,init_normal,init_center);
		Real score(0),best_score(999999),accepted_score(99999);
		Real temperature=2.0;
		Vector trial_normal(init_normal);
		Vector trial_center(init_center);
		Vector orig_trial_center(init_center);
		Vector orig_trial_normal(init_normal);
		Vector best_normal(init_normal);
		Vector best_center(init_center);
		Vector accepted_center(init_center);
		Vector accepted_normal(init_normal);
		score_normal_center(pose,trial_normal,trial_center,best_score);
		accepted_score=best_score;
		bool spanning(check_spanning(pose,trial_normal,trial_center));
		bool best_spanning(spanning);
		bool debug=option[ OptionKeys::membrane::debug ]();
		Real normal_mag=membrane_normal_magnitude_; //bw tmp
		Real center_mag=membrane_center_magnitude_; //bw tmp
		int max_delta_center=membrane_center_max_delta_; //vmyy default 5 A
		Size alpha_start=membrane_normal_start_angle_; //vmyy default 10 degrees
		Size delta_alpha=membrane_normal_delta_angle_; //vmyy default 10 degrees
		Size max_alpha=membrane_normal_max_angle_; //vmyy default 40 degrees
		Size nres=pose.total_residue();
		Size counter=0;
		Size accepted=0;
		Size thermally_accepted=0;

		if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
			for ( Size i = 1; i <= nres; ++i ) {
				Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
				membrane_embed.depth(i)=dot(xyz-best_center,best_normal)+30; //bw check that the values used are shifted 30.
			}
			membrane_embed.set_normal(best_normal);
			membrane_embed.set_center(best_center);
			membrane_embed.spanning()=true;
			membrane_embed.calculated()=true;
			return;
		}

		if ( memb_center_search_ || core::pose::symmetry::is_symmetric( pose ) ) {
			for ( int delta_center = -max_delta_center; delta_center <= max_delta_center; ++delta_center ) {

				if ( Membed_init_ ) break; //pba no mb embed optimization; just intial guess

				trial_center=orig_trial_center;
				search_memb_center (trial_center,trial_normal,delta_center);

				if ( !check_spanning(pose,trial_normal,trial_center) ) {
					if ( debug ) {
						TR << "delta_center= " << delta_center << " not spanning" << std::endl ;
					}
					trial_center=accepted_center;
					trial_normal=accepted_normal;
					continue;
				}

				score_normal_center(pose,trial_normal,trial_center,score);

				if ( score<accepted_score ) {
					if ( score<best_score ) {
						best_score=score;
						best_center=trial_center;
						best_normal=trial_normal;
						best_spanning=true;
					}
					accepted_score=score;
					accepted_center=trial_center;
					accepted_normal=trial_normal;
				}
			}
			//save best projection...
			for ( Size i = 1; i <= nres; ++i ) {
				Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
				membrane_embed.depth(i)=dot(xyz-best_center,best_normal)+30; //bw check that the values used are shifted 30.
			}
			membrane_embed.set_normal(best_normal);
			membrane_embed.set_center(best_center);
			membrane_embed.spanning()=best_spanning;
			membrane_embed.calculated()=true;
		}

		if ( memb_normal_search_ ) {
			for ( Size alpha=alpha_start; alpha<= max_alpha; alpha+=delta_alpha ) {
				for ( Size theta=0; theta<360; theta+=60 ) {
					if ( Membed_init_ ) break; //pba no mb embed optimization; just intial guess
					trial_normal=orig_trial_normal;
					search_memb_normal(trial_normal,alpha,theta);
					if ( !check_spanning(pose,trial_normal,trial_center) ) {
						if ( debug ) {
							TR << "alpha = " << alpha << " not spanning" << std::endl ;
						}
						trial_center=accepted_center;
						trial_normal=accepted_normal;
						continue;
					}

					score_normal_center(pose,trial_normal,trial_center,score);

					if ( score<accepted_score ) {
						if ( score<best_score ) {
							best_score=score;
							best_center=trial_center;
							best_normal=trial_normal;
							best_spanning=true;
						}
						accepted_score=score;
						accepted_center=trial_center;
						accepted_normal=trial_normal;
					}
				}
			}
			//save best projection...
			for ( Size i = 1; i <= nres; ++i ) {
				Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
				membrane_embed.depth(i)=dot(xyz-best_center,best_normal)+30; //bw check that the values used are shifted 30.
			}
			membrane_embed.set_normal(best_normal);
			membrane_embed.set_center(best_center);
			membrane_embed.spanning()=best_spanning;
			membrane_embed.calculated()=true;
		}

		if ( !memb_center_search_ && !memb_normal_search_ ) {
			for ( Size cycles=1; cycles<=membrane_normal_cycles_; ++cycles ) {
				if ( Membed_init_ ) break; //pba no mb embed optimization; just intial guess
				temperature = 2.0/cycles;
				if ( numeric::random::rg().uniform()<0.5 ) { // change center
					rigid_perturb_vector(trial_center,center_mag);
				} else { // change normal
					rot_perturb_vector(trial_normal,normal_mag);
				}
				if ( !check_spanning(pose,trial_normal,trial_center) ) {
					if ( debug ) {
						TR << "Cycle " << cycles << " not spanning" << std::endl ;
					}
					trial_center=accepted_center;
					trial_normal=accepted_normal;
					continue;
				}
				score_normal_center(pose,trial_normal,trial_center,score);
				if ( score<accepted_score ) {
					if ( score<best_score ) {
						best_score=score;
						best_center=trial_center;
						best_normal=trial_normal;
						best_spanning=true; //bw if you are here it is spanning....
					}
					accepted_score=score;
					accepted_center=trial_center;
					accepted_normal=trial_normal;
					++accepted;
				} else {
					++counter;
					Real const boltz_factor=(accepted_score-score)/temperature;
					Real const probability = std::exp( std::min ((core::Real)40.0, std::max((core::Real)-40.0,boltz_factor)) );
					if ( numeric::random::rg().uniform()<probability ) {
						accepted_score=score;
						accepted_center=trial_center;
						accepted_normal=trial_normal;
						++thermally_accepted;
						++accepted;
					} else {
						trial_center=accepted_center;
						trial_normal=accepted_normal;
					}
				}
			}
			//save best projection...
			for ( Size i = 1; i <= nres; ++i ) {
				Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
				membrane_embed.depth(i)=dot(xyz-best_center,best_normal)+30; //bw check that the values used are shifted 30.
			}
			membrane_embed.set_normal(best_normal);
			membrane_embed.set_center(best_center);
			membrane_embed.spanning()=best_spanning;
			membrane_embed.calculated()=true;
		}

		using namespace ObjexxFCL::format;
		TR << "MembraneCenter " << F(8,3,best_center.x())<< F(8,3,best_center.y())<< F(8,3,best_center.z()) << std::endl;
		TR << "MembraneNormal " << F(8,3,best_normal.x())<< F(8,3,best_normal.y())<< F(8,3,best_normal.z()) << std::endl;
		TR << "ATOM   9999  X   MEM A 999    " << F(8,3,best_center.x())<< F(8,3,best_center.y())<< F(8,3,best_center.z()) << std::endl;
		TR << "ATOM   9999  Y   MEM A 999    " << F(8,3,best_center.x()+15.*best_normal.x())<< F(8,3,best_center.y()+15.*best_normal.y())<< F(8,3,best_center.z()+15.*best_normal.z()) << std::endl;
		TR << "ATOM   9999  Z   MEM A 999    " << F(8,3,best_center.x()-15.*best_normal.x())<< F(8,3,best_center.y()-15.*best_normal.y())<< F(8,3,best_center.z()-15.*best_normal.z()) << std::endl;
	}
}

// duplicated in Embedding Factory and Embedding Residues
/// @brief Initialize Membrane Center and Normal parameters (from starting params)
void MembranePotential::init_membrane_center_normal(
	pose::Pose const & pose,
	Vector & normal,
	Vector & center
) const
{
	/// Load from

	if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
		center.zero();
		normal.zero();
		normal.z() = 1.0;

		if ( basic::options::option[basic::options::OptionKeys::membrane::membrane_center].user() ) {
			center.x() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[1];
			center.y() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[2];
			center.z() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[3];
		}
		if ( basic::options::option[basic::options::OptionKeys::membrane::membrane_normal].user() ) {
			normal.x() = basic::options::option[basic::options::OptionKeys::membrane::membrane_normal]()[1];
			normal.y() = basic::options::option[basic::options::OptionKeys::membrane::membrane_normal]()[2];
			normal.z() = basic::options::option[basic::options::OptionKeys::membrane::membrane_normal]()[3];
		}
		return;
	}

	MembraneTopology const & topology( MembraneTopology_from_pose(pose) );
	//Define vectors for inside and outside cap residue
	Vector inside(0);
	Vector outside(0);
	for ( Size i=1; i<=topology.tmhelix(); ++i ) {
		if ( !topology.allow_tmh_scoring(i) ) continue;
		Vector const & start( pose.residue( topology.span_begin(i) ).atom( 2 ).xyz());
		Vector const & end( pose.residue( topology.span_end(i) ).atom( 2 ).xyz());
		// all odd helices goes from outside in (from c++)
		if ( topology.helix_id(i) % 2 == 0 ) {
			inside+=start;
			outside+=end;
		} else {
			outside+=start;
			inside+=end;
		}
	}
	normal=outside-inside;
	normal.normalize();
	center=0.5*(outside+inside)/topology.tmh_inserted();
}

/////////////////////////////// Methods for Monte Carlo Search for Embeddings /////////////////////////////////////////////

// Helper function, will go into Embedding Factory
/// @brief Score Pose using MP Normal and Center
void MembranePotential::score_normal_center(
	pose::Pose const & pose,
	Vector const & normal,
	Vector const & center,
	Real & score
) const
{

	// Compute Starting Conditions
	Size const nres=pose.total_residue();
	MembraneTopology const & topology( MembraneTopology_from_pose(pose) ); // @ra should be SpanningTopology total

	// Initialize Scoring Parameters
	score = 0;
	Real residue_score(0);
	Real tm_projection(0);
	Real non_helix_pen(0);
	Real termini_pen(0);

	// For every residue in the pose, evaluate membrane environment based on center and normal given. Makes corrections
	// for symmetry
	for ( Size i = 1; i <= nres; ++i ) {

		Size rsdSeq(i);

		// Symmetry
		if ( core::pose::symmetry::is_symmetric( pose ) ) {

			using namespace core::conformation::symmetry;

			SymmetricConformation const & symm_conf ( dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
			SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
			if ( !symm_info->bb_is_independent(pose.residue(i).seqpos()) ) {
				rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
			}
			if ( symm_info->is_virtual(i) ) {
				rsdSeq = 0;
			}
		}

		// If any of these conditions apply, skip to the end of the loop
		if ( rsdSeq == 0 ) continue;
		if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
		if ( !topology.allow_scoring(rsdSeq) ) continue;

		// Grab CA coords, compute depth, score based on env and append score
		Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
		core::Real depth = dot( xyz-center, normal ) + 30;
		evaluate_env( pose, pose.residue(i), depth, residue_score );
		score+=residue_score;

	}

	// if the user specified to apply mp penalties, append the score
	if ( Menv_penalties_ ) {
		tm_projection_penalty( pose, normal, center, tm_projection );
		non_helix_in_membrane_penalty( pose, normal, center, non_helix_pen );
		termini_penalty( pose, normal, center, termini_pen );
		score+=tm_projection+non_helix_pen+termini_pen;
	}
}

// Helper function, will go into Embedding Factory
/// @brief Search for the Membrane Normal Parameter
void MembranePotential::search_memb_normal(
	Vector & n,
	Real const & alpha,
	Real const & theta
) const
{
	Real r_alpha = numeric::conversions::radians(alpha);
	Real r_theta = numeric::conversions::radians(theta);
	Vector u(std::sin(r_alpha) * std::cos(r_theta), std::sin(r_alpha) * std::sin(r_theta), std::cos(r_alpha));
	n=rotation_matrix_degrees(u,alpha)*n;
}

// Helper function, will go into Embedding Factory
/// @brief Search for the Membrane Center
void
MembranePotential::search_memb_center(
	Vector & c,
	Vector & n,
	Real const & delta
) const
{
	c = c + delta*n;
}

// Helper function, will go into Embedding Factory
/// @brief Rotatate Vector (Helper method - randomnly perturb for normal search)
void
MembranePotential::rot_perturb_vector(
	Vector & v,
	Real const & std_dev
) const
{
	Vector u( numeric::random::gaussian(), numeric::random::gaussian(), numeric::random::gaussian() ); //bw rotation_matrix will normalize.
	Real alpha( numeric::random::gaussian() * std_dev );
	v = rotation_matrix(u, alpha) * v;
}

void
MembranePotential::rigid_perturb_vector(Vector & v,
	Real const & std_dev) const
{
	Vector u(numeric::random::gaussian(),numeric::random::gaussian(),numeric::random::gaussian()); // there is a weird thing here???
	u.normalize();
	v=v+std_dev*u;
}

// Helper function, will go into Embedding Factory
/// @brief Should get replaced with the refactored check spanning in metrics
bool
MembranePotential::check_spanning(pose::Pose const & pose, Vector const & normal,Vector const & center) const
{
	MembraneTopology const & topology( MembraneTopology_from_pose(pose) );

	for ( Size i=1; i<=topology.tmhelix()-1; ++i ) {
		if ( !topology.allow_tmh_scoring(i) ) continue;
		Vector const & start_i( pose.residue( topology.span_begin(i) ).atom( 2 ).xyz());
		bool start_i_side=(dot(start_i-center,normal) > 0);
		bool span_check=false;
		for ( Size j=i+1; j<=topology.tmhelix(); ++j ) {
			if ( !topology.allow_tmh_scoring(j) || span_check ) continue;
			span_check=true;
			Vector const & start_j( pose.residue( topology.span_begin(j) ).atom( 2 ).xyz());
			bool start_j_side=(dot(start_j-center,normal) > 0);
			bool coord_para=(start_i_side==start_j_side);
			if ( topology.helix_id(i)-topology.helix_id(j) % 2 == 0 ) { // both should be on the same side (parallel)
				if ( !(coord_para) ) {
					return false;
				}
			} else {
				// should be on opposite sides.
				if ( coord_para ) {
					return false;
				}
			}
		}
	}
	return true;
}


////////////////////////// Membrane Penalties //////////////////////////////////////////////
void
MembranePotential::tm_projection_penalty(pose::Pose const & pose,Real & tm_proj) const
{
	tm_proj=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	Vector const normal(MembraneEmbed_from_pose( pose ).normal());
	Vector const center(MembraneEmbed_from_pose( pose ).center());
	tm_projection_penalty(pose,normal,center,tm_proj);
}

void
MembranePotential::tm_projection_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & tm_proj) const
{
	tm_proj=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	MembraneTopology const & topology( MembraneTopology_from_pose(pose) );

	//Define vectors for inside and outside cap residue
	Vector inside(0);
	Vector outside(0);
	tm_proj=0;

	for ( Size i=1; i<=topology.tmhelix(); ++i ) {
		if ( !topology.allow_tmh_scoring(i) ) continue;
		Vector const & start( pose.residue( topology.span_begin(i) ).atom( 2 ).xyz());
		Vector const & end( pose.residue( topology.span_end(i) ).atom( 2 ).xyz());
		Real tm_length=std::abs(dot(start-center,normal)-dot(end-center,normal));
		Real ratio=tm_length/(topology.span_end(i)-topology.span_begin(i)+1);
		if ( tm_length<15 ) {
			tm_proj++;
		}
		if ( ratio<1 || ratio > 1.5 ) {
			tm_proj++;
		}
	}
	tm_proj*=50; //total_embed weight is 0.5 in membrane_score_quick.cc
}

void
MembranePotential::non_helix_in_membrane_penalty(pose::Pose const & pose, Real & non_helix_pen) const
{
	non_helix_pen=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	Vector const normal(MembraneEmbed_from_pose( pose ).normal());
	Vector const center(MembraneEmbed_from_pose( pose ).center());
	non_helix_in_membrane_penalty(pose,normal,center,non_helix_pen);
}

void
MembranePotential::non_helix_in_membrane_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & non_helix_pen) const
{
	non_helix_pen=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	MembraneTopology const & topology( MembraneTopology_from_pose(pose) );
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		Size rsdSeq(i);
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			using namespace core::conformation::symmetry;
			SymmetricConformation const & symm_conf (
				dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
			SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
			if ( !symm_info->bb_is_independent(pose.residue(i).seqpos()) ) {
				rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
			}
			if ( symm_info->is_virtual(i) ) {
				rsdSeq = 0;
			}
		}
		if ( rsdSeq ==0 ) continue; // skip virtual residue

		if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
		if ( !topology.allow_scoring(rsdSeq) ) continue;
		if ( topology.tmregion(rsdSeq) && pose.conformation().secstruct(i)!='H' ) {
			Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
			Real depth=dot(xyz-center,normal)+30;
			if ( depth>18 &&
					depth<42 ) {
				non_helix_pen++;
			}
		}
	}
	non_helix_pen*=10; //total_embed weight is 0.5 in c++
}

void
MembranePotential::termini_penalty(pose::Pose const & pose, Real & termini_pen) const
{
	termini_pen=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	Vector const normal(MembraneEmbed_from_pose( pose ).normal());
	Vector const center(MembraneEmbed_from_pose( pose ).center());
	termini_penalty(pose,normal,center,termini_pen);
}

void
MembranePotential::termini_penalty(pose::Pose const & pose, Vector const & normal,Vector const & center,Real & termini_pen) const
{
	termini_pen=0.0;
	if ( !Menv_penalties_ ) {
		return;
	}
	MembraneTopology const & topology( MembraneTopology_from_pose(pose) );

	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_terminus() ) continue;

		Size rsdSeq(i);
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			using namespace core::conformation::symmetry;
			SymmetricConformation const & symm_conf (
				dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
			SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
			if ( !symm_info->bb_is_independent(pose.residue(i).seqpos()) ) {
				rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
			}
			if ( symm_info->is_virtual(i) ) {
				rsdSeq = 0;
			}
		}
		if ( rsdSeq ==0 ) continue;

		if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
		if ( topology.allow_scoring(rsdSeq) ) {
			Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
			Real depth=dot(xyz-center,normal);
			if ( depth>-12 &&
					depth<12 ) {
				termini_pen++;
			}
		}
	}
	termini_pen*=50;
}

/// @brief Add Const Membrane Embedding to the pose cache
/// @details Pose must already contain a cenlist object or this method will fail.
MembraneEmbed const &
MembraneEmbed_from_pose( pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) );
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
}

/// @brief Add a non const membrane embedding object to the pose cache
/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
MembraneEmbed & nonconst_MembraneEmbed_from_pose( pose::Pose & pose )
{

	if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ) {
		return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
	}

	MembraneEmbedOP membrane_embed( new MembraneEmbed );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED, membrane_embed );
	return *membrane_embed;
}

} // scoring
} // core

#endif // INCLUDED_core_scoring_MembranePotential_cc
