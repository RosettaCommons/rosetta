// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/data/RNA_DMS_LowResolutionPotential.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.hh>
#include <core/pose/rna/RNA_DataInfo.hh>
#include <core/scoring/rna/data/util.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/pose/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

using utility::tools::make_vector1;
using utility::vector1;
using Matrix = numeric::xyzMatrix<core::Real>;
using namespace basic::options;

//////////////////////////////////////////////////////
//
//              O
//              |
// DMS is CH3-O-S-O-CH3   (dimethyl sulfate)
//              |
//              O
//
// puts a methyl group at N1 of adenosine.
//
// This is a fast potential for scoring RNAs during
//  fragment-assembly -- the only feature is
//  whether an A has its Watson-Crick edge paired to
//  another base (or near an O2' -- potential pairing
//  partner)
//
// TO DO: Could move some of the base pair classification stuff to its
//  own class (RNA_BasePairInfo in pose/rna/ ??)
//
// -- rhiju, 2014
//
//////////////////////////////////////////////////////
static basic::Tracer TR( "core.scoring.rna.data.RNA_DMS_LowResolutionPotential" );

namespace core {
namespace scoring {
namespace rna {
namespace data {

//Constructor
RNA_DMS_LowResolutionPotential::RNA_DMS_LowResolutionPotential():
	careful_base_pair_classifier_( option[ OptionKeys::score::rna::DMS_careful_base_pair_classifier ]() )
{
}

//Destructor
RNA_DMS_LowResolutionPotential::~RNA_DMS_LowResolutionPotential() = default;

//////////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_LowResolutionPotential::initialize_DMS_low_resolution_potential() {

	//Read in data file, and fill in private data.
	auto WC_unprotected = read_DMS_low_resolution_stats_file( "scoring/rna/chem_map_lores/dms/ade_lores_not_WC_protected_logstats.txt" );
	auto WC_protected = read_DMS_low_resolution_stats_file( "scoring/rna/chem_map_lores/dms/ade_lores_WC_protected_logstats.txt" );
	is_protected_values_ = make_vector1( false, true );
	numeric::MathMatrix< Real > all_DMS_stats( 2, WC_unprotected.size() );
	all_DMS_stats.replace_col( 0 , WC_unprotected );
	all_DMS_stats.replace_col( 1, WC_protected );

	figure_out_low_resolution_potential( all_DMS_stats );
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_LowResolutionPotential::initialize( core::pose::Pose & pose, bool const rna_base_pair_computed /* = false */ ) {
#ifdef MULTI_THREADED
	utility_exit_with_message("The RNA_DMS_LowResolutionPotential is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.");
#endif
	get_rna_base_pairing_status( pose, wc_edge_paired_, hoog_edge_paired_, sugar_edge_paired_, is_bulged_, rna_base_pair_computed );
}

//////////////////////////////////////////////////////////////////////////////////
numeric::MathVector< core::Real >
RNA_DMS_LowResolutionPotential::read_DMS_low_resolution_stats_file( std::string const & potential_file ) {

	utility::io::izstream stream;
	basic::database::open( stream, potential_file );
	if ( !stream.good() ) utility_exit_with_message( "Unable to open "+potential_file );

	std::string line;

	// get dim
	getline( stream, line );
	std::istringstream l1( line );
	Real dim;
	l1 >> dim;

	// check labels
	getline( stream, line );
	std::istringstream l( line );
	utility::vector1< std::string > labels( 4, "" );
	l >> labels[1] >> labels[2] >> labels[3] >> labels[4];
	runtime_assert( labels[1] == "DMS" );
	runtime_assert( labels[2] == "log-stats" );

	Real  DMS, stats_value, log_stats_value;
	numeric::MathVector< Real > DMS_stats( 80 );
	Size DMS_idx = 0;
	while ( getline( stream, line ) ) {
		std::istringstream l( line );
		l  >> DMS >> log_stats_value;
		stats_value = exp( log_stats_value );
		DMS_values_.insert( DMS );
		DMS_stats( DMS_idx ) = stats_value;
		++DMS_idx;
	}

	return DMS_stats;
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_LowResolutionPotential::figure_out_low_resolution_potential( numeric::MathMatrix< Real > & DMS_stats ) {

	// first of all, need to normalize everything to total.
	Real DMS_stats_total( 0.0 );
	for ( Size h = 1; h <= is_protected_values_.size(); h++ ) {
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {
			DMS_stats_total += DMS_stats( h, k );
		}
	}
	for ( Size h = 1; h <= is_protected_values_.size(); h++ ) {
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {
			DMS_stats( h, k ) /= DMS_stats_total;
		}
	}

	// initialize projections to 0.0.
	numeric::MathVector< Real > p_DMS( DMS_values_.size() );
	numeric::MathVector< Real > p_model( is_protected_values_.size() );

	// fill projections, which give denominator of log-odds score.
	for ( Size h = 1; h <= is_protected_values_.size(); h++ ) {
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {

			p_DMS( k )   += DMS_stats( h, k );
			p_model( h ) += DMS_stats( h, k );

		}
	}

	DMS_low_resolution_potential_ = DMS_stats; // values will be replaced
	for ( Size h = 1; h <= is_protected_values_.size(); h++ ) {
		for ( Size k = 1; k <= DMS_values_.size(); k++ ) {
			DMS_low_resolution_potential_( h, k ) =
				-1.0 * log( DMS_stats( h, k ) / ( p_model( h ) * p_DMS( k ) ) );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
Real
RNA_DMS_LowResolutionPotential::evaluate( core::pose::Pose const & pose,
	pose::rna::RNA_Reactivity const & rna_reactivity ) {

#ifdef MULTI_THREADED
	utility_exit_with_message("The RNA_DMS_LowResolutionPotential is fundamentally non-threadsafe, and cannot be used in a multi-threaded context.");
#endif

	using namespace core::pose::full_model_info;

	if ( DMS_low_resolution_potential_.size() == 0 ) initialize_DMS_low_resolution_potential();

	runtime_assert( rna_reactivity.type() == pose::rna::DMS );
	Size const & pos = rna_reactivity.position();
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	if ( full_model_info.full_sequence()[ pos - 1 ] != 'a' ) return 0.0;
	if ( !full_model_info.working_res().has_value( pos ) ) return 0.0;

	// initialize feature values to what would be appropriate for a bulged A
	bool ade_wc_protected( false );
	Size i( 0 );
	if ( full_model_info.res_list().has_value( pos ) ) {
		i = full_model_info.full_to_sub( pos );
		runtime_assert ( pose.residue( i ).aa() == core::chemical::na_rad );
		// Skip virtual, repulsive
		if ( pose.residue( i ).is_virtual_residue() ) return 0.0;
		if ( pose.residue( i ).has_variant_type( chemical::REPLONLY ) ) return 0.0;

		// no syn-adenosines (some of these are exposed but not DMS-reactive -- chemical understanding
		// is currently incomplete.
		if ( pose.chi( i ) < 0.0 ) return 0.0;

		ade_wc_protected = wc_edge_paired_[ i ] || get_wc_near_o2prime( pose, i );
	}

	Size const protected_idx = get_bool_idx( ade_wc_protected, is_protected_values_ );
	Size const DMS_idx = get_idx( rna_reactivity.value(), DMS_values_ );

	Real score( 0.0 );
	score = DMS_low_resolution_potential_[ protected_idx ][ DMS_idx  ];

	//  TR << pos << " "  <<  pose.pdb_info()->number(i) << " ade_wc_protected " << ade_wc_protected << " (" << protected_idx << ")" << "  value " << rna_reactivity.value() << " (" << DMS_idx << ")" << " SCORE " << score << std::endl;

	return score;
}

///////////////////////////////////////////////////////////////////////////////
// following could move to protocols/, along with rna_base_pair classification code?
void
RNA_DMS_LowResolutionPotential::update_edge_paired(  Size const i, Size const k,
	vector1< bool > & wc_edge_paired,
	vector1< bool > & hoogsteen_edge_paired,
	vector1< bool > & sugar_edge_paired ) {

	if ( k == chemical::rna::WATSON_CRICK ) wc_edge_paired[i] = true;
	else if ( k == chemical::rna::HOOGSTEEN ) hoogsteen_edge_paired[i] = true;
	else {
		runtime_assert( k == chemical::rna::SUGAR );
		sugar_edge_paired[i] = true;
	}
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_DMS_LowResolutionPotential::get_rna_base_pairing_status( core::pose::Pose & pose,
	vector1< bool > & wc_edge_paired,
	vector1< bool > & hoogsteen_edge_paired,
	vector1< bool > & sugar_edge_paired,
	vector1< bool > & is_bulged,
	bool const already_scored ){

	using namespace core::scoring;

	//initialize
	wc_edge_paired        = vector1< bool >( pose.size(), false );
	hoogsteen_edge_paired = vector1< bool >( pose.size(), false );
	sugar_edge_paired     = vector1< bool >( pose.size(), false );
	is_bulged             = vector1< bool >( pose.size(), false );

	if ( careful_base_pair_classifier_ ) {
		ScoreFunctionOP scorefxn( new ScoreFunction );
		scorefxn->set_weight( rna_base_pair, 1.0 );
		scorefxn->set_weight( hbond_sc, 1.0 ); // used for HBondSet
		(*scorefxn)( pose );
		// scorefxn->show( std::cout, pose );

		// score pose with low res scorefunction -- will determine base pairs
		vector1< core::pose::rna::BasePair > base_pair_list;

		// following makes use of atomic-resolution hbond information -- might be more computationally efficient
		// to not compute that information -- and may contribute to smoother landscape at low-resolution stage.
		pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged ); // could also get 'energies' for each base pair...

		for ( pose::rna::BasePair const base_pair : base_pair_list ) {
			update_edge_paired(  base_pair.res1(), base_pair.edge1(), wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
			update_edge_paired(  base_pair.res2(), base_pair.edge2(), wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
		}
	} else {

		if ( !already_scored ) {
			ScoreFunctionOP scorefxn( new ScoreFunction );
			scorefxn->set_weight( rna_base_pair, 1.0 ); // should fill RNA_FilteredBaseBaseInfo
			(*scorefxn)( pose );
		}
		RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
		pose::rna::RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() ); // assume this has been filled.

		pose::rna::EnergyBasePairList const & scored_base_pair_list = rna_filtered_base_base_info.scored_base_pair_list();
		for ( auto const & elem : scored_base_pair_list ) {
			pose::rna::BasePair const & base_pair = elem.second;
			update_edge_paired( base_pair.res1(), base_pair.edge1(), wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
			update_edge_paired( base_pair.res2(), base_pair.edge2(), wc_edge_paired, hoogsteen_edge_paired, sugar_edge_paired );
		}

		pose::rna::EnergyBaseStackList const & scored_base_stack_list = rna_filtered_base_base_info.scored_base_stack_list();
		vector1< bool > is_stacked( pose.size(), false );
		for ( auto const & elem : scored_base_stack_list ) {
			pose::rna::BaseStack const & base_stack = elem.second;
			is_stacked[ base_stack.res1() ] = true;
			is_stacked[ base_stack.res2() ] = true;
		}

		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( wc_edge_paired[ i ] )        continue;
			if ( hoogsteen_edge_paired[ i ] ) continue;
			if ( sugar_edge_paired[ i ] )     continue;
			if ( is_stacked[ i ] )            continue;
			is_bulged[ i ] = true;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// copied code from RNA_LowResolutionPotential.cc -- eval_atom_derivative_rna_base_backbone
//   -- should *not* copy this code!
bool
RNA_DMS_LowResolutionPotential::get_wc_near_o2prime( core::pose::Pose const & pose, Size const i ) {
	using namespace core::scoring;
	using namespace core::scoring::rna;
	RNA_LowResolutionPotential const & rna_low_res_potential = ScoringManager::get_instance()->get_RNA_LowResolutionPotential();
	RNA_ScoringInfo  const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );
	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();
	for ( Size j = 1; j <= pose.size(); j++ ) {
		if ( !pose.residue_type( j ).is_RNA() ) continue;
		if ( i == j ) continue;
		Vector const & o2prime_xyz = pose.residue( j ).xyz( " O2'" );
		Vector const d_ij = o2prime_xyz - centroid_i;
		Real const dist_ij = d_ij.length();
		if ( dist_ij >= rna_low_res_potential.base_backbone_distance_cutoff() ) continue;
		Real const dist_x = dot_product( d_ij, x_i );
		Real const dist_y = dot_product( d_ij, y_i );
		Real const dist_z = dot_product( d_ij, z_i );
		Real const rho = std::sqrt( dist_x * dist_x + dist_y * dist_y );
		if ( std::abs( dist_z ) > rna_low_res_potential.base_backbone_z_cutoff() ) continue; // Look for atoms in the base plane
		if ( rho > rna_low_res_potential.base_backbone_rho_cutoff() ) continue; // Look for atoms in the base plane
		Real atom_cutoff_weight( 1.0 );
		if ( !rna_low_res_potential.check_for_base_neighbor( pose.residue(i), o2prime_xyz, atom_cutoff_weight ) ) continue;
		Real zeta_hoogsteen_cutoff( 60.0 ), zeta_sugar_cutoff( -60.0 );
		rna_low_res_potential.get_zeta_cutoff( pose.residue(i), zeta_hoogsteen_cutoff, zeta_sugar_cutoff );
		Real const zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x ) );
		if ( zeta < zeta_hoogsteen_cutoff && zeta > zeta_sugar_cutoff ) {
			//edge_bin = core::chemical::rna::WATSON_CRICK;  //Watson-Crick edge
			return true;
		}
	}
	return false;
}

} //data
} //rna
} //scoring
} //core
