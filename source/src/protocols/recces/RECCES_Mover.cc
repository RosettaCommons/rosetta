// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/RECCES_Mover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/RECCES_Mover.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.hh>
#include <protocols/recces/util.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.RECCES_Mover" );

using namespace core;
using namespace utility;
using namespace utility::tools;
using namespace protocols::recces;
using namespace protocols::toolbox;

namespace protocols {
namespace recces {

//Constructor
RECCES_Mover::RECCES_Mover( utility::vector1< Real > const & temps,
	utility::vector1< Real > const & st_weights ):
	temps_( temps ),
	n_cycle_( 0 ),
	n_dump_( 0 ),
	dump_pdb_( false ),
	save_scores_( false ),
	a_form_range_( 0.0 )
{
	initialize( st_weights );
}

//Destructor
RECCES_Mover::~RECCES_Mover()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	if ( scorefxn_ == 0 ) scorefxn_ = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	initialize_sampler( pose );
	sampler_->apply( pose );

	run_sampler( pose );

	save_data_to_disk();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::run_sampler( pose::Pose & pose )
{
	using namespace core::pose;

	clock_t const time_start( clock() );

	// Simulated Tempering setup
	moves::SimulatedTempering tempering( pose, scorefxn_, temps_, weights_ );
	// tempering.set_rep_cutoff( 100 );
	set_sampler_gaussian_stdev( tempering.temperature(), pose );

	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1< float > scores;
	update_scores( scores, pose, scorefxn_ );
	utility::vector1< float > const null_arr_;
	data_ = utility::vector1<utility::vector1< float > >( temps_.size(), null_arr_ );

	// need to fix these magic numbers:
	Real const min( -100.05 ), max( 800.05 ), spacing( 0.1 );
	Histogram null_hist( min, max, spacing);
	hist_list_.clear();
	for ( Size n = 1; n <= temps_.size(); n++ ) hist_list_.push_back( null_hist );

	// Useful counters and variables during the loop
	Size n_accept_total( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	Size const n_t_jumps( n_cycle_ / t_jump_interval );
	Size temp_id( tempering.temp_id() );

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	std::cout << "Start the main sampling loop." << std::endl;
	if ( dump_pdb_ ) pose.dump_pdb( "init.pdb" );

	Size curr_dump = 1;

	// Main sampling cycle
	Size pct = 0;
	for ( Size n = 1; n <= n_cycle_; ++n ) {

		if ( n % ( n_cycle_/100 ) == 0 ) TR << ++pct << "% complete." << std::endl;

		// note that sampler holds information on all moving DOFs, and controls next conformation.
		++( *sampler_ );
		sampler_->apply( pose );

		// Check metropolis criterion
		if ( tempering.check_boltzmann( pose ) || n == n_cycle_ ) {
			// before accepting new pose, save history at this pose so far.
			if ( save_scores_ ) fill_data( data_[ temp_id ], curr_counts, scores );
			hist_list_[ temp_id ].add( scores[ 1 ], curr_counts );
			curr_counts = 1;

			// accept new pose.
			++n_accept_total;
			update_scores( scores, pose, scorefxn_ );
			if ( n == n_cycle_ ) break;

			// critical: tell sampler that new DOFs have been accepted.
			sampler_->update();
			update_dump_pdb( n, pose, scores, min_pose, min_score, curr_dump );
		} else {
			++curr_counts;
		}

		// Check temperature jump.
		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;

			// before jumping temperature, save history at this pose so far.
			if ( save_scores_ ) fill_data( data_[ temp_id ], curr_counts, scores );
			hist_list_[ temp_id ].add( scores[ 1 ], curr_counts );
			curr_counts = 1;

			set_sampler_gaussian_stdev( tempering.temperature(), pose);
			temp_id = tempering.temp_id();
		}
	}

	if ( dump_pdb_ ) {
		pose.dump_pdb( "end.pdb" );
		min_pose.dump_pdb( "min.pdb" );
		scorefxn_->show( min_pose );
	}

	TR << "n_cycles: " << n_cycle_ << std::endl;
	TR << "Accept rate: " << double( n_accept_total ) / n_cycle_ << std::endl;
	TR << "T_jump accept rate: " << double( n_t_jumps_accept ) / n_t_jumps << std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	TR << "Time in sampler: " <<  time_in_test << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::initialize( vector1< Real > const & orig_weights ) {
	runtime_assert( temps_.size() != 0 );
	if ( temps_.size() != orig_weights.size() ) weights_.push_back( 0 );
	weights_.append( orig_weights );
	runtime_assert( temps_.size() == weights_.size() );
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Sampler setup
void
RECCES_Mover::initialize_sampler( pose::Pose const & pose )
{
	using namespace protocols::recces::sampler::rna;

	std::string const & rna_secstruct_legacy = farna::secstruct::get_rna_secstruct_legacy_from_const_pose( pose );

	sampler_ = MC_RNA_MultiSuiteOP( new MC_RNA_MultiSuite );
	Size total_len( pose.size() );
	runtime_assert( total_len == rna_secstruct_legacy.size() );
	Size len1 = std::count( rna_secstruct_legacy.begin(), rna_secstruct_legacy.end(), 'H' ) / 2;
	for ( Size i = 1; i <= total_len; ++i ) {
		bool const sample_near_a_form( rna_secstruct_legacy[ i - 1 ] == 'H' );
		if ( i == 1 || ( i > len1 && i != total_len ) ) {
			MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i ) );
			suite_sampler->set_sample_bb( i != 1 );
			suite_sampler->set_sample_lower_nucleoside( true );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range_ );
			sampler_->add_external_loop_rotamer( suite_sampler );
		} else {
			MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i - 1 ) );
			suite_sampler->set_sample_bb( len1 == total_len || i != total_len );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( true );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range_ );
			sampler_->add_external_loop_rotamer( suite_sampler );
		}
	}

	if ( pose.size() == 2 /*  sample_rigid_body_ */ ) {
		MC_RNA_OneJumpOP jump_sampler( new MC_RNA_OneJump( pose, 1 /*jump*/ ) );
		jump_sampler->set_translation_mag( 0.01 ); //Real const translation_mag_( option[ translation_mag]() /* 0.01 */);
		jump_sampler->set_rotation_mag( 1.0 ); // Real const rotation_mag_( option[ rotation_mag ]() /* 1.0 degrees */ );
		sampler::MC_CombOP sampler_comb( sampler_ );
		sampler_comb->add_external_loop_rotamer( jump_sampler );
	}

	sampler_->init();
	sampler_->show( TR, 0 );
}

//////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::set_sampler_gaussian_stdev( Real const & temperature, pose::Pose const & pose )
{
	using namespace core::id;
	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	using namespace protocols::recces::sampler;

	std::string const & rna_secstruct_legacy = farna::secstruct::get_rna_secstruct_legacy_from_const_pose( pose );
	Size const n_rsd( pose.total_residue() );
	Real const bp_stdev(       gaussian_stdev( n_rsd, temperature, true ) );
	Real const dangling_stdev( gaussian_stdev( n_rsd, temperature, false ) );

	// update gaussian stdev for all chi's.
	for ( Size i = 1; i <= n_rsd; ++ i ) {
		Real const stdev = ( rna_secstruct_legacy[ i - 1 ] == 'H' ) ? bp_stdev : dangling_stdev;
		MC_SamplerOP torsion_sampler = sampler_->find( TorsionID( i, TorsionType::CHI, 1 ) );
		runtime_assert( torsion_sampler != 0 ); // these all move in RECCES
		std::dynamic_pointer_cast< MC_OneTorsion >(torsion_sampler)->set_gaussian_stdev( stdev );
	}

	// no need to update sugars -- they do not have gaussian ranges.

	// update gaussian stdev for all suites
	for ( Size i = 1; i < n_rsd; ++ i ) { // watch out: may need to get last residue if we cyclize
		if ( pose.fold_tree().is_cutpoint( i ) ) continue; // watch out: later generalize to cutpoint_closed
		Real const stdev = ( rna_secstruct_legacy[ i - 1 ] == 'H' &&
			rna_secstruct_legacy[ i ] == 'H' )      ? bp_stdev : dangling_stdev;
		vector1< TorsionID > suite_torsion_ids = get_suite_torsion_ids( i );
		for ( auto bb_torsion_id : suite_torsion_ids ) {
			MC_SamplerOP torsion_sampler = sampler_->find( bb_torsion_id );
			runtime_assert( torsion_sampler != 0 ); // these all move in RECCES
			std::dynamic_pointer_cast< MC_OneTorsion >(torsion_sampler)->set_gaussian_stdev( stdev );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::update_dump_pdb( Size const & n,
	pose::Pose const & pose,
	utility::vector1< Real > const & scores,
	pose::Pose & min_pose,
	Real & min_score,
	Size & curr_dump ) const
{
	if ( dump_pdb_ && scores[ 1 ] < min_score ) {
		min_score = scores[ 1 ];
		min_pose = pose;
	}
	if ( n_dump_ != 0 && n * (n_dump_ + 1) / double(n_cycle_) >= curr_dump ) {
		std::ostringstream oss;
		oss << "intermediate" << '_' << curr_dump << ".pdb";
		pose.dump_pdb(oss.str());
		++curr_dump;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::save_data_to_disk() const
{
	for ( Size i = 1; i <= temps_.size(); ++i ) {
		if ( save_scores_ ) {
			std::ostringstream oss;
			oss << out_prefix_ << '_' << std::fixed << std::setprecision( 2 )
				<< temps_[ i ] << ".bin.gz";
			Size const data_dim2( data_dim() );
			Size const data_dim1( data_[ i ].size() / data_dim2 );
			vector2disk_in2d( oss.str(), data_dim1, data_dim2, data_[ i ] );
		}
		std::ostringstream oss;
		oss << out_prefix_ << '_' << std::fixed << std::setprecision( 2 )
			<< temps_[ i ] << ".hist.gz";
		utility::vector1< Size > const & hist( hist_list_[ i ].get_hist() );
		utility::vector1< Real > const & scores( hist_list_[ i ].get_scores() );
		vector2disk_in1d( oss.str(), hist );

		std::ostringstream oss1;
		oss1 << out_prefix_ << "_hist_scores.gz";
		vector2disk_in1d( oss1.str(), scores );
	}
}

} //recces
} //protocols
