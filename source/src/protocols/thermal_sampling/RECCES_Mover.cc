// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/thermal_sampling/RECCES_Mover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/thermal_sampling/RECCES_Mover.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>
#include <protocols/thermal_sampling/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.thermal_sampling.RECCES_Mover" );

using namespace core;
using namespace utility;
using namespace protocols::stepwise;

namespace protocols {
namespace thermal_sampling {

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

	save_data();
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
	set_gaussian_stdev( tempering.temperature() );


	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1< float > scores;
	update_scores( scores, pose, scorefxn_ );
	utility::vector1< float > const null_arr_;
	data_ = utility::vector1<utility::vector1< float > >( temps_.size(), null_arr_ );


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
		if ( n % ( n_cycle_/100 ) == 0 ) {
			++pct;
			std::cout << pct << "% complete." << std::endl;
		}
		++( *sampler_ );
		sampler_->apply( pose );
		if ( tempering.boltzmann( pose ) || n == n_cycle_ ) {
			if ( save_scores_ ) fill_data( data_[ temp_id ], curr_counts, scores );
			++n_accept_total;
			hist_list_[ temp_id ].add( scores[ 1 ], curr_counts );
			update_scores( scores, pose, scorefxn_ );
			if ( n == n_cycle_ ) break;
			sampler_->update();
			curr_counts = 1;
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
		} else {
			++curr_counts;
		}

		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			if ( save_scores_ ) fill_data( data_[ temp_id ], curr_counts, scores );
			hist_list_[ temp_id ].add( scores[ 1 ], curr_counts );
			curr_counts = 1;
			set_gaussian_stdev( tempering.temperature() );
			temp_id = tempering.temp_id();
		}
	}
	if ( dump_pdb_ ) {
		pose.dump_pdb( "end.pdb" );
		min_pose.dump_pdb( "min.pdb" );
		scorefxn_->show( min_pose );
	}

	TR << "n_cycles: " << n_cycle_ << std::endl;
	TR << "Accept rate: " << double( n_accept_total ) / n_cycle_
		<< std::endl;
	TR << "T_jump accept rate: " << double( n_t_jumps_accept ) / n_t_jumps
		<< std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start )
		/ CLOCKS_PER_SEC;
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
	using namespace protocols::stepwise::sampler::rna;

	sampler_ = RNA_MC_MultiSuiteOP( new RNA_MC_MultiSuite );
	Size total_len( pose.size() );
	runtime_assert( total_len == ( bp_rsd_.size() + dangling_rsd_.size() ) );
	Size len1 = bp_rsd_.size() / 2;
	for ( Size i = 1; i <= total_len; ++i ) {
		bool const sample_near_a_form( bp_rsd_.has_value( i ) );
		if ( i == 1 || ( i > len1 && i != total_len ) ) {
			RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( i ) );
			suite_sampler->set_sample_bb( i != 1 );
			suite_sampler->set_sample_lower_nucleoside( true );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range_ );
			sampler_->add_external_loop_rotamer( suite_sampler );
		} else {
			RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( i - 1 ) );
			suite_sampler->set_sample_bb( len1 == total_len || i != total_len );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( true );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range_ );
			sampler_->add_external_loop_rotamer( suite_sampler );
		}
	}
	sampler_->init();
}

//////////////////////////////////////////////////////////////////////////////
// move this out? and do not make this Mover dependent on bp_rsd & dangling_rsd
void
RECCES_Mover::set_gaussian_stdev( Real const & temperature )
{
	Size const n_rsd( bp_rsd_.size() + dangling_rsd_.size() );
	Real const bp_stdev( gaussian_stdev( n_rsd, temperature, true ) );
	Real const dangling_stdev( gaussian_stdev( n_rsd, temperature, false ) );
	for ( Size i = 1; i <= bp_rsd_.size(); ++ i ) {
		sampler_->set_gaussian_stdev( bp_stdev, bp_rsd_[i] );
	}
	for ( Size i = 1; i <= dangling_rsd_.size(); ++ i ) {
		sampler_->set_gaussian_stdev( dangling_stdev, dangling_rsd_[i] );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::save_data() const
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

} //thermal_sampling
} //protocols
