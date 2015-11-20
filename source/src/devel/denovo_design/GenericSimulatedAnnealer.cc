// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/GenericSimulatedAnnealer.cc
/// @brief perform a given mover and sample structures by MonteCarlo with gradual simulated annealing
/// @details The score evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Tom Linsky (tlinsky@uw.edu)


// Unit Headers
#include <devel/denovo_design/GenericSimulatedAnnealer.hh>
#include <devel/denovo_design/GenericSimulatedAnnealerCreator.hh>

// Package Headers

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// C/C++ headers
#include <fstream>
#include <boost/foreach.hpp>

static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.GenericSimulatedAnnealer" );

namespace devel {
namespace denovo_design {


std::string
GenericSimulatedAnnealerCreator::keyname() const
{
	return GenericSimulatedAnnealerCreator::mover_name();
}

protocols::moves::MoverOP
GenericSimulatedAnnealerCreator::create_mover() const {
	return protocols::moves::MoverOP( new GenericSimulatedAnnealer );
}

std::string
GenericSimulatedAnnealerCreator::mover_name()
{
	return "GenericSimulatedAnnealer";
}

/// @brief default constructor
GenericSimulatedAnnealer::GenericSimulatedAnnealer():
	GenericMonteCarloMover(),
	history_( 10 ),
	periodic_mover_( /* NULL */ ),
	eval_period_( 0 ),
	checkpoint_file_( "" ),
	keep_checkpoint_file_( false ),
	anneal_step_( 1 ),
	temp_step_( 1 ),
	current_trial_( 1 )
{
	accepted_scores_.clear();
	start_temperatures_.clear();
}

/// @brief destructor
GenericSimulatedAnnealer::~GenericSimulatedAnnealer(){}

/// @brief clone this object
protocols::moves::MoverOP
GenericSimulatedAnnealer::clone() const
{
	return protocols::moves::MoverOP( new GenericSimulatedAnnealer( *this ) );
}

/// @brief create this type of object
protocols::moves::MoverOP
GenericSimulatedAnnealer::fresh_instance() const
{
	return protocols::moves::MoverOP( new GenericSimulatedAnnealer() );
}

core::Real
GenericSimulatedAnnealer::calc_boltz_score( utility::vector1< core::Real > const & scores ) const
{
	runtime_assert( scores.size() == filters().size() );
	runtime_assert( scores.size() == temperatures().size() );
	utility::vector1< core::Real > const & temps( temperatures() );
	core::Real ranking_score( 0.0 );
	for ( core::Size index=1; index<=filters().size(); ++index ) {
		ranking_score += scores[ index ] / temps[ index ];
	}
	return ranking_score;
}

void
GenericSimulatedAnnealer::scale_temperatures( core::Real const temp_factor )
{
	utility::vector1< core::Real > newtemps( start_temperatures_ );
	for ( core::Size i=1; i<=newtemps.size(); ++i ) {
		newtemps[i] *= temp_factor;
	}
	TR << "Factor=" << temp_factor << "; Setting temperatures to:";
	for ( core::Size j=1; j<=newtemps.size(); ++j ) TR << " " << newtemps[j];
	TR << std::endl;
	temperatures( newtemps );

	// recalculate current_score and lowest_score, which may be different due to different temperatures if boltz_rank option is used
	if ( boltz_rank() ) {
		last_accepted_score( calc_boltz_score( last_accepted_scores() ) );
		lowest_score( calc_boltz_score( lowest_scores() ) );
	}

}

utility::vector1< core::Real >
GenericSimulatedAnnealer::calculate_standardized_scores( core::Size const filterid ) const
{
	runtime_assert( filterid <= filters().size() );
	runtime_assert( filterid != 0 );
	core::Real mean( 0.0 );
	core::Real stdev( 0.0 );

	core::Size const start( accepted_scores_.size() );
	core::Size const end( ( accepted_scores_.size() < history_-1 ) ? 1 : accepted_scores_.size()-history_+1 );
	core::Size const count( start - end + 1 );

	// find mean
	for ( core::Size i = start; i >= end; --i ) {
		mean += accepted_scores_[i][filterid];
	}
	mean /= count;

	// find std dev
	for ( core::Size i =start; i >= end; --i ) {
		stdev += std::pow( accepted_scores_[i][filterid] - mean, 2.0 );
	}
	stdev = std::sqrt( stdev / count );

	utility::vector1< core::Real > retval;
	for ( core::Size i = start; i >= end; --i ) {
		retval.push_back( ( accepted_scores_[i][filterid] - mean ) / stdev );
	}
	return retval;
}

// linear regression helper functions
core::Real calc_sum( utility::vector1< core::Real > const & vec );
core::Real calc_mean( utility::vector1< core::Real > const & vec );
core::Real calc_sum_squares( utility::vector1< core::Real > const & vec );
core::Real calc_sum_xy( utility::vector1< core::Real > const & x,
	utility::vector1< core::Real > const & y );
utility::vector1< core::Real >
linear_regression( utility::vector1< core::Real > const & x,
	utility::vector1< core::Real > const & y );

core::Real
calc_sum( utility::vector1< core::Real > const & vec )
{
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		sum += vec[i];
	}
	return sum;
}

core::Real
calc_mean( utility::vector1< core::Real > const & vec )
{
	return calc_sum( vec ) / vec.size();
}

core::Real
calc_sum_squares( utility::vector1< core::Real > const & vec )
{
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=vec.size(); ++i ) {
		sum += std::pow( vec[i], 2.0 );
	}
	return sum;
}

core::Real
calc_sum_xy( utility::vector1< core::Real > const & x,
	utility::vector1< core::Real > const & y )
{
	core::Real sum( 0.0 );
	runtime_assert( x.size() == y.size() );
	for ( core::Size i=1; i<=x.size(); ++i ) {
		sum += x[i] * y[i];
	}
	return sum;
}

/// @brief Given a set of x and y points, returns the linear regression
/// @brief the first element of the vector is slope
/// @brief the second is intercept
/// @brief the third is the sum of squares error
utility::vector1< core::Real >
linear_regression( utility::vector1< core::Real > const & x,
	utility::vector1< core::Real > const & y )
{
	runtime_assert( x.size() == y.size() );
	core::Size const n( x.size() );
	core::Real const sum_x( calc_sum( x ) );
	core::Real const sum_y( calc_sum( y ) );
	core::Real const sum_xy( calc_sum_xy( x, y ) );
	core::Real const sum_x2( calc_sum_squares( x ) );
	core::Real const sum_y2( calc_sum_squares( y ) );
	utility::vector1< core::Real > retval;
	// m = (n*(Exy) - Ex*Ey) / (n*(E(x^2)) - (Ex)^2)
	retval.push_back( ( n*sum_xy - sum_x*sum_y ) / ( n*sum_x2 - sum_x*sum_x ) );
	// b = ( Ey - m*Ex ) / n
	retval.push_back( ( sum_y - retval[1]*sum_x ) / n );
	// Standard error of estimate Syx
	// Syx = sqrt[ ( E(y2) - m*Exy - b*Ey ) / ( n - 2 ) ]
	core::Real const Syx( std::sqrt( ( sum_y2 - retval[1]*sum_xy - retval[2]*sum_y ) / ( n - 2 ) ) );
	// Standard error of slope Sb
	// Sb = Syx / sqrt( E(x^2) - n*Xbar^2 )
	retval.push_back( Syx / std::sqrt( sum_x2 - sum_x*sum_x / n ) );
	TR << "Regrission returned: ";
	for ( core::Size i=1; i<=retval.size(); ++i ) {
		TR << " " << retval[i];
	}
	TR << std::endl;
	return retval;
}

core::Real
GenericSimulatedAnnealer::calc_temp_factor() const
{
	core::Real const low_temp( 0.01 );
	core::Real const high_temp( 1.0 );
	return (high_temp - low_temp)*std::exp(-(core::Real)(anneal_step_)/filters().size()) + low_temp;
}

void
GenericSimulatedAnnealer::calculate_temps()
{
	// this uses the same temperature scheduler as the packer
	TR << "accepted_scores_=";
	for ( core::Size j=1; j<=accepted_scores_.size(); ++j ) {
		for ( core::Size k=1; k<=accepted_scores_[j].size(); ++k ) TR << " " << accepted_scores_[j][k];
		TR << std::endl;
	}
	TR << std::endl;
	core::Real temp_factor( -1.0 );
	if ( temp_step_ >= history_ ) {
		core::Size const start( accepted_scores_.size() );
		core::Size const end( ( accepted_scores_.size() < history_-1 ) ? 1 : accepted_scores_.size()-history_+1);

		utility::vector1< utility::vector1< core::Real > > lines;
		for ( core::Size i=1; i<=filters().size(); ++i ) {
			utility::vector1< core::Real > x,y;
			core::Size count( 1 );
			for ( core::Size j=start; j>=end; --j ) {
				x.push_back( count );
				++count;
				y.push_back( accepted_scores_[j][i] );
			}
			// now calculate slope
			lines.push_back( linear_regression( x, y ) );
		}
		// if the slope of each parameter is within error of zero, cool
		bool cool( true );
		for ( core::Size i=1; i<=filters().size(); ++i ) {
			TR << "slope=" << lines[i][1] << ", error=" << lines[i][3] << std::endl;
			// if any parameter's slope is not within error of zero, don't cool yet
			if ( ( lines[i][1] - lines[i][3] ) > 0 ) {
				cool = false;
				TR << "NOT COOLING." << std::endl;
			}
		}
		if ( cool ) {
			temp_factor = calc_temp_factor();
			TR << "Cool, factor=" << temp_factor << std::endl;
			temp_step_ = 1;
			++anneal_step_;
			scale_temperatures( temp_factor );
		}
	}
	if ( temp_factor < 0 ) {
		++temp_step_;
		TR << "Not scaling temperatures. Going to step " << temp_step_ << " at anneal_step " << anneal_step_ << std::endl;
	}
}

std::string
GenericSimulatedAnnealer::create_tag( std::string const & suffix ) const
{
	std::string cp = checkpoint_file_ + '_' + suffix;
	for ( std::string::iterator c=cp.begin(); c!=cp.end(); ++c ) {
		if ( ( *c == '/' )	 ||  ( *c == ' ' ) || ( *c == '\t' ) ) {
			*c = '_';
		}
	}
	return cp;
}

void
GenericSimulatedAnnealer::load_checkpoint_file( core::pose::Pose & pose )
{
	runtime_assert( checkpoint_file_ != "" );
	std::string const sf_name = checkpoint_file_ + ".out";
	core::io::silent::SilentFileData sfd( sf_name, false, false, "binary" );
	sfd.read_file( sf_name );
	std::string const best_tag = create_tag( "best" );
	if ( sfd.has_tag( best_tag ) ) {
		core::pose::PoseOP best_pose( new core::pose::Pose() );
		sfd.get_structure( best_tag ).fill_pose( *best_pose );
		lowest_score_pose( best_pose );
		TR << "Loaded lowest score pose " << best_tag << " from " << sf_name << std::endl;
	}
	std::string const last_accepted_tag = create_tag( "last" );
	if ( sfd.has_tag( last_accepted_tag ) ) {
		core::pose::PoseOP last_pose( new core::pose::Pose() );
		sfd.get_structure( last_accepted_tag ).fill_pose( *last_pose );
		last_accepted_pose( last_pose );
		TR << "Loaded last accepted pose " << last_accepted_tag << " from " << sf_name << std::endl;
	}
	runtime_assert( last_accepted_pose() );
	runtime_assert( lowest_score_pose() );
	// read the file with stats
	std::ifstream saved_file( checkpoint_file_.c_str() );
	if ( !saved_file.is_open() ) {
		utility_exit_with_message( "Error: The GenericSimulatedAnnealerMover could not open " + checkpoint_file_ + " for reading." );
	}
	// load trial number
	core::Size num_accepts( 0 );
	saved_file >> current_trial_ >> anneal_step_ >> temp_step_ >> num_accepts;
	TR << "current trial=" << current_trial_ << "; anneal_step=" << anneal_step_ << "; temp_step=" << temp_step_ << "; num_accepts=" << num_accepts << std::endl;
	// safety: clear accepted_scores_
	accepted_scores_.clear();
	for ( core::Size i=1; i<=num_accepts; ++i ) {
		utility::vector1< core::Real > filt_scores;
		for ( core::Size j=1; j<=filters().size(); ++j ) {
			core::Real score( 0.0 );
			saved_file >> score;
			filt_scores.push_back( score );
			TR << "score=" << score << " ";
		}
		TR << std::endl;
		accepted_scores_.push_back( filt_scores );
	}
	// load best scores
	utility::vector1< core::Real > best;
	for ( core::Size i=1; i<= filters().size(); ++i ) {
		core::Real score( 0.0 );
		saved_file >> score;
		best.push_back( score );
		TR << "best score : " << score << " ";
	}
	TR << std::endl;
	lowest_scores( best );
	// load last scores
	utility::vector1< core::Real > last;
	for ( core::Size i=1; i<= filters().size(); ++i ) {
		core::Real score( 0.0 );
		saved_file >> score;
		last.push_back( score );
		TR << "last score : " << score << " ";
	}
	TR << std::endl;
	last_accepted_scores( last );

	// calculate temperature factor and scale temperature
	core::Real const temp_factor( calc_temp_factor() );
	TR << "NEW temp factor is " << temp_factor << std::endl;
	// this will take care of changing up the last/best scores
	scale_temperatures( temp_factor );
	if ( saved_file.bad() ) {
		utility_exit_with_message( "The GenericSimulatedAnnealerMover encountered an error writing to " + checkpoint_file_ );
	}
	saved_file.close();
	pose = *last_accepted_pose();
}

void
GenericSimulatedAnnealer::save_checkpoint_file() const
{
	runtime_assert( checkpoint_file_ != "" );
	std::string const sf_name = checkpoint_file_ + ".out~";
	// save best pose
	if ( lowest_score_pose() && last_accepted_pose() ) {
		core::io::silent::SilentFileData sfd;
		core::io::silent::SilentStructOP best_score = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "binary" );
		best_score->fill_struct( *lowest_score_pose() );
		best_score->set_decoy_tag( create_tag( "best" ) );
		core::io::silent::SilentStructOP last_accepted = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "binary" );
		last_accepted->fill_struct( *last_accepted_pose() );
		last_accepted->set_decoy_tag( create_tag( "last" ) );
		sfd.add_structure_replace_tag_if_necessary( best_score );
		sfd.add_structure_replace_tag_if_necessary( last_accepted );
		sfd.write_all( sf_name );
	}
	// Write out a file with stats
	std::string const tmp_checkpoint( checkpoint_file_ + "~" );
	std::ofstream saved_file( std::string( tmp_checkpoint ).c_str() );
	if ( !saved_file.is_open() ) {
		utility_exit_with_message( "Error: The GenericSimulatedAnnealerMover could not open " + tmp_checkpoint + " for writing." );
	}
	// save trial number
	saved_file << boost::lexical_cast< std::string >( current_trial_ );
	// save annealing step
	saved_file << " " << boost::lexical_cast< std::string >( anneal_step_ );
	// save # acceptances at this annealing step
	saved_file << " " << boost::lexical_cast< std::string >( temp_step_ );
	// save total # of accepted scores
	saved_file << " " << boost::lexical_cast< std::string >( accepted_scores_.size() ) << std::endl;
	// save accepted score list
	for ( core::Size i=1; i<=accepted_scores_.size(); ++i ) {
		for ( core::Size j=1; j<=filters().size(); ++j ) {
			saved_file << accepted_scores_[i][j] << " ";
		}
		saved_file << std::endl;
	}
	// save best scores
	for ( core::Size i=1; i<=lowest_scores().size(); ++i ) {
		saved_file << boost::lexical_cast< std::string >( lowest_scores()[i] ) << " ";
	}
	saved_file << std::endl;
	// save last accepted scores
	for ( core::Size i=1; i<=last_accepted_scores().size(); ++i ) {
		saved_file << boost::lexical_cast< std::string >( last_accepted_scores()[i] ) << " ";
	}
	saved_file << std::endl;

	if ( saved_file.bad() ) {
		utility_exit_with_message( "The GenericSimulatedAnnealerMover encountered an error writing to " + tmp_checkpoint );
	} else {
		saved_file.close();
		// Writing was successful -- therefore, remove the tmp file and overwrite it
		replace_file( tmp_checkpoint, checkpoint_file_ );
		replace_file( sf_name, checkpoint_file_ + ".out" );
	}
}

/// @brief safely replaces a file with another
void replace_file( std::string const & origfile, std::string const & newfile )
{
	// see if newfile already exists
	std::ifstream f( newfile.c_str() );
	std::ifstream f2( origfile.c_str() );
	bool origfile_exists( false );
	bool newfile_exists( false );
	if ( f.good() ) {
		newfile_exists = true;
	}
	f.close();
	if ( f2.good() ) {
		origfile_exists = true;
	}
	f.close();

	if ( newfile_exists &&
			( rename( newfile.c_str(), std::string( newfile + "_tmp" ).c_str() ) != 0 ) ) {
		utility_exit_with_message( "Error renaming " + newfile + " to " + newfile + "_tmp" );
	}
	if ( ( ! origfile_exists ) ||
			( rename( std::string( origfile ).c_str(), newfile.c_str() ) != 0 ) ) {
		TR.Error << "Error moving " << origfile << " to " << newfile << std::endl;
		utility_exit();
	}
	if ( newfile_exists &&
			( remove( std::string( newfile + "_tmp" ).c_str() ) != 0 ) ) {
		utility_exit_with_message( "Error removing " + newfile + "_tmp" );
	}
}

void
GenericSimulatedAnnealer::remove_checkpoint_file() const
{
	runtime_assert( checkpoint_file_ != "" );
	utility::vector1< std::string > files;
	files.push_back( checkpoint_file_ + ".out" );
	files.push_back( checkpoint_file_ );
	for ( core::Size i=1; i<=files.size(); ++i ) {
		if ( remove( files[i].c_str() ) != 0 ) {
			TR.Error << "Error deleting " << files[i] << "!!" << std::endl;
		} else {
			TR.Info << "Deleted checkpoint file " << files[i] << std::endl;
		}
	}
}

bool
GenericSimulatedAnnealer::checkpoint_exists() const
{
	runtime_assert( checkpoint_file_ != "" );
	utility::vector1< std::string > files;
	files.push_back( checkpoint_file_ + ".out" );
	files.push_back( checkpoint_file_ );
	for ( core::Size i=1; i<=files.size(); ++i ) {
		std::ifstream f( files[i].c_str() );
		if ( f.good() ) {
			f.close();
		} else {
			f.close();
			return false;
		}
	}
	return true;
}

void
GenericSimulatedAnnealer::reset( Pose & pose )
{
	GenericMonteCarloMover::reset( pose );
	accepted_scores_.clear();
	anneal_step_ = 0;
	temp_step_ = 1;
	current_trial_ = 0;
	accepted_scores_.push_back( last_accepted_scores() );
	start_temperatures_ = temperatures();
}


/// @Brief applies the mover
void
GenericSimulatedAnnealer::apply( Pose & pose )
{
	if ( !mover() ) {
		TR.Warning << "Mover is empty ! " << std::endl;
		return;
	}
	if ( !filters().size() && !scorefxn() ) {
		TR.Warning << "Both ScorefunctionOP and FilterOP are empty ! " << std::endl;
		return;
	}

	bool const stop_at_start( ( mover_stopping_condition() && mover_stopping_condition()->obj )
		|| ( stopping_condition() && stopping_condition()->apply( pose ) ) );
	if ( stop_at_start ) {
		TR << "MC stopping condition met at the start, so failing without retrying " << std::endl;
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		return;
	}

	//fpd
	if ( mover()->get_additional_output() ) {
		utility_exit_with_message("Movers returning multiple poses are unsupported by GenericSimulatedAnnealer.");
	}

	core::Size trials = maxtrials();
	core::pack::task::TaskFactoryOP tf( task_factory() );
	if ( tf ) {
		core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
		trials = task_scaling() * num_designable( pose, task );
		task_factory( NULL );
	}

	//re-initialize MC statistics
	reset( pose );
	if ( checkpoint_file_ != "" && checkpoint_exists() ) {
		// if the checkpoint information exists, load it up to resume the run,
		TR.Info << "Resuming from checkpoint!" << std::endl;
		load_checkpoint_file( pose );
	} else {
		// otherwise, reset to the beginning
		TR.Info << "Checkpoint does not exist! Starting from the beginning" << std::endl;
		if ( checkpoint_file_ != "" ) {
			TR.Info << "Saving initial state entering the MC trajectory, for use in checkpoint recovery. Checkpointing filename: " << checkpoint_file_ << std::endl;
			save_checkpoint_file();
		}
	}

	core::Size accept_count( 0 );
	core::Size reject_count( 0 );
	++current_trial_;
	for ( ; current_trial_<=trials; ++current_trial_ ) {
		TR.Info << "Starting trial " << current_trial_ << std::endl;
		// apply the mover and test the pose
		TrialResult res( apply_mover( pose ) );
		if ( res == FINISHED ) {
			TR << "MC stopping condition met at trial " << current_trial_ << std::endl;
			break;
		} else if ( res == REJECTED ) {
			++reject_count;
		} else if ( res == ACCEPTED ) {
			++accept_count;
		} else if ( res == FAILED ) {
			break;
		}

		if ( checkpoint_file_ != "" ) {
			save_checkpoint_file();
		}

		// if we are ready for the periodic mover, run it
		if ( periodic_mover_ && eval_period_ && ( current_trial_ % eval_period_ == 0 ) ) {
			TR << "Running periodic mover..." << std::endl;
			periodic_mover_->apply( pose );
			// always accept periodic mover results
			accept( pose, score_pose( pose ), protocols::moves::MCA_accepted_thermally );
		}
	}

	// Output final diagnositics, for potential tuning
	show_scores(TR);
	show_counters(TR);
	TR << "Finished MC. Out of " << trials << ", " << accept_count << " accepted " << " and " << reject_count << " rejected." << std::endl;

	// Recover pose that have the lowest score, or the last accepted pose, as appropriate
	if ( recover_low() ) {
		recover_low( pose );
	} else {
		if ( last_accepted_pose() ) {
			pose = *last_accepted_pose();
		}
	}

	// reset the genericmontecarlo to use the task factory again next time
	if ( tf ) {
		task_factory( tf );
	}

	// clear out the checkpoint files
	if ( !keep_checkpoint_file_ && checkpoint_file_ != "" ) {
		remove_checkpoint_file();
	}

	// reset temperatures to starting values
	temperatures( start_temperatures_ );
}// apply


/// @brief given a pose, score the result
utility::vector1< core::Real >
GenericSimulatedAnnealer::score_pose( core::pose::Pose const & pose ) const
{
	utility::vector1< core::Real > provisional_scores;
	for ( core::Size index( 1 ); index <= filters().size(); ++index ) {
		TR.Debug <<"Filter #"<<index<<std::endl;
		protocols::filters::FilterCOP filter( filters()[ index ] );
		Real const flip( sample_types()[ index ] == "high" ? -1 : 1 );
		core::Real const filter_val( filter->report_sm( pose ) );
		TR<<"Filter "<<index<<" reports "<<filter_val<<" ( best="<<lowest_scores()[index]<<"; last="<<last_accepted_scores()[index]<<" )"<<std::endl;
		provisional_scores.push_back( flip * filter_val );
	}
	return provisional_scores;
}

/// @brief given a modified pose, determines whether we should accept or not, and updates internal class data accordingly
/// @brief uses randomly generated numbers to assess acceptance of scores with temperatures
TrialResult
GenericSimulatedAnnealer::boltzmann_result( core::pose::Pose & pose )
{
	utility::vector1< core::Real > randoms;
	for ( core::Size i=1; i<=filters().size(); ++i ) {
		randoms.push_back( numeric::random::rg().uniform() );
	}
	return boltzmann_result( pose, randoms );
}

/// @brief given a modified pose, determines whether we should accept or not, and updates internal class data accordingly
TrialResult
GenericSimulatedAnnealer::boltzmann_result( core::pose::Pose & pose,
	utility::vector1< core::Real > const & random_nums )
{
	TR << "Randoms=";
	for ( core::Size i=1; i<=random_nums.size(); ++i ) {
		TR << " " << random_nums[i];
	}
	TR << std::endl;
	// check to see if we should accept or not
	if ( ! boltzmann( pose, random_nums ) ) {
		TR << "Rejected" << std::endl;
		pose = *last_accepted_pose();
		return REJECTED;
	} else {
		// save the results as accepted
		accepted_scores_.push_back( last_accepted_scores() );
		// check to see if the stopping condition is satisfied.
		bool const stop( ( mover_stopping_condition() && mover_stopping_condition()->obj )
			|| ( stopping_condition() && stopping_condition()->apply( pose ) ) );
		if ( stop ) {
			return FINISHED;
		}
		calculate_temps();
		return ACCEPTED;
	}
}

TrialResult
GenericSimulatedAnnealer::apply_mover( core::pose::Pose & pose ) {
	// copy the pose so that we can restore it if a rejection occurs
	core::pose::Pose store_pose( pose );
	mover()->apply( pose );
	protocols::moves::MoverStatus ms( mover()->get_last_move_status() );
	if ( ms == protocols::moves::FAIL_RETRY ) {
		TR.Warning << "Mover failed. The mover, " << mover()->get_name() << ", will be performed again." << std::endl;
		pose = store_pose;
		return REJECTED;
	} else if ( ms == protocols::moves::FAIL_DO_NOT_RETRY || ms == protocols::moves::FAIL_BAD_INPUT ) {
		TR.Error << "The mover, " << mover()->get_name() << " failed. Exit from GenericSimulatedAnnealer." << std::endl;
		return FAILED;
	}
	return boltzmann_result( pose );

}

std::string
GenericSimulatedAnnealer::get_name() const {
	return GenericSimulatedAnnealerCreator::mover_name();
}

/// @brief parse xml file
void
GenericSimulatedAnnealer::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	GenericMonteCarloMover::parse_my_tag( tag, data, filters, movers, pose );
	history_ = tag->getOption< core::Size >( "history", history_ );
	eval_period_ = tag->getOption< core::Size >( "eval_period", eval_period_ );
	std::string const mover_name( tag->getOption< std::string >( "periodic_mover", "" ) );
	checkpoint_file_ = tag->getOption< std::string >( "checkpoint_file", checkpoint_file_ );
	keep_checkpoint_file_ = tag->getOption< bool >( "keep_checkpoint_file", keep_checkpoint_file_ );
	if ( mover_name != "" ) {
		Movers_map::const_iterator find_mover( movers.find( mover_name ) );
		if ( find_mover == movers.end() ) {
			TR.Error << "Error! Mover not found in map: " << mover_name << std::endl;
		}
		runtime_assert( find_mover != movers.end() );
		periodic_mover_ = find_mover->second;
	}
}


} // ns denovo_design
} // ns devel
