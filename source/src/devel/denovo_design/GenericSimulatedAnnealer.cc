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
/// @detailed The score evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Tom Linsky (tlinsky@uw.edu)


// Unit Headers
#include <devel/denovo_design/GenericSimulatedAnnealer.hh>
#include <devel/denovo_design/GenericSimulatedAnnealerCreator.hh>

// Package Headers

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// C/C++ headers
#include <boost/foreach.hpp>

static basic::Tracer TR("devel.denovo_design.GenericSimulatedAnnealer");

namespace devel {
namespace denovo_design {

std::string
GenericSimulatedAnnealerCreator::keyname() const
{
  return GenericSimulatedAnnealerCreator::mover_name();
}

protocols::moves::MoverOP
GenericSimulatedAnnealerCreator::create_mover() const {
  return new GenericSimulatedAnnealer;
}

std::string
GenericSimulatedAnnealerCreator::mover_name()
{
  return "GenericSimulatedAnnealer";
}

/// @brief default constructor
GenericSimulatedAnnealer::GenericSimulatedAnnealer():
	GenericMonteCarloMover(),
	history_( 5 ),
	anneal_step_( 1 ),
	temp_step_( 1 )
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
  return new GenericSimulatedAnnealer( *this );
}

/// @brief create this type of object
protocols::moves::MoverOP
GenericSimulatedAnnealer::fresh_instance() const
{
  return new GenericSimulatedAnnealer();
}

void
GenericSimulatedAnnealer::scale_temperatures( core::Real const temp_factor )
{
	utility::vector1< core::Real > newtemps( start_temperatures_ );
	for ( core::Size i=1; i<=newtemps.size(); ++i ) {
		newtemps[i] *= temp_factor;
	}
	TR << "Factor=" << temp_factor << "; Setting temperatures to:";
	for( core::Size j=1; j<=newtemps.size(); ++j ) TR << " " << newtemps[j];
	TR << std::endl;
	temperatures( newtemps );

	// recalculate current_score and lowest_score, which may be different due to different temperatures if boltz_rank option is used
	if ( boltz_rank() ) {
		core::Real last_ranking_score( 0.0 );
		core::Real best_ranking_score( 0.0 );
		utility::vector1< core::Real > const & last_accepted( last_accepted_scores() );
		utility::vector1< core::Real > const & best_accepted( lowest_scores() );
		for ( core::Size index=1; index<=filters().size(); ++index ) {
			last_ranking_score += last_accepted[ index ] / newtemps[ index ];
			best_ranking_score += best_accepted[ index ] / newtemps[ index ];
		}
		last_accepted_score( last_ranking_score );
		lowest_score( best_ranking_score );
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

void
GenericSimulatedAnnealer::calculate_temps()
{
	core::Real const low_temp( 0.01 );
	core::Real const high_temp( 1.0 );
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
			// if any parameter's slope is not within error of zero within one-half SEM, don't cool yet
			if ( ( lines[i][1] - lines[i][3]/2 ) > 0 ) {
				cool = false;
				TR << "NOT COOLING." << std::endl;
			}
		}
		if ( cool ) {
			temp_factor = (high_temp - low_temp)*std::exp(-(core::Real)(anneal_step_)/filters().size()) + low_temp;
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

void
GenericSimulatedAnnealer::reset( Pose & pose )
{
	GenericMonteCarloMover::reset( pose );
	accepted_scores_.clear();
	start_temperatures_.clear();
	anneal_step_ = 1;
	accepted_scores_.push_back( last_accepted_scores() );
	// save initial temperatures
	start_temperatures_ = temperatures();
}


/// @Brief applies the mover
void
GenericSimulatedAnnealer::apply( Pose & pose )
{
  if( !mover() ){
    TR.Warning << "Mover is empty ! " << std::endl;
    return;
  }
  if( !filters().size() && !scorefxn() ){
    TR.Warning << "Both ScorefunctionOP and FilterOP are empty ! " << std::endl;
    return;
  }

	bool const stop_at_start( ( mover_stopping_condition() && mover_stopping_condition()->obj )
														|| stopping_condition()->apply( pose ) );
	if( stop_at_start ){
		TR<<"MC stopping condition met at the start, so failing without retrying "<<std::endl;
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		return;
	}

  //fpd
  if (mover()->get_additional_output())
      utility_exit_with_message("Movers returning multiple poses are unsupported by GenericSimulatedAnnealer.");

	core::Size trials = maxtrials();
	core::pack::task::TaskFactoryOP tf( task_factory() );
	if ( tf ) {
		core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( pose );
		trials = task_scaling() * num_designable( pose, task );
		task_factory( NULL );
	}

	if( saved_accept_file_name() != "" ){
		TR<<"Saving initial pose entering the MC trajectory, for use in checkpoint recovery. Checkpointing filename: "<<saved_accept_file_name()<<std::endl;
		pose.dump_pdb( saved_accept_file_name() );
	}

  PoseOP initial_pose = new Pose( pose );
	reset( pose ); //(re)initialize MC statistics
	core::Size accept( 0 );
	core::Size reject( 0 );

	for ( core::Size i=1; i<=trials; ++i ) {
		TR << "Starting trial " << i << std::endl;
		// apply the mover and test the pose
		TrialResult res( apply_mover( pose ) );
		if ( res == FINISHED ) {
			TR<<"MC stopping condition met at trial "<<i<<std::endl;
			break;
		} else if ( res == REJECTED ) {
			++reject;
		} else if ( res == ACCEPTED ) {
			++accept;
			calculate_temps();
		} else if ( res == FAILED ) {
			break;
		}
	}

	// Output final diagnositics, for potential tuning
	show_scores(TR);
	show_counters(TR);
	TR<<"Finished MC. Out of "<<trials<<" "<<accept<<" accepted "<<" and "<<reject<<" rejected."<<std::endl;

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

}// apply

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

	// check to see if we should accept or not
	if ( ! boltzmann( pose ) ) {
		pose = *last_accepted_pose();
		return REJECTED;
	} else {
		// save the results as accepted
		accepted_scores_.push_back( last_accepted_scores() );
		// check to see if the stopping condition is satisfied.
		bool const stop( ( mover_stopping_condition()() != NULL && mover_stopping_condition()->obj ) || stopping_condition()->apply( pose ) );
		if ( stop ) {
			return FINISHED;
		}
		return ACCEPTED;
	}
}

std::string
GenericSimulatedAnnealer::get_name() const {
  return GenericSimulatedAnnealerCreator::mover_name();
}

/// @brief parse xml file
void
GenericSimulatedAnnealer::parse_my_tag( TagPtr const tag,
																				 protocols::moves::DataMap & data,
																				 protocols::filters::Filters_map const & filters,
																				 protocols::moves::Movers_map const & movers,
																				 core::pose::Pose const & pose )
{
	GenericMonteCarloMover::parse_my_tag( tag, data, filters, movers, pose );
	history_ = tag->getOption< core::Size >( "history", history_ );
}


} // ns denovo_design
} // ns devel
