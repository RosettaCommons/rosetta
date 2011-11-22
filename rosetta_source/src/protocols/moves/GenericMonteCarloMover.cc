// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/GenericMonteCarloMover.cc
/// @brief perform a given mover and sample structures by MonteCarlo
/// @detailed The score evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


// Unit Headers
#include <protocols/moves/GenericMonteCarloMover.hh>
#include <protocols/moves/GenericMonteCarloMoverCreator.hh>

// C/C++ headers
#include <iostream>
// AUTO-REMOVED #include <iterator>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#define foreach BOOST_FOREACH

// Utility headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>

// Package Headers
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("protocols.moves.GenericMonteCarloMover");
static basic::Tracer TR_energies("protocols.moves.GenericMonteCarloMover.individual_energies");
static numeric::random::RandomGenerator mc_RG(61452); // <- Magic number, do not change it!!!

using namespace core;

namespace protocols {
namespace moves {

using namespace ObjexxFCL::fmt;

std::string
GenericMonteCarloMoverCreator::keyname() const
{
  return GenericMonteCarloMoverCreator::mover_name();
}

protocols::moves::MoverOP
GenericMonteCarloMoverCreator::create_mover() const {
  return new GenericMonteCarloMover;
}

std::string
GenericMonteCarloMoverCreator::mover_name()
{
  return "GenericMonteCarlo";
}

/// @brief default constructor
GenericMonteCarloMover::GenericMonteCarloMover():
  Mover("GenericMonteCarlo"),
  maxtrials_( 10 ),
  mover_( NULL ),
  scorefxn_( NULL ),
  temperature_( 0.0 ),
  sample_type_( "low" ),
  drift_( true ),
	preapply_( true ),
  recover_low_( true ),
	rank_by_filter_( 1 ),
	boltz_rank_( false ),
  last_accepted_pose_( NULL ),
  lowest_score_pose_( NULL ),
	stopping_condition_( NULL )
{
  initialize();
}

/// @brief value constructor
GenericMonteCarloMover::GenericMonteCarloMover(
  Size const maxtrials,
  MoverOP const & mover,
  Real const temperature,
  String const sample_type,
  bool const drift ) :
  Super("GenericMonteCarlo"),
  maxtrials_( maxtrials ),
  mover_( mover ),
  scorefxn_( NULL ),
  temperature_( temperature ),
  sample_type_( sample_type ),
  drift_( drift ),
	preapply_( true ),
  recover_low_( true ),
	rank_by_filter_(1),
	boltz_rank_( false ),
  last_accepted_pose_( NULL ),
  lowest_score_pose_( NULL )
{
  initialize();
}

/// @brief value constructor
GenericMonteCarloMover::GenericMonteCarloMover(
  Size const maxtrials,
  MoverOP const & mover,
  ScoreFunctionOP const & sfxn,
  Real const temperature,
  String const sample_type,
  bool const drift ) :
  Super("GenericMonteCarlo"),
  maxtrials_( maxtrials ),
  mover_( mover ),
  scorefxn_( sfxn ),
  temperature_( temperature ),
  sample_type_( sample_type ),
  drift_( drift ),
	preapply_( true ),
  recover_low_( true ),
	rank_by_filter_(1),
	boltz_rank_( false )
{
  initialize();
}

/// @brief destructor
GenericMonteCarloMover::~GenericMonteCarloMover(){}

/// @brief clone this object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::clone() const
{
  return new GenericMonteCarloMover( *this );
}

/// @brief create this type of object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::fresh_instance() const
{
  return new GenericMonteCarloMover();
}

/// @brief initialize
void
GenericMonteCarloMover::initialize()
{
  last_accepted_scores_.clear();
  if( sample_type_ == "high" ){
    flip_sign_ = -1;
  }else if( sample_type_ == "low" ){
    flip_sign_ = 1;
  }else{
    TR << "WARNING: the sample type, " << sample_type_ << ", is not defined." << std::endl;
    runtime_assert( false );
  }
  trial_counter_ = 0;
  accept_counter_ = 0;
  energy_gap_counter_ = 0.0;

  // trigger initialization
  next_trigger_id_ = 1;
}

/// @brief return the last accepted pose
GenericMonteCarloMover::PoseOP
GenericMonteCarloMover::last_accepted_pose() const
{
  return last_accepted_pose_;
}

/// @brief return the last accepted score
GenericMonteCarloMover::Real
GenericMonteCarloMover::last_accepted_score() const
{
  return last_accepted_score_;
}

/// @brief return the lowest score pose
GenericMonteCarloMover::PoseOP
GenericMonteCarloMover::lowest_score_pose() const
{
  return lowest_score_pose_;
}

/// @brief return the lowest score
GenericMonteCarloMover::Real
GenericMonteCarloMover::lowest_score() const
{
  return lowest_score_;
}

/// @brief return the lowest score
GenericMonteCarloMover::Real
GenericMonteCarloMover::current_score() const
{
  return current_score_;
}

/// @brief return mc_accepted
MCA
GenericMonteCarloMover::mc_accpeted() const
{
  return mc_accepted_;
}

/// @brief set max trials of monte carlo iterations
void
GenericMonteCarloMover::set_maxtrials( Size const ntrial )
{
  maxtrials_ = ntrial;
}

/// @brief set mover
void
GenericMonteCarloMover::set_mover( MoverOP const & mover )
{
  mover_ = mover;
}

/// @brief set filter
/// Pose is evaluated by FilterOP which can do report_sm() or ScoreFunctionOP during MC trials
/// You can choose either way FilterOP or ScoreFunction.
void
GenericMonteCarloMover::add_filter( FilterOP filter, bool const adaptive, Real const temp, String const sample_type, bool rank_by)
{
  filters_.push_back( filter );
	if(rank_by) {
		rank_by_filter_ = filters_.size();
	}
  adaptive_.push_back( adaptive );
  temperatures_.push_back( temp );
  sample_types_.push_back( sample_type );
  last_accepted_scores_.assign( filters_.size(), 100000 );
	num_rejections_.push_back(0);
  scorefxn_ = NULL;
}

/// @brief set scorefxn
/// Pose is evaluated by FilterOP which can do report_sm() or ScoreFunctionOP during MC trials
/// You can choose either way FilterOP or ScoreFunction.
void
GenericMonteCarloMover::set_scorefxn( ScoreFunctionOP const & sfxn )
{
  scorefxn_ = sfxn;
  filters_.clear();
}

/// @brief set temperatrue
void
GenericMonteCarloMover::set_temperature( Real const temp )
{
  temperature_ = temp;
}

/// @brief set sample type, high or low
/// when sample_type == high, sample pose which have higher value of scorey
/// when sample_type == low, sample pose which have lower value of score
void
GenericMonteCarloMover::set_sampletype( String const & type )
{
  if( sample_type_ != "high" && sample_type_ != "low" ){
    TR << "WARNING !! the sample type, " << type << ", is not defined." << std::endl;
    runtime_assert( false );
  }
  sample_type_ = type;
}

/// @brief if drift=false, the pose is set back to the initial pose
/// Of course, this is not MC sampling.
void
GenericMonteCarloMover::set_drift( bool const drift ){
  drift_ = drift;
}

/// @brief if preapply=true, auto-accept the first application of the submover,
/// ignoring boltzman criteria.
void
GenericMonteCarloMover::set_preapply( bool const preapply ) {
	preapply_ = preapply;
}

/// @brief if recover_low=true, after apply() the returned
/// is the lowest energy structure, rather than the last accepted structure.
void
GenericMonteCarloMover::set_recover_low( bool const recover_low ){
  recover_low_ = recover_low;
}

/// @brief if boltz_rank=true, rank structures by the temperature-weighted
/// sum of scores, rather than a single filter
void
GenericMonteCarloMover::set_boltz_rank( bool const boltz_rank ){
  boltz_rank_ = boltz_rank;
}

/// @brief show scores of last_accepted_score and lowest_score
void
GenericMonteCarloMover::show_scores( std::ostream & out ) const
{
  if( sample_type_ == "high" ){
    out << "Higher score sampled: trial=" << I( 5, trial_counter_ ) << ", score(current/last_accepted/best)=";
  }else{
    out << "Lower score sampled: trial=" << I( 5, trial_counter_ ) << ", score(current/last_accepted/best)=";
  }
  out << F( 9, 3, flip_sign_*current_score() ) << '/'
      << F( 9, 3, flip_sign_*last_accepted_score() )<< '/'
      << F( 9, 3, flip_sign_*lowest_score() ) << std::endl;
}

/// @brief show counters
void
GenericMonteCarloMover::show_counters( std::ostream & out ) const
{
  String const & mover( mover_->get_name() );
  String evaluation;
  if( scorefxn_ ){
    evaluation = "ScoreXXX"; // should we have name for ScoreFunction ?
  }
  int  const ntrials( trial_counter_ );
  int  const accepts( accept_counter_ );
  Real const energy_gap( energy_gap_counter_ );
  if( accepts > 0 ){
    out << "mover=" << LJ( 16, mover ) << " Score_eval=" << LJ( 16, evaluation ) <<
      " trials= " << I( 5, ntrials ) << "; " <<
      " accepts= " << F( 6, 3, Real( accepts )/ntrials ) << "; " <<
      " energy_gap/trial= " << F( 8, 3, Real( flip_sign_ * energy_gap ) / ntrials ) << std::endl;
  }else{
    out << "mover=" << A( 16, mover ) << " Score_eval=" << A( 16, evaluation )
        << " trials= " << I( 6, ntrials ) <<  " NO ACCEPTS." << std::endl;
  }
	if(num_rejections_.size()) {
		out << "Number of rejections per filter: ";
		for( core::Size ii(1); ii <= num_rejections_.size(); ++ii) {
			out << num_rejections_[ii] << "   ";
		}
		out << std::endl;
	}
}

/// @brief return the simulation state to the lowest energy structure we've seen
void
GenericMonteCarloMover::recover_low( Pose & pose )
{
  if( lowest_score_pose_ ) {
    pose = *lowest_score_pose_;
    *last_accepted_pose_ = *lowest_score_pose_;
  }else{
    // Case of that lowest_score_pose_ was never updated in MonteCarlo
    TR << " No lowest_score_pose_ exists, probably because all movers or filters failed. " << std::endl;
  }
}


/// @brief reset this GenericMonteCarloMover
void
GenericMonteCarloMover::reset( Pose & pose )
{
  if( filters_.size() == 0 ) {
    lowest_score_ = scoring( pose );
	} else {
    last_accepted_scores_.clear();
    for( Size index = 1; index <= filters_.size(); ++index ){
      protocols::filters::FilterCOP filter( filters_[ index ] );
      Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
      last_accepted_scores_.push_back( flip * filter->report_sm( pose ) );
    }
  }// fi filters_.size()

  lowest_score_pose_ = new Pose( pose );
  last_accepted_pose_ = new Pose( pose );
  last_accepted_score_ = lowest_score_;

  trial_counter_ = 0;
  accept_counter_ = 0;
  energy_gap_counter_ = 0;
	if ( num_rejections_.size() ) {
		num_rejections_.assign( num_rejections_.size(), 0 );
	}
  TR << "Initialization done " << std::endl;
}

/// @brief score pose based on filter or scorefxn
GenericMonteCarloMover::Real
GenericMonteCarloMover::scoring( Pose & pose )
{
  Real score( 0.0 );
  if( scorefxn_ ){
    score = flip_sign_ * (*scorefxn_)( pose );
  }
  return score;
}

// @brief
bool
GenericMonteCarloMover::boltzmann( Pose & pose )
{
  ++trial_counter_;
  TR.Debug <<"filters.size() "<<filters_.size()<<std::endl;
  if( filters_.size() ){
    runtime_assert( filters_.size() == adaptive_.size() && filters_.size() == temperatures_.size() && filters_.size() == sample_types_.size() &&  filters_.size() == num_rejections_.size());
    bool accept( false );
    utility::vector1< Real > provisional_scores;
    provisional_scores.clear();
    Real ranking_score( 0.0 );
    for( core::Size index( 1 ); index <= filters_.size(); ++index ){
      TR.Debug <<"Filter #"<<index<<std::endl;
      protocols::filters::FilterCOP filter( filters_[ index ] );
      bool const adaptive( adaptive_[ index ] );
      Real const temp( temperatures_[ index ] );
      Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
			core::Real const filter_val( filter->report_sm( pose ));
			TR<<"Filter "<<index<<" reports "<<filter_val<<std::endl;

      provisional_scores.push_back( flip * filter_val );
      if( index == rank_by_filter_ ) {
        ranking_score = provisional_scores[ rank_by_filter_ ];
			}
      Real const boltz_factor = ( last_accepted_scores_[ index ] - provisional_scores[ index ] ) / temp;
      TR_energies.Debug <<"energy index, last_accepted_score, current_score "<<index<<" "<< last_accepted_scores_[ index ]<<" "<<provisional_scores[ index ]<<std::endl;
      TR.Debug <<"Current, best, boltz "<<provisional_scores[ index ]<<" "<<last_accepted_scores_[ index ]<<" "<<boltz_factor<<std::endl;
      if( !adaptive ) { // return the starting score
        provisional_scores[ index ] = last_accepted_scores_[ index ];
			}
      Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
      Real const random_num( mc_RG.uniform() );
      bool const reject_filter( provisional_scores[ index ] > last_accepted_scores_[ index ] && random_num >= probability );
      if( reject_filter ){
        accept = false;
				++num_rejections_[index];
        break;
      }
      if( !reject_filter ) {
        accept = true;
			}
    }//for index
    if( accept ){
      TR.Debug <<"Accept"<<std::endl;
      mc_accepted_ = MCA_accepted_thermally;
      copy( provisional_scores.begin(), provisional_scores.end(), last_accepted_scores_.begin() );
      ++accept_counter_;
			if(boltz_rank_) {
				ranking_score = 0.0;
				for(core::Size ii(1); ii <= filters_.size(); ++ii ){
					ranking_score += provisional_scores[ii] / temperatures_[ii];
				}
			}
      energy_gap_counter_ += ranking_score - last_accepted_score();
      last_accepted_score_ = ranking_score;
      *last_accepted_pose_ = pose;
      if( ranking_score <= lowest_score() ){
        *lowest_score_pose_ = pose;
        lowest_score_ = ranking_score;
        mc_accepted_ = MCA_accepted_score_beat_low; //3;
      }
      return( true );
    }// fi accept
    else{
      TR.Debug <<"Reject"<<std::endl;
      mc_accepted_ = MCA_rejected;
      return( false );
    }
  }//fi filters_.size()
  else{
    Real score = scoring( pose );
    current_score_ = score; // for debugging
    show_scores( TR.Debug );
    if ( score > last_accepted_score() ) {
      if( temperature_ < 1e-8 ){
        mc_accepted_ = MCA_rejected; // rejected
      }else{
        Real const boltz_factor = ( last_accepted_score() - score ) / temperature_;
        Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
        if ( mc_RG.uniform() >= probability ) {
          mc_accepted_ = MCA_rejected; // rejected
        }else{
          mc_accepted_ = MCA_accepted_thermally; // accepted thermally
        }
      }
    }else{
      mc_accepted_ = MCA_accepted_score_beat_last; // accepted: energy is lower than last_accepted
    }

    if( mc_accepted_ >= 1 ){ // accepted
      ++accept_counter_;
      energy_gap_counter_ += score - last_accepted_score();
      *last_accepted_pose_ = pose;
      last_accepted_score_ = score;
      if ( score < lowest_score() ){  // energy is lower than last_accepted
        lowest_score_ = score;
        *lowest_score_pose_ = pose;
        mc_accepted_ = MCA_accepted_score_beat_low; //3;
      }
      return true;
    }else{ //rejected
      return false;
    }
  }
} // boltzmann

/// @Brief
void
GenericMonteCarloMover::apply( Pose & pose )
{
  using protocols::moves::FAIL_DO_NOT_RETRY;
  using protocols::moves::FAIL_BAD_INPUT;
  using protocols::moves::FAIL_RETRY;

  if( !mover_ ){
    TR.Warning << "Mover is empty ! " << std::endl;
    return;
  }
  if( !filters_.size() && !scorefxn_ ){
    TR.Warning << "Both ScorefunctionOP and FilterOP are empty ! " << std::endl;
    return;
  }

  //fpd
  if (mover_->get_additional_output())
      utility_exit_with_message("Movers returning multiple poses are unsupported by GenericMontoCarloMover.");

  PoseOP initial_pose = new Pose( pose );
	reset( pose ); //(re)initialize MC statistics
  MoverStatus ms( FAIL_RETRY );
	core::Size accept( 0 ), reject( 0 );
  for( Size i=1; i<=maxtrials_; i++ ){
    TR.Debug <<"Trial number: "<<i<<std::endl;
		bool const stop( stopping_condition()->apply( pose ) );
		if( stop ){
			TR<<"MC stopping condition met at trial "<<i<<std::endl;
			break;
		}
    Pose store_pose( pose );
    // Mover apply
    mover_->apply( pose );
    ms = mover_->get_last_move_status();
    if( ms == FAIL_RETRY ){
      TR.Warning << "Mover failed. The mover, " << mover_->get_name() << ", is performed again. " << std::endl;
      pose = store_pose;
      continue;
    }else if( ms == FAIL_DO_NOT_RETRY || ms == FAIL_BAD_INPUT ){
      TR.Error << "Mover failed. Exit from GenericMonteCarloMover." << std::endl;
      break;
    }
    // MonteCarlo
    if( preapply_ && i==1 ){ // Auto-accept first application in order to deal with movers that e.g. change the length of the pose.
      reset( pose );
    }else{
      if( ! boltzmann( pose ) ){ // evaluate pose by scorefxn_ or filter_.report_sm()
        pose = *last_accepted_pose();
				reject++;
      }
			else accept++;
    }
    if( !drift_ ){
      pose = (*initial_pose); // set back pose to initial one, of course this is not way of Monte Carlo
    }

    // Iterate over the list of triggers, executing each of them in turn
    fire_all_triggers(i, maxtrials_, pose, score_function());
  } // i<=maxtrials_

	// Output final diagnositics, for potential tuning
	show_scores(TR);
	show_counters(TR);
	TR<<"Finished MC. Out of "<<maxtrials_<<" "<<accept<<" accepted "<<" and "<<reject<<" rejected."<<std::endl;
  // Recover pose that have the lowest score, or the last accepted pose, as appropriate
	if(recover_low_) {
		recover_low( pose );
	} else {
		if (last_accepted_pose()) {
			pose = *last_accepted_pose();
		}
	}

}// apply

std::string
GenericMonteCarloMover::get_name() const {
  return GenericMonteCarloMoverCreator::mover_name();
}

/// @brief parse xml file
void
GenericMonteCarloMover::parse_my_tag( TagPtr const tag, DataMap & data, Filters_map const &filters, Movers_map const &movers, Pose const & )
{
	maxtrials_ = tag->getOption< core::Size >( "trials", 10 );
	temperature_ = tag->getOption< Real >( "temperature", 0.0 );

	String const  mover_name( tag->getOption< String >( "mover_name" ,""));
	String const filter_name( tag->getOption< String >( "filter_name", "true_filter" ) );
	Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
	Filters_map::const_iterator find_filter( filters.find( filter_name ));
	if( find_mover == movers.end() && mover_name != "" ) {
		TR.Error << "ERROR !! mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if( find_filter == filters.end() ) {
		TR.Error << "ERROR !! filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}
	if( mover_name != "" )
		mover_ = find_mover->second;
	bool const adaptive( tag->getOption< bool >( "adaptive", true ) );
	add_filter( find_filter->second->clone(), adaptive, temperature_, sample_type_ );
	String const sfxn ( tag->getOption< String >( "scorefxn_name", "" ) );
	if( sfxn != "" ){
		//scorefxn_ = new ScoreFunction( *data.get< ScoreFunction * >( "scorefxns", sfxn ));
		scorefxn_ = data.get< ScoreFunction * >( "scorefxns", sfxn )->clone();   //fpd use clone
		TR << "Score evaluation during MC is done by" << sfxn << ", ";
		TR << filter_name << " is ignored." << std::endl;
		filters_.clear();
	}else{
		scorefxn_ = NULL;
	}

	if( filter_name == "true_filter" && !scorefxn_ ){
		TR.Error << "You need to set filter_name or scorefxn_name for MC criteria." << std::endl;
		runtime_assert( false );
	}

	if( filters_.size() == 0 ){
		TR << "Apply mover of " << mover_name << ", and evaluate score by " << sfxn
		<< " at Temperature=" << temperature_ << ", ntrails= " << maxtrials_ << std::endl;
	}

	stopping_condition( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "stopping_condition", "false_filter" ), filters ) );
	TR<<"Generic MC using stopping condition "<< stopping_condition()->get_user_defined_name()<<std::endl;
	drift_ = tag->getOption< bool >( "drift", 1 );
	preapply_ = tag->getOption< bool >( "preapply", 1 ); // Default true for historical reasons
	recover_low_ = tag->getOption< bool >( "recover_low", 1 );
	boltz_rank_ = tag->getOption< bool >( "bolz_rank", 0 );
	sample_type_ = tag->getOption< String >( "sample_type", "low" );

	utility::vector1< TagPtr > const branch_tags( tag->getTags() );
	foreach( TagPtr const btag, branch_tags ){
		if( btag->getName() == "Filters" ){
			utility::vector1< TagPtr > const filters_tags( btag->getTags() );
			foreach( TagPtr const ftag, filters_tags ){
				String const filter_name( ftag->getOption< String >( "filter_name" ) );
				Filters_map::const_iterator find_filt( filters.find( filter_name ));
				if( find_filt == filters.end() ) {
					TR.Error << "Error !! filter not found in map: \n" << tag << std::endl;
					runtime_assert( find_filt != filters.end() );
				}
				Real const temp( ftag->getOption< Real >( "temperature", 1 ) );
				bool const adap( ftag->getOption< bool >( "adaptive", true ));
				String const samp_type( ftag->getOption< String >( "sample_type", "low" ));
				bool const rank( ftag->getOption< bool >( "rank", false ) );
				if( rank ) {
					if( rank_by_filter_ != 1 ) {
						TR.Warning << "WARNING Multiple filters set to rank! Using most recent (currently "<< filter_name << ")." << std::endl;
					}
					if( boltz_rank_ ) {
						TR.Warning << "WARNING Setting of rank on sub-filter "<< filter_name << " will be ignored as parent is set to Boltzmann rank all." << std::endl;
					}
				}
				add_filter( find_filt->second, adap, temp, samp_type, rank );
			} //foreach ftag
		}// fi Filters
		else
			utility_exit_with_message( "tag name " + btag->getName() + " unrecognized." );
	}//foreach btag

  if( !drift_ ){
    TR << "Pose is set back to initial pose every after applying mover and score evaluation." << std::endl;
  }

  if( sample_type_ == "high" ){
    TR << "Pose that have higher score is sampled." << sfxn << std::endl;
  }

  initialize();
}

void GenericMonteCarloMover::fire_all_triggers(
	Size cycle,
	Size num_cycles,
	const Pose& pose,
	ScoreFunctionOP scoring)
{
	boost::unordered_map<Size, GenericMonteCarloMoverTrigger>::iterator i;
	for (i = triggers_.begin(); i != triggers_.end(); ++i) {
		GenericMonteCarloMoverTrigger& t = i->second;
		bool rescore = t(cycle, num_cycles, pose, scoring);

		if (rescore) {
			last_accepted_score_ = scoring->score(*last_accepted_pose_);
			lowest_score_ = scoring->score(*lowest_score_pose_);
		}
	}
}

Size GenericMonteCarloMover::add_trigger(const GenericMonteCarloMoverTrigger& trigger) {
  Size tid = next_trigger_id_++;
  triggers_[tid] = trigger;
  return tid;
}

void GenericMonteCarloMover::remove_trigger(Size trigger_id) {
  boost::unordered_map<Size, GenericMonteCarloMoverTrigger>::iterator i = triggers_.find(trigger_id);
  if (i == triggers_.end()) {
    TR.Warning << "Attempt to remove invalid trigger_id => " << trigger_id << std::endl;
    return;
  }
  triggers_.erase(i);
}

Size GenericMonteCarloMover::num_triggers() const {
  return triggers_.size();
}

void
GenericMonteCarloMover::stopping_condition( protocols::filters::FilterOP f ){
	stopping_condition_ = f;
}

protocols::filters::FilterOP
GenericMonteCarloMover::stopping_condition() const{
	return stopping_condition_;
}

} // ns moves
} // ns protocols
