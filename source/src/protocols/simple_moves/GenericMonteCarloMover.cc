// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/GenericMonteCarloMover.cc
/// @brief perform a given mover and sample structures by MonteCarlo
/// @details The score evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


// Unit Headers
#include <protocols/simple_moves/GenericMonteCarloMover.hh>
#include <protocols/simple_moves/GenericMonteCarloMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>

// C/C++ headers
#include <iostream>
#include <iterator>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <algorithm>

// Utility headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/filters/Filter.hh>

// Package Headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>
#include <core/pose/util.hh>
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/filters/BasicFilters.hh>
//////////////////////////////////////////////////
// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
///////////////////////////////////////////////////
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.GenericMonteCarloMover" );
static THREAD_LOCAL basic::Tracer TR_energies( "protocols.simple_moves.GenericMonteCarloMover.individual_energies" );

using namespace core;
using namespace protocols::moves;

namespace protocols {
namespace simple_moves {

using namespace ObjexxFCL::format;

/// @brief default constructor
GenericMonteCarloMover::GenericMonteCarloMover():
	Super( "GenericMonteCarlo" )
{
	initialize();
}


/// @brief value constructor without a score function

GenericMonteCarloMover::GenericMonteCarloMover(
	Size const maxtrials,
	Size const max_accepted_trials,
	Size const task_scaling,
	MoverOP const & mover,
	Real const temperature,
	String  sample_type,
	bool const drift ) :
	Super("GenericMonteCarlo"),
	maxtrials_( maxtrials ),
	max_accepted_trials_( max_accepted_trials ),
	task_scaling_( task_scaling ),
	mover_(std::move( mover )),
	temperature_( temperature ),
	sample_type_(std::move( sample_type )),
	drift_( drift )
{
	initialize();
}


/// @brief value constructor with a TaskFactory
GenericMonteCarloMover::GenericMonteCarloMover(
	Size const maxtrials,
	Size const max_accepted_trials,
	Size const task_scaling,
	MoverOP const & mover,
	TaskFactoryOP factory_in,
	Real const temperature,
	String  sample_type,
	bool const drift ) :
	Super("GenericMonteCarlo"),
	maxtrials_( maxtrials ),
	max_accepted_trials_( max_accepted_trials ),
	task_scaling_( task_scaling ),
	mover_(std::move( mover )),
	task_( /* NULL */ ),
	factory_ (std::move(factory_in)),
	temperature_( temperature ),
	sample_type_(std::move( sample_type )),
	drift_( drift )
{
	initialize();
}


/// @brief destructor
GenericMonteCarloMover::~GenericMonteCarloMover()= default;

/// @brief clone this object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::clone() const
{
	return GenericMonteCarloMover::MoverOP( new GenericMonteCarloMover( *this ) );
}

void GenericMonteCarloMover::task_factory( core::pack::task::TaskFactoryOP tf ) { factory_ = tf; }

/// @brief create this type of object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::fresh_instance() const
{
	return GenericMonteCarloMover::MoverOP( new GenericMonteCarloMover() );
}

/// @brief initialize
void
GenericMonteCarloMover::initialize()
{
	last_accepted_scores_.clear();
	lowest_scores_.clear();

	runtime_assert_string_msg( sample_type_ == "high" || sample_type_ == "low", "Error!  The sample type, " + sample_type_ + ", is not defined." );

	if ( sample_type_ == "high" ) {
		flip_sign_ = -1;
	} else if ( sample_type_ == "low" ) {
		flip_sign_ = 1;
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

/// @brief Set the lowest scoring pose
void
GenericMonteCarloMover::lowest_score_pose( core::pose::PoseOP pose )
{
	lowest_score_pose_ = pose;
}

/// @brief Set the last pose that was accepted
void
GenericMonteCarloMover::last_accepted_pose( core::pose::PoseOP pose )
{
	last_accepted_pose_ = pose;
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

/// @brief set the maximum number of trials that should be accepted
void
GenericMonteCarloMover::set_max_accepted_trials( Size const n_max_accepted_trials )
{
	max_accepted_trials_ = n_max_accepted_trials;
}

/// @brief set task multiplier to calculate trials from task
void
GenericMonteCarloMover::set_task_scaling( Size const scaling )
{
	task_scaling_ = scaling;
}

/// @brief set mover
void
GenericMonteCarloMover::set_mover( MoverOP mover )
{
	mover_ = mover;
}

/// @brief set filter
/// Pose is evaluated by FilterOP which can do report_sm() or ScoreFunctionOP during MC trials
/// You can choose either way FilterOP or ScoreFunction.
void
GenericMonteCarloMover::add_filter( FilterOP filter, bool const adaptive, Real const temp, String const & sample_type, bool rank_by)
{
	filters_.push_back( filter );
	if ( rank_by ) {
		rank_by_filter_ = filters_.size();
	}
	adaptive_.push_back( adaptive );
	temperatures_.push_back( temp );
	sample_types_.push_back( sample_type );
	last_accepted_scores_.assign( filters_.size(), 100000 );
	num_rejections_.push_back(0);
	scorefxn_ = nullptr;
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
	runtime_assert_string_msg( type == "high" || type == "low", "Error!  The sample type, " + type + ", is not defined." );
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
	if ( sample_type_ == "high" ) {
		out << "Higher score sampled: trial=" << I( 5, trial_counter_ ) << ", score(current/last_accepted/best)=";
	} else {
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
	if ( scorefxn_ ) {
		evaluation = "ScoreXXX"; // should we have name for ScoreFunction ?
	}
	int  const ntrials( trial_counter_ );
	int  const accepts( accept_counter_ );
	Real const energy_gap( energy_gap_counter_ );
	if ( accepts > 0 ) {
		out << "mover=" << LJ( 16, mover ) << " Score_eval=" << LJ( 16, evaluation ) <<
			" trials= " << I( 5, ntrials ) << "; " <<
			" accepts= " << F( 6, 3, Real( accepts )/ntrials ) << "; " <<
			" energy_gap/trial= " << F( 8, 3, Real( flip_sign_ * energy_gap ) / ntrials ) << std::endl;
	} else {
		out << "mover=" << A( 16, mover ) << " Score_eval=" << A( 16, evaluation )
			<< " trials= " << I( 6, ntrials ) <<  " NO ACCEPTS." << std::endl;
	}
	if ( num_rejections_.size() ) {
		out << "Number of rejections per filter: ";
		for ( core::Size ii(1); ii <= num_rejections_.size(); ++ii ) {
			out << num_rejections_[ii] << "   ";
		}
		out << std::endl;
	}
}

/// @brief return the simulation state to the lowest energy structure we've seen
void
GenericMonteCarloMover::recover_low( Pose & pose )
{
	if ( lowest_score_pose_ ) {
		pose = *lowest_score_pose_;
		*last_accepted_pose_ = *lowest_score_pose_;
	} else {
		// Case of that lowest_score_pose_ was never updated in MonteCarlo
		TR << " No lowest_score_pose_ exists, probably because all movers or filters failed. " << std::endl;
	}
}

/// @brief reset this GenericMonteCarloMover
void
GenericMonteCarloMover::reset( Pose & pose )
{
	last_accepted_scores_.clear();
	lowest_scores_.clear();
	if ( filters_.size() == 0 ) {
		if ( pose.size() > 0 ) {
			lowest_score_ = scoring( pose );
		} else {
			lowest_score_ = 999999999.99;
		}
		last_accepted_scores_.push_back( lowest_score_ );
		lowest_scores_.push_back( lowest_score_ );
	} else {
		core::Real ranking_score( 0.0 );
		for ( Size index = 1; index <= filters_.size(); ++index ) {
			Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
			Real score = flip * 99999999.99;
			if ( pose.size() > 0 ) {
				protocols::filters::FilterCOP filter( filters_[ index ] );
				score = flip * filter->report_sm( pose );
			}

			last_accepted_scores_.push_back( score );
			lowest_scores_.push_back( score );
			if ( index == rank_by_filter_ ) {
				ranking_score = last_accepted_scores_[ rank_by_filter_ ];
			}
		}
		if ( boltz_rank_ ) {
			ranking_score = 0.0;
			for ( Size index = 1; index <= filters_.size(); ++index ) {
				ranking_score += last_accepted_scores_[ index ] / temperatures_[ index ];
			}
		}
		lowest_score_ = ranking_score;
	}// fi filters_.size()

	lowest_score_pose_ = PoseOP( new Pose( pose ) );
	last_accepted_pose_ = PoseOP( new Pose( pose ) );
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
	if ( scorefxn_ ) {
		score = flip_sign_ * (*scorefxn_)( pose );
	}
	return score;
}
GenericMonteCarloMover::Size
GenericMonteCarloMover::num_designable( Pose & pose, PackerTaskOP & task )
{
	Size number_designable = 0;
	TR << "Designable Residues: ";
	for ( Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( task->design_residue(i) ) {
			++number_designable;
			TR << i << ", ";
		}
	}
	TR << std::endl << "Calculated number designable residues" << ": ";
	TR << number_designable << std::endl;
	return number_designable;
}

/// @brief
bool GenericMonteCarloMover::boltzmann( Pose & pose )
{
	return boltzmann( pose, generate_random() );
}

/// @brief Does what the mover needs to do when a pose is accepted, given a pose and scores
void
GenericMonteCarloMover::accept( Pose & pose,
	utility::vector1< core::Real > const & provisional_scores,
	MCA const mca_status )
{
	mc_accepted_ = mca_status;
	copy( provisional_scores.begin(), provisional_scores.end(), last_accepted_scores_.begin() );
	++accept_counter_;
	core::Real ranking_score( provisional_scores[rank_by_filter_] );
	if ( boltz_rank_ && filters_.size() > 0 ) {
		ranking_score = 0.0;
		for ( core::Size ii(1); ii <= filters_.size(); ++ii ) {
			ranking_score += provisional_scores[ii] / temperatures_[ii];
		}
	}
	energy_gap_counter_ += ranking_score - last_accepted_score();
	*last_accepted_pose_ = pose;
	last_accepted_score_ = ranking_score;
	if ( ranking_score < lowest_score() ) {
		*lowest_score_pose_ = pose;
		lowest_score_ = ranking_score;
		copy( provisional_scores.begin(), provisional_scores.end(), lowest_scores_.begin() );
		mc_accepted_ = MCA_accepted_score_beat_low;
		if ( saved_accept_file_name_ != "" ) {
			if ( mover_tag_ != nullptr ) {
				TR<<"Adding accepted mover tag to pose comments"<<std::endl;
				//      std::ofstream f;
				//      std::string const fname( saved_accept_file_name_ + ".mover_tag" );
				//      f.open( fname.c_str(), std::ios::out );
				//      if( !f.good() )
				//       utility_exit_with_message( "Unable to open MC mover_tag file " + fname );
				//      f<<mover_tag_->obj;
				//      f.close();
				core::pose::add_comment( pose, user_defined_mover_name_, mover_tag_->obj ); /// adding comment to the pose to save the mover's tag since it's accepted
			}
			TR<<"Dumping accepted file to disk as: "<<saved_accept_file_name_<<std::endl;
			if ( basic::options::option[basic::options::OptionKeys::out::file::pdb_comments]() ) { //if pose has comment field then dump it with comments
				//std::ofstream out( saved_accept_file_name_.c_str() );
				pose.dump_pdb(saved_accept_file_name_);
				/*              out << "##Begin comments##" << std::endl;
				using namespace std;
				map< string, string > const comments = core::pose::get_all_comments(pose);
				for( auto const & comment : comments ) {
				out << comment.first<<" "<<comment.second << std::endl;
				}
				out << "##End comments##" << std::endl;
				*/
			} else {
				pose.dump_pdb( saved_accept_file_name_ );
			}
		}
	}
}

/// @brief Uses filters or scorefunction to evaluate a pose and either accept or reject it
/// @details Uses the filters or scorefunction to evaluates a pose and compare it to the
/// last accepted and best poses.
/// @param pose = the pose to be evaluated
/// @param random_nums = a vector of random numbers between 0 and 1, to be used for accepting poses that are not necessarily improved on filter/scorefxn value over the last accepted pose. Values of 0 indicate that only values that are improvements will be accepted, and values of 1 indicate that everything will be accepted regardless of scorefxn/filter value. Must have the same size as filters_

bool
GenericMonteCarloMover::boltzmann( Pose & pose, utility::vector1< core::Real > const & random_nums )
{
	++trial_counter_;
	TR.Debug <<"filters.size() "<<filters_.size()<<std::endl;
	core::Real filter_val(0.0);
	if ( filters_.size() ) {
		runtime_assert( filters_.size() == adaptive_.size() );
		runtime_assert( filters_.size() == temperatures_.size() );
		runtime_assert( filters_.size() == sample_types_.size() );
		runtime_assert( filters_.size() == num_rejections_.size() );
		runtime_assert( filters_.size() == random_nums.size() );
		bool accepted( false );
		utility::vector1< Real > provisional_scores;
		provisional_scores.clear();
		for ( core::Size index( 1 ); index <= filters_.size(); ++index ) {
			runtime_assert( random_nums[ index ] >= 0.0 );
			runtime_assert( random_nums[ index ] <= 1.0 );
			TR.Debug <<"Filter #"<<index<<std::endl;
			protocols::filters::FilterCOP filter( filters_[ index ] );
			bool const adaptive( adaptive_[ index ] );
			Real const temp( temperatures_[ index ] );
			Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
			filter_val = filter->report_sm( pose );
			TR<<"Filter "<<index<<" reports "<<filter_val<<" ( best="<<lowest_scores_[index]<<"; last="<<last_accepted_scores_[index]<<" )"<<std::endl;

			provisional_scores.push_back( flip * filter_val );
			Real const boltz_factor = ( last_accepted_scores_[ index ] - provisional_scores[ index ] ) / temp;
			TR_energies.Debug <<"energy index, last_accepted_score, current_score "<<index<<" "<< last_accepted_scores_[ index ]<<" "<<provisional_scores[ index ]<<std::endl;
			TR.Debug <<"Current, best, boltz "<<provisional_scores[ index ]<<" "<<last_accepted_scores_[ index ]<<" "<<boltz_factor<<std::endl;
			if ( !adaptive ) { // return the starting score
				provisional_scores[ index ] = last_accepted_scores_[ index ];
			}
			Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
			bool const reject_filter( provisional_scores[ index ] > last_accepted_scores_[ index ]
				&& random_nums[ index ] >= probability );
			if ( reject_filter ) {
				accepted = false;
				++num_rejections_[index];
				break;
			}
			if ( !reject_filter ) {
				accepted = true;
			}
		}//for index

		if ( progress_file_ != "" ) { //write progress data to file
			/// write a table to a progress file that has the following structure
			/// Trial# filter_val pose_comments protein_sequence

			std::ofstream data;
			data.open( progress_file().c_str(), std::ios::app );
			if ( !data.good() ) {
				utility_exit_with_message( "Unable to open GenericMonteCarlo progress file for writing: " + progress_file() + "\n" );
			}

			std::string pose_sequence( "" );
			for ( core::Size chaini = 1 ; chaini <= pose.conformation().num_chains(); ++chaini ) {
				pose_sequence += pose.chain_sequence( chaini );
			}

			using namespace std;
			string stringed_comments("");
			map< string, string > const comments = core::pose::get_all_comments(pose);
			for ( auto const & comment : comments ) {
				stringed_comments += comment.first + ":" + comment.second + " ";
			}
			data<<trial_counter_<<" "<<accepted<<" "<<filter_val<<" "<<stringed_comments<<" "<<pose_sequence<<'\n';
			data.flush();
		}

		if ( accepted ) {
			TR<<"Accept"<<std::endl;
			accept( pose, provisional_scores, MCA_accepted_thermally );
			return true;
		} else { // fi accepted
			TR<<"Reject"<<std::endl;
			mc_accepted_ = MCA_rejected;
			return false;
		}
	} else { //fi filters_.size()
		runtime_assert( random_nums.size() >= 1 );
		MCA mc_status( MCA_rejected );
		utility::vector1< core::Real > provisional_score;
		provisional_score.clear();
		provisional_score.push_back( scoring( pose ) );
		current_score_ = provisional_score[1]; // for debugging
		last_tested_scores_.clear();
		last_tested_scores_ = provisional_score;
		show_scores( TR.Debug );
		if ( provisional_score[1] > last_accepted_score() ) {
			if ( temperature_ >= 1e-8 ) {
				Real const boltz_factor = ( last_accepted_score() - provisional_score[1] ) / temperature_;
				Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
				if ( random_nums[1] < probability ) {
					//if ( numeric::random::rg().uniform() < probability ) {
					mc_status = MCA_accepted_thermally; // accepted thermally
				}
			}
		} else {
			mc_status = MCA_accepted_score_beat_last; // accepted: energy is lower than last_accepted
		}

		if ( mc_status >= 1 ) { // accepted
			accept( pose, provisional_score, mc_status );
			return true;
		} else {
			// rejected
			mc_accepted_ = mc_status;
			return false;
		}
	}
} // boltzmann

core::Size
GenericMonteCarloMover::load_trial_number_from_checkpoint( core::pose::Pose & pose ){
	using namespace std;
	core::Size trial( 1 );
	if ( saved_trial_number_file_ != "" ) {
		ifstream f( saved_trial_number_file_.c_str(), ios::in );
		if ( f.good() ) {
			core::Size const begin = f.tellg();
			f.seekg( 0, ios::end );
			core::Size const end = f.tellg();
			if ( end - begin != 0 ) { //file size != 0
				f.seekg( 0, ios::beg );// return to the beginning

				TR<<"Loading trial number from checkpoint"<<std::endl;
				std::string line;
				getline( f, line );
				std::istringstream line_stream( line );
				line_stream >> trial;
				TR<<"Loaded trial number: "<<trial<<std::endl;
			}//fi end-begin
			f.close();
		}//fi f.good()
	}// fi saved_trial_number_file_ != ""
	if ( mover_tag_ != nullptr ) {
		std::string const fname( saved_trial_number_file_ + ".mover_tag" );
		ifstream f_mover_tag( fname.c_str(), ios::in );
		if ( f_mover_tag.good() ) {
			f_mover_tag >> mover_tag_->obj;
			TR<<"Loaded mover_tag from checkpointing file: "<<mover_tag_->obj<<std::endl;
		} else {
			TR<<"File containing mover_tag "<<fname<<" not found. Not loading movertag from checkpoint"<<std::endl;
		}
	}
	if ( reset_baselines() ) {
		/// see if any subfilters need to be reset
		bool call_reset( false );
		using namespace protocols::filters;
		using namespace protocols::simple_filters;
		for ( FilterOP filter : filters_ ) {
			if ( filter->get_type() == "Operator" ) {
				TR<<"Resetting Operator filter's baseline"<<std::endl;
				OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( filter ) );
				operator_filter->reset_baseline( pose, trial != 1/*if trial>1, attempt to read the baselines from checkpointing files. Otherwise, don't use the checkpointing files*/ );
				call_reset = true;
			} else if ( filter->get_type() == "CompoundStatement" ) { // fi Operator /// User defined filters with confidence!=1 in RosettaScripts are all CompoundFilter, so poke inside...
				CompoundFilterOP comp_filt_op( utility::pointer::dynamic_pointer_cast< protocols::filters::CompoundFilter > ( filter ) );
				runtime_assert( comp_filt_op != nullptr );
				for ( auto & cs_it : *comp_filt_op ) {
					FilterOP filt( cs_it.first );
					if ( filt->get_type() == "Operator" ) {
						TR<<"Resetting Operator filter's baseline"<<std::endl;
						OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( filt ) );
						operator_filter->reset_baseline( pose, trial != 1/*if trial>1, attempt to read the baselines from checkpointing files. Otherwise, don't use the checkpointing files*/ );
						call_reset = true;
					}// fi Operator
				}// for cs_it
			}//elseif CompoundStatement
		} //foreach
		if ( call_reset ) {
			TR<<"Resetting Boltzmann's accepted values, in case filters have changed during loading"<<std::endl;
			reset( pose );// this is called to reset the boltzmann scores with the new baselines
		}
	}/// fi reset_baselines
	return trial;
}

void
GenericMonteCarloMover::save_trial_number_to_checkpoint( core::Size const i ) const{
	if ( saved_trial_number_file_ == "" ) {
		return;
	}
	std::ofstream f;
	f.open( saved_trial_number_file_.c_str(), std::ios::out );
	if ( !f.good() ) {
		utility_exit_with_message( "Unable to open MC checkpointing file " + saved_trial_number_file_ );
	}

	f<<i;
	f.close();
}

/// @Brief
///comment
void
GenericMonteCarloMover::apply( Pose & pose )
{
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_BAD_INPUT;
	using protocols::moves::FAIL_RETRY;

	if ( !mover_ ) {
		TR.Warning << "Mover is empty ! " << std::endl;
		return;
	}
	if ( !filters_.size() && !scorefxn_ ) {
		TR.Warning << "Both ScorefunctionOP and FilterOP are empty ! " << std::endl;
		return;
	}

	core::pack::task::PackerTaskOP task;

	if ( factory_ ) {
		task = factory_->create_task_and_apply_taskoperations( pose );
	} else {
		TR << "No task inputted" << std::endl; //LGN
	}

	if ( factory_ ) {
		number_designable_ = num_designable (pose, task);
		TR << "Input number of trials is " << maxtrials_ << std::endl;
		TR << "The task_scaling is: " << task_scaling_ << std::endl;
		maxtrials_ = task_scaling_ * number_designable_;
		TR << "Resetting number of trials based on your input task to: " << maxtrials_ << std::endl;
	}
	TR << "The number of trials for this run is: " << maxtrials_ << std::endl;

	bool const stop_at_start( ( mover_stopping_condition_ && mover_stopping_condition_->obj ) ||
		( stopping_condition() && stopping_condition()->apply( pose ) ) );
	if ( stop_at_start ) {
		TR<<"MC stopping condition met at the start, so failing without retrying "<<std::endl;
		set_last_move_status( FAIL_DO_NOT_RETRY );
		return;
	}

	//fpd
	if ( mover_->get_additional_output() ) {
		utility_exit_with_message("Movers returning multiple poses are unsupported by GenericMontoCarloMover.");
	}
	if ( saved_accept_file_name_ != "" ) {
		TR<<"Saving initial pose entering the MC trajectory, for use in checkpoint recovery. Checkpointing filename: "<<saved_accept_file_name_<<std::endl;
		if ( basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments]() ) { //if pose has comment field then dump it with comments
			//std::ofstream out( saved_accept_file_name_.c_str() );
			pose.dump_pdb(saved_accept_file_name_);
			/*              out << "##Begin comments##" << std::endl;
			TR<<"SJF DEbug 1"<<std::endl;
			using namespace std;
			map< string, string > const comments = core::pose::get_all_comments(pose);
			for( auto const & comment : comments ) {
			out << comment.first<<" "<<comment.second << std::endl;
			TR<<"writing comment "<<i->first<<" "<<i->second<<std::endl;
			}
			TR<<"SJF debug 2"<<std::endl;
			out << "##End comments##" << std::endl;
			*/
		} else {
			pose.dump_pdb( saved_accept_file_name_ );
		}

	}

	PoseOP initial_pose( new Pose( pose ) );
	//(re)initialize MC statistics
	if ( preapply_ ) {
		core::pose::Pose tmppose;
		reset( tmppose );
	} else {
		reset( pose );
	}
	protocols::moves::MoverStatus ms( FAIL_RETRY );
	core::Size accept( 0 ), reject( 0 );
	using namespace protocols::rosetta_scripts;

	ParsedProtocolOP mover_pp( utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol > ( mover_ ) );
	// if the mover isn't a parsed protocol, mover_pp will be NULL and adaptive_movers basically does nothing -- we throw an error
	if ( adaptive_movers() && !mover_pp ) {
		throw utility::excn::EXCN_BadInput( "In GenericMonteCarlo, adaptive_movers is true, but the user-specified mover is not a ParsedProtocol.  The user-specified mover must be a ParsedProtocol in order to use adaptive_movers." );
	}

	if ( adaptive_movers() ) {
		bool is_single_random = ( mover_pp->mode() == "single_random" );
		if ( !is_single_random ) { // dig in one level (at most) to find the correct ParsedProtocol; if this becomes more generally useful then it would make sense to generatlize this to look for all parsedprotocols of type single_random that are being called by the MC mover. A simple recursion could do it, but I'm not sure how useful this would be
			for ( ParsedProtocol::MoverFilterPair const & mfp : *mover_pp ) {
				ParsedProtocolOP tmp( utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol > ( mfp.first.first ) );
				if ( tmp && tmp->mode() == "single_random" ) { /// the parsedprotocol mover must be run in mode single_random for the apply_probabilities to be modified
					mover_pp = tmp;
					is_single_random = true;
				}
			}//foreach mfp
		}//fi mover_pp && is_single_random
		runtime_assert( is_single_random ); /// yes, we found a single-random parsedprotocol mover somewhere; notice that this part boils down to finding the first parsed protocol of single-random mode up to a depth level of 2
	}//fi adaptive_movers()
	utility::vector1< core::Size > mover_accepts;// count how many accepts each mover in the parsedprotocol made
	mover_accepts.clear();
	if ( adaptive_movers() ) {
		runtime_assert( mover_pp != nullptr );
		runtime_assert( mover_pp->mode() == "single_random" );
		mover_accepts = utility::vector1< core::Size >( mover_pp->size(), 1 ); /// each mover gets a pseudocount of 1. This ensures that the running probability of each mover never goes to 0
	}
	for ( Size i=load_trial_number_from_checkpoint( pose ); i<=maxtrials_; i++ ) {
		TR<<"Trial number: "<<i<<std::endl;
		if ( i > 1 && adaptive_movers() && i % adaptation_period() == 0 ) {
			/// The probability for each mover within a single-random parsedprotocol is determined by the number of accepts it had during the previous adaptation period:
			/// each mover is assigned a pseducount of 1, and then any additional accept favors it over others. At the adaptation stage, the total number of accepts (including pseudocounts) is used to normalize the individual movers' number of accepts and the probability is the mover's accepts / by the total accepts
			core::Size sum_prev_accepts( 0 );
			for ( core::Size const a : mover_accepts ) {
				sum_prev_accepts += a;
			}
			utility::vector1< core::Real > new_probabilities;
			new_probabilities = utility::vector1< core::Real >( mover_accepts.size(), 0.0 );
			TR<<"Adapting running probabilities: old/new probabilities per mover\n";
			for ( core::Size ma = 1; ma<= mover_accepts.size(); ma++ ) {
				new_probabilities[ ma ] = ( core::Real )mover_accepts[ ma ] / ( core::Real )sum_prev_accepts;
				TR<<mover_pp->apply_probability()[ ma ]<<"/"<<new_probabilities[ ma ]<<'\n';
			}
			TR.flush();
			mover_pp->apply_probability( new_probabilities );
			mover_accepts = utility::vector1< core::Size >( mover_accepts.size(), 1 );
		}
		Pose store_pose( pose );
		// Mover apply
		mover_->apply( pose );
		ms = mover_->get_last_move_status();
		if ( ms == FAIL_RETRY ) {
			TR.Warning << "Mover failed. The mover, " << mover_->get_name() << ", is performed again. " << std::endl;
			pose = store_pose;
			continue;
		} else if ( ms == FAIL_DO_NOT_RETRY || ms == FAIL_BAD_INPUT ) {
			TR.Error << "Mover failed. Exit from GenericMonteCarloMover." << std::endl;
			break;
		}

		pose.energies().clear();
		// check stopping condition
		bool const stop( ( mover_stopping_condition_ && mover_stopping_condition_->obj ) ||
			( stopping_condition() && stopping_condition()->apply( pose ) ) ||
			( max_accepted_trials() != 0 && max_accepted_trials() == accept ) );
		if ( stop ) {
			TR<<"MC stopping condition met at trial "<<i;
			if ( max_accepted_trials() != 0 && max_accepted_trials() == accept ) {
				TR << " because maximum number of accepted moves was achieved";
			}
			TR << std::endl;
			break;
		}

		// MonteCarlo
		if ( preapply_ && i==1 ) { // Auto-accept first application in order to deal with movers that e.g. change the length of the pose.
			reset( pose );
		} else {
			if ( ! boltzmann( pose ) ) { // evaluate pose by scorefxn_ or filter_.report_sm()
				pose = *last_accepted_pose();
				reject++;
			} else {
				accept++;
				if ( adaptive_movers() ) {
					mover_accepts[ mover_pp->last_attempted_mover_idx() ]++;
				}
			}
		}

		if ( !drift_ ) {
			pose = (*initial_pose); // set back pose to initial one, of course this is not way of Monte Carlo
		}

		// Iterate over the list of triggers, executing each of them in turn
		fire_all_triggers(i, maxtrials_, pose, score_function());
		save_trial_number_to_checkpoint( i+1 );
	} // i<=maxtrials_

	// Output final diagnositics, for potential tuning
	show_scores(TR);
	show_counters(TR);
	TR<<"Finished MC. Out of "<<maxtrials_<<" "<<accept<<" accepted "<<" and "<<reject<<" rejected."<<std::endl;
	// Recover pose that have the lowest score, or the last accepted pose, as appropriate
	if ( recover_low_ ) {
		recover_low( pose );
	} else {
		if ( last_accepted_pose() ) {
			pose = *last_accepted_pose();
		}
	}
	pose.energies().clear();

}// apply

/// @brief generate a set of random numbers to accompany the filters or scorefxn
utility::vector1< core::Real >
GenericMonteCarloMover::generate_random() const {
	utility::vector1< core::Real > randoms;
	core::Size vec_size( filters_.size() );
	if ( !vec_size ) {
		vec_size = 1;
	}
	for ( core::Size i=1; i<=vec_size; ++i ) {
		randoms.push_back( numeric::random::rg().uniform() );
	}
	return randoms;
}

/// @brief parse xml file
void
GenericMonteCarloMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &filters, Movers_map const &movers, Pose const & )
{
	//using core::pack::task::operation::TaskOperation;
	//using core::pack::task::TaskFactoryOP;
	//using core::pack::task::TaskFactory;

	maxtrials_ = tag->getOption< core::Size >( "trials", 10 );
	max_accepted_trials_ = tag->getOption< core::Size >( "max_accepted_trials", maxtrials_ );
	temperature_ = tag->getOption< Real >( "temperature", 0.0 );
	task_scaling_ = tag->getOption< core::Size >( "task_scaling", 5 );
	adaptive_movers( tag->getOption< bool >( "adaptive_movers", false ) );
	if ( adaptive_movers() ) {
		adaptation_period( tag->getOption< core::Size >( "adaptation_period", std::max( (int) maxtrials_ / 10, 10 ) ) );
	}

	String const  user_defined_mover_name_( tag->getOption< String >( "mover_name" ,""));
	if ( data.has( "stopping_condition", user_defined_mover_name_ ) ) {
		TR<<user_defined_mover_name_<<" defines its own stopping condition, and GenericMC will respect this stopping condition"<<std::endl;
		mover_stopping_condition_ = data.get_ptr< basic::datacache::DataMapObj< bool > >( "stopping_condition", user_defined_mover_name_ );
	}

	String const filter_name( tag->getOption< String >( "filter_name", "true_filter" ) );
	auto  find_mover ( movers.find( user_defined_mover_name_ ));
	auto find_filter( filters.find( filter_name ));
	if ( find_mover == movers.end() && user_defined_mover_name_ != "" ) {
		TR.Error << "mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if ( find_filter == filters.end() ) {
		TR.Error << "filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}
	if ( user_defined_mover_name_ != "" ) {
		mover_ = find_mover->second;
	}

	if ( adaptive_movers() ) { /// adaptive movers only works if the mover being called is of type parsedprotocol
		runtime_assert( utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::ParsedProtocol >( mover_ ) != nullptr );
	}

	bool const adaptive( tag->getOption< bool >( "adaptive", true ) );
	sample_type_ = tag->getOption< String >( "sample_type", "low" );
	add_filter( find_filter->second->clone(), adaptive, temperature_, sample_type_ );
	String const sfxn ( tag->getOption< String >( "scorefxn_name", "" ) );

	if ( sfxn != "" ) {
		//scorefxn_ = new ScoreFunction( *data.get< ScoreFunction * >( "scorefxns", sfxn ));
		scorefxn_ = data.get< ScoreFunction * >( "scorefxns", sfxn )->clone();   //fpd use clone
		TR << "Score evaluation during MC is done by" << sfxn << ", ";
		if ( !keep_filters_ ) {
			TR << filter_name << " is ignored." << std::endl;
			filters_.clear();
		}
	} else {
		scorefxn_ = nullptr;
	}

	parse_task_operations( tag, data, filters, movers );

	runtime_assert_string_msg( filter_name != "true_filter" || scorefxn_, "You need to set filter_nam2e or scorefxn_name for MC criteria." );

	if ( filters_.size() == 0 ) {
		TR << "Apply mover of " << user_defined_mover_name_ << ", and evaluate score by " << sfxn
			<< " at Temperature=" << temperature_ << ", ntrails= " << maxtrials_ << std::endl;
	}

	stopping_condition( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "stopping_condition", "false_filter" ), filters ) );
	if ( tag->hasOption( "stopping_condition" ) ) {
		TR<<"Generic MC using stopping condition "<< stopping_condition()->get_user_defined_name()<<std::endl;
	}
	drift_ = tag->getOption< bool >( "drift", 1 );
	preapply_ = tag->getOption< bool >( "preapply", 1 ); // Default true for historical reasons
	recover_low_ = tag->getOption< bool >( "recover_low", 1 );
	boltz_rank_ = tag->getOption< bool >( "bolz_rank", 0 );

	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	for ( TagCOP const btag : branch_tags ) {
		if ( btag->getName() == "Filters" ) {
			utility::vector1< TagCOP > const filters_tags( btag->getTags() );
			for ( TagCOP const ftag : filters_tags ) {
				String const filter_name( ftag->getOption< String >( "filter_name" ) );
				auto find_filt( filters.find( filter_name ));
				if ( find_filt == filters.end() ) {
					TR.Error << "filter not found in map: \n" << tag << std::endl;
					runtime_assert( find_filt != filters.end() );
				}
				Real const temp( ftag->getOption< Real >( "temperature", 1 ) );
				bool const adap( ftag->getOption< bool >( "adaptive", true ));
				String const samp_type( ftag->getOption< String >( "sample_type", "low" ));
				bool const rank( ftag->getOption< bool >( "rank", false ) );
				if ( rank ) {
					if ( rank_by_filter_ != 1 ) {
						TR.Warning << "Multiple filters set to rank! Using most recent (currently "<< filter_name << ")." << std::endl;
					}
					if ( boltz_rank_ ) {
						TR.Warning << "Setting of rank on sub-filter "<< filter_name << " will be ignored as parent is set to Boltzmann rank all." << std::endl;
					}
				}
				add_filter( find_filt->second, adap, temp, samp_type, rank );
			} //foreach ftag
		} else { // fi Filters
			utility_exit_with_message( "tag name " + btag->getName() + " unrecognized." );
		}
	}//foreach btag

	if ( !drift_ ) {
		TR << "Pose is set back to initial pose every after applying mover and score evaluation." << std::endl;
	}

	if ( sample_type_ == "high" ) {
		TR << "Pose that have higher score is sampled." << sfxn << std::endl;
	}

	saved_accept_file_name_ = tag->getOption< std::string >( "saved_accept_file_name", "" );
	saved_trial_number_file_ = tag->getOption< std::string >( "saved_trial_number_file", "" );
	if ( tag->hasOption( "mover_tag" ) ) { // mover_tag isn't used. It's good to have this if needed in future to communicate splice -> monte carlo
		mover_tag_ = basic::datacache::get_set_from_datamap< basic::datacache::DataMapObj< std::string > >( "tags", tag->getOption< std::string >( "mover_tag" ), data );
	}
	reset_baselines( tag->getOption< bool >( "reset_baselines", true ) );
	progress_file( tag->getOption< std::string >( "progress_file", "" ) );
	initialize();
}

/// @brief parse "task_operations" XML option
void GenericMonteCarloMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	Movers_map const &
)
{
	if ( ( tag->hasOption("task_operations") ) ) {
		TR << "Found a task operation" << std::endl;
		TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
		if ( new_task_factory == nullptr ) return;
		task_factory( new_task_factory );
	}
}


void GenericMonteCarloMover::fire_all_triggers(
	Size cycle,
	Size num_cycles,
	const Pose& pose,
	ScoreFunctionOP scoring)
{
	boost::unordered_map<Size, GenericMonteCarloMoverTrigger>::iterator i;
	for ( i = triggers_.begin(); i != triggers_.end(); ++i ) {
		GenericMonteCarloMoverTrigger& t = i->second;
		bool rescore = t(cycle, num_cycles, pose, scoring);

		if ( rescore ) {
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
	if ( i == triggers_.end() ) {
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

std::string
GenericMonteCarloMover::saved_accept_file_name() const{
	return saved_accept_file_name_;
}

void
GenericMonteCarloMover::saved_accept_file_name( std::string const s ){
	saved_accept_file_name_ = s;
}

std::string
GenericMonteCarloMover::saved_trial_number_file() const{
	return saved_trial_number_file_;
}

void
GenericMonteCarloMover::saved_trial_number_file( std::string const s ){
	saved_trial_number_file_ = s;
}

std::string GenericMonteCarloMover::get_name() const {
	return mover_name();
}

std::string GenericMonteCarloMover::mover_name() {
	return "GenericMonteCarlo";
}

std::string gen_mc_ct_namer( std::string const & element_name ) {
	return "generic_monte_carlo_" + element_name + "_type";
}


utility::tag::XMLSchemaComplexTypeGeneratorOP
GenericMonteCarloMover::define_composition_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	typedef XMLSchemaAttribute Attr;
	AttributeList attlist;

	attlist + Attr::attribute_w_default(
		"trials", xsct_non_negative_integer,
		"The number of trials in a GenericMonteCarlo run. "
		"Will be automatically calculated when a task is included. "
		"See task_scaling for more information. ",
		"10" )
		+ Attr::attribute_w_default(
		"temperature", xsct_real,
		"Temperature value which affects the acceptance of new solutions "
		"according to the MC criterion",
		"0.0" )
		+ Attr::attribute_w_default(
		"task_scaling", xsct_non_negative_integer,
		"If a task is included, the number of designable residues will "
		"be calculated and the number of trials will be automatically "
		"set as task_scaling * (number designable residues). "
		"For example, if there are 10 designable residues and "
		"task_scaling is 5 (the default) the number of trials will be 50.",
		"5" )
		+ Attr::attribute_w_default(
		"adaptive_movers", xsct_rosetta_bool,
		"If the mover you call or a submover of that mover is of type "
		"ParsedProtocol with mode single-random, then GenericMonteCarlo "
		"can 'learn' the best sampling strategy by adapting the apply "
		"probabilities of individual movers within that ParsedProtocol. "
		"For each adaptation period (say 20 mover applies) the number "
		"of accepts of each submover is recorded (with pseudocounts of "
		"1 for each mover) and during the next adaptation period the apply "
		"probabilities of the submovers in the ParsedProtocol are adjusted "
		"according in proportion to the acceptance probabilities of the "
		"previous stage. Due to the pseudocounts, all movers have at "
		"least some chance of being called.",
		"false" )
		+ Attr(
		"adaptation_period", xsct_non_negative_integer,
		"goes together with adaptive_movers, defined above. Defaults to "
		"max( max_trials/10, 10 ) but can be set to any integer." )
		+ Attr(
		"mover_name", xs_string,
		"Mover to be used for sampling." )
		+ Attr::attribute_w_default(
		"filter_name", xs_string,
		"Single filter case. See Filters subcategory for applying multiple filters.",
		"true_filter" )
		+ Attr::attribute_w_default(
		"adaptive", xsct_rosetta_bool,
		"If the mover you call or a submover of that mover is of type "
		"ParsedProtocol with mode single-random, then GenericMonteCarlo can "
		"'learn' the best sampling strategy by adapting the apply "
		"probabilities of individual movers within that ParsedProtocol. "
		"For each adaptation period (say 20 mover applies) the number of "
		"'accepts of each submover is recorded (with pseudocounts of 1 for "
		"each mover) and during the next adaptation period the apply "
		"probabilities of the submovers in the ParsedProtocol are "
		"adjusted according in proportion to the acceptance probabilities "
		"of the previous stage. Due to the pseudocounts, all movers have "
		"at least some chance of being called.",
		"true" )
		+ Attr::attribute_w_default(
		"sample_type", xs_string,
		"low - MC samples by minimizing the structure score "
		"(REU: the more negative the better); "
		"high - MC samples by maximizing the structure score",
		"low" )
		+ Attr(
		"scorefxn_name", xs_string,
		"As alternative to scoring by filters (See Filters subcategory), "
		"structures may be scored by a given scoring function." );

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + Attr::attribute_w_default(
		"stopping_condition", xs_string,
		"stops before trials are done if a filter evaluates to true",
		"false_filter" )
		+ Attr::attribute_w_default(
		"drift", xsct_rosetta_bool,
		"true - the state of the pose at the end of the previous iteration "
		"will be the starting state for the next iteration; false - the state "
		"of the pose at the start of each iteration will be reset to the "
		"state when the mover is first called ( Of course, this is not MC ).",
		"true" )
		+ Attr::attribute_w_default(
		"preapply", xsct_rosetta_bool,
		"true - Automatically accept the first application of the sub-mover, "
		"ignoring the Boltzmann criterion. false - apply Boltzmann accept/reject "
		"to all applications of the mover. Though defaulting to true for "
		"historical reasons, it is highly recommended to set this to false "
		"unless you know you need it to be true.",
		"true" )
		+ Attr::attribute_w_default(
		"recover_low", xsct_rosetta_bool,
		"true - at the end of application, the pose is set to the lowest "
		"(or highest if sample_type=\"high\") scoring pose; false - the "
		"pose after apply completes is the last accepted pose",
		"true" )
		+ Attr::attribute_w_default(
		"bolz_rank", xsct_rosetta_bool,
		"For use with multiple filters. If no sub-filters are set with "
		"rank=1, the first filter is used for ranking. As a special case, "
		"if boltz_rank is set to true, the ranking score is a temperature-"
		"weighted sum of all filter values. (This value is equivalent to the "
		"effective value optimized by the MC protocol.) This boltz_rank score "
		"is computed by the equation SUM( multiplier * filter_value / filter_"
		"temperature ) over all filter values, where filter_value is the value "
		"returned by the filter and multiplier is 1 if the filter sample_type "
		"is low and -1 if the filter sample_type is high.",
		"false" )
		+ Attr(
		"saved_accept_file_name", xs_string,
		"save the most recent accepted structure in a temporary PDB? This allows "
		"recovery by checkpointing. Note that different processes would need to "
		"work from different directories or somehow control the checkpointing "
		"file name, else confusion will reign." )
		+ Attr(
		"saved_trial_number_file", xs_string,
		"Checkpointing file for the current trial number. "
		"Allows the mover to recover after failure." )
		+ Attr(
		"mover_tag", xs_string,
		"this is used by the called movers to set a certain tag. If saved_accept_file_name_ is set, then at exit the tag coming from the chosen mover is written to disk as, {saved_accept_file_name}.mover_tag. To work, mover_tag_ must be exposed to the movers being called.")
		+ Attr::attribute_w_default(
		"reset_baselines", xsct_rosetta_bool,
		"If the filter is of type Sigmoid/Operator/CompoundStatement, look for "
		"all subfilters of type Sigmoid and reset their baseline to the pose's "
		"current filter evaluation at trial=1. Useful in cases where you want "
		"to set the thresholds relative to the pose's evaluation at the start "
		"of the MC trajectory.",
		"true" )
		+ Attr(
		"progress_file", xs_string,
		"If specified opens a file in which each trial's outcome is reported "
		"(trial number, accept/reject, filter value, and pose sequence). "
		"Useful to monitor progress in MC" );

	// Subelements!
	// Filters
	AttributeList filter_attributes;
	filter_attributes
		+ Attr( "filter_name", xs_string, "name of the filter declaration" )
		+ Attr::attribute_w_default(
		"temperature", xsct_real,
		"Filter specific temperature",
		"1" )
		+ Attr::attribute_w_default(
		"adaptive", xsct_rosetta_bool,
		"Incorporate this filter when using adaptive scores.",
		"true" )
		+ Attr::attribute_w_default(
		"sample_type", xs_string,
		"high - the filter evalutas good structures with positive scores; "
		"low - the filter evaulates good structures with negative scores",
		"low" )
		+ Attr::attribute_w_default(
		"rank", xsct_rosetta_bool,
		"Give this filter a rank of 1. Can only be applied for one filter. "
		"A filter with rank 1 will score the structure solitarily.",
		"false" );
	XMLSchemaSimpleSubelementList filters_subelement;
	filters_subelement.add_simple_subelement(
		"AND", filter_attributes,
		"Filter that may be used as ensemble to score the strucutres." );

	XMLSchemaComplexTypeGenerator filters_ct_gen;
	// no attributes, not even a name attribute
	filters_ct_gen.element_name( "Filters" )
		.complex_type_naming_func( & gen_mc_ct_namer )
		.set_subelements_repeatable( filters_subelement )
		.description(
		"Multiple filters can be defined for an MC mover. These filters are then applied "
		"sequentially in the order listed and only if the pose passes the Metropolis "
		"criterion for all filters is it accepted. This allows the extension of MC to a "
		"multicriterion framework where more than one criterion is optimized, say the total "
		"score and the binding energy. See demos/rosetta_scripts/experimental/computational_"
		"affinity_maturation_strategy2 for an example. "
		"It's recommended to list the computationally expensive filters last, "
		"as later filters will only be calculated if the earlier filters all pass." );

	filters_ct_gen.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList genmc_subelements;
	genmc_subelements.add_already_defined_subelement( "Filters", & gen_mc_ct_namer );
	XMLSchemaComplexTypeGeneratorOP genmc_ct_gen( new XMLSchemaComplexTypeGenerator);
	genmc_ct_gen->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.description(
		"Allows sampling structures by MonteCarlo with a mover. "
		"The score evaluation of pose during MC are done by Filters that can do "
		"report_sm(), not only ScoreFunctions")
		.set_subelements_repeatable( genmc_subelements );

	return genmc_ct_gen;
}


void GenericMonteCarloMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( mover_name() )
		.write_complex_type_to_schema( xsd );
}

std::string GenericMonteCarloMoverCreator::keyname() const {
	return GenericMonteCarloMover::mover_name();
}

protocols::moves::MoverOP
GenericMonteCarloMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new GenericMonteCarloMover );
}

void GenericMonteCarloMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GenericMonteCarloMover::provide_xml_schema( xsd );
}

} // ns simple_moves
} // ns protocols
