// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rosetta_scripts/ParsedProtocol.cc
/// @brief  A mover that applies a protocol parsed from a ROSETTASCRIPTS script
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Luki Goldschmidt (lugo@uw.edu)

// Unit Headers
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/ParsedProtocolCreator.hh>

// Project Headers
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pose/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/moves/ResId.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/moves/MoverFactory.hh>

#include <protocols/rosetta_scripts/MultiplePoseMover.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

// JD2 headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

//Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// C++ headers
#include <map>
#include <string>
#include <algorithm>

#define foreach BOOST_FOREACH

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.ParsedProtocol" );
static basic::Tracer TR_call_order( "protocols.rosetta_scripts.ParsedProtocol_call_order" );
static basic::Tracer TR_report( "protocols.rosetta_scripts.ParsedProtocol.REPORT" );
static numeric::random::RandomGenerator RG(48569);

typedef core::Real Real;
typedef core::pose::Pose Pose;

using namespace core;
using namespace std;

std::string
ParsedProtocolCreator::keyname() const
{
	return ParsedProtocolCreator::mover_name();
}

protocols::moves::MoverOP
ParsedProtocolCreator::create_mover() const {
	return new ParsedProtocol;
}

std::string
ParsedProtocolCreator::mover_name()
{
	return "ParsedProtocol";
}

ParsedProtocol::ParsedProtocol() :
	protocols::moves::Mover( "ParsedProtocol" ),
	final_scorefxn_( 0 ), // By default, don't rescore with any scorefunction.
	mode_("sequence"),
	last_attempted_mover_idx_( 0 ),
	report_call_order_( false ),
	last_mover_(NULL),
	resume_support_(false)
{
}

ParsedProtocol::~ParsedProtocol() {
}
/// @detailed Takes care of the docking, design and filtering moves. pre_cycle and pose_cycle can
/// be setup in derived classes to setup variables before and after these cycles.
void
ParsedProtocol::apply( Pose & pose )
{
	TR.Debug << "ParsedProtocol::apply BEGIN 0x" << this << std::endl;
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobOutputterOP job_outputter( protocols::jd2::JobDistributor::get_instance()->job_outputter() );
		job_outputter->add_output_observer( this );
	}

	protocols::moves::Mover::set_last_move_status( protocols::moves::FAIL_RETRY );
	pose.update_residue_neighbors();

	//fpd search the mover-filter pairs backwards for movers that have remaining poses
	MoverFilterVector::const_reverse_iterator rmover_it = movers_.rbegin();
  const MoverFilterVector::const_reverse_iterator movers_crend = movers_.rend();

  if(resume_support_ && mode_ == "sequence") {
		for ( rmover_it=movers_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
			core::pose::PoseOP checkpoint = (*rmover_it).first.first->get_additional_output();

			// otherwise continue where we left off
			if (checkpoint) {
				std::string const mover_name( rmover_it->first.first->get_name() );

				// if mode_ is not 'sequence' then checkpointing is unsupported
				if (checkpoint && mode_ != "sequence")
					utility_exit_with_message("Mover "+mover_name+" returned multiple poses in a ParsedProtocol with mode!=sequence");

				TR<<"=======================RESUMING FROM "<<mover_name<<"======================="<<std::endl;
				pose = *checkpoint;

				if( ! apply_filter( pose, *rmover_it) ) {
//					final_score(pose);
					if ( protocols::jd2::jd2_used() ) {
						protocols::jd2::JobOutputterOP job_outputter( protocols::jd2::JobDistributor::get_instance()->job_outputter() );
						job_outputter->remove_output_observer( this );
					}
					return;
				} else {
					break;
				}
			}
		}
	} else {
    rmover_it = movers_.rend();
  }

	if(mode_ == "sequence"){
		sequence_protocol(pose, rmover_it.base());
	}else if(mode_ =="random_order"){
		random_order_protocol(pose);
	}else if(mode_ =="single_random"){
		random_single_protocol(pose);
	}else
	{
		TR <<"WARNING: mode is " << mode_ << " .This is not a valid ParsedProtocol Mode, your pose is being ignored" <<std::endl;
	}

	if ( get_last_move_status() == protocols::moves::MS_SUCCESS ) {// no point scoring a failed trajectory (and sometimes you get etable vs. pose atomset mismatches
		final_score(pose);
	}

	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobOutputterOP job_outputter( protocols::jd2::JobDistributor::get_instance()->job_outputter() );
		job_outputter->remove_output_observer( this );
	}

	TR.Debug << "ParsedProtocol::apply   END 0x" << this << std::endl;
}

std::string
ParsedProtocol::get_name() const {
	return ParsedProtocolCreator::mover_name();
}

void ParsedProtocol::final_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	final_scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionCOP ParsedProtocol::final_scorefxn() const
{
	return final_scorefxn_;
}

void
ParsedProtocol::final_score(core::pose::Pose & pose) const {
	core::scoring::ScoreFunctionCOP scorefxn( final_scorefxn() );

	if( ! scorefxn ) { return; }

	if (core::pose::symmetry::is_symmetric(pose)) {
		scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn );
	}
	(*scorefxn)(pose);
}

void
ParsedProtocol::report_all( Pose const & pose ) const {
  foreach ( MoverFilterPair mover_pair, movers_ ) {
    if ( mover_pair.report_filter_at_end_ ) {
      TR_report<<"============Begin report for "<<mover_pair.second->get_user_defined_name()<<"=================="<<std::endl;
      mover_pair.second->report( TR_report, pose );
      TR_report<<"============End report for "<<mover_pair.second->get_user_defined_name()<<"=================="<<std::endl;
    }
	}
	TR_report.flush();
}

void
ParsedProtocol::add_values_to_job( Pose const & pose, protocols::jd2::JobOP job ) const {
	for ( MoverFilterVector::const_iterator mover_it = movers_.begin();
       mover_it != movers_.end(); ++mover_it ) {
    if ( mover_it->report_filter_at_end_ ) {
      core::Real const filter_value( (*mover_it).second->report_sm( pose ) );
      if ( filter_value > -9999 ) {
        job->add_string_real_pair((*mover_it).second->get_user_defined_name(), filter_value);
      }
    }
  }
}

void
ParsedProtocol::report_filters_to_job( Pose const & pose) const {
	using protocols::jd2::JobDistributor;
	protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
	add_values_to_job( pose, job_me );
}

void
ParsedProtocol::report_filters_to_pose( Pose & pose ) {
	for( utility::vector1< MoverFilterPair >::const_iterator mover_it = movers_.begin();
			mover_it!=movers_.end(); ++mover_it ) {
		core::Real const filter_value( (*mover_it).second->report_sm( pose ) );
		if( filter_value > -9999 ) {
			setPoseExtraScores(pose, (*mover_it).second->get_user_defined_name(), (float)filter_value);
		}
	}
}

//void
//ParsedProtocol::report_all_sm( std::map< std::string, core::Real > & score_map, Pose const & pose ) const {
//	for( utility::vector1< mover_filter_pair >::const_iterator mover_it = movers_.begin();
//			mover_it!=movers_.end(); ++mover_it ) {
//		core::Real const filter_value( (*mover_it).second->report_sm( pose ) );
//		if( filter_value >= -9999 ) {
//			score_map[ (*mover_it).second->get_user_defined_name() ] = filter_value;
//		}
//	}
//}

ParsedProtocol::iterator
ParsedProtocol::begin(){
	return movers_.begin();
}

ParsedProtocol::const_iterator
ParsedProtocol::begin() const{
	return movers_.begin();
}

ParsedProtocol::iterator
ParsedProtocol::end(){
	return movers_.end();
}

ParsedProtocol::const_iterator
ParsedProtocol::end() const{
	return movers_.end();
}

/// @details sets resid for the constituent filters and movers
void
ParsedProtocol::set_resid( core::Size const resid ){
	for( iterator it( movers_.begin() ); it!=movers_.end(); ++it ){
		using namespace protocols::moves;
		modify_ResId_based_object( it->first.first, resid );
		modify_ResId_based_object( it->second, resid );
	}
}

protocols::moves::MoverOP ParsedProtocol::clone() const
{
	return protocols::moves::MoverOP( new protocols::rosetta_scripts::ParsedProtocol( *this ) );
}

std::pair< moves::MoverOP, std::string >
parse_mover_subtag( utility::tag::TagCOP const tag_ptr,
                    basic::datacache::DataMap& data,
                    protocols::filters::Filters_map const& filters,
                    protocols::moves::Movers_map const& movers,
                    core::pose::Pose const& pose ) {
	using namespace protocols::moves;

  MoverOP mover_to_add = NULL;
  std::string mover_name; // user must specify a mover name. there is no valid default.

  runtime_assert( !( tag_ptr->hasOption("mover_name") && tag_ptr->hasOption("mover") ) );
  if( tag_ptr->hasOption( "mover_name" ) ){
    mover_name = tag_ptr->getOption<string>( "mover_name" );
    Movers_map::const_iterator find_mover( movers.find( mover_name ) );
    if( find_mover == movers.end() ) {
      throw utility::excn::EXCN_RosettaScriptsOption("Mover " + mover_name + " not found in map");
    }
    mover_to_add = find_mover->second;
  }
  else if( tag_ptr->hasOption( "mover" ) ){
    mover_name = tag_ptr->getOption<string>( "mover" );
    Movers_map::const_iterator find_mover( movers.find( mover_name ) );
    if( find_mover == movers.end() ) {
      throw utility::excn::EXCN_RosettaScriptsOption("Mover " + mover_name + " not found in map");
    }
    mover_to_add = find_mover->second;
  } else if( tag_ptr->getName() != "Add" ){
    MoverOP new_mover( MoverFactory::get_instance()->newMover( tag_ptr, data, filters, movers, pose ) );
    assert( new_mover );

    if( tag_ptr->hasOption("name") ){
      TR.Warning << "The mover named '" << tag_ptr->getOption<std::string>( "name" )
      << " will not be accessible, since it was defined on-the-fly in the PROTOCOLS"
      << " section." << std::endl;
    }

    TR.Info << "Defined on-the-fly '" << new_mover->name() << "' mover of type "
    << tag_ptr->getName() << std::endl;
    mover_to_add = new_mover;
    mover_name = "OnTheFly("+tag_ptr->getName()+")";
  } else {
    mover_to_add = new NullMover();
    mover_name = "NULL_MOVER";
  }

  if( data.has( "stopping_condition", mover_name ) && tag_ptr->getParent()->hasOption( "name" ) ){
    TR.Info << "ParsedProtocol's mover " << mover_name
            << " requests its own stopping condition. This ParsedProtocol's stopping_condition will point at the mover's"
            << std::endl;
    data.add( "stopping_condition", tag_ptr->getParent()->getOption< std::string >( "name" ), data.get< basic::datacache::DataMapObj< bool > * >( "stopping_condition", mover_name ) );
  }
  
  return std::make_pair( mover_to_add, mover_name );
}

void
ParsedProtocol::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & pose )
{
	using namespace protocols::moves;
	using namespace utility::tag;

	TR<<"ParsedProtocol mover with the following movers and filters\n";

	mode_=tag->getOption<string>("mode", "sequence");
	if(mode_ != "sequence" && mode_ != "random_order" && mode_ != "single_random"){
		throw utility::excn::EXCN_RosettaScriptsOption("Error: mode must be sequence, random_order, or single_random");
	}

	utility::vector0< TagCOP > const & dd_tags( tag->getTags() );
	utility::vector1< core::Real > a_probability( dd_tags.size(), 1.0/dd_tags.size() );
	core::Size count( 1 );
	for( utility::vector0< TagCOP >::const_iterator dd_it=dd_tags.begin(); dd_it!=dd_tags.end(); ++dd_it ) {
		TagCOP const tag_ptr( *dd_it );

    std::pair< MoverOP, std::string > mover_add_pair = parse_mover_subtag( tag_ptr, data, filters, movers, pose );
    std::string const& mover_name( mover_add_pair.second );
    MoverOP mover_to_add( mover_add_pair.first );

		protocols::filters::FilterOP filter_to_add;
		bool filter_defined( false );

		std::string filter_name="true_filter"; // used in case user does not specify a filter name.
		runtime_assert( !( tag_ptr->hasOption("filter_name") && tag_ptr->hasOption( "filter" ) ) );
		if( tag_ptr->hasOption( "filter_name" ) ){
			filter_name = tag_ptr->getOption<string>( "filter_name", "true_filter" );
			filter_defined = true;
		} else if( tag_ptr->hasOption( "filter" ) ){
			filter_name = tag_ptr->getOption<string>( "filter", "true_filter" );
			filter_defined = true;
		}

		if( filter_defined ){
			protocols::filters::Filters_map::const_iterator find_filter( filters.find( filter_name ));
			if( find_filter == filters.end() ) {
				throw utility::excn::EXCN_RosettaScriptsOption("Filter " + filter_name + " not found in map");
			}
			filter_to_add = find_filter->second;
		} else {
			filter_to_add = new protocols::filters::TrueFilter;
		}

		TR << "added mover \"" << mover_name << "\" with filter \"" << filter_name << "\"\n";
		if( mode_ == "single_random" ){
			a_probability[ count ] = tag_ptr->getOption< core::Real >( "apply_probability", 1.0/dd_tags.size() );
			TR<<"and execution probability of "<<a_probability[ count ]<<'\n';
		}
		count++;

    bool const report_at_end( tag_ptr->getOption< bool >( "report_at_end", true ) );
    movers_.push_back( MoverFilterPair( mover_to_add->clone(), mover_name, filter_to_add, report_at_end ) );
	}
	if( mode_ == "single_random" )
		apply_probability( a_probability );
	report_call_order( tag->getOption< bool >( "report_call_order", false ) );

	resume_support_ = tag->getOption<bool>("resume_support", false);
	if(resume_support_)
		TR << "Legacy protocol resume support enabled." << std::endl;

	TR.flush();
}

/// @detailed Looks for any submovers that have additional output poses to process.
/// If any are found, run remainder of parsed protocol.
core::pose::PoseOP
ParsedProtocol::get_additional_output( )
{
	core::pose::PoseOP pose=NULL;
	MoverFilterVector::const_reverse_iterator rmover_it;
  const MoverFilterVector::const_reverse_iterator movers_crend = movers_.rend();

	if(!resume_support_) {
		// Get output from last specified mover
		const utility::vector1< MoverFilterPair >::const_reverse_iterator movers_crend = movers_.rend();
		for ( rmover_it=movers_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
			protocols::moves::MoverOP mover = (*rmover_it).first.first;
			if(mover && mover->get_name() != "NullMover") {
				return mover->get_additional_output();
			}
		}
		return NULL;
	}

	// Legacy Protocol resume support:

	//fpd search the mover-filter pairs backwards; look for movers that have remaining poses
	for ( rmover_it=movers_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
		core::pose::PoseOP checkpoint = (*rmover_it).first.first->get_additional_output();
		if (checkpoint) {
			std::string const mover_name( rmover_it->first.first->get_name() );
			pose = checkpoint;

			if ( ! apply_filter( *pose, *rmover_it) ) {
				return pose;
			} else {
				break;
			}
		}
	}

	// no saved poses?  return now
	if ( !pose ) return NULL;

	// if mode_ is not 'sequence' then checkpointing is unsupported
	if ( mode_ != "sequence" )
		utility_exit_with_message("ParsedProtocol returned multiple poses in a ParsedProtocol with mode!=sequence");

	// otherwise pick up from the checkpoint
	for ( MoverFilterVector::const_iterator mover_it = rmover_it.base();
       mover_it!=movers_.end(); ++mover_it ) {
    apply_mover( *pose, *mover_it );
    if ( !apply_filter( *pose, *mover_it ) ) {
			return pose;
		}
	}
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
	TR<<"setting status to success"<<std::endl;

	// report filter values to the job object as string_real_pair
	//report_filters_to_job( *pose );
	// report filter values to pose DataCache
	report_filters_to_pose( *pose );
	// report filter values to tracer output
	report_all( *pose );

	return pose;
}

void
ParsedProtocol::apply_mover( Pose & pose, MoverFilterPair const & mover_pair ) {
  std::cout << mover_pair.first.first << std::endl;
	std::string const mover_name( mover_pair.first.first->get_name() );
	std::string const mover_user_name( mover_pair.first.second);

	// If the mover about to be applied is a MultiplePoserMover,
	// tell it about the previous mover where it should pull poses from
	TR.Debug << "apply_mover_filter_pair: Last mover: " << ( last_mover_ ? last_mover_->get_name() : "(None)" ) << std::endl;
	TR.Debug << "apply_mover_filter_pair: This mover: " << ( mover_pair.first.first ? mover_pair.first.first->get_name() : "(None)" ) << std::endl;

	if(last_mover_ &&
		mover_pair.first.first &&
		mover_pair.first.first->get_name() == "MultiplePoseMover"
	) {
		// Is there a better way of handling this cast?
		protocols::moves::Mover & mover ( *mover_pair.first.first );
		protocols::rosetta_scripts::MultiplePoseMoverOP mp_mover ( dynamic_cast< protocols::rosetta_scripts::MultiplePoseMover * >( &mover ) );
		if(mp_mover)
			mp_mover->set_previous_mover(last_mover_);
		else
			throw utility::excn::EXCN_Msg_Exception("Dynamic cast of MultiplePoseMover failed; this should not happen");
	}

	mover_pair.first.first->set_native_pose( get_native_pose() );
	TR<<"=======================BEGIN MOVER "<<mover_name<<" - "<<mover_user_name<<"======================="<<std::endl;
	mover_pair.first.first->apply( pose );

	if(mover_pair.first.first && mover_pair.first.first->get_name()	!= "NullMover") {
		last_mover_ = mover_pair.first.first;
	}
}

bool
ParsedProtocol::apply_filter( Pose & pose, MoverFilterPair const & mover_pair) {
	std::string const filter_name( mover_pair.second->get_user_defined_name() );

  TR << "=======================BEGIN FILTER " << filter_name << "=======================" << std::endl;
	info().insert( info().end(), mover_pair.first.first->info().begin(), mover_pair.first.first->info().end() );
	pose.update_residue_neighbors();
	moves::MoverStatus status( mover_pair.first.first->get_last_move_status() );
  bool const pass( status==protocols::moves::MS_SUCCESS  && mover_pair.filter().apply( pose ) );
  if ( !mover_pair.report_filter_at_end_ ) { //report filter now
    using namespace protocols::jd2;
    protocols::jd2::JobOP job( JobDistributor::get_instance()->current_job() );
    core::Real const filter_value(  mover_pair.filter().report_sm( pose ) );
    job->add_string_real_pair( mover_pair.filter().get_user_defined_name(), filter_value );
    TR_report << "============Begin report for " << mover_pair.second->get_user_defined_name() << "==================" << std::endl;
    mover_pair.filter().report( TR_report, pose );
    TR_report << "============End report for " << mover_pair.second->get_user_defined_name() << "==================" << std::endl;
  }
  TR << "=======================END FILTER " << filter_name << "=======================" << std::endl;

  //filter failed -- set status in mover and return false
  if ( !pass ) {
    if ( status != protocols::moves::MS_SUCCESS ) {
			TR << "Mover " << mover_pair.first.first->get_name() << " reports failure!" << std::endl;
			protocols::moves::Mover::set_last_move_status( status );
		} else {
			TR << "Filter " << filter_name << " reports failure!" << std::endl;
      protocols::moves::Mover::set_last_move_status( protocols::moves::FAIL_RETRY );
		}
		return false;
	}
	return true;
}

void ParsedProtocol::finish_protocol(Pose & pose) {
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
  TR.Info << "setting status to success" << std::endl; // JRP/OFL that sounds like a debug statement to me...

	// report filter values to the job object as string_real_pair
	//report_filters_to_job( pose );
	// report filter values to pose DataCache
	report_filters_to_pose( pose );
	// report filter values to tracer output
	report_all( pose );

  using namespace protocols::jd2;
	protocols::jd2::JobOP job2( jd2::JobDistributor::get_instance()->current_job() );
	std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job2 ) );
	if( report_call_order() ){
		TR_call_order << job_name<<" ";
		foreach( MoverFilterPair const p, movers_ ) {
			TR_call_order<<p.first.second<<" ";
		}
		TR_call_order<<std::endl;
	}
}

void
ParsedProtocol::sequence_protocol( Pose & pose,
                                   MoverFilterVector::const_iterator mover_it_in ) {
  for ( MoverFilterVector::const_iterator mover_it = mover_it_in;
       mover_it!=movers_.end(); ++mover_it ) {
    apply_mover( pose, *mover_it );
    if ( !apply_filter( pose, *mover_it ) ) {
      return;
    }
  }
  // we're done! mark as success
  finish_protocol( pose );
}

void ParsedProtocol::random_order_protocol(Pose & pose){
	numeric::random::random_permutation(movers_.begin(),movers_.end(),RG);
  for ( MoverFilterVector::const_iterator it = movers_.begin(); it != movers_.end(); ++it ) {
    apply_mover( pose, *it );
    if ( !apply_filter( pose, *it ) ) {
      return;
    }
  }
	// we're done! mark as success
	finish_protocol( pose );
}

utility::vector1< core::Real >
ParsedProtocol::apply_probability() {
	core::Real sum( 0 );
	foreach( core::Real const prob, apply_probability_ )
		sum += prob;
	runtime_assert( sum >= 0.999 && sum <= 1.001 );
	return apply_probability_;
}

void
ParsedProtocol::apply_probability( utility::vector1< core::Real > a ){
	apply_probability_ = a;
	runtime_assert( apply_probability_.size() == movers_.size() );
	core::Real sum( 0 );
	foreach( core::Real const prob, apply_probability_ )
		sum += prob;
	runtime_assert( sum >= 0.999 && sum <= 1.001 );
}

void ParsedProtocol::random_single_protocol(Pose & pose){
	core::Real const random_num( RG.uniform() );
	core::Real sum( 0.0 );
	core::Size mover_index( 0 );
	foreach( core::Real const probability, apply_probability() ){
		sum += probability; mover_index++;
		if( sum >= random_num )
			break;
	}
	last_attempted_mover_idx( mover_index );
  apply_mover( pose, movers_[ mover_index ] );
  if ( !apply_filter( pose, movers_[ mover_index ] ) ) {
    return;
  }

	// we're done! mark as success
	finish_protocol( pose );
}

} //rosetta_scripts
} //protocols

