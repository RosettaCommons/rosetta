// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/ParsedProtocol.cc
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit Headers
#include <protocols/moves/ParsedProtocol.hh>
#include <protocols/moves/ParsedProtocolCreator.hh>
#include <protocols/moves/NullMover.hh>

#include <protocols/viewer/viewers.hh>
// Project Headers
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/MoverStatus.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <basic/Tracer.hh>


#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/ResId.hh>

// JD2 headers
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

//Numeric Headers
#include <numeric/random/random.hh>

// C++ headers
#include <map>
#include <string>
#include <set>
#include <algorithm>

//Auto Headers
#include <utility/options/keys/BooleanOptionKey.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace protocols {
namespace moves {

static basic::Tracer TR( "protocols.moves.ParsedProtocol" );
static basic::Tracer TR_report( "protocols.moves.ParsedProtocol.REPORT" );
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

/// @detailed Takes care of the docking, design and filtering moves. pre_cycle and pose_cycle can
/// be setup in derived classes to setup variables before and after these cycles.
void
ParsedProtocol::apply( Pose & pose )
{
//	runtime_assert( movers_.size() );

	protocols::moves::Mover::set_last_move_status( protocols::moves::FAIL_RETRY );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );

	pose.update_residue_neighbors();

	//fpd search the mover-filter pairs backwards for movers that have remaining poses
	utility::vector1< mover_filter_pair >::const_reverse_iterator rmover_it;
	const utility::vector1< mover_filter_pair >::const_reverse_iterator movers_crend = movers_.rend();
	if(mode_ == "sequence")
	{
		for ( rmover_it=movers_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
			core::pose::PoseOP checkpoint = (*rmover_it).first->get_additional_output();

			// otherwise continue where we left off
			if (checkpoint) {
				std::string const mover_name( rmover_it->first->get_name() );

				// if mode_ is not 'sequence' then checkpointing is unsupported
				if (checkpoint && mode_ != "sequence")
					utility_exit_with_message("Mover "+mover_name+" returned multiple poses in a ParsedProtocol with mode!=sequence");

				std::string const filter_name( rmover_it->second->get_user_defined_name() );
				TR<<"=======================RESUMING FROM "<<mover_name<<"======================="<<std::endl;
				pose = *checkpoint;
				TR<<"=======================BEGIN FILTER "<<filter_name<<"=======================\n{"<<std::endl;
				info().insert( info().end(), rmover_it->first->info().begin(), rmover_it->first->info().end() );
				pose.update_residue_neighbors();
				moves::MoverStatus status( (*rmover_it).first->get_last_move_status() );
				bool const pass( status==protocols::moves::MS_SUCCESS  && (*rmover_it).second->apply( pose ) );
				TR<<"\n}\n=======================END FILTER "<<filter_name<<"======================="<<std::endl;
				if( !pass ) {
					if( status != protocols::moves::MS_SUCCESS )
						protocols::moves::Mover::set_last_move_status( status );
					return;
				}
				break;
			}
		}
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

}

std::string
ParsedProtocol::get_name() const {
	return ParsedProtocolCreator::mover_name();
}

void
ParsedProtocol::report_all( Pose const & pose ) const {
	TR_report<<"=============Starting final report================"<<std::endl;
	foreach(mover_filter_pair mover_pair, movers_){
		TR_report<<"============Begin report for "<<mover_pair.second->get_user_defined_name()<<"=================="<<std::endl;
		mover_pair.second->report( TR_report, pose );
		TR_report<<"============End report for "<<mover_pair.second->get_user_defined_name()<<"=================="<<std::endl;
	}
	TR_report.flush();
}

void
ParsedProtocol::report_filters_to_job( Pose const & pose) const {
	using protocols::jd2::JobDistributor;
	protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
	for( utility::vector1< mover_filter_pair >::const_iterator mover_it = movers_.begin();
			 mover_it!=movers_.end(); ++mover_it ) {
		core::Real const filter_value( (*mover_it).second->report_sm( pose ) );
		if( filter_value > -9999 )
			job_me->add_string_real_pair((*mover_it).second->get_user_defined_name(), filter_value);
	}
}

void
ParsedProtocol::report_all_sm( std::map< std::string, core::Real > & score_map, Pose const & pose ) const {
	for( utility::vector1< mover_filter_pair >::const_iterator mover_it = movers_.begin();
		 mover_it!=movers_.end(); ++mover_it ) {
		 core::Real const filter_value( (*mover_it).second->report_sm( pose ) );
		 if( filter_value >= -9999 )
			score_map[ (*mover_it).second->get_user_defined_name() ] = filter_value;
	}
}

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
		modify_ResId_based_object( it->first, resid );
		modify_ResId_based_object( it->second, resid );
	}
}

void
ParsedProtocol::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap &,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & )
{
	using namespace protocols::moves;
	using namespace utility::tag;

	TR<<"ParsedProtocol mover with the following movers and filters\n";

	mode_=tag->getOption<string>("mode", "sequence");

	if(mode_ != "sequence" && mode_ != "random_order" && mode_ != "single_random"){
		utility_exit_with_message("Error: mode must be sequence, random_order, or single_random");
	}

	utility::vector0< TagPtr > const dd_tags( tag->getTags() );
	for( utility::vector0< TagPtr >::const_iterator dd_it=dd_tags.begin(); dd_it!=dd_tags.end(); ++dd_it ) {
		TagPtr const tag_ptr = *dd_it;

		MoverOP mover_to_add;
		protocols::filters::FilterOP filter_to_add;

		bool mover_defined( false ), filter_defined( false );

		std::string mover_name="null"; // user must specify a mover name. there is no valid default.
		runtime_assert( !( tag_ptr->hasOption("mover_name") && tag_ptr->hasOption("mover") ) );
		if( tag_ptr->hasOption( "mover_name" ) ){
			mover_name = tag_ptr->getOption<string>( "mover_name", "null" );
			mover_defined = true;
		}
		else if( tag_ptr->hasOption( "mover" ) ){
			mover_name = tag_ptr->getOption<string>( "mover", "null" );
			mover_defined = true;
		}
		//runtime_assert( mover_name ); // redundant with mover find below

		std::string filter_name="true_filter"; // used in case user does not specify a filter name.
		runtime_assert( !( tag_ptr->hasOption("filter_name") && tag_ptr->hasOption( "filter" ) ) );
		if( tag_ptr->hasOption( "filter_name" ) ){
			filter_name = tag_ptr->getOption<string>( "filter_name", "true_filter" );
			filter_defined = true;
		} else if( tag_ptr->hasOption( "filter" ) ){
			filter_name = tag_ptr->getOption<string>( "filter", "true_filter" );
			filter_defined = true;
		}

		if( mover_defined ){
			Movers_map::const_iterator find_mover( movers.find( mover_name ) );
			if( find_mover == movers.end() ) {
				TR.Error<<"mover not found in map. skipping:\n"<<tag_ptr<<std::endl;
				runtime_assert( find_mover != movers.end() );
				continue;
			}
			mover_to_add = find_mover->second;
		}	else {
			mover_to_add = new NullMover;
		}
		if( filter_defined ){
			protocols::filters::Filters_map::const_iterator find_filter( filters.find( filter_name ));
			if( find_filter == filters.end() ) {
				TR.Error<<"filter not found in map. skipping:\n"<<tag_ptr<<std::endl;
				runtime_assert( find_filter != filters.end() );
				continue;
			}
			filter_to_add = find_filter->second;
		} else {
			filter_to_add = new protocols::filters::TrueFilter;
		}
		add_mover( mover_to_add, filter_to_add );
		TR << "added mover \"" << mover_name << "\" with filter \"" << filter_name << "\"\n";
	}
	TR.flush();
}

/// @detailed Looks for any submovers that have additional output poses to process.
/// If any are found, run remainder of parsed protocol.
core::pose::PoseOP
ParsedProtocol::get_additional_output( )
{
	//fpd search the mover-filter pairs backwards; look for movers that have remaining poses
	core::pose::PoseOP pose=NULL;
	utility::vector1< mover_filter_pair >::const_reverse_iterator rmover_it;
	const utility::vector1< mover_filter_pair >::const_reverse_iterator movers_crend = movers_.rend();
	for ( rmover_it=movers_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
		core::pose::PoseOP checkpoint = (*rmover_it).first->get_additional_output();
		if (checkpoint) {
			std::string const mover_name( rmover_it->first->get_name() );
			std::string const filter_name( rmover_it->second->get_user_defined_name() );
			TR<<"=======================RESUMING FROM "<<mover_name<<"======================="<<std::endl;
			pose = checkpoint;
			TR<<"=======================BEGIN FILTER "<<filter_name<<"=======================\n{"<<std::endl;
			info().insert( info().end(), rmover_it->first->info().begin(), rmover_it->first->info().end() );
			pose->update_residue_neighbors();
			moves::MoverStatus status( (*rmover_it).first->get_last_move_status() );
			bool const pass( status==protocols::moves::MS_SUCCESS  && (*rmover_it).second->apply( *pose ) );
			TR<<"\n}\n=======================END FILTER "<<filter_name<<"======================="<<std::endl;
			if( !pass ) {
				if( status != protocols::moves::MS_SUCCESS )
					protocols::moves::Mover::set_last_move_status( status );
				return pose;
			}
			break;
		}
	}

	// no saved poses?  return now
	if (!pose) return NULL;

	// if mode_ is not 'sequence' then checkpointing is unsupported
	if (mode_ != "sequence")
		utility_exit_with_message("ParsedProtocol returned multiple poses in a ParsedProtocol with mode!=sequence");

	// otherwise pick up from the checkpoint
	for( utility::vector1< mover_filter_pair >::const_iterator mover_it = rmover_it.base();
		 mover_it!=movers_.end(); ++mover_it ) {
		std::string const mover_name( mover_it->first->get_name() );
		std::string const filter_name( mover_it->second->get_user_defined_name() );

		(*mover_it).first->set_native_pose( get_native_pose() );
		TR<<"=======================BEGIN MOVER "<<mover_name<<"=======================\n{"<<std::endl;
		(*mover_it).first->apply( *pose );
		TR<<"\n}\n=======================END MOVER "<<mover_name<<"======================="<<std::endl;
		TR<<"=======================BEGIN FILTER "<<filter_name<<"=======================\n{"<<std::endl;
		info().insert( info().end(), mover_it->first->info().begin(), mover_it->first->info().end() );
		pose->update_residue_neighbors();
		moves::MoverStatus status( (*mover_it).first->get_last_move_status() );
		bool const pass( status==protocols::moves::MS_SUCCESS  && (*mover_it).second->apply( *pose ) );
		TR<<"\n}\n=======================END FILTER "<<filter_name<<"======================="<<std::endl;
		if( !pass ) {
			if( status != protocols::moves::MS_SUCCESS )
				protocols::moves::Mover::set_last_move_status( status );
			return pose;
		}
	}
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
	TR<<"setting status to success"<<std::endl;

	// report filter values to the job object as string_real_pair
	report_filters_to_job( *pose );
	// report filter values to tracer output
	report_all( *pose );
	// rescore the pose with either score12 or a user-specified scorefunction. this ensures that all output files end up with scores.
//	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
//	(*scorefxn)(*pose);

	return pose;
}

bool ParsedProtocol::apply_mover_filter_pair(Pose & pose, mover_filter_pair const & mover_pair)
{
	std::string const mover_name( mover_pair.first->get_name() );
	std::string const filter_name( mover_pair.second->get_user_defined_name() );

	mover_pair.first->set_native_pose( get_native_pose() );
	TR<<"=======================BEGIN MOVER "<<mover_name<<"=======================\n{"<<std::endl;
	mover_pair.first->apply( pose );
	TR<<"\n}\n=======================END MOVER "<<mover_name<<"======================="<<std::endl;
	// collect Mover info: jd2 JobDistributor passes this info to Job,
	// and JobOutputters may then write this info to output files
	TR<<"=======================BEGIN FILTER "<<filter_name<<"=======================\n{"<<std::endl;
	info().insert( info().end(), mover_pair.first->info().begin(), mover_pair.first->info().end() );
	pose.update_residue_neighbors();
	moves::MoverStatus status( mover_pair.first->get_last_move_status() );
	bool const pass( status==protocols::moves::MS_SUCCESS  && mover_pair.second->apply( pose ) );
	TR<<"\n}\n=======================END FILTER "<<filter_name<<"======================="<<std::endl;
	if( !pass ) {
		if( status != protocols::moves::MS_SUCCESS )
			protocols::moves::Mover::set_last_move_status( status );
		return false;
	}
	return true;
}


void ParsedProtocol::finish_protocol(Pose & pose) {
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
	TR<<"setting status to success"<<std::endl;

	// report filter values to the job object as string_real_pair
	report_filters_to_job( pose );
	// report filter values to tracer output
	report_all( pose );
	// rescore the pose with either score12 or a user-specified scorefunction. this ensures that all output files end up with scores.
//	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
//	(*scorefxn)(pose);
}


void ParsedProtocol::sequence_protocol(Pose & pose, utility::vector1< mover_filter_pair >::const_iterator mover_it_in)
{
	//bool last_mover=false;
	for( utility::vector1< mover_filter_pair >::const_iterator mover_it = mover_it_in;
			mover_it!=movers_.end(); ++mover_it ) {

		if(!apply_mover_filter_pair(pose, *mover_it)) {
			return;
		}
	}

	// we're done! mark as success
	finish_protocol( pose );
}


void ParsedProtocol::random_order_protocol(Pose & pose){
	std::random_shuffle(movers_.begin(),movers_.end());
	for(utility::vector1<mover_filter_pair>::const_iterator it = movers_.begin(); it != movers_.end();++it)
	{
		if(!apply_mover_filter_pair(pose, *it))
		{
			return;
		}
	}
	// we're done! mark as success
	finish_protocol( pose );
}


void ParsedProtocol::random_single_protocol(Pose & pose){
	core::Size index=RG.random_range(1,movers_.size());
	if(!apply_mover_filter_pair(pose, movers_[index])) {
		return;
	}

	// we're done! mark as success
	finish_protocol( pose );
}

} //protein_interface_design
} //protocols

