// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/StepWiseBasicWorkingParameters.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/working_parameters/StepWiseBasicWorkingParameters.hh>
#include <core/pose/Pose.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.working_parameters.StepWiseBasicWorkingParameters" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace working_parameters {

//Constructor
StepWiseBasicWorkingParameters::StepWiseBasicWorkingParameters():
	full_sequence_( "" ),
	working_sequence_( "" ),
	moving_res_( 0 ),
	working_moving_res_( 0 ),
	working_moving_suite_( 0 ),
	is_prepend_( false ),
	is_internal_( false ),
	gap_size_( GAP_SIZE_DUMMY )
{}

//Destructor
StepWiseBasicWorkingParameters::~StepWiseBasicWorkingParameters()
{}

//////////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseBasicWorkingParameters::actually_moving_res() const{
	return ( is_prepend_ ?
		working_moving_suite_list_[1] :
		working_moving_suite_list_[1]+1 );
}

//////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
StepWiseBasicWorkingParameters::working_res_list() const{
	if ( sub_to_full_.size() == 0 ) utility_exit_with_message( "sub_to_full_.size() == 0. Cannot output working_res_list" );
	utility::vector1< Size > working_res_list;
	for ( std::map< core::Size, core::Size > ::const_iterator it = sub_to_full_.begin(), end = sub_to_full_.end(); it != end; ++it ) {
		working_res_list.push_back( it->second );
	}
	return working_res_list;
}

//////////////////////////////////////////////////////////////////////////////////////////
void StepWiseBasicWorkingParameters::set_working_moving_res_list( utility::vector1< Size > const & setting ){
	working_moving_res_list_ = setting;
	working_moving_res_ = ( working_moving_res_list_.size() > 0 ) ? working_moving_res_list_[1] : 0;
	update_working_moving_suite(); //Be careful..currently when working_moving_res_list is first initialize in WP_Setup, is_Prepend is not setup yet.
}
//////////////////////////////////////////////////////////////////////////////////////////
void StepWiseBasicWorkingParameters::set_full_to_sub( std::map< core::Size, core::Size > const & setting ){
	full_to_sub_ = setting;
	sub_to_full_ = create_sub_to_full_map( full_to_sub_ ); //Parin Jan 18, 2009
}

///////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
StepWiseBasicWorkingParameters::apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector){

	utility::vector1< core::Size > working_res_vector;
	for ( Size n = 1; n <= res_vector.size(); n++ ) {
		if ( !is_working_res_[ res_vector[ n ] ] ) continue;

		working_res_vector.push_back( full_to_sub_[ res_vector[ n ] ]);
	}

	return working_res_vector;
}

///////////////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseBasicWorkingParameters::apply_full_to_sub_mapping( Size const res ) const {
	std::map<Size,Size>::const_iterator iter = full_to_sub_.find( res );
	if ( iter == full_to_sub_.end() ) return 0;
	return iter->second;
}

//////////////////////////////////////////////////////////////////////////////////////////
void StepWiseBasicWorkingParameters::set_partition_definition( utility::vector1< Size > const & partition_definition_vector ){
	partition_definition_.dimension( partition_definition_vector.size() );
	for ( Size n = 1; n <= partition_definition_vector.size(); n++ ) {
		partition_definition_( n ) = partition_definition_vector[ n ];
	}
}
//////////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseCOP
StepWiseBasicWorkingParameters::working_native_pose() const{
	return working_native_pose_;
}

//////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicWorkingParameters::set_working_native_pose( core::pose::PoseOP & pose ){
	working_native_pose_ = pose;
}

//////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicWorkingParameters::set_working_native_pose( core::pose::PoseCOP pose ){
	working_native_pose_ = pose;
}

//////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseBasicWorkingParameters::update_working_moving_suite(){

	runtime_assert( working_moving_res_list_.size() > 0  || working_moving_res_ == 0 );
	runtime_assert( working_moving_res_list_.size() == 0 || working_moving_res_ != 0 );
	if ( working_moving_res_list_.size() == 0 ) return;
	runtime_assert( working_moving_res_ > 0 );

	working_moving_suite_list_.clear();
	if ( is_prepend_ ) {
		working_moving_suite_ = working_moving_res_;
		working_moving_suite_list_ = working_moving_res_list_;
	} else {
		working_moving_suite_ = working_moving_res_ - 1;
		for ( Size n = 1; n <= working_moving_res_list_.size(); n++ ) {
			working_moving_suite_list_.push_back( working_moving_res_list_[n] - 1 );
		}
	}

	//check
	if ( working_moving_suite_ != working_moving_suite_list_[1] ) utility_exit_with_message( "working_moving_suite_ != working_moving_suite_list_[1]" );
}

//////////////////////////////////////////////////////////////////////////////////////////
std::map< core::Size, core::Size >  //Parin Jan 18, 2009
StepWiseBasicWorkingParameters::create_sub_to_full_map( std::map< core::Size, core::Size > const & full_to_sub ) const{
	std::map< core::Size, core::Size > ::const_iterator it, end ;
	std::map< core::Size, core::Size > sub_to_full;
	sub_to_full.clear();
	for ( it = full_to_sub.begin(), end = full_to_sub.end(); it != end; ++it ) {
		sub_to_full[it->second] = it->first;
	}
	return sub_to_full;
}

//////////////////////////////////////////////////////////////////////////////////////////
void StepWiseBasicWorkingParameters::update_working_sequence(){

	if ( full_sequence_.size() == 0 ) return;
	if ( is_working_res_.size() == 0 ) return;

	working_sequence_ = "";

	for ( Size full_seq_num = 1, string_index = 0; string_index < full_sequence_.size(); ++full_seq_num, ++string_index ) {
		if ( !is_working_res_[ full_seq_num ] ) continue;

		if ( full_sequence_[ string_index ] == 'X' ) {
			working_sequence_ += full_sequence_[ string_index++ ];
			working_sequence_ += full_sequence_[ string_index++ ];
			working_sequence_ += full_sequence_[ string_index++ ];
			working_sequence_ += full_sequence_[ string_index++ ];
			working_sequence_ += full_sequence_[ string_index++ ];
		}
		working_sequence_ += full_sequence_[ string_index ]; //i-1 because std::string elements starts at 0...
	}
}


} //working_parameters
} //modeler
} //stepwise
} //protocols
