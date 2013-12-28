// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseJobParameters
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseJobParameters.hh>
#include <core/pose/Pose.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;
using core::Real;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseJobParameters::StepWiseJobParameters():
		sequence_( "" ),
		working_sequence_( "" ),
		which_chain_has_moving_res_( 0 ),
		gap_size_( 0 ),
		first_chain_break_res_( 0 ),
		is_prepend_( false ),
		is_internal_( false )
	{
	}

	StepWiseJobParameters::~StepWiseJobParameters(){}

	//////////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseJobParameters::sequence() const{
		return sequence_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseJobParameters::working_sequence() const{
		return working_sequence_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_res_list() const{
		return working_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_moving_res_list() const{
		return working_moving_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_bridge_res() const{
		return working_bridge_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_moving_suite_list() const{
		return working_moving_suite_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< Size > const & StepWiseJobParameters::is_working_res() const{
		return is_working_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< Size > const & StepWiseJobParameters::is_moving_res() const{
		return is_moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > & StepWiseJobParameters::full_to_sub(){
		return full_to_sub_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > & StepWiseJobParameters::sub_to_full(){
		return sub_to_full_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< std::pair< core::Size, core::Size > > const & StepWiseJobParameters::chain_boundaries() const{
		return chain_boundaries_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseJobParameters::which_chain_has_moving_res() const{
		return which_chain_has_moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseJobParameters::gap_size() const{
		return gap_size_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseJobParameters::first_chain_break_res() const{
		return first_chain_break_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseJobParameters::is_prepend() const{
		return is_prepend_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseJobParameters::is_internal() const{
		return is_internal_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< bool > const & StepWiseJobParameters::partition_definition() const{
		return partition_definition_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_fixed_res() const{
		return working_fixed_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_terminal_res() const{
		return working_terminal_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_superimpose_res() const{
		return working_superimpose_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::working_calc_rms_res() const{
		return working_calc_rms_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseJobParameters::moving_pos() const{
		return moving_pos_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_sequence( std::string const & setting ){
		sequence_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_working_sequence( std::string const & setting ){
		working_sequence_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_working_res_list( utility::vector1< Size > const & setting ){
		working_res_list_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_working_moving_res_list( utility::vector1< Size > const & setting ){
		working_moving_res_list_ = setting; //Could just store moving_res internally and full_to_sub map + moving_res to calculate and return working_moving_res when working_moving_res() is called
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_working_moving_suite_list( utility::vector1< Size > const & setting ){
		working_moving_suite_list_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_is_working_res( ObjexxFCL::FArray1D< Size > const & setting ){
		is_working_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_is_moving_res( ObjexxFCL::FArray1D< Size > const & setting ){
		is_moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_full_to_sub( std::map< core::Size, core::Size > const & setting ){
		full_to_sub_ = setting;
		sub_to_full_=create_sub_to_full_map(full_to_sub_); //Parin Jan 18, 2009
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size >  //Parin Jan 18, 2009
	StepWiseJobParameters::create_sub_to_full_map(std::map< core::Size, core::Size > const & full_to_sub) const{

		std::map< core::Size, core::Size > ::const_iterator it;
		std::map< core::Size, core::Size > sub_to_full;
		sub_to_full.clear();

		for (it=full_to_sub.begin(); it!=full_to_sub.end(); it++ ){
			sub_to_full[it->second]=it->first;
		}

		//output for debug
		//		std::cout << "full_to_sub" << std::endl;
		for (it=full_to_sub.begin(); it!=full_to_sub.end(); it++ ){
			//			std::cout << it->first << " => " << it->second << std::endl;
		}

		//		std::cout << "sub_to_full" << std::endl;
		for (it=sub_to_full.begin(); it!=sub_to_full.end(); it++ ){
			//		std::cout << it->first << " => " << it->second << std::endl;
		}

		return sub_to_full;

	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< bool > const
	StepWiseJobParameters::is_pre_proline() const {

		utility::vector1< bool > is_pre_proline;

		for ( Size i = 1; i <= working_res_list_.size(); i++ ){

			Size const & full_seq_pos = working_res_list_[ i ];

			if ( full_seq_pos == sequence_.size() ) {
				is_pre_proline.push_back( false );
			} else {
				is_pre_proline.push_back( sequence_[ full_seq_pos ] == 'P' ); /*note offset by one -- this is the *next* sequence position*/
			}
		}

		return is_pre_proline;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting ){
		chain_boundaries_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_which_chain_has_moving_res( Size const & setting ){
		which_chain_has_moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_gap_size( Size const & setting ){
		gap_size_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_first_chain_break_res( Size const & setting ){
		first_chain_break_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_is_prepend( bool const & setting ){
		is_prepend_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_is_internal( bool const & setting ){
		is_internal_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_partition_definition( ObjexxFCL::FArray1D< bool > const & setting ){
		partition_definition_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	core::pose::PoseCOP
	StepWiseJobParameters::working_native_pose() const{
		return working_native_pose_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_working_native_pose( core::pose::PoseOP & pose ){
		working_native_pose_ = pose;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseJobParameters::actually_moving_res() const{ //Should rewrite code so that there is moving_suite and moving_res......having actually_moving_res is confusing..Parin Jan 30, 2010
		return ( is_prepend_ ?
						 working_moving_suite_list_[1] :
						 working_moving_suite_list_[1]+1 );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res){
		working_fixed_res_ = working_fixed_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res){
		working_terminal_res_ = working_terminal_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_working_superimpose_res(	utility::vector1< core::Size > const & working_superimpose_res){
		working_superimpose_res_ = working_superimpose_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_working_calc_rms_res(	utility::vector1< core::Size > const & working_calc_rms_res){
		working_calc_rms_res_ = working_calc_rms_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseJobParameters::set_working_bridge_res( utility::vector1< Size > const & setting ){
		working_bridge_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseJobParameters::set_moving_pos(	utility::vector1< core::Size > const & moving_pos){
		moving_pos_ = moving_pos;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	StepWiseJobParameters::apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector){

		utility::vector1< core::Size > working_res_vector;
		for ( Size n = 1; n <= res_vector.size(); n++ ) {
			if ( !is_working_res_( res_vector[ n ] ) ) continue;
			working_res_vector.push_back( full_to_sub_[ res_vector[ n ] ]);
		}

		return working_res_vector;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseJobParameters::apply_full_to_sub_mapping( Size const res ) const {

		std::map<Size,Size>::const_iterator iter = full_to_sub_.find( res );
		if ( iter == full_to_sub_.end() ) return 0;
		return iter->second;

	}


} //protein
} //enumerate
} //stepwise
} //protocols

