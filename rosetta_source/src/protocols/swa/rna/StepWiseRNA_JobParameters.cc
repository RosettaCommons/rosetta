// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_JobParameters
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>

#include <core/pose/Pose.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <string>

#include <utility/vector1.hh>


using namespace core;
using core::Real;

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_JobParameters::StepWiseRNA_JobParameters():
		sequence_( "" ),
		working_sequence_( "" ),
		moving_res_( 0 ),
		working_moving_res_( 0 ),
		working_moving_suite_( 0 ),
		which_chain_has_moving_res_( 0 ),
		gap_size_( 0 ),
		five_prime_chain_break_res_( 0 ),
		Is_prepend_( false ),
		Is_internal_( false )
	{
	}

	StepWiseRNA_JobParameters::~StepWiseRNA_JobParameters(){}

	//////////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseRNA_JobParameters::sequence() const{
		return sequence_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseRNA_JobParameters::working_sequence() const{
		return working_sequence_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::moving_res() const{
		return moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::working_moving_res() const{
		return working_moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::working_moving_suite() const{
		return working_moving_suite_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< Size > const & StepWiseRNA_JobParameters::is_working_res() const{
		return is_working_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > & StepWiseRNA_JobParameters::full_to_sub(){
		return full_to_sub_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< std::pair< core::Size, core::Size > > const & StepWiseRNA_JobParameters::chain_boundaries() const{
		return chain_boundaries_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::which_chain_has_moving_res() const{
		return which_chain_has_moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::gap_size() const{
		return gap_size_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::five_prime_chain_break_res() const{
		return five_prime_chain_break_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseRNA_JobParameters::Is_prepend() const{
		return Is_prepend_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseRNA_JobParameters::Is_internal() const{
		return Is_internal_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< bool > const & StepWiseRNA_JobParameters::partition_definition() const{
		return partition_definition_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_fixed_res() const{
		return working_fixed_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_terminal_res() const{
		return working_terminal_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::moving_pos() const{
		return moving_pos_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_sequence( std::string const & setting ){
		sequence_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_working_sequence( std::string const & setting ){
		working_sequence_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_moving_res( Size const & setting ){
		moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_working_moving_res( Size const & setting ){
		working_moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_working_moving_suite( Size const & setting ){
		working_moving_suite_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_is_working_res( ObjexxFCL::FArray1D< Size > const & setting ){
		is_working_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_full_to_sub( std::map< core::Size, core::Size > & setting ){
		full_to_sub_ = setting;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting ){
		chain_boundaries_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_which_chain_has_moving_res( Size const & setting ){
		which_chain_has_moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_gap_size( Size const & setting ){
		gap_size_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_five_prime_chain_break_res( Size const & setting ){
		five_prime_chain_break_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_Is_prepend( bool const & setting ){
		Is_prepend_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_Is_internal( bool const & setting ){
		Is_internal_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_partition_definition( ObjexxFCL::FArray1D< bool > const & setting ){
		partition_definition_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	core::pose::PoseCOP
	StepWiseRNA_JobParameters::working_native_pose() const{
		return working_native_pose_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_native_pose( core::pose::PoseOP & pose ){
		working_native_pose_ = pose;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseRNA_JobParameters::actually_moving_res() const{
		return ( Is_prepend_ ?
						 working_moving_suite_ :
						 working_moving_suite_+1 );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res){
		working_fixed_res_ = working_fixed_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res){
		working_terminal_res_ = working_terminal_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_moving_pos(	utility::vector1< core::Size > const & moving_pos){
		moving_pos_ = moving_pos;
	}



}
}
}
