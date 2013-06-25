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
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using core::Real;

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_JobParameters::StepWiseRNA_JobParameters():
		output_extra_RMSDs_(false),
		Is_simple_full_length_job_params_(false),
		full_sequence_( "" ),
		working_sequence_( "" ),
		moving_res_( 0 ),
		working_moving_res_( 0 ),
		working_moving_suite_( 0 ),
		//which_chain_has_moving_res_( 0 ),
		gap_size_( 0 ),
		five_prime_chain_break_res_( 0 ),
		Is_prepend_( false ),
		Is_internal_( false ),
		add_virt_res_as_root_( false )
	{

		//These vectors and map should be empty to begin with, but not harm to ensure this.
		is_working_res_.clear();
		working_moving_res_list_.clear();
		working_moving_suite_list_.clear();
		full_to_sub_.clear();
		sub_to_full_.clear();
		Is_prepend_map_.clear();
		chain_boundaries_.clear();
		partition_definition_.clear();
		working_fixed_res_.clear();
		rmsd_res_list_.clear();
		working_terminal_res_.clear();
		working_moving_partition_pos_.clear();
		fold_tree_.clear();
		input_res_vectors_.clear();
		cutpoint_closed_list_.clear();
		working_best_alignment_.clear();
		native_alignment_.clear();
		working_native_alignment_.clear();
		global_sample_res_list_.clear();
		working_global_sample_res_list_.clear();
		force_syn_chi_res_list_.clear();
		working_force_syn_chi_res_list_.clear();
		force_north_ribose_list_.clear();
		working_force_north_ribose_list_.clear();
		force_south_ribose_list_.clear();
		working_force_south_ribose_list_.clear();
		protonated_H1_adenosine_list_.clear();
		working_protonated_H1_adenosine_list_.clear();
	}

	StepWiseRNA_JobParameters::~StepWiseRNA_JobParameters(){}

	//////////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseRNA_JobParameters::output_extra_RMSDs() const{
		return output_extra_RMSDs_;
	}
	///////////////////////////////////////////////////////////////////////////////////////
	bool const & StepWiseRNA_JobParameters::Is_simple_full_length_job_params() const{
		return Is_simple_full_length_job_params_;
	}
	///////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseRNA_JobParameters::full_sequence() const{

		if(full_sequence_.size()==0) utility_exit_with_message( "full_sequence_.size()==0" );
		return full_sequence_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::string const & StepWiseRNA_JobParameters::working_sequence() const{

		if(working_sequence_.size()==0) utility_exit_with_message( "working_sequence_.size()==0" );
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
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_moving_res_list() const{
		if(working_moving_res_list_.size()==0) utility_exit_with_message("working_moving_res_list_.size()==0");
		return working_moving_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size StepWiseRNA_JobParameters::working_reference_res() const{ //the last static_residues that this attach to the moving residues

		Size const num_nucleotides=working_moving_res_list_.size();

		//check that moving_res_ and working_moving_res_list list are intialized (note this is not foolproof)
		if(num_nucleotides==0) utility_exit_with_message( "num_building_nucleotides==0!!" );
		if(working_moving_res_==0) utility_exit_with_message( "working_moving_res_==0!!" );


		Size const working_reference_res_ = (Is_prepend_) ? working_moving_res_ + num_nucleotides : working_moving_res_ - num_nucleotides;

		return working_reference_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size const & StepWiseRNA_JobParameters::working_moving_suite() const{
		return working_moving_suite_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_moving_suite_list() const{
		if(working_moving_suite_list_.size()==0) utility_exit_with_message("working_moving_suite_list_.size()==0");
		return working_moving_suite_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::is_working_res() const{
		if(is_working_res_.size()==0) utility_exit_with_message("is_working_res_.size()==0");
		return is_working_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > & StepWiseRNA_JobParameters::full_to_sub(){
		if(full_to_sub_.size()==0) utility_exit_with_message("full_to_sub_.size()==0");
		return full_to_sub_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > & StepWiseRNA_JobParameters::sub_to_full(){
		if(sub_to_full_.size()==0) utility_exit_with_message("sub_to_full_.size()==0");
		return sub_to_full_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > const & StepWiseRNA_JobParameters::const_full_to_sub() const{
		if(full_to_sub_.size()==0) utility_exit_with_message("full_to_sub_.size()==0");
		return full_to_sub_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size > const & StepWiseRNA_JobParameters::const_sub_to_full() const{
		if(sub_to_full_.size()==0) utility_exit_with_message("sub_to_full_.size()==0");
		return sub_to_full_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	core::kinematics::FoldTree const & StepWiseRNA_JobParameters::fold_tree() const{
		if(fold_tree_.size()==0) utility_exit_with_message("fold_tree_.size()==0"); //Thie is the number of edge. simple_tree have 1 edge!
		return fold_tree_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, bool > const & StepWiseRNA_JobParameters::Is_prepend_map() const {
		if(Is_prepend_map_.size()==0) utility_exit_with_message("Is_prepend_map_.size()==0");
		return Is_prepend_map_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< std::pair< core::Size, core::Size > > const & StepWiseRNA_JobParameters::chain_boundaries() const{
		return chain_boundaries_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	//Size const & StepWiseRNA_JobParameters::which_chain_has_moving_res() const{
	//	return which_chain_has_moving_res_;
	//}
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
	bool const & StepWiseRNA_JobParameters::add_virt_res_as_root() const{
		return add_virt_res_as_root_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D< bool > const & StepWiseRNA_JobParameters::partition_definition() const{
		if(partition_definition_.size()==0) utility_exit_with_message("partition_definition_.size()==0!");
		return partition_definition_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_fixed_res() const{
		return working_fixed_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::rmsd_res_list() const{
		if(rmsd_res_list_.size()==0) utility_exit_with_message("rmsd_res_list_.size()==0");
		return rmsd_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_terminal_res() const{
		return working_terminal_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_moving_partition_pos() const{
		return working_moving_partition_pos_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > > const & StepWiseRNA_JobParameters::input_res_vectors() const {
		if(input_res_vectors_.size()==0) utility_exit_with_message("input_res_vectors_.size()==0");
		return input_res_vectors_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const & StepWiseRNA_JobParameters::cutpoint_closed_list() const {
		return cutpoint_closed_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_best_alignment() const{

		if(working_best_alignment_.size()==0) utility_exit_with_message("working_best_alignment_.size()==0!");

		return working_best_alignment_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::native_alignment() const{
		return native_alignment_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_native_alignment() const{
		return working_native_alignment_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::global_sample_res_list() const{
		return global_sample_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_global_sample_res_list() const{
		return working_global_sample_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::force_syn_chi_res_list() const{
		return force_syn_chi_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_force_syn_chi_res_list() const{
		return working_force_syn_chi_res_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::force_north_ribose_list() const{
		return force_north_ribose_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_force_north_ribose_list() const{
		return working_force_north_ribose_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::force_south_ribose_list() const{
		return force_south_ribose_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_force_south_ribose_list() const{
		return working_force_south_ribose_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::protonated_H1_adenosine_list() const{
		return protonated_H1_adenosine_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const & StepWiseRNA_JobParameters::working_protonated_H1_adenosine_list() const{
		return working_protonated_H1_adenosine_list_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_output_extra_RMSDs( bool const & setting ){
		output_extra_RMSDs_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_Is_simple_full_length_job_params( bool const & setting ){ //Oct 31, 2011
		Is_simple_full_length_job_params_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_full_sequence( std::string const & setting ){
		full_sequence_ = setting;

		update_working_sequence();

		std::cout << "full_sequence_= " << full_sequence_ << std::endl;
		std::cout << "working_sequence_= " <<  working_sequence_ << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_moving_res( Size const & setting ){
		moving_res_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_working_moving_res_list( utility::vector1< Size > const & setting ){
		working_moving_res_list_ = setting;
		working_moving_res_ = working_moving_res_list_[1];

		update_working_moving_suite(); //Be careful..currently when working_moving_res_list is first initialize in JP_Setup, Is_Prepend is not setup yet.
	}



	//////////////////////////////////////////////////////////////////////////////////////////
	//void StepWiseRNA_JobParameters::set_working_moving_res( Size const & setting ){
	//	working_moving_res_ = setting;
	//}
	//////////////////////////////////////////////////////////////////////////////////////////
	//void StepWiseRNA_JobParameters::set_working_moving_suite_list( utility::vector1< Size > const & setting ){
	//	working_moving_suite_list_ = setting;
	//}
	////////////////////////////////////////////////////////////////////////////////////////////
	//void StepWiseRNA_JobParameters::set_working_moving_suite( Size const & setting ){
	//	working_moving_suite_ = setting;
	//}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_is_working_res( utility::vector1< core::Size > const & setting ){
		is_working_res_ = setting;

		update_working_sequence();
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_full_to_sub( std::map< core::Size, core::Size > const & setting ){
		full_to_sub_ = setting;
		sub_to_full_=create_sub_to_full_map(full_to_sub_); //Parin Jan 18, 2009
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_fold_tree(core::kinematics::FoldTree const & setting){
		fold_tree_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_Is_prepend_map( std::map< core::Size, bool > const & setting ){
		Is_prepend_map_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseRNA_JobParameters::set_chain_boundaries( utility::vector1< std::pair< core::Size, core::Size > > const & setting ){
		chain_boundaries_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	//void StepWiseRNA_JobParameters::set_which_chain_has_moving_res( Size const & setting ){
	//	which_chain_has_moving_res_ = setting;
	//}
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
		update_working_moving_suite();
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
	void
	StepWiseRNA_JobParameters::set_working_native_pose( core::pose::PoseCOP pose ){
		working_native_pose_ = pose;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseRNA_JobParameters::actually_moving_res() const{
			return working_moving_res_;

	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_fixed_res(	utility::vector1< core::Size > const & working_fixed_res){
		working_fixed_res_ = working_fixed_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_rmsd_res_list(	utility::vector1< core::Size > const & rmsd_res_list){
		rmsd_res_list_ = rmsd_res_list;
		sort_seq_num_list(rmsd_res_list_);
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_terminal_res(	utility::vector1< core::Size > const & working_terminal_res){
		working_terminal_res_ = working_terminal_res;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_moving_partition_pos(	utility::vector1< core::Size > const & working_moving_partition_pos){
		working_moving_partition_pos_ = working_moving_partition_pos;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_input_res_vectors(	utility::vector1< utility::vector1< Size > > const & setting){
		input_res_vectors_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_cutpoint_closed_list( utility::vector1< Size >  const & setting){
		cutpoint_closed_list_= setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_best_alignment( utility::vector1< core::Size > const & setting){
		working_best_alignment_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_native_alignment( utility::vector1< core::Size > const & setting){
		native_alignment_=setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_working_native_alignment( utility::vector1< core::Size > const & setting){
		working_native_alignment_=setting;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_global_sample_res_list( utility::vector1< core::Size > const & setting){
		global_sample_res_list_=setting;
		working_global_sample_res_list_=apply_full_to_sub_mapping( global_sample_res_list_, is_working_res_, full_to_sub_ );

	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting){
		force_syn_chi_res_list_=setting;
		working_force_syn_chi_res_list_=apply_full_to_sub_mapping( force_syn_chi_res_list_, is_working_res_, full_to_sub_ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_force_north_ribose_list( utility::vector1< core::Size > const & setting){
		force_north_ribose_list_=setting;
		working_force_north_ribose_list_=apply_full_to_sub_mapping( force_north_ribose_list_, is_working_res_, full_to_sub_ );

		for(Size n=1; n<=force_north_ribose_list_.size(); n++){
			if(Contain_seq_num(force_north_ribose_list_[n], force_south_ribose_list_)){
				utility_exit_with_message("seq_num= " + ObjexxFCL::string_of(force_north_ribose_list_[n]) + " is in both force_north_ribose_list_ and force_south_ribose_list_! " );
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_force_south_ribose_list( utility::vector1< core::Size > const & setting){
		force_south_ribose_list_=setting;
		working_force_south_ribose_list_=apply_full_to_sub_mapping( force_south_ribose_list_, is_working_res_, full_to_sub_ );

		for(Size n=1; n<=force_north_ribose_list_.size(); n++){
			if(Contain_seq_num(force_north_ribose_list_[n], force_south_ribose_list_)){
				utility_exit_with_message("seq_num= " + ObjexxFCL::string_of(force_north_ribose_list_[n]) + " is in both force_north_ribose_list_ and force_south_ribose_list_! " );
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::set_protonated_H1_adenosine_list( utility::vector1< core::Size > const & setting){
		protonated_H1_adenosine_list_=setting;
		working_protonated_H1_adenosine_list_=apply_full_to_sub_mapping( protonated_H1_adenosine_list_, is_working_res_, full_to_sub_ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParameters::update_working_moving_suite(){

		//std::cout << "update_working_moving_suite " << std::endl;
		if( (working_moving_res_list_.size()==0) && (working_moving_res_!=0) ){
			utility_exit_with_message( "working_moving_res_list_.size()==0) && (working_moving_res_!=0");
		}

		if( (working_moving_res_list_.size()!=0) && (working_moving_res_==0) ){
			utility_exit_with_message( "working_moving_res_list_.size()!=0) && (working_moving_res_==0");
		}

		if(working_moving_res_list_.size()==0) return;

		if(working_moving_res_==0) utility_exit_with_message( "working_moving_res_list_.size()!=0) && (working_moving_res_==0");



		working_moving_suite_list_.clear();

		if(Is_prepend_ ){
			working_moving_suite_ = working_moving_res_;
			working_moving_suite_list_=working_moving_res_list_;
		}else{
			working_moving_suite_ = working_moving_res_ - 1;

			for(Size n=1; n<=working_moving_res_list_.size(); n++){
				working_moving_suite_list_.push_back(working_moving_res_list_[n] - 1);
			}

		}

		//check
		if(working_moving_suite_!=working_moving_suite_list_[1]) utility_exit_with_message( "working_moving_suite_!=working_moving_suite_list_[1]");

	}

	//////////////////////////////////////////////////////////////////////////////////////////
	std::map< core::Size, core::Size >  //Parin Jan 18, 2009
	StepWiseRNA_JobParameters::create_sub_to_full_map(std::map< core::Size, core::Size > const & full_to_sub) const{

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

	void StepWiseRNA_JobParameters::update_working_sequence(){

		if(full_sequence_.size()==0) return;
		if(is_working_res_.size()==0) return;

		working_sequence_="";

		for ( Size full_seq_num = 1; full_seq_num <= full_sequence_.size(); full_seq_num++ ) {

			if ( is_working_res_[ full_seq_num ] ) {
				working_sequence_ += full_sequence_[ full_seq_num-1 ]; //i-1 because std::string elements starts at 0...
			}

		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////


}
}
}
