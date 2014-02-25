// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiator.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiator.hh>
#include <protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarSampler.hh>
#include <protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/RNA_Util.hh> // for SYN
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

///////////////////////////////////////////////////////////////////////////////////
//////////////////////Build previously virtualized sugar/////////////////////
// A virtualized sugar occurs when a previous move was a 'floating base'
// step, which only samples euler angles of a base, but virtualized the
// attached sugar and any residues connecting that nucleotide to the
//'instantiated' body of the RNA.
//
// There are potentially four different virtualized sugar positions,
//  and here is a diagram of the most extreme case, merging of
//  two helices whose edge base pairs were all somehow created through
//  bulge skipping moves, leaving 4 virtual riboses (marked *):
//
//             A1 -- U16                 |
//             C2 -- G15                 | POSE 1
// bulge -> (U3)       (U14) <- bulge    |
//            *G4 -- C13*                |
//    moving-> |      x <- chain break to close
//            *G5 -- C12*                |
// bulge -> (U6)       (U11) <- bulge    | POSE2
//             C7 -- G10                 |
//             C8 -- G9                  |
//
// (1) Anchor_sugar (most common case) :
//      This is the sugar of the nucleotide that was built immediately
//      before the current/moving nucleotide along the chain. Corresponds
//      to sugar of residue (moving_res - 1) if current step is built in
//      the forward (3') direction. Likewise, corresponds to sugar of
//      residue (moving_res + 1) if current step built in the backward
//      (5') direction.
//
// (2) Current_sugar (rare) :
//      The sugar of the current/moving nucleotide. This sugar can
//      be virtual in the situation where the current step is a step to
//      combine two moving_elements that were previously built with SWA. It is
//      possible that the current nucleotide (in the moving moving_element) was
//      previously built with a 'floating base' step and therefore will
//      contain a virtual sugar.
//
// (3) Five_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 5' of the chain-closure phosphate group.
//      Five_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.
//
// (4) Three_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 3' of the chain-closure phosphate group.
//      Three_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_VirtualSugarJustInTimeInstantiator" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace sugar {

//Constructor
StepWiseRNA_VirtualSugarJustInTimeInstantiator::StepWiseRNA_VirtualSugarJustInTimeInstantiator( StepWiseRNA_JobParametersCOP & job_parameters ):
	job_parameters_( job_parameters ),
	is_prepend_(  job_parameters_->is_prepend() ),
	moving_res_(  job_parameters_->working_moving_res() ),
	five_prime_chain_break_res_( job_parameters_->five_prime_chain_break_res() ),
	three_prime_chain_break_res_( five_prime_chain_break_res_ + 1 ),
	num_nucleotides_(  job_parameters_->working_moving_res_list().size() ),
	gap_size_( job_parameters_->gap_size() ),
	rebuild_bulge_mode_( job_parameters_->rebuild_bulge_mode() ),
	is_anchor_sugar_virt_( false ), // will be updated below, based on pose
	is_current_sugar_virt_( false ), // will be updated below, based on pose
	is_five_prime_chain_break_sugar_virt_( false ), // will be updated below, based on pose
	is_three_prime_chain_break_sugar_virt_( false ), // will be updated below, based on pose
	num_virtual_sugar_( 0 ),  // will be updated below, based on pose
	keep_base_fixed_( false )
{}

//Destructor
StepWiseRNA_VirtualSugarJustInTimeInstantiator::~StepWiseRNA_VirtualSugarJustInTimeInstantiator()
{}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::apply( core::pose::Pose & pose ){
	success_ = do_the_sampling( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::do_the_sampling( core::pose::Pose & pose ){

	output_title_text( "Build previously virtualize sugar", TR.Debug );
	if ( rebuild_bulge_mode_ ) {
		TR << TR.Red << "In: REBUILD_BULGE_MODE " << TR.Reset << std::endl;
		apply_virtual_rna_residue_variant_type( pose,  moving_res_ );
	}

	// infer anchor/reference res for each virtual sugar based on pose fold tree and variants.
	reference_res_for_each_virtual_sugar_ = get_reference_res_for_each_virtual_sugar_based_on_fold_tree( pose );

	if ( !initialize_parameters( pose ) ) return false;

	if ( is_anchor_sugar_virt_                  && !do_sugar_sampling( pose, anchor_sugar_modeling_, "anchor" ) ) return false;
	if ( is_current_sugar_virt_                 && !do_sugar_sampling( pose, current_sugar_modeling_, "current" ) ) return false;
	if ( is_five_prime_chain_break_sugar_virt_  && !do_sugar_sampling( pose, five_prime_chain_break_sugar_modeling_, "five_prime_CB" ) ) return false;
	if ( is_three_prime_chain_break_sugar_virt_ && !do_sugar_sampling( pose, three_prime_chain_break_sugar_modeling_, "three_prime_CB" ) ) return false;

	/////////////Sort the pose_data_list by score..should be determined according to torsional potential score.	/////////////////////////
	std::sort( anchor_sugar_modeling_.pose_list.begin(),                  anchor_sugar_modeling_.pose_list.end(),         sort_pose_by_score );
	std::sort( current_sugar_modeling_.pose_list.begin(),                 current_sugar_modeling_.pose_list.end(),           sort_pose_by_score );
	std::sort( five_prime_chain_break_sugar_modeling_.pose_list.begin(),  five_prime_chain_break_sugar_modeling_.pose_list.end(),  sort_pose_by_score );
	std::sort( three_prime_chain_break_sugar_modeling_.pose_list.begin(), three_prime_chain_break_sugar_modeling_.pose_list.end(), sort_pose_by_score );
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::do_sugar_sampling( pose::Pose & viewer_pose, SugarModeling & sugar_modeling, std::string const name ){

	StepWiseRNA_VirtualSugarSampler virtual_sugar_sampler( job_parameters_, sugar_modeling );
	virtual_sugar_sampler.set_tag( name );
	virtual_sugar_sampler.set_scorefxn( scorefxn_ );
	virtual_sugar_sampler.set_integration_test_mode( options_->integration_test_mode() );
	virtual_sugar_sampler.set_use_phenix_geo( options_->use_phenix_geo() );
	virtual_sugar_sampler.set_legacy_mode( options_->virtual_sugar_legacy_mode() );
	virtual_sugar_sampler.set_keep_base_fixed( options_->virtual_sugar_keep_base_fixed() );
	virtual_sugar_sampler.set_choose_random( options_->choose_random() );

	virtual_sugar_sampler.apply( viewer_pose );

	if ( sugar_modeling.pose_list.size() == 0 ){
		TR.Debug << name + " sugar modeling unsuccessful!" << std::endl;
		return false;
	}
	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_sugar_virtual( pose::Pose const & pose, Size const & n ){
	return ( pose.residue( n ).has_variant_type( "VIRTUAL_RIBOSE" ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::setup_sugar_modeling( pose::Pose const & pose, Size const moving_res, SugarModeling & sugar_modeling ){

	using namespace core::chemical::rna;

	if ( !is_sugar_virtual( pose, moving_res ) ) return false;

	sugar_modeling = SugarModeling( moving_res, reference_res_for_each_virtual_sugar_[ moving_res ] );
	if ( job_parameters_->working_force_syn_chi_res_list().has_value( moving_res ) ) sugar_modeling.moving_res_base_state = SYN;

	// model bulge?
	// this is assumed to be the residue immediately adjacent to the moving_residue,
	// in the same direction as moving_res. This would be the case in Parin's
	// usual dinucleotide move.  But we want to be able to totally leave out that
	// filler base, in which case there's no bulge to rebuild. -- rhiju
	Size const bulge_res_ = sugar_modeling.bulge_res;
	if ( bulge_res_ == sugar_modeling.reference_res ||
			 bulge_res_ == 0 ||
			 !pose.residue_type( bulge_res_ ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
		TR.Debug <<  "CHECKING BULGE AT " << bulge_res_ << " is it bulge? " <<  pose.residue_type( bulge_res_ ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) << "  moving_res: " << sugar_modeling.moving_res << "  reference_res " << sugar_modeling.reference_res << std::endl;
		sugar_modeling.bulge_res   = 0;
		sugar_modeling.bulge_suite = 0;
	}

	return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::is_anchor_sugar_virtual( core::pose::Pose const & pose ) {
	//Check if anchor sugar is virtual, if virtual then need to sample it.
	Size const anchor_moving_res = job_parameters_->working_reference_res();
	return setup_sugar_modeling( pose, anchor_moving_res, anchor_sugar_modeling_ );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////New June 12, 2011/////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::is_current_sugar_virtual( core::pose::Pose const & pose ) {
	// Check if curr sugar is virtual. If virtual then need to sample it.
	// This occur when combining two moving_element and the moving_res
	// in the moving_element was built with a dinucleotide move.
	Size const virtual_sugar_res = moving_res_;
	return setup_sugar_modeling( pose, virtual_sugar_res, current_sugar_modeling_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::is_five_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) {
	if ( gap_size_ != 0 ) return false;
	// Make sure to not over count number of virtual_sugar-virtual_bulge
	// pairs to be build!
	if ( moving_res_ == five_prime_chain_break_res_ ) return false;
	return setup_sugar_modeling( pose, five_prime_chain_break_res_, five_prime_chain_break_sugar_modeling_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::is_three_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) {
	if ( gap_size_ != 0 ) return false;
	// Make sure to not over count number of virtual_sugar-virtual_bulge
	// pairs to be build!
	if ( moving_res_ == three_prime_chain_break_res_ ) return false;
	return setup_sugar_modeling( pose, three_prime_chain_break_res_, three_prime_chain_break_sugar_modeling_ );
}

////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::sampling_sugar() const{
	return ( anchor_sugar_modeling_.sample_sugar ||
					 current_sugar_modeling_.sample_sugar ||
					 five_prime_chain_break_sugar_modeling_.sample_sugar ||
					 three_prime_chain_break_sugar_modeling_.sample_sugar );
}

////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::sampling_sugar_at_chain_break() const{
	return ( five_prime_chain_break_sugar_modeling_.sample_sugar ||
					 three_prime_chain_break_sugar_modeling_.sample_sugar );
}


////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
																																							 utility::vector1< PoseOP > & pose_data_list ) {

	if ( !sampling_sugar() ) utility_exit_with_message( "pose_data_list is empty for all 4 possible virtual sugars!" );

	for ( Size anchor_sugar_ID = 1;
				anchor_sugar_ID <= anchor_sugar_modeling_.pose_list.size() || anchor_sugar_ID == 1;
				anchor_sugar_ID++ ){
		for ( Size current_sugar_ID = 1;
					current_sugar_ID <= current_sugar_modeling_.pose_list.size() || current_sugar_ID == 1;
					current_sugar_ID++ ){
			for ( Size five_prime_chain_break_sugar_ID = 1;
						five_prime_chain_break_sugar_ID <= five_prime_chain_break_sugar_modeling_.pose_list.size() || five_prime_chain_break_sugar_ID == 1;
						five_prime_chain_break_sugar_ID++ ){
				for ( Size three_prime_chain_break_sugar_ID = 1;
							three_prime_chain_break_sugar_ID <= three_prime_chain_break_sugar_modeling_.pose_list.size() || three_prime_chain_break_sugar_ID == 1; three_prime_chain_break_sugar_ID++ ){

					PoseOP start_pose = pose.clone();
					tag_into_pose( *start_pose, "" );
					instantiate_sugar( *start_pose, anchor_sugar_modeling_, anchor_sugar_ID );
					instantiate_sugar( *start_pose, current_sugar_modeling_, current_sugar_ID );
					instantiate_sugar( *start_pose, five_prime_chain_break_sugar_modeling_, five_prime_chain_break_sugar_ID );
					instantiate_sugar( *start_pose, three_prime_chain_break_sugar_modeling_, three_prime_chain_break_sugar_ID );
					pose_data_list.push_back( start_pose );

				}
			}
		}
	}

	///////Ok, finally have to remove clashes that may arise due to the fact that the floating base sugar sampling and minimization were done individually of each other///
	// NO! Floating bases are bulged. No reason for a virtual sugar sampler to move them. That can be another Sampler's job. --rd2013.
	if ( options_->virtual_sugar_legacy_mode() ){
		utility::vector1< SugarModeling > sampled_sugar_modeling_list;
		if ( anchor_sugar_modeling_.pose_list.size() > 0 )                  sampled_sugar_modeling_list.push_back( anchor_sugar_modeling_ );
		if ( current_sugar_modeling_.pose_list.size() > 0 )                 sampled_sugar_modeling_list.push_back( current_sugar_modeling_ );
		if ( five_prime_chain_break_sugar_modeling_.pose_list.size() > 0 )  sampled_sugar_modeling_list.push_back( five_prime_chain_break_sugar_modeling_ );
		if ( three_prime_chain_break_sugar_modeling_.pose_list.size() > 0 ) sampled_sugar_modeling_list.push_back( three_prime_chain_break_sugar_modeling_ );
		Pose pose_copy = pose;
		minimize_all_sampled_floating_bases( pose_copy, sampled_sugar_modeling_list, pose_data_list, scorefxn_, job_parameters_, true /*virtual_sugar_is_from_prior_step*/ );
	}

	TR.Debug << "pose_data_list.size() = " << pose_data_list.size() << std::endl;

}


////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::prepare_from_prior_sampled_sugar_jobs_for_chain_break( pose::Pose const & pose,
																																															utility::vector1< PoseOP > & pose_data_list ) {

	if ( !sampling_sugar_at_chain_break() ) utility_exit_with_message( "pose_data_list is empty for all 4 possible virtual sugars!" );

	for ( Size five_prime_chain_break_sugar_ID = 1;
				five_prime_chain_break_sugar_ID <= five_prime_chain_break_sugar_modeling_.pose_list.size() || five_prime_chain_break_sugar_ID == 1;
				five_prime_chain_break_sugar_ID++ ){
		for ( Size three_prime_chain_break_sugar_ID = 1;
					three_prime_chain_break_sugar_ID <= three_prime_chain_break_sugar_modeling_.pose_list.size() || three_prime_chain_break_sugar_ID == 1; three_prime_chain_break_sugar_ID++ ){

			PoseOP start_pose = pose.clone();
			tag_into_pose( *start_pose, "" );
			instantiate_sugar( *start_pose, five_prime_chain_break_sugar_modeling_, five_prime_chain_break_sugar_ID );
			instantiate_sugar( *start_pose, three_prime_chain_break_sugar_modeling_, three_prime_chain_break_sugar_ID );
			pose_data_list.push_back( start_pose );

		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::instantiate_sugar( pose::Pose & pose,
																																	 SugarModeling const & sugar_modeling,
																																	 Size const sugar_ID ){
	if ( sugar_modeling.pose_list.size() > 0 ) {
		tag_into_pose( pose,   tag_from_pose(pose) + tag_from_pose( *sugar_modeling.pose_list[ sugar_ID ] ) );
		copy_bulge_res_and_sugar_torsion( sugar_modeling, pose, ( *sugar_modeling.pose_list[ sugar_ID ] ), true /*instantiate_sugar*/ );
	} else{
		tag_into_pose( pose, tag_from_pose( pose ) + "_null" );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::initialize_parameters( pose::Pose const & pose ){

	// These parameters will be used for sampling virtual sugars and associated bulge. Reset them here.
	anchor_sugar_modeling_ = SugarModeling();
	current_sugar_modeling_ = SugarModeling();
	five_prime_chain_break_sugar_modeling_ = SugarModeling();
	three_prime_chain_break_sugar_modeling_ = SugarModeling();

	is_anchor_sugar_virt_         = is_anchor_sugar_virtual( pose );
	is_current_sugar_virt_           = !job_parameters_->floating_base() && is_current_sugar_virtual( pose ); // occurs less frequently -- typically only when connecting two pieces.
	is_five_prime_chain_break_sugar_virt_  = is_five_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop
	is_three_prime_chain_break_sugar_virt_ = is_three_prime_chain_break_sugar_virtual( pose ); // may occur when closing loop

	output_boolean( " is_anchor_sugar_virt_ = ", is_anchor_sugar_virt_, TR.Debug );
	output_boolean( " is_current_sugar_virt_ = ", is_current_sugar_virt_, TR.Debug );
	output_boolean( " is_five_prime_chain_break_sugar_virt_ = ", is_five_prime_chain_break_sugar_virt_, TR.Debug );
	output_boolean( " is_three_prime_chain_break_sugar_virt_ = ", is_three_prime_chain_break_sugar_virt_, TR.Debug );
	TR.Debug << std::endl;

	num_virtual_sugar_ = 0;
	if ( is_anchor_sugar_virt_ )                  num_virtual_sugar_++;
	if ( is_current_sugar_virt_ )                 num_virtual_sugar_++;
	if ( is_five_prime_chain_break_sugar_virt_ )  num_virtual_sugar_++;
	if ( is_three_prime_chain_break_sugar_virt_ ) num_virtual_sugar_++;
	TR.Debug << "num_virtual_sugar_ = " << num_virtual_sugar_ << std::endl;

	return do_consistency_checks();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseRNA_VirtualSugarJustInTimeInstantiator::do_consistency_checks(){

	bool const floating_base = job_parameters_->floating_base();
	if ( options_->sampler_assert_no_virt_sugar_sampling() /*run checks*/ ){
		if ( floating_base && num_nucleotides_ == 2 ){ //Hacky..ok the only acception right now is in floating_base + dinucleotide mode.
			runtime_assert( num_virtual_sugar_ <= 1 );
			if ( num_virtual_sugar_ == 1 ) runtime_assert( is_anchor_sugar_virt_ );
		} else{
			runtime_assert( num_virtual_sugar_ == 0 );
		}
	}

	if ( is_current_sugar_virt_ ){ //Consistency test.
		// runtime_assert( gap_size_ == 0 ); // not true for internal suites.
	}
	// not true after some kinds of merges!
	//	runtime_assert( gap_size_ == 0 || num_virtual_sugar_ <= 1 );

	if ( floating_base ){
		if ( is_five_prime_chain_break_sugar_virt_ ){ //This is rare since floating_base sampling is not often used at chain-closure step!
			TR.Debug << "WARNING: floating_base and is_five_prime_chain_break_sugar_virt_ case. Code not implemented yet, early return!" << std::endl;
			return false;
		}
		if ( is_three_prime_chain_break_sugar_virt_ ){ //This is rare since floating_base sampling is not often used at chain-closure step!
			TR.Debug << "WARNING: floating_base and is_three_prime_chain_break_sugar_virt_ case. Code not implemented yet, early return!" << std::endl;
			return false;
		}
		runtime_assert( !is_current_sugar_virt_ );
	}

	runtime_assert( !options_->do_not_sample_multiple_virtual_sugar() ||  !options_->sample_ONLY_multiple_virtual_sugar() );

	// this was just for testing.
	if ( options_->do_not_sample_multiple_virtual_sugar() &&  num_virtual_sugar_ > 1 ) return false;
	if ( options_->sample_ONLY_multiple_virtual_sugar() ){
		runtime_assert( gap_size_ == 0 );
		if ( num_virtual_sugar_ <= 1 ) return false;
	}

	return true;

}

/////////////////////
SugarModeling const &
StepWiseRNA_VirtualSugarJustInTimeInstantiator::anchor_sugar_modeling(){ return anchor_sugar_modeling_; }

/////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){ scorefxn_ = scorefxn; }

/////////////////////
std::string
StepWiseRNA_VirtualSugarJustInTimeInstantiator::get_name() const {
	return "StepWiseRNA_VirtualSugarJustInTimeInstantiator";
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_VirtualSugarJustInTimeInstantiator::set_options( StepWiseRNA_ModelerOptionsCOP options ){
	options_ = options;
}

} //sugar
} //rna
} //sampling
} //stepwise
} //protocols
