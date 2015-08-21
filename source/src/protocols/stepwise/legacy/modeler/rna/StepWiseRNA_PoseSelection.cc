// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/StepWiseRNA_PoseSelection.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSelection.hh>
#include <protocols/stepwise/legacy/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <basic/Tracer.hh>

using core::pose::sort_pose_by_score;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.StepWiseRNA_PoseSelection" );

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

StepWiseRNA_CountStruct local_count_data;

//Constructor
StepWiseRNA_PoseSelection::StepWiseRNA_PoseSelection( working_parameters::StepWiseWorkingParametersCOP & working_parameters,
	core::scoring::ScoreFunctionCOP scorefxn,
	StepWiseRNA_CountStruct & count_data ):
	working_parameters_( working_parameters ),
	num_pose_kept_( 108 ),
	multiplier_( 2 ), //Sort and cluster poses when the number of pose is pose_list exceed multiplier*num_pose_kept,
	cluster_rmsd_( 0.5001 ),
	PBP_clustering_at_chain_closure_( false ), //New option Aug 15 2010
	distinguish_pucker_( true ),
	current_score_cutoff_( 999999.9 ), //Feb 02, 2012
	verbose_( false ),
	count_data_( count_data )
{
	initialize_modeler_scorefxn( scorefxn );
}
//Destructor
StepWiseRNA_PoseSelection::~StepWiseRNA_PoseSelection()
{}

////////////////////Setup modeler scoring//////////////////////////////////////////////////////////////////////////////
//1. Want to increase fa_rep during the minimization phase but want to keep it at 0.12 during the sample phase
//2. Sugar scoring is always turned off during modeler stage since it screw up pose selection. (TURN IT BACK ON: RD 01/31/2010)
//3. Harmonic and Linear RNA_Chain_break scoring is always turned off during modeler stage
void
StepWiseRNA_PoseSelection::initialize_modeler_scorefxn( core::scoring::ScoreFunctionCOP & scorefxn ){
	using namespace core::scoring;
	modeler_scorefxn_ = get_modeler_scorefxn( scorefxn );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseRNA_PoseSelection::pose_selection_by_full_score( pose::Pose & current_pose, std::string const & tag ){

	using namespace core::scoring;

	count_data_.full_score_count++;

	Real const current_score = ( *modeler_scorefxn_ )( current_pose );

	update_pose_list( tag, current_pose, current_score );

	if ( ( pose_list_.size() == ( num_pose_kept_ * multiplier_ ) ) ) {
		std::sort( pose_list_.begin(), pose_list_.end(), sort_pose_by_score );
		cluster_pose_list();
		if ( pose_list_.size() > num_pose_kept_ ) {
			pose_list_.erase( pose_list_.begin() + num_pose_kept_, pose_list_.end() );
			TR.Debug << "after erasing.. pose_list.size() = " << pose_list_.size() << std::endl;
		} else {
			TR.Debug << "pose_list.size() = " << pose_list_.size() << std::endl;
		}
	}

	return current_score;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dec 18, 2009...took off alot of optimization from this code since it is very fast (not rate limiting) anyways.
void
StepWiseRNA_PoseSelection::cluster_pose_list(){

	bool const is_prepend(  working_parameters_->is_prepend() );
	Size const actually_moving_res = working_parameters_->actually_moving_res();

	utility::vector1< bool > pose_state_list( pose_list_.size(), true );

	Size num_clustered_pose = 0;

	for ( Size i = 1; i <= pose_list_.size(); i++ ) {

		if ( pose_state_list[i] == true ) {
			num_clustered_pose++;
			for ( Size j = i + 1; j <= pose_list_.size(); j++ ) {

				Real rmsd;
				if ( PBP_clustering_at_chain_closure_ && working_parameters_->gap_size() == 0 ) { //new option Aug 15, 2010..include both phosphates in rmsd calculation at chain_break
					rmsd = phosphate_base_phosphate_rmsd( ( *pose_list_[i] ), ( *pose_list_[j] ), actually_moving_res,  false /*ignore_virtual_atom*/ );
				} else {
					rmsd = suite_rmsd( ( *pose_list_[i] ), ( *pose_list_[j] ), actually_moving_res, is_prepend, false /*ignore_virtual_atom*/ );
				}

				bool const same_pucker = is_same_sugar_pucker( ( *pose_list_[i] ), ( *pose_list_[j] ), actually_moving_res );

				if ( rmsd < cluster_rmsd_ && ( same_pucker || !distinguish_pucker_ ) ) {
					pose_state_list[j] = false;
					if ( verbose_ ) {
						TR.Debug << "rmsd = " << rmsd << "  pose " << tag_from_pose( *pose_list_[j] ) << " is a neighbor of pose " << tag_from_pose( *pose_list_[i] );
						TR.Debug << " same_pucker = "; output_boolean( same_pucker, TR.Debug );
						print_sugar_pucker_state( " center_pucker = ", core::pose::rna::get_residue_pucker_state( ( *pose_list_[i] ), actually_moving_res ), TR.Debug );
						print_sugar_pucker_state( " curr_pucker = ", core::pose::rna::get_residue_pucker_state( ( *pose_list_[j] ), actually_moving_res ), TR.Debug );
						TR.Debug << std::endl;
					}
				}
			}
		}
	}

	utility::vector1< pose::PoseOP > clustered_pose_list;

	for ( Size i = 1; i <= pose_list_.size(); i++ ) {
		if ( pose_state_list[i] == true ) {
			clustered_pose_list.push_back( pose_list_[i] );
		}
	}

	pose_list_ = clustered_pose_list;

	//check if pose_list size is equal to or exceed num_pose_kept_. Important to get score_cutoff here right after clustering.
	if ( pose_list_.size() >= num_pose_kept_ ) {
		current_score_cutoff_ = total_energy_from_pose( *pose_list_[num_pose_kept_] );
	} else {
		//keep on adding pose to list if there are still not enough clusters
		current_score_cutoff_ = 999999.9; //Feb 02, 2012
	}
	////////////////////////////////////////////


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSelection::finalize( bool const do_clustering /* = true */ ){

	std::sort( pose_list_.begin(), pose_list_.end(), sort_pose_by_score );

	if ( do_clustering ) cluster_pose_list();

	if ( pose_list_.size() > num_pose_kept_ ) pose_list_.erase( pose_list_.begin() + num_pose_kept_, pose_list_.end() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSelection::update_pose_list(
	std::string const & tag,
	pose::Pose const & current_pose,
	Real const & current_score ) {

	bool add_pose_to_list = ( current_score < current_score_cutoff_ );

	//The order of evaluation of the two expression in the if statement is important!
	if ( add_pose_to_list ) {

		if ( verbose_ ) {
			TR.Debug << "tag = " << tag << " current_score_cutoff_ " << current_score_cutoff_ << " score = " << current_score;
		}
		pose_list_.push_back( current_pose.clone() );
		if ( verbose_ ) TR.Debug << " pose_list.size = " << pose_list_.size() << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSelection::set_num_pose_kept( core::Size const & num_pose_kept ){
	num_pose_kept_ = num_pose_kept ;
}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSelection::set_cluster_rmsd( Real const & setting ){
	cluster_rmsd_ = setting;
	TR.Debug << "Set cluster_rmsd to " << cluster_rmsd_ << std::endl;
}

//////////////////////////////////////////////////////////////////
utility::vector1< pose::PoseOP >
StepWiseRNA_PoseSelection::pose_list(){ return pose_list_;}

//////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSelection::set_pose_list( utility::vector1< pose::PoseOP > & pose_list ){ pose_list_ = pose_list; }

} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
