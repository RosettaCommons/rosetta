// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/StepWiseModeler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/StepWiseConnectionSampler.hh>
#include <protocols/stepwise/modeler/StepWiseMinimizer.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/scoring_util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/packer/StepWiseMasterPacker.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/farna/FARNA_Optimizer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.StepWiseModeler" );

using namespace core;
using namespace core::pose;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler::protein;
using utility::tools::make_vector1;

namespace protocols {
namespace stepwise {
namespace modeler {

//Constructor
StepWiseModeler::StepWiseModeler( Size const moving_res, core::scoring::ScoreFunctionCOP scorefxn ):
	scorefxn_( scorefxn ),
	figure_out_prepack_res_( false ),
	prepack_res_was_inputted_( false )
{
	set_moving_res_and_reset( moving_res );
}

//Constructor
StepWiseModeler::StepWiseModeler( core::scoring::ScoreFunctionCOP scorefxn ):
	moving_res_( 0 ),
	scorefxn_( scorefxn ),
	figure_out_prepack_res_( false ),
	prepack_res_was_inputted_( false )
{}

//Destructor
StepWiseModeler::~StepWiseModeler()
{}

//////////////////////////////////////////////////////////////////////////////
StepWiseModeler::StepWiseModeler( StepWiseModeler const & src ):
	Mover( src )
{
	*this = src;
}

//////////////////////////////////////////////////////////////////////////////
StepWiseModelerOP
StepWiseModeler::clone_modeler() const {
	return StepWiseModelerOP( new StepWiseModeler( *this ) );
}

//////////////////////////////////////////////////////////////////////////////
StepWiseModeler &
StepWiseModeler::operator=( StepWiseModeler const & src ){

	moving_res_ = src.moving_res_;
	moving_res_list_ = src.moving_res_list_;

	options_ = src.options_;
	scorefxn_ = src.scorefxn_;
	working_prepack_res_ = src.working_prepack_res_;
	working_minimize_res_ = src.working_minimize_res_;
	input_streams_ = src.input_streams_;
	figure_out_prepack_res_ = src.figure_out_prepack_res_;
	prepack_res_was_inputted_ = src.prepack_res_was_inputted_;
	precomputed_library_mover_ = src.precomputed_library_mover_;

	return *this;
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::apply( pose::Pose & pose ){

	initialize( pose );

	do_prepacking( pose );
	do_sampling( pose );
	if ( sampling_successful() ) do_minimizing( pose );

	reinitialize( pose );
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::initialize( pose::Pose & pose ){
	initialize_working_parameters_and_root( pose );
	runtime_assert( moving_res_list_.size() > 0 || prepack_res_was_inputted_ || figure_out_prepack_res_ ); // otherwise this is a no op.
	initialize_scorefunctions( pose );
	pose_list_.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::do_prepacking( core::pose::Pose & pose ) {

	master_packer_ = packer::StepWiseMasterPackerOP( new packer::StepWiseMasterPacker( working_parameters_, options_->get_sampler_options() ) );
	master_packer_->set_scorefxn( pack_scorefxn_ );
	master_packer_->initialize( pose );

	if ( working_prepack_res_.size() == 0 ) return;

	master_packer_->set_working_pack_res( working_prepack_res_ );
	master_packer_->do_prepack( pose ); // will split at moving_res_list

}

////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::do_sampling( core::pose::Pose & pose ) {

	StepWiseConnectionSampler stepwise_sampler( working_parameters_ );
	stepwise_sampler.set_options( options_->get_sampler_options() ); // careful!
	stepwise_sampler.set_scorefxn( sample_scorefxn_ );
	stepwise_sampler.set_input_streams( input_streams_ );
	stepwise_sampler.set_master_packer( master_packer_ );

	stepwise_sampler.apply( pose );
	pose_list_                 = stepwise_sampler.get_pose_list();

}

////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::do_minimizing( core::pose::Pose & pose ) {
	moves::MoverOP optimizer;
	if ( options_->lores() ) {
		optimizer = moves::MoverOP( new protocols::farna::FARNA_Optimizer( pose_list_, scorefxn_, 100 /* cycles */ ) );
	} else {
		StepWiseMinimizerOP stepwise_minimizer( new StepWiseMinimizer( pose_list_,
			working_parameters_,
			options_,
			scorefxn_ ) );
		if ( master_packer_->packer()->working_pack_res_was_inputted() ) {
			stepwise_minimizer->set_working_pack_res( master_packer_->packer()->previous_working_pack_res() );
		}
		optimizer = stepwise_minimizer;
	}

	if ( optimizer == 0 ) return;

	optimizer->apply( pose ); // will save work in pose_list_
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::reinitialize( pose::Pose & pose ){

	// if using native constraints, those got added to pose during align_pose_and_add_rmsd_constraints()
	pose.constraint_set( cst_set_ );
	( *scorefxn_ )( pose ); // updating constraints clears score.

	// Important: make sure that the next time this is used, job parameters is set explicitly -- or it will be reset.
	working_parameters_.reset();
	moving_res_list_.clear();
	working_prepack_res_.clear();
	working_minimize_res_.clear();
	figure_out_prepack_res_ = false;
	prepack_res_was_inputted_ = false;
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::set_moving_res_and_reset( Size const moving_res ){
	moving_res_ = moving_res;
	working_parameters_.reset();
}

//////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::initialize_working_parameters_and_root( pose::Pose & pose ){
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
	look_for_precompute_move( pose );
	figure_out_moving_res_list( pose );

	if ( working_parameters_ != 0 ) return;

	revise_root_and_moving_res_list( pose, moving_res_list_ ); // specify reference_res_? [i.e. anchor_res?]

	Real const rmsd_screen = options_->rmsd_screen();
	cst_set_ = pose.constraint_set()->clone();
	if ( !options_->disallow_realign() ) align::align_pose_and_add_rmsd_constraints( pose, get_native_pose(), moving_res_list_, rmsd_screen );

	working_parameters_ = working_parameters::setup_working_parameters_for_swa( moving_res_list_, pose, get_native_pose(),
		options_->bridge_res(), working_minimize_res_  );

	if ( figure_out_prepack_res_ ) {
		( *scorefxn_ )( pose ); // ideally would just be able to assemble a nbr list without going through scoring...
		working_prepack_res_ = packer::figure_out_working_interface_res( pose, get_all_working_moving_res( working_parameters_ ) );
	}
}

////////////////////////////////////////////////////////////////////
void
StepWiseModeler::figure_out_moving_res_list( pose::Pose const & pose ){
	if ( moving_res_list_.size() > 0 ) return; // in principle, SWA can explicitly give moving_res_list

	// otherwise, figure it out from moving_res_ (single residue!)
	figure_out_moving_res_list_from_most_distal_res( pose, moving_res_ );
}


////////////////////////////////////////////////////////////////////
void
StepWiseModeler::figure_out_moving_res_list_from_most_distal_res( pose::Pose const & pose, Size const moving_res ) {

	moving_res_list_.clear();
	if ( moving_res == 0 ) return;

	moving_res_list_.push_back( moving_res );

	// this goes back one residue (default behavior in protein case -- two amino acids)
	figure_out_protein_modeling_info( pose, moving_res, moving_res_list_ );

}

//////////////////////////////////////////////////////////////////////////////
bool
StepWiseModeler::sampling_successful() {
	Size const num_sampled = pose_list_.size();
	if ( num_sampled == 0 ) {
		if ( options_->num_random_samples() > 0 ) TR << TR.Red << "WARNING! WARNING! WARNING! pose_list_.size() == 0! " << TR.Reset << std::endl;
		if ( options_ && !options_->output_minimized_pose_list() ) return false; // don't do a minimize...
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::initialize_scorefunctions( pose::Pose const & pose ){
	sample_scorefxn_ = initialize_sample_scorefxn( scorefxn_, pose, options_ );
	pack_scorefxn_   = initialize_pack_scorefxn( sample_scorefxn_, pose );
}


////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseModeler::look_for_precompute_move( pose::Pose & pose ) {

	//////////////////////////////////////////////////////////////
	// This should now be deprecated in favor of SubMotifLibrary!
	//              -- rhiju, feb 2015
	//////////////////////////////////////////////////////////////
	if ( moving_res_ == 0 ) return; // a deletion -- don't trigger any kind of structural change.
	if ( precomputed_library_mover_ == 0 ) return;
	if ( !precomputed_library_mover_->has_precomputed_move( pose ) ) return;

	precomputed_library_mover_->apply( pose );

	// taken from DeleteMover -- some code duplication...
	// no further backbone resampling:
	set_moving_res_and_reset( 0 );
	// just pack and minimize:
	set_working_prepack_res( const_full_model_info( pose ).res_list() );
	set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
}

/////////////////////////////////////////////////////////////////////
void
StepWiseModeler::set_working_parameters( working_parameters::StepWiseWorkingParametersCOP working_parameters ){
	working_parameters_ = working_parameters;
}

/////////////////////////////////////////////////////////////////////
working_parameters::StepWiseWorkingParametersCOP
StepWiseModeler::working_parameters() { return working_parameters_; }

/////////////////////////////////////////////////////////////////////
void
StepWiseModeler::set_input_streams( utility::vector1< protein::InputStreamWithResidueInfoOP > const & input_streams ){
	input_streams_ = input_streams;
}

} //modeler
} //stepwise
} //protocols
