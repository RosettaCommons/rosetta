// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <protocols/simple_moves/GreenPacker.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.rna.O2PrimePacker" ) ;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace o2prime {

	//Constructor
	O2PrimePacker::O2PrimePacker( pose::Pose const & pose,
																core::scoring::ScoreFunctionCOP const & scorefxn,
																utility::vector1< core::Size > moving_res ):
		pose_with_original_HO2prime_torsion_( pose ),
		moving_res_( moving_res ),
		o2prime_pack_pose_( pose ),
		use_green_packer_( false )
	{
		if ( use_green_packer_ ){
			initialize_o2prime_green_packer();
		} else {
			initialize_o2prime_packer_task();
		}
		initialize_o2prime_pack_scorefxn( scorefxn );
	}

	//Destructor
	O2PrimePacker::~O2PrimePacker()
	{}


	////////////////////////////////////////////////////////////////////////
	void
	O2PrimePacker::initialize_o2prime_pack_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){

		using namespace core::scoring;

		o2prime_pack_scorefxn_ = new ScoreFunction;
		// Each of the following terms have been pretty optimized for the packer (trie, etc.)
		o2prime_pack_scorefxn_->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
		o2prime_pack_scorefxn_->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
		o2prime_pack_scorefxn_->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
		o2prime_pack_scorefxn_->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
		///Warning, don't include hbond_intra, since hbond_intra HAS NOT been been optimized for packing!
		o2prime_pack_scorefxn_->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) ); //// note that geom_sol is not optimized well --> replace with lk_sol for now. //IS THIS A MISTAKE???
		o2prime_pack_scorefxn_->set_energy_method_options( scorefxn->energy_method_options() ); //This set NO_HB_ENV_DEP, INCLUDE_INTRA_RES_RNA_HB and etcs.
	}

	////////////////////////////////////////////////////////////////////////
	void
	O2PrimePacker::initialize_o2prime_packer_task(){
		utility::vector1< core::Size > const O2prime_pack_seq_num = get_surrounding_O2prime_hydrogen( o2prime_pack_pose_, moving_res_, false /*verbose*/ );
		o2prime_pack_task_ = create_standard_o2prime_pack_task( o2prime_pack_pose_, O2prime_pack_seq_num );
	}

////////////////////////////////////////////////////////////////////////////////
void
O2PrimePacker::initialize_o2prime_green_packer()
{
	using namespace protocols::simple_moves;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	o2prime_green_packer_ = new protocols::simple_moves::GreenPacker;

	if ( partition_definition_.size() == 0 ) utility_exit_with_message( "To use green packer, make sure to set partition definition." );
	bool const root_partition = partition_definition_( o2prime_pack_pose_.fold_tree().root() );

	Size const nres = o2prime_pack_pose_.total_residue();
	UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new UserDefinedGroupDiscriminator );
	utility::vector1< Size > group_ids;

	Size current_group = 0;
	Size spectator_group = 1;

	for ( Size i = 1; i <= nres; i++ ) {
		if ( partition_definition_( i ) != root_partition ) {
			current_group = 0;
			TR.Debug << "GREENPACKER SAMPLER " << i << std::endl;
		} else {
			TR.Debug << "GREENPACKER SPECTATOR   " << i <<  " --> group " << spectator_group << std::endl;
		}
		group_ids.push_back( current_group );
	}

	user_defined_group_discriminator->set_group_ids( group_ids );
	o2prime_green_packer_->set_scorefunction( *o2prime_pack_scorefxn_ );
	o2prime_green_packer_->set_group_discriminator( user_defined_group_discriminator );

	TaskFactoryOP task_factory( new TaskFactory );
	task_factory->push_back( new InitializeFromCommandline );
	task_factory->push_back( new RestrictToRepacking );
	task_factory->push_back( new IncludeCurrent );
	for ( Size i = 1; i <= nres; i++ ) {
		if ( !o2prime_pack_pose_.residue( i ).is_RNA() ) continue;
		task_factory->push_back( new ExtraChiCutoff( i, 0 ) );
		task_factory->push_back( new ExtraRotamers( i, 4 /*ex4*/ ) );
	}

	o2prime_green_packer_->set_task_factory( task_factory );
	o2prime_green_packer_->set_reference_round_task_factory( task_factory );

	// This should also initialize rotamers, etc...
	o2prime_green_packer_->apply( o2prime_pack_pose_ );
}


////////////////////////////////////////////////////////////////////////
void
O2PrimePacker::sample_o2prime_hydrogen(){

	using namespace core::id;
	using namespace core::conformation;

	//reset the HO2prime torsion to its starting value as to prevent randomness due to conformations sampling order...
	stepwise::sampling::rna::copy_all_o2prime_torsions( o2prime_pack_pose_, pose_with_original_HO2prime_torsion_ );

	//TR.Debug << "Packing 2'-OH ... ";
	if ( use_green_packer_ ) {
		o2prime_green_packer_->apply( o2prime_pack_pose_ );
	} else {
		//problem with bulge variant -- need to initialize_o2prime_packer_task each time.
		initialize_o2prime_packer_task();
		pack::rotamer_trials( o2prime_pack_pose_,
													*o2prime_pack_scorefxn_,
													o2prime_pack_task_ );
	}
}

////////////////////////////////////////////////////////////////////////
void
O2PrimePacker::copy_all_o2prime_torsions( core::pose::Pose & mod_pose ) const {
	stepwise::sampling::rna::copy_all_o2prime_torsions( mod_pose, o2prime_pack_pose_ );
}

////////////////////////////////////////////////////////////////////////
pose::Pose &
O2PrimePacker::pose(){ return o2prime_pack_pose_; }

} //o2prime
} //rna
} //sampling
} //stepwise
} //protocols
