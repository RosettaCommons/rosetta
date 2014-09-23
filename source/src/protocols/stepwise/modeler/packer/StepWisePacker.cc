// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWisePacker
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/packer/SideChainCopier.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/util.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <basic/Tracer.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <utility/exit.hh>

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using namespace core;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise modeler of proteins and RNA.
//
//  Main input is obligate_pack_res -- pack these side chains as well
//  as first-level neighbors.
//
//  Will update soon to accept actual residues to pack [?], which can
//  be set *externally*.
//
//   -- Rhiju, 2014.
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.packer.StepWisePacker" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWisePacker::StepWisePacker( utility::vector1< Size > const & working_moving_res_list ):
		working_moving_res_list_( working_moving_res_list ),
		use_packer_instead_of_rotamer_trials_( false ),
		allow_virtual_side_chains_( false ),
		allow_virtual_o2prime_hydrogens_( false ),
		pack_o2prime_hydrogens_( true ),
		working_pack_res_was_inputted_( false )
  {
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWisePacker::~StepWisePacker()
  {}

	/////////////////////
	std::string
	StepWisePacker::get_name() const {
		return "StepWisePacker";
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::apply( core::pose::Pose & pose ) {
		if ( pose_init_ ) reinstate_side_chain_angles( pose, *pose_init_ );
		figure_out_neighbors( pose );
		setup_pack_task( pose );
		do_packing( pose );
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::do_packing( core::pose::Pose & pose ) {
		//		TR << TR.Blue << "before pack " << (*scorefxn_)( pose ) << " with working_pack_res " << make_tag_with_dashes( previous_working_pack_res_ ) << TR.Reset << std::endl;
		if ( use_packer_instead_of_rotamer_trials_ ) {
			pack::pack_rotamers(  pose, *scorefxn_, pack_task_ );
		} else {
			pack::rotamer_trials( pose, *scorefxn_, pack_task_ );
		}
		//		TR << TR.Red << "after pack " << (*scorefxn_)( pose ) << " with working_pack_res " << make_tag_with_dashes( previous_working_pack_res_ ) << TR.Reset << std::endl;
		// reset.
		previous_working_pack_res_ = working_pack_res_;
		working_pack_res_.clear();
	}


  //////////////////////////////////////////////////////////////////////////////////////
	// obligate_working_pack_res defines a set of residues whose neighbors *must* be packed.
	// I will change this soon to smartly look at moving residues & partitions -- likely
	// make this an external function.
	void
	StepWisePacker::figure_out_neighbors( core::pose::Pose & pose ) {
		using namespace core::scoring;

		if ( pack_all_side_chains_ ) {
			working_pack_res_was_inputted_ = false;
			working_pack_res_ = get_all_residues( pose );
			return;
		}

		if ( working_pack_res_.size() > 0 ){
			working_pack_res_was_inputted_ = true;
			return; // user inputted
		}

		working_pack_res_was_inputted_ = false;
		( *scorefxn_ )( pose ); // currently needs to occur before interface_res determination.
		working_pack_res_ = figure_out_working_interface_res( pose, working_moving_res_list_ );

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::setup_pack_task( pose::Pose const & pose ) {

		pack_task_ = pack::task::TaskFactory::create_packer_task( pose ); // create form scratch?
		pack_task_->restrict_to_repacking();

		for (Size i = 1; i <= pose.total_residue(); i++) {

			if ( working_pack_res_.has_value( i ) )  {

					pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
					pack_task_->nonconst_residue_task(i).or_include_current( true );

					if ( pose.residue(i).is_protein() ) {
						pack_task_->nonconst_residue_task(i).or_ex1( true );
						pack_task_->nonconst_residue_task(i).or_ex2( true );
						pack_task_->nonconst_residue_task(i).or_include_virtual_side_chain( allow_virtual_side_chains_ );
					} else if ( pose.residue(i).is_RNA() ) {
						if ( pack_o2prime_hydrogens_ ) {
							pack_task_->nonconst_residue_task(i).or_ex4( true );
							pack_task_->nonconst_residue_task(i).or_include_virtual_side_chain( allow_virtual_o2prime_hydrogens_ );
						}
					}

			} else {
				pack_task_->nonconst_residue_task(i).prevent_repacking();
			}
		}

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::reinstate_side_chain_angles( pose::Pose & pose, pose::Pose const & src_pose ){

		SideChainCopier copier( src_pose, previous_working_pack_res_ /*set during previous pack*/,
														pack_o2prime_hydrogens_ );
		copier.apply( pose );

	}

  //////////////////////////////////////////////////////////////////////////
	// Splits into far-apart partitions and packs.
	//  Packing can include virtual side chains (if -allow_virtual_side_chains) is on.
	//  Which side-chains will be packed? Will depend on working_pack_res, which needs to be set separately.
	//
	void
	StepWisePacker::do_prepack( pose::Pose & pose ){

		pose::Pose pose_to_split = pose;
		pose_to_split.remove_constraints(); // floating point errors if coordinate constraints are in there.
		split_pose( pose_to_split, working_moving_res_list_ );

		pose_init_.reset(); //don't copy any side chains into pose.
		runtime_assert( working_pack_res_.size() > 0 );
		apply( pose_to_split );

		SideChainCopier copier( pose_to_split, pack_o2prime_hydrogens_ );
		copier.apply( pose );
		reset( pose );
	}

	////////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::reset( core::pose::Pose const & pose ){
		// save this pose for later.
		pose_init_ = pose.clone();
		previous_working_pack_res_  = get_all_residues( pose ); // for resetting side chains.
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){
		scorefxn_ = scorefxn;
	}


} //packer
} //modeler
} //stepwise
} //protocols
