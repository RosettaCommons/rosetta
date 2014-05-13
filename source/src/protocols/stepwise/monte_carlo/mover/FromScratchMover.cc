// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/FromScratchMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>
#include <protocols/stepwise/sampling/StepWiseModeler.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <basic/Tracer.hh>

//Req'd on WIN32
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>

static basic::Tracer TR( "protocols.stepwise.monte_carlo.rna.FromScratchMover" );
using namespace protocols::stepwise::sampling;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	//Constructor
	FromScratchMover::FromScratchMover()
	{}

	//Destructor
	FromScratchMover::~FromScratchMover()
	{}

  //////////////////////////////////////////////////////////////////////////
  void
  FromScratchMover::apply( core::pose::Pose &  )
	{
		std::cout << "not defined" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
  void
  FromScratchMover::apply( core::pose::Pose & pose,
															 utility::vector1< Size > const & residues_to_instantiate_in_full_model_numbering ) const
	{
    using namespace core::chemical;
		using namespace core::pose;
		using namespace core::pose::full_model_info;

		// an alias:
		utility::vector1< Size > const & resnum = residues_to_instantiate_in_full_model_numbering;

		// only do dinucleotides for now.
		runtime_assert( resnum.size() == 2 );

		std::string new_sequence;
		std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
		for ( Size n = 1; n <= resnum.size(); n++ ) {
			char newrestype = full_sequence[ resnum[n]-1 ];
			sampling::rna::choose_random_if_unspecified_nucleotide( newrestype );
			new_sequence += newrestype;
		}

		ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		Pose new_pose;
		make_pose_from_sequence( new_pose, new_sequence, *rsd_set );

		update_full_model_info_and_switch_focus_to_new_pose( pose, new_pose, resnum );
		fix_up_jump_atoms_and_residue_type_variants( pose );

    sample_by_swa( pose, 2 );

	}


	//////////////////////////////////////////////////////////////////////////////
	void
	FromScratchMover::update_full_model_info_and_switch_focus_to_new_pose( pose::Pose & pose, pose::Pose & new_pose, utility::vector1< Size > const & resnum ) const {
		// prepare full_model_info for this new pose
		FullModelInfoOP new_full_model_info = nonconst_full_model_info( pose ).clone_info();
		FullModelInfoOP full_model_info     = nonconst_full_model_info( pose ).clone_info();

		// relieve original pose of holding information on other poses.
		full_model_info->clear_other_pose_list();
		set_full_model_info( pose, full_model_info );

		new_full_model_info->set_res_list( resnum );
		if ( pose.total_residue() > 1 ) new_full_model_info->add_other_pose( pose.clone() );
		set_full_model_info( new_pose, new_full_model_info );
		update_pdb_info_from_full_model_info( new_pose ); // for output pdb or silent file -- residue numbering.

		pose = new_pose; // switch focus to new pose.
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	FromScratchMover::sample_by_swa( pose::Pose & pose, Size const sample_res ) const {
		stepwise_modeler_->set_moving_res_and_reset( sample_res );
		stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
		stepwise_modeler_->apply( pose );
	}

	///////////////////////////////////////////////////////////////////
	void
	FromScratchMover::set_stepwise_modeler( protocols::stepwise::sampling::StepWiseModelerOP stepwise_modeler ){
		stepwise_modeler_ = stepwise_modeler;
	}


	///////////////////////////////////////////////////////////////////////////////
	std::string
	FromScratchMover::get_name() const {
		return "FromScratchMover";
	}

} //rna
} //monte_carlo
} //stepwise
} //protocols
