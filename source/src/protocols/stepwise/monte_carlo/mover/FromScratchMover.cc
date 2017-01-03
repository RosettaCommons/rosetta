// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/FromScratchMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <basic/Tracer.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.FromScratchMover" );
using namespace protocols::stepwise::modeler;

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

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
	utility::vector1< Size > resnum = residues_to_instantiate_in_full_model_numbering;

	// only do dinucleotides for now.
	runtime_assert( resnum.size() == 2 );

	std::string new_sequence;
	std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
	std::map< Size, std::string > const & nc_res_map = const_full_model_info( pose ).full_model_parameters()->non_standard_residue_map();
	for ( Size n = 1; n <= resnum.size(); n++ ) {
		char newrestype = full_sequence[ resnum[n]-1 ];
		modeler::rna::choose_random_if_unspecified_nucleotide( newrestype );
		new_sequence += newrestype;
		if ( nc_res_map.find( resnum[n] ) != nc_res_map.end() ) {
			new_sequence += "[" + nc_res_map.at( resnum[n] ) + "]";
		}
	}

	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	PoseOP new_pose( new Pose );
	make_pose_from_sequence( *new_pose, new_sequence, *rsd_set );

	// calebgeniesse: add virt root for density scoring, only when instantiating initial pose, for now
	if ( pose.total_residue() < 1 && full_sequence[full_sequence.size()-1] == 'X' ) {
		TR << "Adding virtual residue as root" << std::endl;
		core::pose::addVirtualResAsRoot( *new_pose );
		resnum.push_back( full_sequence.size() );
	}

	update_full_model_info_and_switch_focus_to_new_pose( pose, *new_pose, resnum );
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
	new_full_model_info->update_submotif_info_list();
	if ( pose.size() > 1 ) new_full_model_info->add_other_pose( pose.clone() );
	set_full_model_info( new_pose, new_full_model_info );
	update_pose_objects_from_full_model_info( new_pose ); // for output pdb or silent file -- residue numbering.

	protocols::scoring::set_vdw_cached_rep_screen_info_from_pose( new_pose, pose );

	pose = new_pose; // switch focus to new pose.
}

//////////////////////////////////////////////////////////////////////////////
void
FromScratchMover::sample_by_swa( pose::Pose & pose, Size const sample_res ) const {
	if ( stepwise_modeler_ == 0 ) return;
	stepwise_modeler_->set_moving_res_and_reset( sample_res );
	stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
	// currently, this will actually look for a precomputed move, without doing anything.
	stepwise_modeler_->apply( pose );
}

///////////////////////////////////////////////////////////////////
void
FromScratchMover::set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler ){
	stepwise_modeler_ = stepwise_modeler;
}


///////////////////////////////////////////////////////////////////////////////
std::string
FromScratchMover::get_name() const {
	return "FromScratchMover";
}

} //mover
} //monte_carlo
} //stepwise
} //protocols
