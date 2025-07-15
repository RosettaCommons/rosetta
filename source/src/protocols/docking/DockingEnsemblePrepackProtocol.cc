// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingEnsemblePrepackProtocol.cc
/// @brief Prepacking of the bound structure before docking with ensembles
/// @author Monica Berrondo
/// @author Modified by Jeliazko Jeliazkov -- check_ensemble_member_compatibility()
/// @author Modified by Ameya Harmalkar and Rituparna Samanta --

// Unit Headers
#include <protocols/docking/DockingEnsemblePrepackProtocol.hh>

// Package headers
#include <protocols/docking/util.hh>
#include <protocols/docking/DockingEnsemble.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/docking/SidechainMinMover.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <protocols/membrane/util.hh>
#include <core/conformation/Conformation.hh>


#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/membrane/MPFastRelaxMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
// Utility Headers
#include <basic/Tracer.hh>


#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>


using namespace protocols::moves;
using namespace core;
using namespace pack::task;
using namespace protocols::membrane;
using namespace core::conformation::membrane;

static basic::Tracer TR( "protocols.docking.DockingEnsemblePrepackProtocol" );

namespace protocols {
namespace docking {

DockingEnsemblePrepackProtocol::DockingEnsemblePrepackProtocol(): DockingHighRes()
{
	Mover::type( "DockingEnsemblePrepackProtocol" );
	setup_defaults();
	register_options();
	init_from_options();
}


void DockingEnsemblePrepackProtocol::setup_defaults()
{
	trans_magnitude_ = 1000.0;
	pack_operations_ = utility::pointer::make_shared< SequenceMover >();

	ensemble1_ = nullptr;
	ensemble2_ = nullptr;
	ensemble1_filename_ = "";
	ensemble2_filename_ = "";
	membrane_ = false;
	movers_setup_ = false;
}

DockingEnsemblePrepackProtocol::~DockingEnsemblePrepackProtocol()= default;

void DockingEnsemblePrepackProtocol::init_from_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::docking::dock_rtmin ].user() ) {
		set_rt_min(option[ OptionKeys::docking::dock_rtmin ]());
	}

	if ( option[ OptionKeys::docking::sc_min ].user() ) {
		set_sc_min(option[ OptionKeys::docking::sc_min ]());
	}

	if ( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	if ( option[ OptionKeys::docking::ensemble1 ].user() ) {
		set_ensemble1(option[ OptionKeys::docking::ensemble1 ]());
	}

	if ( option[ OptionKeys::docking::ensemble2 ].user() ) {
		set_ensemble2(option[ OptionKeys::docking::ensemble2 ]());
	}

	if ( option[ OptionKeys::mp::setup::spanfiles ].user() ) {
		membrane_ = true;
	}

	if ( option[ OptionKeys::mp::setup::span1 ].user() ) {
		set_spanfile1(option[ OptionKeys::mp::setup::span1 ]());
	}

	if ( option[ OptionKeys::mp::setup::span2 ].user() ) {
		set_spanfile2(option[ OptionKeys::mp::setup::span2 ]());
	}

}

void DockingEnsemblePrepackProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::dock_rtmin );
	option.add_relevant( OptionKeys::docking::sc_min );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::docking::ensemble1 );
	option.add_relevant( OptionKeys::docking::ensemble2 );
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
	option.add_relevant( OptionKeys::mp::setup::span1 );
	option.add_relevant( OptionKeys::mp::setup::span2 );
}

void DockingEnsemblePrepackProtocol::setup_pack_operation_movers()
{
	prepack_full_repack_ = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >();
	prepack_full_repack_->score_function( scorefxn_pack() );
	prepack_full_repack_->task_factory( task_factory() );
	pack_operations_->add_mover(prepack_full_repack_);

	if ( rt_min() ) {
		rtmin_mover_ = utility::pointer::make_shared< protocols::minimization_packing::RotamerTrialsMinMover >( );
		rtmin_mover_->score_function( scorefxn_pack() );
		rtmin_mover_->task_factory( task_factory() );
		pack_operations_->add_mover( rtmin_mover_ );
	}
	if ( sc_min() ) {
		scmin_mover_ = utility::pointer::make_shared< SidechainMinMover >();
		scmin_mover_->set_scorefxn( scorefxn_pack() );
		scmin_mover_->set_task_factory( task_factory() );
		pack_operations_->add_mover( scmin_mover_ );
	}
	movers_setup_ = true;
}

void DockingEnsemblePrepackProtocol::finalize_setup( pose::Pose & pose ) {

	// pose.dump_pdb( "pose_before_setupmem.pdb" );
	// create a membrane protein from the pose
	if ( membrane_ ) {
		membrane::AddMembraneMoverOP add_mem( new membrane::AddMembraneMover( ) );
		add_mem->apply( pose );

		membrane::TransformIntoMembraneMoverOP transform_into_memb( utility::pointer::make_shared< TransformIntoMembraneMover >() );
		transform_into_memb->apply( pose );

		// pose.dump_pdb( "pose_in_finalize_setupmem.pdb" );

		//Set the foldtree for membrane proteins
		core::Size dock_jump = create_membrane_docking_foldtree_from_partners( pose, partners() );

		// set DockJumps in DockingHighres protocols
		DockJumps dock_jumps;
		dock_jumps.push_back( static_cast< int >( dock_jump ) );
		set_movable_jumps( dock_jumps );

		//a default score function was set up here.
		//core::scoring::ScoreFunctionOP mem_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
		//set_scorefxn( mem_sfxn );
		if ( span1_filename_ == "" || span2_filename_ == "" ) utility_exit_with_message( "Must define span file for both partners");

	} else {

		docking::setup_foldtree( pose, partners(), movable_jumps() );
	}

	tf2()->set_prepack_only(true);
	tf2()->create_and_attach_task_factory( this, pose );
	if ( !movers_setup_ ) {

		setup_pack_operation_movers();
	}

	core::scoring::ScoreFunctionOP scorefxn_low = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );

	if ( ensemble1_filename_ == "" || ensemble2_filename_ == "" ) utility_exit_with_message( "Must define ensemble file for both partners");
	runtime_assert( movable_jumps().size() == 1 ); // ensemble mode is only allowed for traditional docking

	// TR << "movable_jumps: " << movable_jumps()[1] << std::endl;

	core::Size const rb_jump = movable_jumps()[1];
	core::Size start_res(1), end_res(1), cutpoint(pose.fold_tree().cutpoint_by_jump( rb_jump ));

	TR << "Ensemble 1: " << ensemble1_filename_ << std::endl;
	start_res = 1;
	end_res = cutpoint;

	ensemble1_ = utility::pointer::make_shared< DockingEnsemble >( start_res, end_res, rb_jump, ensemble1_filename_, "dock_ens_conf1", scorefxn_low, scorefxn() );
	//ensemble1_ is the pointer described in DockingEnsemble.hh
	TR << "Ensemble 2: " << ensemble2_filename_ << std::endl;
	start_res = cutpoint + 1;
	if ( membrane_ ) {
		end_res = pose.size()-1;
	} else {
		end_res = pose.size();
	}
	ensemble2_ = utility::pointer::make_shared< DockingEnsemble >( start_res, end_res, rb_jump, ensemble2_filename_, "dock_ens_conf2", scorefxn_low, scorefxn() );

}

utility::vector1< std::string > DockingEnsemblePrepackProtocol::get_pose_chains( core::pose::Pose & pose ) {
	// returns chains in order of appearance in pose.pdb_info().chain()
	utility::vector1< std::string > chain_list;
	std::string current_chain;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( i==1 ) {
			current_chain = pose.pdb_info()->chain(i);
			chain_list.push_back(current_chain);
		} else {
			if ( pose.pdb_info()->chain(i) != current_chain ) {
				chain_list.push_back(pose.pdb_info()->chain(i));
				current_chain = pose.pdb_info()->chain(i);
			}
		}
	}
	return chain_list;
}

void DockingEnsemblePrepackProtocol::check_ensemble_member_compatibility() {

	// use the -partners flag to get vectors of chains
	// compare to chains of each ensemble
	utility::vector1< std::string > partner1_chains;
	utility::vector1< std::string > partner2_chains;
	utility::vector1< std::string > * current_partner = &partner1_chains;

	// get the string from the -partners flag
	std::string partners( get_partners() );

	// loop over the string, character by character and split it on the underscore
	// assume here that ensemble1 is reported first (in my experience this is always the case)
	// (Note that due to input format limitations, only single letter chains are supported)
	for ( char & partner : partners ) {

		if ( partner == '_' ) {
			current_partner = &partner2_chains;
			continue;
		}

		current_partner->push_back(std::string{partner});

	}
	// ensemble1_/ensemble2_ must be present i.e. not NULL
	// note ensembles index at 1
	if ( !ensemble1_ || !ensemble2_ ) utility_exit_with_message( "Ensembles must be loaded, otherwise comparison is nonsensical!" );

	TR.Debug << "Ensemble 1 length is: " << ensemble1_->size() << std::endl;
	TR.Debug << "Ensemble 2 length is: " << ensemble2_->size() << std::endl;

	// check if partners flag has at least two partners "A_B" before doing partners flag comparisons
	if ( partners.size() > 2 ) {

		// chain based checks only: ensemble 1
		for ( core::Size i=1; i<=ensemble1_->size(); ++i ) { // outer loop
			core::pose::Pose c1 = ensemble1_->get_conformer(i);
			utility::vector1< std::string > chains = get_pose_chains( c1 );

			// if there is a different number of chains in any ensemble member vs. the partners flag, error!
			if ( chains.size() != partner1_chains.size() ) {
				std::string exit_message = "Ensemble 1 member differs in number of chains from partners flag!\n";
				exit_message = exit_message + "Partner flag has " + std::to_string(partner2_chains.size()) + " chains.\n";
				exit_message = exit_message + "Member " + std::to_string(i) + " has " + std::to_string(chains.size()) + " chains!\n";
				utility_exit_with_message( exit_message );
			}

			// if the chain identities are not equivalent, error!
			// assuming ensemble 1 chains are first reported in parterns flag
			for ( core::Size k=1; k<chains.size(); ++k ) {
				if ( chains[k] != partner1_chains[k] ) {
					std::string exit_message = "Ensemble 1 member differs in chain identity from partners flag!\n";
					exit_message = exit_message + "Member " + std::to_string(i) + ": " + chains[k] + " vs. " + partner1_chains[k] + "\n";
					utility_exit_with_message( exit_message );
				}
			}
		}

		// chain based checks only: ensemble 2
		for ( core::Size i=1; i<=ensemble2_->size(); ++i ) { // outer loop
			core::pose::Pose c1 = ensemble2_->get_conformer(i);
			utility::vector1< std::string > chains = get_pose_chains( c1 );

			// if there is a different number of chains in any ensemble member vs. the partners flag, error!
			if ( chains.size() != partner2_chains.size() ) {
				std::string exit_message = "Ensemble 2 member differs in number of chains from partners flag!\n";
				exit_message = exit_message + "Partner flag has " + std::to_string(partner2_chains.size()) + " chains.\n";
				exit_message = exit_message + "Member " + std::to_string(i) + " has " + std::to_string(chains.size()) + " chains!\n";
				utility_exit_with_message( exit_message );
			}

			// if the chain identities are not equivalent, error!
			// assuming ensemble 2 chains are second reported in parterns flag
			for ( core::Size k=1; k<chains.size(); ++k ) {
				if ( chains[k] != partner2_chains[k] ) {
					std::string exit_message = "Ensemble 2 member differs in chain identity from partners flag!\n";
					exit_message = exit_message + "Member " + std::to_string(i) + ": " + chains[k] + " vs. " + partner2_chains[k] + "\n";
					utility_exit_with_message( exit_message );
				}
			}
		}

	} else {
		TR.Warning << "-partners is not specified and prepack partially cannot check compatibility of ensembles. EnsembleDock may fail." << std::endl;
	}

	// loop over all member pairs and compare sequences
	// ensures ensemble docking doesn't break when swapping
	// Ensemble 1
	for ( core::Size i=1; i<=ensemble1_->size(); ++i ) { // outer loop
		core::pose::Pose c1 = ensemble1_->get_conformer(i);

		for ( core::Size j=1; j<i; ++j ) { // inner loop
			core::pose::Pose c2 = ensemble1_->get_conformer(j);

			TR.Debug << "Comparing Ensemble 1 members to each other..." << std::endl;
			TR.Debug << "Ensemble 1, member " << std::to_string(i) << ": " << c1.sequence() << std::endl;
			TR.Debug << "Ensemble 1, member " << std::to_string(j) << ": " << c2.sequence() << std::endl;

			// if sequence of conformer i, does not match j, then error!
			if ( c1.sequence() != c2.sequence() ) {
				std::string exit_message = "Ensemble 1 members are unequal!\n";
				exit_message = exit_message + "Member " + std::to_string(i) + ": " + c1.sequence() + "\n";
				exit_message = exit_message + "Member " + std::to_string(j) + ": " + c2.sequence() + "\n";
				utility_exit_with_message( exit_message );
			}

		} // end inner loop
	} // end outer looop


	// Ensemble 2
	for ( core::Size i=1; i<=ensemble2_->size(); ++i ) { // outer loop
		core::pose::Pose c1 = ensemble2_->get_conformer(i);

		for ( core::Size j=1; j<i; ++j ) { // inner loop
			core::pose::Pose c2 = ensemble2_->get_conformer(j);

			TR.Debug << "Comparing Ensemble 2 members to each other..." << std::endl;
			TR.Debug << "Ensemble 2, member " << std::to_string(i) << ": " << c1.sequence() << std::endl;
			TR.Debug << "Ensemble 2, member " << std::to_string(j) << ": " << c2.sequence() << std::endl;

			// if sequence of conformer i, does not match j, then error!
			if ( c1.sequence() != c2.sequence() ) {
				std::string exit_message = "Ensemble 2 members are unequal!\n";
				exit_message = exit_message + "Member " + std::to_string(i) + ": " + c1.sequence() + "\n";
				exit_message = exit_message + "Member " + std::to_string(j) + ": " + c2.sequence() + "\n";
				utility_exit_with_message( exit_message );
			}

		} // end inner loop
	} // end outer looop
}

void DockingEnsemblePrepackProtocol::move_away( core::pose::Pose & pose )
{
	for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {

		if ( membrane_ ) {
			// get membrane axis
			core::Vector trans_axis( membrane_axis( pose, *jump ) );
			TR << "trans axis: " << trans_axis.x() <<"i+ " << trans_axis.y() << "j +" << trans_axis.z() << "k" << std::endl;
			//this axis is calculated in a hacky way.

			// create new translation mover
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(trans_axis, *jump) );

			// do actual translation
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);

		} else {

			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );

			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
	}
}

void DockingEnsemblePrepackProtocol::move_back( core::pose::Pose & pose )
{
	for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {

		if ( membrane_ ) {

			core::Vector trans_axis( protocols::membrane::membrane_axis( pose, *jump ) );
			TR << "membrane trans axis: " << trans_axis.x() <<"i +" << trans_axis.y() << "j +" << trans_axis.z() << "k" << std::endl;

			rigid::RigidBodyTransMoverOP translate_back ( new rigid::RigidBodyTransMover(trans_axis, *jump) );
			translate_back->step_size( trans_magnitude_ );

			// translate_back->trans_axis().negate();
			translate_back->apply(pose);

		} else {
			rigid::RigidBodyTransMoverOP translate_back ( new rigid::RigidBodyTransMover(pose, *jump) );

			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);
		}
	}
}

void DockingEnsemblePrepackProtocol::apply( core::pose::Pose & pose )
{   using namespace core::conformation::membrane;
	using core::conformation::membrane::ImplicitLipidInfo;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;

	finalize_setup(pose);
	check_ensemble_member_compatibility();
	protocols::docking::ConformerSwitchMoverOP switch_mover;
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom( core::chemical::FA_STANDARD );
	core::pose::Pose starting_pose;
	core::pose::PoseOP partner1( new core::pose::Pose() );
	core::pose::PoseOP partner1_mem( new core::pose::Pose() );
	core::pose::PoseOP partner2( new core::pose::Pose() );
	core::pose::PoseOP partner2_mem( new core::pose::Pose() );
	core::pose::PoseOP pose_move( new core::pose::Pose() );

	starting_pose = pose;
	// TR << "Fold_tree of mem protein: "<< starting_pose.fold_tree() << std::endl;
	move_away( pose );

	/*TR << "move_away applied"<<std::endl;
	pose.dump_pdb("pose_moved_away.pdb");*/
	ensemble1_->set_packer( pack_operations_ );
	ensemble2_->set_packer( pack_operations_ );

	// This splits the pose along the movable jump and adds it as the first member of the ensemble.
	core::pose::partition_pose_by_jump( pose, movable_jumps()[1], *partner1, *partner2 );

	//TR <<  "partioning passed" << std::endl;
	// SSRB: This is an inelegant solution to the problem that partition_pose_by_jump messes up the residue numbering.
	// To get scoring to work, we dump the pose and read it back in.
	partner1->dump_pdb( "partner1_temp.pdb" );
	partner2->dump_pdb( "partner2_temp.pdb" );

	core::import_pose::pose_from_file( *partner1, "partner1_temp.pdb" , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( *partner2, "partner2_temp.pdb" , core::import_pose::PDB_file);

	//remove( "partner2_temp.pdb" );
	if ( membrane_ ) {

		// Let's check if the protein chains have membranes
		TR.Debug << "Does Partner 1 have a mem atom? : " << partner1->conformation().is_membrane() << std::endl;
		TR.Debug << "Does Partner 2 have a mem atom? : " << partner2->conformation().is_membrane() << std::endl;

		std::string spanfile_partner1 = span1_filename_;
		std::string spanfile_partner2 = span2_filename_;

		TR << "gathering pose membrane info " <<std::endl;

		// pose.membrane_info()->implicit_lipids()->show();//(partner2->conformation());
		std::string lipid_composition( pose.membrane_info()->implicit_lipids()->lipid_composition_name() );
		core::Real membrane_temperature( pose.membrane_info()->implicit_lipids()->temperature() );
		core::Size mem_rsd_p1( partner1->size() );
		//partner2 is the one without membrane
		membrane::AddMembraneMoverOP add_mem2( new membrane::AddMembraneMover( spanfile_partner2 ) );
		add_mem2->apply( (*partner2) );

		/* //Set the foldtree for membrane proteins
		core::Size dock_jump2 = create_membrane_foldtree_anchor_com( *partner2 );
		// set DockJumps in DockingHighres protocols
		DockJumps dock_jumps2;
		dock_jumps2.push_back( static_cast< int >( dock_jump2 ) );
		set_movable_jumps( dock_jumps2 );
		TR << "Fold_tree of mem protein partner2: " << partner2->fold_tree() << std::endl;
		*/
		membrane::AddMembraneMoverOP add_mem1( new membrane::AddMembraneMover( spanfile_partner1, mem_rsd_p1, lipid_composition, membrane_temperature ) );
		add_mem1->apply( (*partner1) );

	}
	partner2->dump_pdb( "partner2_temp2.pdb" );
	partner1->dump_pdb( "partner1_temp2.pdb");


	//It is necessary to relax the input pose, else the score comparison is really poor
	//protocols::relax::membrane::MPFastRelaxMoverOP start_struc_relax( new protocols::relax::membrane::MPFastRelaxMover( scorefxn_pack() ));//, 1 ) ); //Just one cycle to save time
	protocols::relax::FastRelaxOP start_struc_relax( new protocols::relax::FastRelax( scorefxn_pack(), 1 ) ); //Just one cycle to save time
	start_struc_relax->set_task_factory( task_factory()->clone() );
	start_struc_relax->constrain_relax_to_start_coords( true );

	// TR << "applying structural relax on partner1" <<std::endl;
	start_struc_relax->apply( *partner1 );
	// TR << "applying structural relax on partner2" <<std::endl;
	start_struc_relax->apply( *partner2 );

	ensemble1_->calculate_highres_ref_energy( *partner1, "partner1" );
	ensemble2_->calculate_highres_ref_energy( *partner2, "partner2" );

	// TR << " calculate high res ref energy "<<std::endl;
	to_centroid.apply( pose );

	ensemble1_->calculate_lowres_ref_energy( pose );
	ensemble2_->calculate_lowres_ref_energy( pose );
	to_fullatom.apply( pose );
	// required for score function reset

	/*TR << " calculate low res ref energy "<<std::endl;
	pose.dump_pdb("pose_before_moveback.pdb");*/

	// Abort if the membrane is not fixed. Leads to significant rounding errors
	if ( membrane_ ) {
		TR << "The pose has a mem atom. Check for translations while prepacking" << std::endl;
		if ( is_membrane_fixed( pose ) ) {
			TR << "Membrane orientation is fixed and the pose is moveable" << std::endl;
		} else {
			std::string exit_message = "ERROR: Translations while prepacking are not accurate. Either membrane is movable or pose is fixed\n";
			utility_exit_with_message( exit_message );
		}
	}

	move_back( pose );

	// TR << " Fold tree pose: " << pose.fold_tree () <<std::endl;
	// TR << " move back done "<<std::endl;
	//pose.dump_pdb("pose_after_moveback.pdb");
	starting_pose = pose;

	switch_mover = utility::pointer::make_shared< protocols::docking::ConformerSwitchMover >( ensemble1_ );
	membrane::AddMembraneMoverOP add_mem( new membrane::AddMembraneMover( ) );
	for ( core::Size i=1; i<=ensemble1_->size(); ++i ) {

		to_centroid.apply( pose );
		switch_mover->switch_conformer( pose, i );
		// TR << "switch_mover" << std::endl;
		move_away( pose );
		// TR << "move away done" << std::endl;
		/*pose.dump_pdb("pose_after_moveaway_ensemble.pdb");*/

		ensemble1_->calculate_lowres_ref_energy( pose );
		// TR << i << " inside ensemble loop" << std::endl;
		move_back( pose );
		// TR << "move back done" << std::endl;

		if ( membrane_ ) {
			add_mem->apply( ensemble1_->get_conformer(i) );
		}

		ensemble1_->calculate_highres_ref_energy( i ); // also does the dump_pdb

	}

	TR << "outside the first ensemble loop" << std::endl;
	ensemble1_->update_pdblist_file( "partner1" );

	// reset to starting pose
	pose = starting_pose;
	switch_mover = utility::pointer::make_shared< protocols::docking::ConformerSwitchMover >( ensemble2_ );
	for ( core::Size i=1; i<=ensemble2_->size(); ++i ) {

		to_centroid.apply( pose );
		switch_mover->switch_conformer( pose, i );

		move_away( pose );

		ensemble2_->calculate_lowres_ref_energy( pose );

		move_back( pose );
		if ( membrane_ ) {
			add_mem->apply( ensemble2_->get_conformer(i) );
		}
		//ensemble2_->set_packer( pack_operations_ );
		ensemble2_->calculate_highres_ref_energy( i ); // also does the dump_pdb
	}

	TR << "outside the second ensemble loop" << std::endl;
	ensemble2_->update_pdblist_file( "partner2" );

	// reset to starting pose
	pose = starting_pose;
	move_away( pose );
	pack_operations_->apply( pose );
	move_back( pose );

	// for the sake of naming consistency (JRJ)
	// get prefix, append _prepack.pdb, output
	( *scorefxn_pack() )( pose );
	//to_fullatom.apply( pose ); // go high res
	//pose.dump_pdb( basename + ".prepack.pdb" ); //SSRB:Taking out .prepack.pdb because JD2 already handles the output


}

std::string DockingEnsemblePrepackProtocol::get_name() const {
	return "DockingEnsemblePrepackProtocol";
}

void DockingEnsemblePrepackProtocol::set_ensemble1( std::string const &ensemble1 ) {
	ensemble1_filename_ = copy_ensemble( ensemble1 );
}

void DockingEnsemblePrepackProtocol::set_ensemble2( std::string const &ensemble2 ) {
	ensemble2_filename_ = copy_ensemble( ensemble2 );
}

void DockingEnsemblePrepackProtocol::set_spanfile1( std::string const &spanfile1 ) {
	span1_filename_ = spanfile1 ;
}

void DockingEnsemblePrepackProtocol::set_spanfile2( std::string const &spanfile2 ) {
	span2_filename_ = spanfile2 ;
}

}//docking
}//protocols
