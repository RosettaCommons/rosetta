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

// Unit Headers
#include <protocols/docking/DockingEnsemblePrepackProtocol.hh>

// Package headers
#include <protocols/docking/util.hh>
#include <protocols/docking/DockingEnsemble.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/docking/SidechainMinMover.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>


using namespace protocols::moves;
using namespace core;
using namespace pack::task;

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
	pack_operations_ = SequenceMoverOP( new SequenceMover() );

	ensemble1_ = nullptr;
	ensemble2_ = nullptr;
	ensemble1_filename_ = "";
	ensemble2_filename_ = "";
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
}

void DockingEnsemblePrepackProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::dock_rtmin );
	option.add_relevant( OptionKeys::docking::sc_min );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::docking::ensemble1 );
	option.add_relevant( OptionKeys::docking::ensemble2 );
}

void DockingEnsemblePrepackProtocol::setup_pack_operation_movers()
{
	prepack_full_repack_ = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover() );
	prepack_full_repack_->score_function( scorefxn_pack() );
	prepack_full_repack_->task_factory( task_factory() );
	pack_operations_->add_mover(prepack_full_repack_);

	if ( rt_min() ) {
		rtmin_mover_ = protocols::minimization_packing::RotamerTrialsMinMoverOP( new protocols::minimization_packing::RotamerTrialsMinMover( ) );
		rtmin_mover_->score_function( scorefxn_pack() );
		rtmin_mover_->task_factory( task_factory() );
		pack_operations_->add_mover( rtmin_mover_ );
	}
	if ( sc_min() ) {
		scmin_mover_ = SidechainMinMoverOP( new SidechainMinMover() );
		scmin_mover_->set_scorefxn( scorefxn_pack() );
		scmin_mover_->set_task_factory( task_factory() );
		pack_operations_->add_mover( scmin_mover_ );
	}
	movers_setup_ = true;
}

void DockingEnsemblePrepackProtocol::finalize_setup( pose::Pose & pose ) {
	setup_foldtree( pose, partners(), movable_jumps() );
	tf2()->set_prepack_only(true);
	tf2()->create_and_attach_task_factory( this, pose );
	if ( !movers_setup_ ) {
		setup_pack_operation_movers();
	}

	core::scoring::ScoreFunctionOP scorefxn_low = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );

	if ( ensemble1_filename_ == "" || ensemble2_filename_ == "" ) utility_exit_with_message( "Must define ensemble file for both partners");
	runtime_assert( movable_jumps().size() == 1 ); // ensemble mode is only allowed for traditional docking
	core::Size const rb_jump = movable_jumps()[1];
	Size start_res(1), end_res(1), cutpoint(pose.fold_tree().cutpoint_by_jump( rb_jump ));

	TR << "Ensemble 1: " << ensemble1_filename_ << std::endl;
	start_res = 1;
	end_res = cutpoint;

	ensemble1_ = DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble1_filename_, "dock_ens_conf1", scorefxn_low, scorefxn() ) );

	TR << "Ensemble 2: " << ensemble2_filename_ << std::endl;
	start_res = cutpoint + 1;
	end_res = pose.size();

	ensemble2_ = DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble2_filename_, "dock_ens_conf2", scorefxn_low, scorefxn() ) );
}

utility::vector1< char > DockingEnsemblePrepackProtocol::get_pose_chains( core::pose::Pose & pose ) {
	// returns chains in order of appearance in pose.pdb_info().chain()
	utility::vector1< char > chain_list;
	char current_chain;
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
	utility::vector1< char > partner1_chains;
	utility::vector1< char > partner2_chains;
	utility::vector1< char > * current_partner = &partner1_chains;

	// get the string from the -partners flag
	std::string partners( get_partners() );

	// loop over the string, character by character and split it on the underscore
	// assume here that ensemble1 is reported first (in my experience this is always the case)
	for ( char & partner : partners ) {

		if ( partner == '_' ) {
			current_partner = &partner2_chains;
			continue;
		}

		current_partner->push_back(partner);

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
			utility::vector1< char > chains = get_pose_chains( c1 );

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
			utility::vector1< char > chains = get_pose_chains( c1 );

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
			if ( c1.sequence().compare(c2.sequence()) != 0 ) {
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
			if ( c1.sequence().compare(c2.sequence()) != 0 ) {
				std::string exit_message = "Ensemble 2 members are unequal!\n";
				exit_message = exit_message + "Member " + std::to_string(i) + ": " + c1.sequence() + "\n";
				exit_message = exit_message + "Member " + std::to_string(j) + ": " + c2.sequence() + "\n";
				utility_exit_with_message( exit_message );
			}

		} // end inner loop
	} // end outer looop

}

void DockingEnsemblePrepackProtocol::apply( core::pose::Pose & pose )
{
	finalize_setup(pose);
	check_ensemble_member_compatibility();
	protocols::docking::ConformerSwitchMoverOP switch_mover;
	core::pose::Pose starting_pose;

	starting_pose = pose;

	switch_mover = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble1_ ) );
	for ( Size i=1; i<=ensemble1_->size(); ++i ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );
		switch_mover->switch_conformer( pose, i );

		//Move each partners away from the others
		for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
		ensemble1_->calculate_lowres_ref_energy( pose );
		//bringing the packed structures together
		for (  DockJumps::const_iterator jump= movable_jumps().begin() ; jump != movable_jumps().end(); ++jump ) {
			rigid::RigidBodyTransMoverOP translate_back( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);
		}

		ensemble1_->set_packer( pack_operations_ );

		ensemble1_->calculate_highres_ref_energy( i ); // also does the dump_pdb
	}
	ensemble1_->update_pdblist_file();

	// reset to starting pose
	pose = starting_pose;
	switch_mover = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble2_ ) );
	for ( Size i=1; i<=ensemble2_->size(); ++i ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );

		switch_mover->switch_conformer( pose, i );

		//Move each partners away from the others
		for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
		ensemble2_->calculate_lowres_ref_energy( pose );
		//bringing the packed structures together
		for (  DockJumps::const_iterator jump= movable_jumps().begin() ; jump != movable_jumps().end(); ++jump ) {
			rigid::RigidBodyTransMoverOP translate_back( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);
		}

		ensemble2_->set_packer( pack_operations_ );

		ensemble2_->calculate_highres_ref_energy( i ); // also does the dump_pdb
	}
	ensemble2_->update_pdblist_file();

	// for the sake of naming consistency (JRJ)
	// get prefix, append _prepack.pdb, output
	std::string basename = utility::file::file_basename(pose.pdb_info()->name());
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom( core::chemical::FA_STANDARD );
	to_fullatom.apply( pose ); // go high res
	pose.dump_pdb( basename + ".prepack.pdb" );
}

std::string DockingEnsemblePrepackProtocol::get_name() const {
	return "DockingEnsemblePrepackProtocol";
}

}
}
