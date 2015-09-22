// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesFixBBProtocol.cc
///
/// @brief
/// @author Florian Richter


#include <protocols/enzdes/EnzdesFixBBProtocol.hh>
#include <protocols/enzdes/EnzdesBaseProtocol.hh>
#include <protocols/enzdes/EnzdesMovers.hh>
#include <protocols/enzdes/ModifyStoredLigandRBConfsMovers.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/ConsensusDesignMover.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace enzdes {

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.EnzdesFixBBProtocol" );

EnzdesFixBBProtocol::EnzdesFixBBProtocol()
: EnzdesBaseProtocol(),
	start_from_random_rb_conf_( basic::options::option[basic::options::OptionKeys::enzdes::start_from_random_rb_conf] )
{}

EnzdesFixBBProtocol::~EnzdesFixBBProtocol(){}

void
EnzdesFixBBProtocol::apply(
	core::pose::Pose & pose
){


	using namespace protocols::moves;
	using namespace core::pack::task;

	// Scoring function already set up by superclass
	tr.Info << "starting apply function..." << std::endl;

	//set the native pose if requested
	if ( ! basic::options::option[basic::options::OptionKeys::in::file::native].user() ) {

		core::pose::PoseOP natpose( new core::pose::Pose( pose ) );
		(*scorefxn_)( *natpose );
		this->set_native_pose( natpose );
	}

	//set up constraints (read cstfile, do mapping, etc, then add to pose)
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ) {
		enable_constraint_scoreterms();
		setup_enzdes_constraints( pose, false );
	}

	if ( start_from_random_rb_conf_ ) {
		ApplyRandomStoredRBConf ranconf;
		ranconf.apply( pose );
	}

	//create packer task (read resfile, etc)
	PackerTaskOP design_pack_task;

	tr.Info << "Done setting up the task and constraints... " << std::endl;
	//score pose to make sure everything is initialised correctly
	(*scorefxn_)( pose );

	//cst opt stage, if demanded
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_opt] ) {
		design_pack_task =  create_enzdes_pack_task( pose );
		tr.Info << "starting cst_opt minimization..." << std::endl;
		cst_minimize(pose, design_pack_task, true);
		(*scorefxn_)( pose );
		tr.Info << "done cst_opt minimization." << std::endl;
	}


	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_predock] ) {
		//design_pack_task =  create_enzdes_pack_task( pose );
		PredesignPerturbMoverOP predock( new PredesignPerturbMover() );
		predock->set_ligand( get_ligand_id(pose, pose.num_jump()) );
		predock->apply(pose);
		(*scorefxn_)( pose );
	}


	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_design] ) {

		design_pack_task = create_enzdes_pack_task( pose ); //make a new task in case the ligand has moved a lot
		tr.Info << "starting cst_design, " << basic::options::option[basic::options::OptionKeys::enzdes::design_min_cycles] << " cycles of design/minimization ... " << std::endl;

		core::Size design_min_cycles = basic::options::option[basic::options::OptionKeys::enzdes::design_min_cycles];

		//  bool favor_native_res(false);
		//  if( basic::options::option[basic::options::OptionKeys::enzdes::favor_native_res].user() ) favor_native_res = true;

		enzdes_pack( pose, design_pack_task, scorefxn_, design_min_cycles, basic::options::option[basic::options::OptionKeys::enzdes::cst_min], false, true );

		design_pack_task = create_enzdes_pack_task( pose, false );

		if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_min] ) {
			remove_enzdes_constraints( pose, true );
			cst_minimize(pose, design_pack_task);
			add_pregenerated_enzdes_constraints( pose );
		}
		(*scorefxn_)( pose );

	} else if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_min] ) { //if cst_design

		design_pack_task = create_enzdes_pack_task( pose );

		cst_minimize(pose, design_pack_task);

		(*scorefxn_)( pose );
	}

	if ( basic::options::option[basic::options::OptionKeys::enzdes::make_consensus_mutations] ) {

		simple_moves::ConsensusDesignMover consensus_mover( create_enzdes_pack_task( pose, false ), scorefxn_ );
		consensus_mover.set_invert_task( true );
		consensus_mover.set_use_seqprof_constraints( true );
		consensus_mover.set_sasa_cutoff( 1.0 ); //let's prevent totally buried residues from being redesigned
		consensus_mover.apply( pose );
	}

	PackerTaskOP repack_task;


	//do a repack without constraints
	if ( ! basic::options::option[basic::options::OptionKeys::enzdes::no_unconstrained_repack] ) {

		remove_enzdes_constraints( pose, true );
		(*scorefxn_)( pose );
		repack_task = create_enzdes_pack_task( pose, false ); //remake task in case the ligand has moved a lot
		tr.Info << "Starting after design unconstrained repack/minimization... " << std::endl;
		protocols::simple_moves::PackRotamersMoverOP enzdes_repack( new protocols::simple_moves::PackRotamersMover(scorefxn_, repack_task) );
		enzdes_repack->apply( pose );

		if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_min] ) cst_minimize(pose, repack_task);

		//and turn constraints back on for the final scoring
		add_pregenerated_enzdes_constraints( pose );
		tr.Info <<"Finished after design unconstrained repack/minimization... " << std::endl;

		(*scorefxn_)( pose );

	}

	if ( basic::options::option[basic::options::OptionKeys::enzdes::cst_dock] ) {


		//note: this is not really ready to go yet, still to be developed
		//to do: write a wrapper that executes the docking protocol a number
		//of times, either till a maximum number of decoys is produced, or
		//until a decoy is found that has a better energy than the designed pose

		tr.Info << "Starting ligand docking... " << std::endl;
		remove_enzdes_constraints( pose, false );
		//write stuff for ligand docking protocol
		protocols::ligand_docking::LigandDockProtocolOP dock_lig_protocol( new protocols::ligand_docking::LigandDockProtocol() );
		dock_lig_protocol->apply( pose );

		add_pregenerated_enzdes_constraints( pose );
		(*scorefxn_)( pose );

	}


} //apply function


std::string
EnzdesFixBBProtocol::get_name() const {
	return "EnzdesFixBBProtocol";
}

void
EnzdesFixBBProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::enzdes::cstfile );
	option.add_relevant( OptionKeys::enzdes::cst_opt );
	option.add_relevant( OptionKeys::enzdes::cst_predock );

	PredesignPerturbMover::register_options();

	option.add_relevant( OptionKeys::enzdes::cst_design );
	option.add_relevant( OptionKeys::enzdes::design_min_cycles );
	option.add_relevant( OptionKeys::enzdes::no_unconstrained_repack );

	option.add_relevant( OptionKeys::in::file::pssm);
}

} //namespace enzdes
} //namespace protocols

