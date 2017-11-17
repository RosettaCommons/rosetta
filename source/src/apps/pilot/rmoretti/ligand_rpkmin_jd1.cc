// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   demo/ian/jobdist_rpkmin.cc
///
/// @brief
/// @author Ian Davis (ian.w.davis@gmail.com)


#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>


#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//////////////////////////////////////////////////////////////////////////////
class LigandRepackMinimizeProtocol; // fwd declaration
typedef utility::pointer::shared_ptr< LigandRepackMinimizeProtocol > LigandRepackMinimizeProtocolOP;
typedef utility::pointer::shared_ptr< LigandRepackMinimizeProtocol const > LigandRepackMinimizeProtocolCOP;

class LigandRepackMinimizeProtocol : public protocols::ligand_docking::LigandBaseProtocol
{
public:

	inline
	LigandRepackMinimizeProtocol():
		LigandBaseProtocol()
	{
		Mover::type( "LigandRepackMinimizeProtocol" );
	}

	virtual void apply( core::pose::Pose & pose );

}; // class LigandRepackMinimizeProtocol


/// @brief Creates a new hierarchy of Movers for each Pose passed in.
/// @details Some Movers (e.g. repack) require knowledge of the Pose to create,
///  and are only valid for that Pose and other conformations of it
///  (i.e. poses with the same number of residues, jumps, etc).
///  In general, we expect each Pose coming in here to be from a different PDB file.
///  The cost of creating these objects anew should be minimal compared
///  to the processing work they do.
void
LigandRepackMinimizeProtocol::apply( core::pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace core::pack::task;

	//int jump_id = pose.num_jump(); // assume ligand attached by last jump
	//if ( jump_id == 0 ) {
	// utility_exit_with_message("Pose has no jumps!");
	//}

	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	// Scoring function already set up by superclass
	// BUG:  currently need to score pose explicitly to get everything initialized properly
	(*scorefxn_)( pose );

	// Repack all sidechains
	PackerTaskOP pack_task = core::pack::task::TaskFactory::create_packer_task(pose);
	pack_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	pack_task->restrict_to_repacking(); // all residues
	pack_task->or_include_current(true); // may already be in lowest E conf
	pack_task->append_rotamerset_operation( unboundrot_ );
	// Disable packing completely for ligands b/c we want them to stay put
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( !pose.residue(i).is_polymer() ) {
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	//core::pack::rtmin(pose, *scorefxn_, pack_task);
	//return;

	protocols::simple_moves::PackRotamersMoverOP fullRepack( new protocols::simple_moves::PackRotamersMover(scorefxn_, pack_task) );
	fullRepack->apply(pose);

	// Is this necessary?  How well does the repack converge?
	protocols::simple_moves::RotamerTrialsMoverOP rotamerTrials( new protocols::simple_moves::RotamerTrialsMover(scorefxn_, *pack_task) );
	rotamerTrials->apply(pose);

	// Set up move map for minimizing.
	// Only want to minimize protein sc;  keep ligand the same as reference point
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	for ( int i = 1, end_i = pose.size(); i <= end_i; ++i ) {
		if ( pose.residue(i).is_polymer() ) {
			movemap->set_chi(i, true);
		}
	}
	protocols::simple_moves::MinMoverOP dfpMinTightTol( new protocols::simple_moves::MinMover( movemap, scorefxn_, "lbfgs_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/ ) );
	dfpMinTightTol->min_options()->nblist_auto_update(true);
	dfpMinTightTol->apply(pose);

}


//////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		OPT(in::path::database);
		OPT(in::file::extra_res_fa);
		OPT(packing::unboundrot);
		OPT(packing::ex1::ex1);
		OPT(packing::ex1aro::ex1aro);
		OPT(packing::ex2::ex2);
		OPT(packing::extrachi_cutoff);
		OPT(packing::no_optH);
		OPT(packing::flip_HNQ);
		OPT(docking::ligand::soft_rep);
		OPT(docking::ligand::old_estat);

		OPT(in::file::s);
		OPT(out::nstruct);
		OPT(out::path::pdb);

		// Parses command line options and inits RNG.
		// Need to do this before trying to check options.
		// Doesn't seem to hurt to do it again if already done once (?)
		devel::init(argc, argv);

		// Build overall docking protocol Mover
		LigandRepackMinimizeProtocolOP dockingProtocol( new LigandRepackMinimizeProtocol() );

		protocols::jobdist::main_plain_pdb_mover(*dockingProtocol, dockingProtocol->get_scorefxn());

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

