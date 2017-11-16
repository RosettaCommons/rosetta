// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test_bbmc.cc
/// @brief run Monte Carlo for sampling protein conformation
/// @author Yuan Liu
/// @details
/// Modified from Colin's backrub.cc, put all backbone algorithm(backrub, bbg, conrot)
/// and all sidechain algorithm(sc, scmc), and MonteCarlo/ReplicaExchange together

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ReplicaExchangeMC.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <basic/options/option.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <basic/options/option_macros.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

//Backbone Gaussian Mover
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/BBConRotMover.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/exit.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/MultiDimensionalHistogram.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/io/mpistream.hh>
#include <fstream>
#include <string>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/prof.hh>

static basic::Tracer TR( "pilot.wendao.bbmc" );

//params for all
OPT_1GRP_KEY(Integer, mc, ntrials) //how many steps
OPT_1GRP_KEY(Real, mc, sm_prob) //prob of smallmover
OPT_1GRP_KEY(Real, mc, sm_angle_max)
OPT_1GRP_KEY(Real, mc, backrub_prob) //prob of backrubmover
OPT_1GRP_KEY(Real, mc, conrot_prob) //prob of conrotmover
OPT_1GRP_KEY(Real, mc, kt) //temperatrue
OPT_1GRP_KEY(Real, mc, mm_bend_weight) //mm bend energy
OPT_1GRP_KEY(Boolean, mc, detailed_balance) //detailed balance correction
OPT_1GRP_KEY(Boolean, mc, initial_pack) //for packer
OPT_1GRP_KEY(File, mc, movemap) //movemap
OPT_1GRP_KEY(Integer, mc, rmsd_region_start) //beginning of the rmsd region
OPT_1GRP_KEY(Integer, mc, rmsd_region_stop) //end of the rmsd region
OPT_1GRP_KEY(File, mc, minimize_movemap) //movemap
OPT_1GRP_KEY(RealVector, mc, trajectory_tlist)
OPT_1GRP_KEY(Integer, mc, score_stride) //score only silent file
OPT_1GRP_KEY(Integer, mc, trajectory_stride)
OPT_1GRP_KEY(Integer, mc, output_stride)
OPT_1GRP_KEY(Integer, mc, cluster_ndx)
OPT_1GRP_KEY(Boolean, mc, centroid)
OPT_1GRP_KEY(Boolean, mc, noscore)

//replica
OPT_1GRP_KEY(Boolean, mc, replica)
OPT_1GRP_KEY(String, mc, re_pdb_prefix)
OPT_1GRP_KEY(Boolean, mc, re_pdb_suffix)
OPT_1GRP_KEY(Integer, mc, re_ninterval)
OPT_1GRP_KEY(RealVector, mc, re_tlist)

//params for sidechain
OPT_1GRP_KEY(Real, mc, sc_prob)
OPT_1GRP_KEY(Real, mc, sc_prob_uniform)
OPT_1GRP_KEY(Real, mc, sc_prob_withinrot)
OPT_1GRP_KEY(Real, mc, sc_prob_random_pert_current)
OPT_1GRP_KEY(Boolean, mc, fast_sc)
OPT_1GRP_KEY(Real, mc, fast_sc_prob)//two ways for calling fastsc
OPT_1GRP_KEY(Boolean, mc, sc_strategy2)
OPT_1GRP_KEY(Boolean, mc, fast_sc_strategy2)
OPT_1GRP_KEY(Integer, mc, sc_ntrials)
OPT_1GRP_KEY(IntegerVector, mc, sc_statistic)
OPT_1GRP_KEY(Integer, mc, bb_dih_statistic)
OPT_1GRP_KEY(String,mc,restart_from_silent)

OPT_1GRP_KEY(Real, mc, near_native_threshold )
OPT_1GRP_KEY(Boolean, mc, follow_classic_naming_convention )
OPT_1GRP_KEY(String,mc, movable_segment)

void *my_main( void* );

int main( int argc, char * argv [] )
{

	//all
	NEW_OPT(mc::kt, "value of kT for Monte Carlo", 0.56);
	NEW_OPT(mc::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(mc::trajectory_stride, "write out a trajectory frame every N steps", 0);
	//sc
	NEW_OPT(mc::sc_prob, "probability of making a side chain move", 0);
	NEW_OPT(mc::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.0);
	NEW_OPT(mc::sc_prob_withinrot, "probability of sampling within the current rotamer", 0.0);
	NEW_OPT(mc::sc_prob_random_pert_current, "probability of sampling within the current rotamer", 0.1);
	NEW_OPT(mc::sc_ntrials, "fast sidechainmover(scmc)'s internal sc move", 100 );
}


void *
my_main( void* )
{
	using namespace core;
	using namespace core::io::silent;
	using namespace core::io::pdb;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::moves;

	// create a TaskFactory with the resfile
	using namespace core::pack::task;
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	} else {
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( new operation::PreserveCBeta );


	//setup the BBGMover
	protocols::simple_moves::BBG8T3AMover bbgmover;

	//setup Movemap
	if ( option[ mc::movemap ].user() ) {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->init_from_file(option[ mc::movemap ]);
		bbgmover.movemap(movemap);
	}

	//setup switch mover
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");

	//setup native pose for rmsd
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	PoseOP native_pose;

	//setup starting pose
	PoseOP pose = new Pose();
	Pose &p(*pose);

	//if constraints are specified, add them!
	if ( option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		core::scoring::constraints::add_constraints_from_cmdline( p, *score_fxn);
	}


	//setup Monte Carlo
	MonteCarloOP mc;
	core::Real kT = option[ mc::kt ];
	mc = new MonteCarlo(p, *score_fxn, kT);

	//setup SidechainMC mover
	protocols::simple_moves::sidechain_moves::SidechainMCMover scmc;

	// Sidechain
	pack::task::PackerTaskOP pt = core::pack::task::TaskFactory::create_packer_task( p );
	scmc.set_task( pt );
	pt->restrict_to_repacking();
	scmc.init_task( p );
	scmc.set_ntrials( option[mc::sc_ntrials] );
	TR << "sc_ntrials are " << scmc.ntrials() << std::endl;

	scmc.set_prob_uniform( option[ mc::sc_prob_uniform ] );
	scmc.set_prob_withinrot( option[ mc::sc_prob_withinrot ] );
	scmc.set_prob_random_pert_current( option[ mc::sc_prob_random_pert_current ] );
	scmc.set_preserve_detailed_balance( true );
	scmc.set_temperature( kT ); //only for intra mc criteria
	scmc.set_scorefunction( *scfxn );
	scmc.setup( scfxn );

	Pose start_pose(p); //save the starting structure so we can compare by rmsd later ek 1-21-2011

	Size ntrials = option[ mc::ntrials ];

	for ( Size i = 1; i <= ntrials; ++i ) {
		//init
		string move_type( "" );
		Real proposal_density_ratio=1.0;
		//random number
		core::Real prob = numeric::random::rg().uniform();

		//choose one
		if ( prob > 0.75 ) { //bbg
			bbgmover.apply(p);
		} else {
			scmc.apply(p);
		}

		//sync temperature
		//reset sidechainmover's temperature
		//if (option[mc::fast_sc_strategy2]) scmc.set_sampling_temperature(mc->temperature());
		scmc.set_temperature(mc->temperature()); //should be
	}
	return 0;
}


