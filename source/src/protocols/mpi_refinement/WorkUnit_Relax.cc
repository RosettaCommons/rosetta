// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mpi_refinement/WorkUnit_Sampler.cc
/// @brief
/// @author Mike Tyka
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

//#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/mpi_refinement/WorkUnit_Relax.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh> // strand pairings
#include <basic/options/option.hh>

//#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
//
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/LoopHashSampler.hh>
//
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
//
#include <protocols/md/CartesianMD.hh>
#include <protocols/relax/FastRelax.hh>
//
//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>
#include <utility> //for std::pair

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

// Contains mover - bbG, MD, Relax

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("WorkUnit_Sampler.Relax");

////////////////////////////////////////////
//////// WorkUnit bbGauss
WorkUnit_bbGauss::WorkUnit_bbGauss( core::Size const nstruct,
	core::Real const kT,
	bool const centroid,
	bool const on_defined_segment )
//: WorkUnit_SilentStructStore()
{
	set_defaults();
	set_nstruct( nstruct );
	set_kT( kT );

	if ( centroid ) {
		set_centroid( 1 );
	} else {
		set_centroid( 0 );
	}

	if ( on_defined_segment ) {
		set_segdef( 1 );
	} else {
		set_segdef( 0 );
	}
}

void
WorkUnit_bbGauss::set_defaults(){}

void
WorkUnit_bbGauss::run()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		//TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	core::Size const nstep_store( 5000 ); // Hard-coded

	SilentStructCOP start_struct = decoys().get_struct(0);
	core::pose::Pose pose;
	decoys().get_pose( 0, pose );

	// clear the store of structures
	decoys().clear();

	// scmc setup
	core::pack::task::PackerTaskOP pt
		= core::pack::task::TaskFactory::create_packer_task( pose );

	protocols::simple_moves::sidechain_moves::SidechainMCMover scmc;

	scmc.set_task( pt );
	pt->restrict_to_repacking();

	protocols::moves::MoverOP tofa
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	bool centroid = (get_centroid() == 1) ? true : false;

	if ( centroid ) tocen->apply( pose );
	core::scoring::ScoreFunctionCOP sfxn_loc = centroid ?
		core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth" ) :
		core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );

	scmc.init_task( pose );
	scmc.set_ntrials( 100 );
	scmc.set_prob_uniform( 0.0 );
	scmc.set_prob_withinrot( 0.0 );
	scmc.set_prob_random_pert_current( 0.1 );
	scmc.set_preserve_detailed_balance( false );
	scmc.set_temperature( get_kT() );
	scmc.set_scorefunction( *sfxn_loc );
	scmc.setup( sfxn_loc );

	// Setup MC / bbgmover
	protocols::moves::MonteCarlo mc( pose, *sfxn_loc, get_kT() );

	// I don't know why yet, but using bbgmover as local class
	// breaks run; let's initiate it temporarily every time
	protocols::simple_moves::BBG8T3AMover bbgmover;

	core::kinematics::MoveMapOP mm_local( new core::kinematics::MoveMap );
	mm_local->set_jump( true );
	mm_local->set_bb( true );
	mm_local->set_chi( true );
	bbgmover.movemap( mm_local );

	// Minimizer
	core::optimization::MinimizerOptions minoption( "lbfgs_armijo_nonmonotone",
		0.001, true, false, false );
	minoption.max_iter( 50 );
	core::optimization::AtomTreeMinimizer minimizer;

	// Run MC
	core::Size nstruct( get_nstruct() );

	//TR << "Executing WorkUnit_bbGauss_Mover on ssid " << start_struct->get_energy("ssid")
	// << " with report_step " << nstep_store;
	// << " at kT = " << get_kT() << std::endl;
	// << std::endl;

	//core::Size starttime = time(NULL);
	core::Size istep( 0 );
	while ( true ) {
		istep++;
		core::Real prob = numeric::random::rg().uniform();
		core::Real proposal_density_ratio( 1.0 );
		std::string movetype;

		if ( prob > 0.0 ) {
			bbgmover.apply( pose );
			movetype = bbgmover.type();
			proposal_density_ratio = bbgmover.last_proposal_density_ratio();

		} else {
			scmc.apply( pose );
			movetype = "sc";
		}

		mc.boltzmann( pose, movetype, proposal_density_ratio );

		if ( istep%nstep_store == 0 ) {
			// Run short minimization before storing?

			core::pose::Pose pose_tmp( pose );
			if ( pose_tmp.is_centroid() ) tofa->apply( pose_tmp );

			store_to_decoys( start_struct, pose_tmp );
		}

		if ( decoys().size() >= nstruct ) break;
	}

	//core::Size endtime = time(NULL);
	//TR.Debug << "Build " << decoys().size() << " structures in ";
	//TR.Debug << endtime - starttime << " s " << std::endl;

}

////////////////////////////////////////////
//////// WorkUnit MD
WorkUnit_MD::WorkUnit_MD( core::Size const relaxtype,
	core::Size const scoretype,
	core::Size const nstruct,
	core::Real const cstweight,
	bool const looponly )
//: WorkUnit_SilentStructStore()
{
	set_defaults();

	set_relaxtype( relaxtype );
	set_scoretype( scoretype );
	set_nstruct( nstruct );
	set_cstweight( cstweight );
	if ( looponly ) {
		set_options( "looponly" );
	} else {
		set_options( "no" );
	}
}

void
WorkUnit_MD::set_defaults() {}

void
WorkUnit_MD::run()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		//TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	SilentStructCOP start_struct = decoys().get_struct(0);
	core::pose::Pose pose;
	decoys().get_pose( 0, pose );
	decoys().clear();

	core::Size const nstep_store( 500 ); // Hard-coded!

	//TR << "Executing WorkUnit_MD_Mover with relaxtype " << get_relaxtype()
	// << ", scoretype " << get_scoretype();
	// << ", nsteps " << get_nstruct()*nstep_store
	//  << std::endl;
	//core::Size starttime = time(NULL);

	bool partial_sampling = false;
	bool centroid = false;

	//std::string mdoption = get_options();
	std::string mdoption("");
	if ( mdoption.compare("cen") == 0 ) {
		centroid = true;
	} else if ( mdoption.compare( "looponly" ) == 0 ) {
		partial_sampling = true;
	}

	protocols::moves::MoverOP tofa
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	// Starts here
	// 1. Score function setup
	bool softpack = ( get_relaxtype() == 2 )? true: false;

	core::Real cst_weight( 0.0 );
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		cst_weight = get_cstweight();
		//TR << "Applying cst_fa_file given by user with weight " << get_cstweight() << std::endl;
	}

	std::string sfxn_name("");
	if ( centroid ) {
		sfxn_name = "cen_cart";
		softpack = false;
	} else if ( get_scoretype() == 1 ) { // facts
		sfxn_name = "scorefacts_cart";
	} else {
		sfxn_name = "talaris2013_cart";
	}

	core::scoring::ScoreFunctionOP sfxn_sampling = get_energy( sfxn_name,
		softpack, cst_weight );

	core::scoring::ScoreFunctionOP sfxn_pack = get_energy( "talaris2013_cart", softpack );

	// 2. MD setup
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	if ( partial_sampling ) {
		mm = get_movemap( pose, "looponly", true );
	} else {
		mm = get_movemap( pose, "full", true );
	}

	if ( sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) > 0.0 ) {
		TR << "Sampling with elec_dens_fast : " << sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) << std::endl;
	}

	protocols::md::CartesianMD MD( pose, sfxn_sampling, mm );

	core::Size neqstep( 2 );  // first 1ps is for eq and will be removed
	MD.set_store_trj( true );
	MD.set_nstep( get_nstruct()*nstep_store );
	MD.set_reportstep( nstep_store );
	MD.set_temperature( 150.0 );

	// Turn on if no other cst_fa_file provided
	// BTW, MD instance will take care of applying cst_fa_file...
	if ( !basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() &&
			get_cstweight() > 0.0 ) {
		//TR << "Applying uniform coord_constraint with weight " << get_cstweight() << std::endl;
		core::Real stdev = std::sqrt( 1.0/get_cstweight() );
		MD.set_constraint( stdev );
	}

	core::pose::Pose pose_work( pose );
	if ( centroid ) tocen->apply( pose_work );

	sfxn_sampling->score( pose );

	// 3. Run!
	repack( pose, sfxn_pack );  // repack before MD

	MD.apply( pose_work );
	utility::vector1< core::pose::Pose > poses_out = MD.dump_poses( pose );

	for ( Size i = 1+neqstep; i <= poses_out.size(); ++i ) {
		core::pose::Pose pose_out( poses_out[i] );
		if ( pose_out.is_centroid() ) tofa->apply( pose_out );

		store_to_decoys( start_struct, pose_out, "_"+string_of( i ) );

		if ( decoys().store().size() >= get_nstruct() ) break;
	}

	//core::Size endtime = time(NULL);
	//TR.Debug << "Build " << decoys().size() << " structures in ";
	//TR.Debug << endtime - starttime << " s " << std::endl;

	// Revert parameters if necessary
	if ( get_scoretype() == 1 ) revert_facts_params();

} // WorkUnit MD

////////////////////////////////////////////
//////// WorkUnit Relax
WorkUnit_Relax::WorkUnit_Relax( core::Size const relaxtype,
	core::Size const scoretype,
	core::Size const nrepeat,
	core::Real const cstweight
)
//: WorkUnit_SilentStructStore()
{
	set_defaults();
	set_relaxtype( relaxtype );
	set_scoretype( scoretype );
	set_nrepeat( nrepeat );
	set_cstweight( cstweight );
}

void
WorkUnit_Relax::set_defaults(){}

std::vector< std::string >
WorkUnit_Relax::set_relax_schedule() const
{
	// Current setup:

	// 0: Null
	// 1: CartStd, 2 cycles
	// 2: DualStd, 5 cycles
	// 3: CartProb, 2 cycles
	// 4: DualProb, 5 cycles
	// 5: DualProb, 2 cycles
	// 6: TorsStd, 2 cycles
	// 7: CartStd-fixed, 2 cycles
	// 8: CartProb-fixed, 2 cycles
	// 9: DualProb-fixed, 5 cycles

	std::vector< std::string > cmdlines;

	if ( get_relaxtype() ==  0 ) return cmdlines;

	// Special schedule for rerelax-style
	if ( get_relaxtype() == 10 ) { //prv 5
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 0.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 11 ) { // veryshort, prv 6
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 0.0 10");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 12 ) { // prv7
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 1.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 13 ) { //prv 8
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 10.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 14 ) { //prv 9
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001  0.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 15 ) {
		// special relax, ramping but restricting CA, soft pack
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "repeat 1" );
		cmdlines.push_back( "scale:coordinate_constraint 10.0" );
		cmdlines.push_back( "scale:fa_rep 0.01" );
		cmdlines.push_back( "scale:fa_atr 0.0" );
		cmdlines.push_back( "repack" );
		cmdlines.push_back( "scale:fa_rep 1.0" );
		cmdlines.push_back( "scale:fa_atr 1.0" );
		cmdlines.push_back( "min 0.01 50" );
		cmdlines.push_back( "ramp_repack_min 0.250 0.01    10.0  50");
		cmdlines.push_back( "ramp_repack_min 0.550 0.01    10.0 100");
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 10.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 16 ) {
		cmdlines.push_back( "switch:torsion" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.001    1.0 100");
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001  1.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	} else if ( get_relaxtype() == 17 ) {
		cmdlines.push_back( "switch:torsion" );
		cmdlines.push_back( "ramp_repack_min 0.02  0.001    1.0  50");
		cmdlines.push_back( "ramp_repack_min 0.25  0.001    0.5  50");
		cmdlines.push_back( "ramp_repack_min 0.55  0.001    0.0 100");
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001  0.0 200");
		cmdlines.push_back( "switch:cartesian" );
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001  0.0 200");
		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
		return cmdlines;
	}

	// Torsion first
	if ( get_relaxtype() == 2 || get_relaxtype() == 4 || get_relaxtype() == 5
			|| get_relaxtype() == 6 || get_relaxtype() == 9 ) {
		cmdlines.push_back( "switch:torsion" );
		if ( get_relaxtype() == 5 || get_relaxtype() == 6 ) {
			cmdlines.push_back( "repeat 2" );
		} else {
			cmdlines.push_back( "repeat 3" );
		}
		// iter1: stdpack or probpack
		if ( get_relaxtype() == 4 || get_relaxtype() == 5 ) { //probpack
			cmdlines.push_back( "scale:coordinate_constraint 1.0" );
			cmdlines.push_back( "scale:fa_rep 0.01" );
			cmdlines.push_back( "scale:fa_atr 0.0" );
			cmdlines.push_back( "repack" );
			cmdlines.push_back( "scale:fa_rep 1.0" );
			//cmdlines.push_back( "scale:fa_rep 0.02" );
			cmdlines.push_back( "scale:fa_atr 1.0" );
			cmdlines.push_back( "min 0.01 50" );

		} else { //stdpack
			cmdlines.push_back( "ramp_repack_min 0.02  0.01    1.0  50");
		}

		cmdlines.push_back( "accept_to_best" );
		cmdlines.push_back( "endrepeat" );
	}

	if ( get_relaxtype() == 6 ) return cmdlines;

	// Add cartesian
	cmdlines.push_back( "switch:cartesian" );
	if ( get_relaxtype() == 5 ) {
		cmdlines.push_back( "repeat 2" );
	} else {
		cmdlines.push_back( "repeat 1" );
	}

	if ( get_relaxtype() == 1 || get_relaxtype() == 2 || get_relaxtype() == 7 ) { //stdpack
		cmdlines.push_back( "ramp_repack_min 0.02  0.01    1.0  50");

	} else if ( get_relaxtype() == 3 || get_relaxtype() == 4 ||
			get_relaxtype() == 5 || get_relaxtype() == 8 || get_relaxtype() == 9 ) { // probpack
		cmdlines.push_back( "scale:coordinate_constraint 1.0" );
		cmdlines.push_back( "scale:fa_rep 0.01" );
		cmdlines.push_back( "scale:fa_atr 0.0" );
		cmdlines.push_back( "repack" );
		cmdlines.push_back( "scale:fa_rep 1.0" );
		//cmdlines.push_back( "scale:fa_rep 0.02" );
		cmdlines.push_back( "scale:fa_atr 1.0" );
		cmdlines.push_back( "min 0.01 50" );
	}

	// fixed relax
	if ( get_relaxtype() == 7 || get_relaxtype() == 8 || get_relaxtype() == 9 ) {
		cmdlines.push_back( "ramp_repack_min 0.250 0.01    0.5  50");
		cmdlines.push_back( "ramp_repack_min 0.550 0.01    0.1 100");
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 0.1 200");
	} else { // normal relax
		cmdlines.push_back( "ramp_repack_min 0.250 0.01    0.5  50");
		cmdlines.push_back( "ramp_repack_min 0.550 0.01    0.0 100");
		cmdlines.push_back( "ramp_repack_min 1.0   0.00001 0.0 200");
	}

	cmdlines.push_back( "accept_to_best" );
	cmdlines.push_back( "endrepeat" );

	return cmdlines;
}

void
WorkUnit_Relax::run()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		//TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	// get a copy, clear the store of structures
	protocols::wum::SilentStructStore decoys_in = decoys();
	decoys().clear();

	//TR << "Executing WorkUnit_Relax_Mover on ndecoys " << decoys_in.size()
	//  << ", relaxtype " << get_relaxtype()
	//  << ", scoretype " << get_scoretype()
	//  <<", cstweight " << get_cstweight()
	//  << std::endl;

	//core::Size starttime = time(NULL);

	// Start here
	// 1. Score function setup
	std::string sfxn_name("");
	if ( get_scoretype() == 1 ) {
		sfxn_name = "scorefacts_cart";
	} else {
		sfxn_name = "talaris2013_cart";
	}

	core::scoring::ScoreFunctionOP sfxn_sampling = get_energy( sfxn_name, false, get_cstweight() );

	if ( sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) > 0.0 ) {
		TR << "Sampling with elec_dens_fast : " << sfxn_sampling->get_weight( core::scoring::elec_dens_fast ) << std::endl;
	}

	core::pose::Pose pose;

	// Relax schedule setup
	protocols::relax::FastRelax relax( sfxn_sampling );
	//relax.cst_calpha_only( true ); // only calpha; default is all backbones
	//relax.set_movemap( mm_ );

	std::vector< std::string > const script_lines = set_relax_schedule();
	relax.set_script_from_lines( script_lines );
	relax.min_type( "lbfgs_armijo_nonmonotone" );

	// Run
	for ( core::Size i = 0; i < decoys_in.size(); ++i ) {
		SilentStructCOP start_struct = decoys_in.get_struct(0);
		decoys_in.get_pose( i, pose );

		// call user-defined cstfile
		if ( option[ constraints::cst_fa_file ].user() ) {
			//TR << "Applying constraints from cmd" << std::endl;
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
		}

		for ( core::Size j = 1; j <= get_nrepeat(); ++j ) {
			core::pose::Pose pose_work( pose );
			relax.apply( pose_work );

			store_to_decoys( start_struct, pose_work );
		}
	}


	//core::Size endtime = time(NULL);
	//TR.Debug << "Build " << decoys().size() << " structures in ";
	//TR.Debug << endtime - starttime << " s " << std::endl;

	// Revert parameters
	if ( get_scoretype() == 1 ) revert_facts_params();

} // WorkUnit_Relax


}
}

