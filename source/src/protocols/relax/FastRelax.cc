// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/FastRelax.cc
/// @brief The FastRelax Protocol
/// @details
/// @author Mike Tyka
/// @author Roland A. Pache
/// @author Jared Adolf-Bryfogle
/*
Format for relax script file

repeat <loopcount>

starts a section to be repeated <loopcount> times - end the block with "endrepeat".
Nested loop are NOT supported.

endrepeat

is the end marker of a loop. loops may *NOT* be nested! (this is not a programming language. well it is. but a very simple one.)

accept_to_best

compares the energy of the current pose and the best_pose and replaces the best pose
with the current pose if it's energy is lower.

load_best

Load the best_pose into the current working pose, replacing it.

load_start

Load the starting pose into the current working pose, replacing it.

Sets the best_pose to the current pose whatever the energy of the current pose.

dump <number>

DUmps a pdbfile called dump_<number>.pdb

dumpall:<true/false>

Will dump a pdb for every repack/min/ramp_repack_min hereafter with an incrementing number

scale:<scoretype> <scale>

Scales the scoretype's default value by whatever the scale number is.
For ex, scale:fa_rep 0.1 will take the original fa_rep weight and set it
to 0.1 * original weight.

rscale:<scoretype> <lower limit> <upper limit>

Like scale, but picks a random number between the limits to multiply by

weight:<scoretype> <weight>

Sets the weight of the scoretype to whatever the weight number is.
ALSO CHANGES DEFAULT WEIGHT.  This is so in weight, scale, scale routines
the scale will be using the user-defined weight, which seems to make more sense.
For ex, weight:fa_rep 0.2 will take set fa_rep to 0.2

show_weights

Outputs the current weights.  If a parameter is not outputted, then its weight
is 0. Most useful when redirecting stdout to a file, and only one input structure.

coord_cst_weight <scale>

Sets the coordinate_constraint weight to <scale>*start_weight.
This is used when using the commandline options -constrain_relax_to_native_coords and
-constrain_relax_to_start_coords.

switch:<torsion/cartesian>

Switches to the torsion minimizer (AtomTreeMinimzer, dfp minimizer) or the cartesian
minimizer (CartesianMinimizer, lbfgs).  Obviously ignores command line flags.

repack

Triggers a full repack

min <tolerance>

Triggers a dfp minimization with a given tolernace (0.001 is a decent value)

ramp_repack_min <scale:fa_rep> <min_tolerance> [ <coord_cst_weight = 0.0> ]

Causes a typical ramp, then repack then minimize cycle - i.e. bascially combines
the three commands scale:fa_rep, repack and min (and possible coord_cst_weight)

These two command scripts are equivalent:
--------
scale:fa_rep 0.5
repack
min 0.001
--------
and

--------
ramp_repack_min 0.5 0.001
--------

batch_shave <keep_proportion>

Valid for batchrelax only - ignored in normal FastRelax
In Batchrelax it will remove the worst structures from the current pool and leave only
<keep_proportion>. I.e. the command

batch_shave 0.75

Will remove the worst 75% of structures from the current working pool.


exit

will quit with immediate effect


A typical FastRelax command_script is: (thi sin fact is the default command script)

repeat 5
ramp_repack_min 0.02  0.01     1.0
ramp_repack_min 0.250 0.01     0.5
ramp_repack_min 0.550 0.01     0.0
ramp_repack_min 1     0.00001  0.0
accept_to_best
endrepeat


*/
//Project Headers
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/Ramady.hh>
#include <protocols/relax/FastRelaxCreator.hh>
#include <protocols/relax/util.hh>
#include <protocols/relax/cst_util.hh>

//Core Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Protocol Headers
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/md/CartesianMD.hh>
#include <protocols/hydrate/Hydrate.hh>

//Basic Headers
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.relax.FastRelax" );

using namespace core;
using namespace core::io::silent;

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {
////////////////////////////////////////////////////////////////////////////////////////////////////


using namespace ObjexxFCL;


// XRW TEMP std::string
// XRW TEMP FastRelaxCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FastRelax::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FastRelaxCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FastRelax() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FastRelax::mover_name()
// XRW TEMP {
// XRW TEMP  return "FastRelax";
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  FastRelax main code:
///  ---------------------------------------------------------------------------------

FastRelax::FastRelax(
	core::Size                     standard_repeats
) :
	RelaxProtocolBase( std::string("FastRelax" ) ),
	checkpoints_( std::string("FastRelax") )
{
	set_to_default();
	if ( standard_repeats == 0 ) standard_repeats = default_repeats_;
	if ( explicit_ramp_constraints() && ! ramp_down_constraints() ) {
		read_script_file( "NO CST RAMPING", standard_repeats );
	} else {
		read_script_file( "", standard_repeats );
	}
}


FastRelax::FastRelax(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::Size                     standard_repeats
) :
	RelaxProtocolBase("FastRelax", scorefxn_in ),
	checkpoints_("FastRelax")
{
	set_to_default();
	if ( standard_repeats == 0 ) standard_repeats = default_repeats_;
	if ( explicit_ramp_constraints() && ! ramp_down_constraints() ) {
		read_script_file( "NO CST RAMPING", standard_repeats );
	} else {
		read_script_file( "", standard_repeats );
	}
}


FastRelax::FastRelax(
	core::scoring::ScoreFunctionOP scorefxn_in,
	const std::string &          script_file
) :
	RelaxProtocolBase("FastRelax", scorefxn_in ),
	checkpoints_("FastRelax")
{
	set_to_default();
	read_script_file( script_file );
}

FastRelax::FastRelax(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::Size                     standard_repeats,
	const std::string &          script_file
) :
	RelaxProtocolBase("FastRelax", scorefxn_in ),
	checkpoints_("FastRelax")
{
	set_to_default();
	if ( standard_repeats == 0 ) standard_repeats = default_repeats_;
	read_script_file( script_file , standard_repeats );
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FastRelax::~FastRelax() = default;


/// Return a copy of ourselves
protocols::moves::MoverOP
FastRelax::clone() const {
	return protocols::moves::MoverOP( new FastRelax(*this) );
}

void
FastRelax::dualspace(bool val){
	dualspace_ = val;
}

bool
FastRelax::dualspace() {
	return dualspace_;
}

void
FastRelax::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	using namespace basic::options;
	using core::pack::task::operation::TaskOperationCOP;

	bool nonideal = false; // Are we setup for non-ideal Relax?

	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data )->clone() );

	//Make sure we have a taskfactory before we overwrite our null in the base class.
	set_enable_design( !( tag->getOption<bool>("disable_design", !enable_design_ ) ) );

	core::pack::task::TaskFactoryOP tf = protocols::rosetta_scripts::parse_task_operations( tag, data );
	if ( tf->size() > 0 ) {
		if ( !enable_design_ ) {
			tf->push_back(TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ));
		}
		set_task_factory( tf );
	}

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi( true );
	mm->set_bb( true );
	mm->set_jump( true );

	if ( tag->getOption< bool >( "bondangle", false ) ) {
		// minimize_bond_angles( true ); // Don't set here, to avoid warning later
		nonideal = true;
		mm->set( core::id::THETA, true );
	}

	if ( tag->getOption< bool >( "bondlength", false ) ) {
		// minimize_bond_lengths( true ); // Don't set here, to avoid warning later
		nonideal = true;
		mm->set( core::id::D, true );
	}

	set_movemap(mm); // Set the base movemap - The MoveMapFactories will modify this

	core::select::movemap::MoveMapFactoryOP mmf_legacy = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false );
	core::select::movemap::MoveMapFactoryOP mmf = core::select::movemap::parse_movemap_factory( tag, data );

	// If we've used the newer MoveMapFactory syntax, use that - otherwise use the legacy syntax.
	if ( mmf ) {
		set_movemap_factory( mmf );
		if ( mmf_legacy ) {
			TR.Warning << "Unintended behavior may result from using both a MoveMapFactory and a MoveMap in the " << tag->getOption< std::string >( "name", "(unnamed FastRelaxMover)" ) << " FastRelaxMover; the MoveMapFatory has the potential to overwrite the settings provided by MoveMap" << std::endl;
		}
	} else if ( mmf_legacy ) {
		set_movemap_factory( mmf_legacy );
	} // else don't bother setting mmf

	default_repeats_ = tag->getOption< int >( "repeats", option[ OptionKeys::relax::default_repeats ]() );
	std::string script_file("");
	if ( tag->hasOption("relaxscript") ) {
		script_file = tag->getOption< std::string >("relaxscript");
		script_file_specified_ = true;
	}

	// Only support single file; is there a way to support multiple input files in rosetta scripts?
	std::string cstfile = tag->getOption< std::string >("cst_file", "");
	if ( cstfile != "" ) add_cst_files( cstfile );

	bool batch = tag->getOption< bool >( "batch", false );

	cartesian (tag->getOption< bool >( "cartesian", cartesian() ) );
	dualspace (tag->getOption< bool >( "dualspace", dualspace() ) );

	ramp_down_constraints( tag->getOption< bool >( "ramp_down_constraints", ramp_down_constraints() ) );


	if ( tag->hasOption( "min_type" ) ) {
		min_type( tag->getOption< std::string >( "min_type" ) );
	} else {
		//fpd if no minimizer is specified, and we're doing flexible bond minimization
		//fpd   we should use lbfgs (dfpmin is way too slow)
		if ( cartesian() || nonideal ) {
			min_type( "lbfgs_armijo_nonmonotone" );
		}
	}

	if ( batch ) {
		set_script_to_batchrelax_default( default_repeats_ );
	} else if ( script_file == "" ) {
		read_script_file( "", default_repeats_ );
	} else {
		read_script_file( script_file );
	}

	delete_virtual_residues_after_FastRelax_ = tag->getOption<bool>("delete_virtual_residues_after_FastRelax", false);
	/// @brief if one uses constrain_relax_to_native_coords/constrain_relax_to_start_coords for the FastRelax,
	//      virtual residues need to be deleted to calculate rmsd after FastRelax

} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////////////////////////
void FastRelax::set_to_default( )
{
	using namespace basic::options;

	default_repeats_ = basic::options::option[ OptionKeys::relax::default_repeats ]();
	ramady_ = basic::options::option[ OptionKeys::relax::ramady ]();
	repack_ = basic::options::option[ OptionKeys::relax::chi_move]();
	test_cycles_ = basic::options::option[ OptionKeys::run::test_cycles ]();
	script_max_accept_ = basic::options::option[ OptionKeys::relax::script_max_accept ]();
	symmetric_rmsd_ = option[ basic::options::OptionKeys::evaluation::symmetric_rmsd ]();

	force_nonideal_ = false;

	dna_move_ = option[ basic::options::OptionKeys::relax::dna_move]();

	//fpd additional ramady options
	ramady_num_rebuild_ = basic::options::option[ OptionKeys::relax::ramady_max_rebuild ]();
	ramady_cutoff_ = basic::options::option[ OptionKeys::relax::ramady_cutoff ]();
	ramady_force_ = basic::options::option[ OptionKeys::relax::ramady_force ]();
	ramady_rms_limit_ = basic::options::option[ OptionKeys::relax::ramady_rms_limit ]();
	dualspace_ = basic::options::option[ basic::options::OptionKeys::relax::dualspace ]();
	script_file_specified_ = false;

	// cartesian

	// dumpall
	dumpall_ = false;

	// design
	enable_design_ = false;

	// delete_virtual_residues_after_FastRelax_ = 'false' by default in this constructor, since parse_my_tag in rosetta_scripts alone is not enough to specifically designate default value for a user who uses command-line
	delete_virtual_residues_after_FastRelax_ = false;
} //set_to_default

///Make sure correct mintype for cartesian or dualspace during apply
void FastRelax::check_nonideal_mintype() {
	// This should probably look at the movemap too -- We can set variable lengths/angles in the movemap rather than directly in the Mover
	if ( dualspace() || cartesian() || minimize_bond_lengths() || minimize_bond_angles() ) {
		if ( !basic::options::option[ basic::options::OptionKeys::relax::min_type ].user() ) {
			min_type("lbfgs_armijo_nonmonotone");
		}
	}
}

// grab the score and remember the pose if the score is better then ever before.
void FastRelax::cmd_accept_to_best(
	const core::scoring::ScoreFunctionOP local_scorefxn,
	core::pose::Pose &pose,
	core::pose::Pose &best_pose,
	const core::pose::Pose &start_pose,
	core::Real       &best_score,
	core::Size       &accept_count
)
{
	using namespace core::scoring;
	using namespace core::conformation;
	core::Real score = (*local_scorefxn)( pose );
	if ( ( score < best_score) || (accept_count == 0) ) {
		best_score = score;
		best_pose = pose;
	}
#ifdef BOINC_GRAPHICS
	boinc::Boinc::update_graphics_low_energy( best_pose, best_score  );
	boinc::Boinc::update_graphics_last_accepted( pose, score );
	//boinc::Boinc::update_mc_trial_info( total_count , "FastRelax" );  // total_count not defined
#endif
	core::Real rms = 0, irms = 0;
	if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
		rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
		irms = CA_rmsd_symmetric( start_pose , best_pose );
	} else {
		rms = native_CA_rmsd( *get_native_pose() , best_pose );
		irms = native_CA_rmsd( start_pose , best_pose );
	}
	TR << "MRP: " << accept_count << "  " << score << "  " << best_score << "  " << rms << "  " << irms << "  " << std::endl;
} //cmd_accept_to_best


void FastRelax::do_minimize(
	core::pose::Pose &pose,
	core::Real tolerance,
	core::kinematics::MoveMapOP local_movemap,
	core::scoring::ScoreFunctionOP local_scorefxn
){
	using namespace core::scoring;
	using namespace core::conformation;

	//Why create a new min_mover every time we minimize?
	protocols::minimization_packing::MinMoverOP min_mover;
	min_mover = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( local_movemap, local_scorefxn, min_type(), tolerance, true ) );

	min_mover->cartesian( cartesian() );
	if ( max_iter() > 0 ) min_mover->max_iter( max_iter() );

	min_mover->apply( pose );

	// If relax membrane, update structure based embeddings
	bool membrane = pose.conformation().is_membrane();
	if ( membrane ) {
		//pose.membrane_info()->update_structure_based_embeddings();
	}

}

void FastRelax::do_md(
	core::pose::Pose &pose,
	core::Real nstep,
	core::Real temp0,
	core::kinematics::MoveMapOP movemap_in,
	core::scoring::ScoreFunctionOP local_scorefxn
){

	// UM... this looks like a bug to the extent that you are performing a *shallow* copy
	// of the input movemap.
	core::kinematics::MoveMapOP local_movemap = movemap_in;

	protocols::md::CartesianMD MD_mover( pose, local_scorefxn, local_movemap );

	MD_mover.set_nstep( (Size)(nstep) );
	MD_mover.set_temperature( temp0 );

	MD_mover.apply( pose );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void FastRelax::apply( core::pose::Pose & pose ){
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::kinematics;
	using namespace protocols;
	using namespace basic::options;

	TR.Debug   << "================== FastRelax: " << script_.size() << " ===============================" << std::endl;
	if ( pose.size() == 0 ) {
		TR.Warning << "Pose has no residues. Doing a FastRelax would be pointless. Skipping." << std::endl;
		return;
	}

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "TESTING");
#endif
	// One out of 10 times dont bother doing Ramady Relax. You may wonder why - the reason is that occasionally
	// there are phi-pso pairs outside of the commonly allowed ranges which are *correct* for some sort of quirk
	// of nature. In that case using ramady relax screws things up, so to allow for that chance only execture ramady relax
	// 90% of the time.
	bool do_rama_repair = ramady_;
	if ( !ramady_force_ && numeric::random::uniform() <= 0.1 ) do_rama_repair = false;


	// Relax is a fullatom protocol so switch the residue set here. The rotamers
	// wont be properly packed at this stage but it doesnt really matter.
	// They'll get repacked shortly.
	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
	}

	// Make a local copy of movemap, called local movemap. The reason is that
	// the "constrain to coordinates" codes will wanna mess witht the movemap..
	core::kinematics::MoveMapOP local_movemap = get_movemap()->clone();
	initialize_movemap( pose, *local_movemap );

	//set_movemap(local_movemap);  ///fpd NO!  The whole point of the local copy is that subsequent movemap
	///           changes do not then modify the stored movemap
	///         If you absolutely must do this use clone!
	check_nonideal_mintype();

	// Deal with constraint options and add coodrinate constraints for all or parts of the structure.
	set_up_constraints( pose, *local_movemap );

	// make a copy of the energy function too. SInce we're going to be ramping around with weights,
	// we dont want to modify the existing scorefunction
	ScoreFunctionOP local_scorefxn( get_scorefxn()->clone() );
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, local_scorefxn ); //If the pose is symmetric but the scorefunction is not, remedy this.

	// Remember the oroiginal weights - we're gonna be changing these during the ramp ups/downs
	core::scoring::EnergyMap full_weights = local_scorefxn->weights();

	// Make DNA Rigid or setup DNA-specific relax settings.  Use the orbitals scorefunction when relaxing with DNA
	if ( dna_move_ ) {
		setup_for_dna( *local_scorefxn );
	} else {
		make_dna_rigid( pose, *local_movemap );
	}


	// Make sure we only allow symmetrical degrees of freedom to move and convert the local_movemap
	// to a local movemap
	if ( core::pose::symmetry::is_symmetric( pose )  )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *local_movemap );
	}

	//Change behavior of Task to be initialized in PackRotamersMover to allow design directly within FastRelax
	// Jadolfbr 5/2/2013
	TaskFactoryOP local_tf( new TaskFactory() );

	//If a user gives a TaskFactory, completely respect it.

	if ( get_task_factory() ) {
		local_tf = get_task_factory()->clone();
	} else {
		local_tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));
		if ( option[ OptionKeys::relax::respect_resfile]() && option[ OptionKeys::packing::resfile].user() ) {
			local_tf->push_back(TaskOperationCOP( new ReadResfile() ));
			TR << "Using Resfile for packing step. " <<std::endl;
		} else {
			//Keep the same behavior as before if no resfile given for design.
			//Though, as mentioned in the doc, movemap now overrides chi_move as it was supposed to.

			local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
			PreventRepackingOP turn_off_packing( new PreventRepacking() );
			for ( Size pos = 1; pos <= pose.size(); ++pos ) {
				if ( ! local_movemap->get_chi(pos) ) {
					turn_off_packing->include_residue(pos);
				}
			}
			local_tf->push_back(turn_off_packing);
		}
	}
	//Include current rotamer by default - as before.
	local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));

	if ( limit_aroma_chi2() ) {
		local_tf->push_back(TaskOperationCOP( new task_operations::LimitAromaChi2Operation() ));
	}

	protocols::minimization_packing::PackRotamersMoverOP pack_full_repack_( new protocols::minimization_packing::PackRotamersMover( local_scorefxn ) );
	pack_full_repack_->task_factory(local_tf);

	(*local_scorefxn)( pose );


	/// prints something like this ***1***C***1*********2***C********3****C****2********3*****
	///                            **********xxxxxxxxxxxxx************************************
	if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
		kinematics::simple_visualize_fold_tree_and_movemap_bb_chi( pose.fold_tree(),  *local_movemap, TR );
	}


	/// OK, start of actual protocol
	pose::Pose start_pose=pose;
	pose::Pose best_pose=pose;
	core::Real best_score=100000000.0;
	core::Size accept_count = 0;
	core::Size chk_counter = 0;

	int repeat_step=0;
	int repeat_count=-1;
	int total_repeat_count = 0;


	// Obtain the native pose
	if ( !get_native_pose() ) set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new Pose( start_pose ) ) ) );

	// Statistic information
	std::vector< core::Real > best_score_log;
	std::vector< core::Real > curr_score_log;

	// for dry runs, literally just score the pose and quit.
	if ( dry_run() ) {
		(*local_scorefxn)( pose );
		return;
	}

	// Deal with disulphides - i have no idea what this does, ask Rob Vernon, he put this here.
	apply_disulfides(pose);

	int total_count=0;
	int dump_counter=0;
	for ( core::Size ncmd = 0; ncmd < script_.size(); ncmd ++ ) {
		total_count++;
		if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) local_scorefxn->show( TR, pose );
		if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) pose.constraint_set()->show_numbers( TR.Debug );

		// No MC is used, so update graphics manually
#ifdef BOINC_GRAPHICS
		boinc::Boinc::update_graphics_current( pose );
#endif

		RelaxScriptCommand cmd = script_[ncmd];

		TR.Debug << "Command: " << cmd.command << std::endl;

		if ( cmd.command == "repeat" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "ERROR: Syntax: " + cmd.command + "<number_of_repeats> " ); }
			repeat_count = (int) cmd.param1;
			repeat_step = ncmd;
		} else if ( cmd.command == "endrepeat" ) {
			TR.Debug << "CMD:  Repeat: " << repeat_count << std::endl;
			repeat_count -- ;
			total_repeat_count ++ ;
			if ( repeat_count <= 0 ) {}
			else {
				ncmd = repeat_step;
			}
		} else if ( cmd.command == "dump" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "ERROR: Syntax: " + cmd.command + "<number> " ); }
			pose.dump_pdb( "dump_" + right_string_of( (int) cmd.param1, 4, '0' ) );
		} else if ( cmd.command.substr(0,7) == "dumpall" ) {
			if ( cmd.command.substr(8) == "true" ) {
				dumpall_ = true;
			}
			if ( cmd.command.substr(8) == "false" ) {
				dumpall_ = false;
			}
		} else if ( cmd.command == "repack" ) {
			//if( cmd.nparams < 0 ){ utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }
			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if ( !checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true ) ) {
				// hydrate/SPaDES protocol
				if ( option[ OptionKeys::relax::use_explicit_water ]() ) {
					//if ( option[ OptionKeys::relax::enforce_waters_during_relax ]() == true ) {
					// // keep waters enforced on every ramping step UNTIL final repeat count and fa_rep < enforce_until_ramp_step
					// // below code overwrites the enforce_all_waters option
					// if ( cmd.param1 < option[ OptionKeys::relax::enforce_until_ramp_step ]() || repeat_count != 1 ) {
					//  basic::options::option[ basic::options::OptionKeys::hydrate::force_enforce_all_waters ].value( true );
					// } else { // last ramping step
					//  basic::options::option[ basic::options::OptionKeys::hydrate::force_enforce_all_waters ].value( false );
					// }
					//}
					protocols::hydrate::HydrateOP hydrate_protocol( new protocols::hydrate::Hydrate( local_scorefxn ) );
					hydrate_protocol->apply( pose );
				} else {
					pack_full_repack_->apply( pose );
				}
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}
			if ( dumpall_ ) {
				pose.dump_pdb( "dump_" + right_string_of( dump_counter, 4, '0' ) );
				dump_counter++;
			}
		} else if ( cmd.command == "min" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "ERROR: Syntax " + cmd.command + " <min_tolerance>  " ); }

			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if ( !checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true ) ) {
				do_minimize( pose, cmd.param1, local_movemap, local_scorefxn );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}
			if ( dumpall_ ) {
				pose.dump_pdb( "dump_" + right_string_of( dump_counter, 4, '0' ) );
				dump_counter++;
			}
		} else if ( cmd.command == "md" ) {
			if ( cmd.nparams < 2 ) { utility_exit_with_message( "ERROR: Syntax " + cmd.command + " <mdstep> <temperature> " ); }

			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if ( !checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true ) ) {
				// param1: nstep; param2: temperature
				//if ( input_csts ) pose.constraint_set( input_csts );
				do_md( pose, cmd.param1, cmd.param2, local_movemap, local_scorefxn );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}
		} else if ( cmd.command.substr(0,5) == "scale" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(6));
			local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * cmd.param1 );
		}   else if ( cmd.command.substr(0,6) == "rscale" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
			local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * ((cmd.param2 - cmd.param1 ) * numeric::random::uniform() + cmd.param1 ));

		} else if ( cmd.command.substr(0,9) == "reference" ) {
			// added for design applications
			if ( cmd.nparams < 20 ) {
				utility_exit_with_message( "ERROR: Syntax " + cmd.command + " <20 reference weights in ACDEF order>" );
			} else {
				//local_scorefxn->energy_method_options().set_method_weights(
				// scoring::ScoreTypeManager::score_type_from_name( "ref" ), cmd.params_vec );
				methods::EnergyMethodOptions eopts( local_scorefxn->energy_method_options() );

				eopts.set_method_weights(
					scoring::ScoreTypeManager::score_type_from_name( "ref" ), cmd.params_vec );

				local_scorefxn->set_energy_method_options( eopts );

				/*
				utility::vector1< Real > refw =
				local_scorefxn->energy_method_options().method_weights(
				scoring::ScoreTypeManager::score_type_from_name( "ref" ) );
				TR << "Check, Refweight: " << refw[1] << " " << refw[2] << " " << refw[3]
				<< std::endl;
				*/
			}

		}   else if ( cmd.command.substr(0,6) == "switch" ) {
			// no input validation as of now, relax will just die
			if ( cmd.command.substr(7) == "torsion" ) {
				TR << "Using AtomTreeMinimizer"  << std::endl;
				cartesian( false );
			} else if ( cmd.command.substr(7) == "cartesian" ) {
				TR << "Using CartesianMinizer with lbfgs"  << std::endl;
				cartesian( true );
			}
		}   else if ( cmd.command.substr(0,6) == "weight" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
			local_scorefxn->set_weight( scale_param, cmd.param1 );
			// I'm not too sure if the changing the default weight makes sense
			full_weights[ scale_param ] = cmd.param1;
		}   else if ( cmd.command == "batch_shave" ) {  // grab the score and remember the pose if the score is better then ever before.
		} else if ( cmd.command == "show_weights" ) {
			local_scorefxn->show(TR, pose);
		}   else if ( cmd.command == "ramp_repack_min" ) {
			if ( cmd.nparams < 2 ) { utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }

			// The first parameter is the relate repulsive weight
			local_scorefxn->set_weight( scoring::fa_rep, full_weights[ scoring::fa_rep ] * cmd.param1 );

			// The third paramter is the coordinate constraint weight
			if ( ( constrain_coords() || ramp_down_constraints() ) && (cmd.nparams >= 3) ) {
				set_constraint_weight( local_scorefxn, full_weights, cmd.param3, pose );
			}

			// The fourth paramter is the minimization
			if ( cmd.nparams >= 4 ) {
				auto const iter_cmd = (Size)(cmd.param4);
				max_iter( iter_cmd );
				//Size const original_iter = max_iter();
				//max_iter( original_iter );
			}

			// decide when to call ramady repair code
			if ( total_repeat_count > 1 && repeat_count > 2 ) {
				if ( cmd.param1 < 0.2 ) {
					if ( do_rama_repair ) {
						fix_worst_bad_ramas( pose, ramady_num_rebuild_, 0.0, ramady_cutoff_, ramady_rms_limit_);
					}
				}
			}

			chk_counter++;
			std::string checkpoint_id = "chk" + string_of( chk_counter );
			if ( !checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true ) ) {
				// hydrate/SPaDES protocol
				if ( option[ OptionKeys::relax::use_explicit_water ]() ) {
					//if ( option[ OptionKeys::relax::enforce_waters_during_relax ]() == true ) {
					// // keep waters enforced on every ramping step UNTIL final repeat count and fa_rep < enforce_until_ramp_step
					// // below code overwrites the enforce_all_waters option
					// if ( cmd.param1 < option[ OptionKeys::relax::enforce_until_ramp_step ]() || repeat_count != 1 ) {
					//  basic::options::option[ basic::options::OptionKeys::hydrate::force_enforce_all_waters ].value( true );
					// } else { // last ramping step
					//  basic::options::option[ basic::options::OptionKeys::hydrate::force_enforce_all_waters ].value( false );
					// }
					//}
					protocols::hydrate::HydrateOP hydrate_protocol( new protocols::hydrate::Hydrate( local_scorefxn ) );
					hydrate_protocol->apply( pose );
				} else {
					pack_full_repack_->apply( pose );
				}
				do_minimize( pose, cmd.param2, local_movemap, local_scorefxn );
				checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
			}


			// print some debug info
			if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
				core::Real imedscore = (*local_scorefxn)( pose );
				core::pose::setPoseExtraScore( pose, "R" + right_string_of( total_count ,5,'0'), imedscore );
			}

			if ( dumpall_ ) {
				pose.dump_pdb( "dump_" + right_string_of( dump_counter, 4, '0' ) );
				dump_counter++;
			}
		} else if ( cmd.command == "accept_to_best" ) {
			// grab the score and remember the pose if the score is better then ever before.
			core::Real score = (*local_scorefxn)( pose );
			if ( ( score < best_score) || (accept_count == 0) ) {
				best_score = score;
				best_pose = pose;
			}
#ifdef BOINC_GRAPHICS
			boinc::Boinc::update_graphics_low_energy( best_pose, best_score  );
			boinc::Boinc::update_graphics_last_accepted( pose, score );
			boinc::Boinc::update_mc_trial_info( total_count , "FastRelax" );
#endif
			core::Real rms = 0, irms = 0;
			if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
				rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
				irms = CA_rmsd_symmetric( start_pose , best_pose );
			} else {
				rms = native_CA_rmsd( *get_native_pose() , best_pose );
				irms = native_CA_rmsd( start_pose , best_pose );
			}
			TR << "MRP: " << accept_count << "  " << score << "  " << best_score << "  " << rms << "  " << irms << "  " << std::endl;
			best_score_log.push_back( best_score );
			curr_score_log.push_back( score );

			accept_count++;
			if ( accept_count > script_max_accept_ ) break;
			if ( test_cycles_ || dry_run() ) break;

		} else if ( cmd.command == "load_best" ) {
			pose = best_pose;
		} else if ( cmd.command == "load_start" ) {
			pose = start_pose;
		} else if ( cmd.command == "exit" ) {
			utility_exit_with_message( "EXIT INVOKED" );
		} else {

			utility_exit_with_message( "Unknown command: " + cmd.command );
		}

		core::Real rms( -1.0 );
		if ( get_native_pose() ) {
			rms =  native_CA_rmsd( *get_native_pose() , pose );
			if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
				rms = CA_rmsd_symmetric( *get_native_pose() , best_pose );
			}
		}
		core::Real irms =  native_CA_rmsd( start_pose , pose );
		if ( core::pose::symmetry::is_symmetric( pose ) && symmetric_rmsd_ ) {
			irms =  CA_rmsd_symmetric( start_pose , best_pose ); //need to make a symmetrical verision?
		}

		TR << "CMD: " <<  cmd.command << "  "
			<< (*local_scorefxn)( pose ) << "  "
			<< rms << "  "
			<< irms << "  "
			<< local_scorefxn->get_weight( scoring::fa_rep )
			<< std::endl;
	}

	pose = best_pose;
	(*local_scorefxn)( pose );

	if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
		for ( Size j = 0; j < best_score_log.size(); j++ ) {
			core::pose::setPoseExtraScore( pose, "B" + right_string_of(j,3,'0'), best_score_log[j]);
		}

		for ( Size j = 0; j < curr_score_log.size(); j ++ ) {
			core::pose::setPoseExtraScore( pose, "S" + right_string_of(j,3,'0'), curr_score_log[j] );
		}
	}

	checkpoints_.clear_checkpoints();

	if ( constrain_coords() ) {
		if ( delete_virtual_residues_after_FastRelax_ ) {
			// remove extra virtual atom cooordinate constraints
			protocols::relax::delete_virtual_residues( pose );
		}
	}
} //apply

// XRW TEMP std::string
// XRW TEMP FastRelax::get_name() const {
// XRW TEMP  return "FastRelax";
// XRW TEMP }

// Override the stored script with the default script for batchrelax
void FastRelax::set_script_to_batchrelax_default( core::Size repeats ) {
	script_.clear();

	std::vector< std::string > filelines;
	std::string line;

	//fpd
	//fpd sets a "reasonable" script for performing batch relax
	//fpd uses 'default_repeats'
	if ( repeats == 0 ) {
		repeats = default_repeats_;
	}
	runtime_assert( repeats > 0 );

	// repeat 1
	filelines.emplace_back("ramp_repack_min 0.02  0.01"  );
	filelines.emplace_back("batch_shave 0.25"  );
	filelines.emplace_back("ramp_repack_min 0.250 0.01"  );
	filelines.emplace_back("batch_shave 0.25"  );
	filelines.emplace_back("ramp_repack_min 0.550 0.01"  );
	filelines.emplace_back("batch_shave 0.25"  );
	filelines.emplace_back("ramp_repack_min 1     0.00001"  );
	filelines.emplace_back("accept_to_best"  );

	// repeats 2->n
	for ( core::Size i=2; i<=repeats; ++i ) {
		filelines.emplace_back("ramp_repack_min 0.02  0.01"  );
		filelines.emplace_back("ramp_repack_min 0.250 0.01"  );
		filelines.emplace_back("batch_shave 0.25"  );
		filelines.emplace_back("ramp_repack_min 0.550 0.01"  );
		filelines.emplace_back("ramp_repack_min 1     0.00001"  );
		filelines.emplace_back("accept_to_best"  );
	}

	set_script_from_lines( filelines );
}

// Override the stored script with the default script for batchrelax
void FastRelax::set_script_from_lines( std::vector< std::string > const & filelines )
{
	std::string line;

	script_.clear();

	core::Size i;
	for ( i =0; i< filelines.size(); i++ ) {
		line = filelines[i];
		TR.Debug << line << std::endl;
		utility::vector1< std::string > tokens ( utility::split( line ) );

		if ( tokens.size() > 20 ) { // large data; now used for reference weights storage only.
			RelaxScriptCommand newcmd;
			newcmd.command = tokens[1];
			for ( Size i = 2; i <= tokens.size(); ++i ) {
				newcmd.params_vec.push_back( atof(tokens[i].c_str()) );
			}
			newcmd.nparams = newcmd.params_vec.size();
			script_.push_back( newcmd );

		} else if ( tokens.size() > 0 ) {
			RelaxScriptCommand newcmd;
			newcmd.command = tokens[1];

			if ( tokens.size() > 1 ) { newcmd.param1 = atof(tokens[2].c_str()); newcmd.nparams = 1;}
			if ( tokens.size() > 2 ) { newcmd.param2 = atof(tokens[3].c_str()); newcmd.nparams = 2;}
			if ( tokens.size() > 3 ) { newcmd.param3 = atof(tokens[4].c_str()); newcmd.nparams = 3;}
			if ( tokens.size() > 4 ) { newcmd.param4 = atof(tokens[5].c_str()); newcmd.nparams = 4;}

			script_.push_back( newcmd );
		}
	}
}

void FastRelax::read_script_file( const std::string &script_file, core::Size standard_repeats ){
	using namespace ObjexxFCL;
	script_.clear();
	std::vector< std::string > filelines;
	std::string line;

	runtime_assert( standard_repeats > 0 );
	if ( script_file == "" && dualspace_ ) {
		TR << "================== Using dualspace script ==================" << std::endl;
		filelines.emplace_back("switch:torsion"           );
		filelines.push_back( "repeat " + string_of( standard_repeats - 1 ) );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     0.5"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     0.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  0.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );

		filelines.emplace_back("switch:cartesian"           );
		filelines.emplace_back("repeat 1"                  );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     0.5"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     0.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  0.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );
	} else if ( script_file == "" ) {
		TR << "================== Using default script ==================" << std::endl;
		filelines.push_back( "repeat " + string_of( standard_repeats )  );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     0.5"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     0.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  0.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );
	} else if ( script_file == "NO CST RAMPING" && dualspace_ ) {
		TR << "================== Using dualspace script - no constraint ramping ==================" << std::endl;
		filelines.emplace_back("switch:torsion"                  );
		filelines.push_back( "repeat " + string_of( standard_repeats - 1 ) );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  1.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );

		filelines.emplace_back("switch:cartesian"                  );
		filelines.emplace_back("repeat 1"                                );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  1.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );
		/*
		}else if (script_file == "NO CST RAMPING" && basic::options::option[ basic::options::OptionKeys::relax::dualfaster ]() ){
		TR << "================== Using faster dualspace script ==================" << std::endl;
		filelines.push_back( "switch:torsion"                  );
		filelines.push_back( "repeat 3"                                );
		filelines.push_back( "ramp_repack_min 0.02  0.01     1.0  50"  );
		filelines.push_back( "ramp_repack_min 0.250 0.01     1.0  50"  );
		filelines.push_back( "ramp_repack_min 0.550 0.01     1.0 100"  );
		filelines.push_back( "ramp_repack_min 1     0.00001  1.0 200"  );
		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );

		filelines.push_back( "switch:cartesian"                  );
		filelines.push_back( "repeat 2"                                );
		filelines.push_back( "ramp_repack_min 0.02  0.01     1.0  50"  );
		filelines.push_back( "ramp_repack_min 0.250 0.01     1.0  50"  );
		filelines.push_back( "ramp_repack_min 0.550 0.01     1.0 100"  );
		filelines.push_back( "ramp_repack_min 1     0.00001  1.0 200"  );
		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );
		*/
	} else if ( script_file == "NO CST RAMPING" ) {
		TR << "================== Using default script - no constraint ramping ==================" << std::endl;
		filelines.push_back( "repeat " + string_of( standard_repeats )  );
		filelines.emplace_back("ramp_repack_min 0.02  0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.250 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 0.550 0.01     1.0"      );
		filelines.emplace_back("ramp_repack_min 1     0.00001  1.0"      );
		filelines.emplace_back("accept_to_best"                  );
		filelines.emplace_back("endrepeat "                      );
	} else {
		std::ifstream infile( script_file.c_str() );
		TR.Debug << "================== Reading script file: ==================" << std::endl;
		if ( !infile.good() ) {
			utility_exit_with_message( "[ERROR] Error opening script file '" + script_file + "'" );
		}
		while ( getline(infile,line) ) {
			filelines.push_back( line );
		}
		infile.close();
	}

	set_script_from_lines( filelines );
}

// Batch Relax stuff

struct SRelaxPose {
	SRelaxPose() = default;

	core::Real initial_score;
	core::Real initial_rms;
	bool active{true};
	SilentStructOP current_struct;
	SilentStructOP start_struct;
	SilentStructOP best_struct;
	core::Real current_score;
	core::Real best_score;
	core::Size accept_count = 0;
	std::vector< core::Real > best_score_log;
	std::vector< core::Real > curr_score_log;

	core::Size mem_footprint(){
		return sizeof ( SRelaxPose ) + current_struct->mem_footprint() + start_struct->mem_footprint() + best_struct->mem_footprint() + (best_score_log.size()+curr_score_log.size())*sizeof( core::Real);
	}
};


void FastRelax::batch_apply(
	std::vector < SilentStructOP > & input_structs,
	core::scoring::constraints::ConstraintSetOP input_csts, // = NULL
	core::Real decay_rate // = 0.5
){
	using namespace basic::options;
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::pack::task;
	using namespace core::kinematics;
	using namespace protocols;


	check_nonideal_mintype();
	PackerTaskOP task_;
	protocols::minimization_packing::PackRotamersMoverOP pack_full_repack_;
	core::kinematics::MoveMapOP local_movemap = get_movemap()->clone();
	core::pose::Pose pose;


	TR.Debug  << "================== RelaxScriptBatchRelax: " << script_.size() << " ===============================" << std::endl;
	TR.Info  << "BatchRelax: Size: " << input_structs.size() << " Scriptlength: " <<  script_.size() << std::endl;

	// 432mb

	if ( input_structs.size()  < 1 ) return;

	ScoreFunctionOP local_scorefxn( get_scorefxn()->clone() );

	core::scoring::EnergyMap full_weights = local_scorefxn->weights();
	if ( dry_run() ) {
		return;
	}
	// Set up rama-repair if requested - one out of 10 times dont bother
	bool do_rama_repair = ramady_;
	if ( !ramady_force_ && numeric::random::uniform() <= 0.1 ) do_rama_repair = false;

	// 432/436mb

	// create a local array of relax_decoys
	std::vector < SRelaxPose > relax_decoys;
	core::Size total_mem = 0;

	for ( core::Size i = 0; i < input_structs.size(); ++i ) {
		TR.Debug << "iClock" << clock() << std::endl;

		SRelaxPose new_relax_decoy;

		input_structs[i]->fill_pose( pose );
		if ( input_csts ) pose.constraint_set( input_csts );
		// 432mb
		if ( !pose.is_fullatom() ) {
			TR.Debug << "Switching struct to fullatom" << std::endl;
			core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
			//core::Real afterscore = (*local_scorefxn)(pose);
		}
		// 453mb


		//pose.dump_pdb("test_" + string_of( i ) + ".pdb" );
		new_relax_decoy.active = true;
		new_relax_decoy.initial_score = (input_structs[i])->get_energy("censcore");

		new_relax_decoy.initial_rms = 0;
		if ( get_native_pose() ) {
			new_relax_decoy.initial_rms = native_CA_rmsd( *get_native_pose() , pose );
		}

		new_relax_decoy.best_score = (*local_scorefxn)(pose);
		new_relax_decoy.current_score = new_relax_decoy.best_score;
		SilentFileOptions opts;
		new_relax_decoy.current_struct = force_nonideal_?
			SilentStructFactory::get_instance()->get_silent_struct("binary", opts) :
			SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		new_relax_decoy.start_struct =  force_nonideal_?
			SilentStructFactory::get_instance()->get_silent_struct("binary", opts ) :
			SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		new_relax_decoy.best_struct =  force_nonideal_?
			SilentStructFactory::get_instance()->get_silent_struct("binary", opts) :
			SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		new_relax_decoy.current_struct->fill_struct( pose );
		new_relax_decoy.start_struct->fill_struct( pose );
		new_relax_decoy.best_struct->fill_struct( pose );
		new_relax_decoy.current_struct->copy_scores( *(input_structs[i]) );
		new_relax_decoy.start_struct->copy_scores( *(input_structs[i]) );
		new_relax_decoy.best_struct->copy_scores( *(input_structs[i]) );

		TR.Trace << "Fillstruct: " <<  new_relax_decoy.best_score << std::endl;

		initialize_movemap( pose, *local_movemap );
		TR.Trace << "SRelaxPose mem: " <<  new_relax_decoy.mem_footprint() << std::endl;
		total_mem += new_relax_decoy.mem_footprint();
		relax_decoys.push_back( new_relax_decoy );
		// 453 mb
		if ( i == 0 ) {
			if ( get_task_factory() ) {
				pack_full_repack_ = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover() );
				pack_full_repack_->score_function(local_scorefxn);
				pack_full_repack_->task_factory(get_task_factory());
			} else {
				task_ = TaskFactory::create_packer_task( pose );

				bool const repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
				utility::vector1<bool> allow_repack( pose.size(), repack);

				if ( !basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ) {
					for ( Size pos = 1; pos <= pose.size(); pos++ ) {
						allow_repack[ pos ] = local_movemap->get_chi( pos );
					}
				}
				task_->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
				task_->or_include_current( true );
				pack_full_repack_ = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover( local_scorefxn, task_ ) );
			}

			// Make sure we only allow symmetrical degrees of freedom to move
			if ( core::pose::symmetry::is_symmetric( pose )  )  {
				core::pose::symmetry::make_symmetric_movemap( pose, *local_movemap );
			}

			if ( basic::options::option[ basic::options::OptionKeys::run::debug ]() ) {
				kinematics::simple_visualize_fold_tree_and_movemap_bb_chi( pose.fold_tree(),  *local_movemap, TR );
			}
		}
		// 453 mb
	}
	relax_decoys.reserve( relax_decoys.size() );

	TR.Debug << "BatchRelax mem: " << total_mem << std::endl;

	/// RUN

	int repeat_step=0;
	int repeat_count=-1;
	int total_repeat_count = 0;

	int total_count=0;
	for ( core::Size ncmd = 0; ncmd < script_.size(); ncmd ++ ) {
		total_count++;

		RelaxScriptCommand cmd = script_[ncmd];
		TR.Debug << "Command: " << cmd.command << std::endl;

		if ( cmd.command == "repeat" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "ERROR: Syntax: " + cmd.command + "<number_of_repeats> " ); }
			repeat_count = (int) cmd.param1;
			repeat_step = ncmd;
		} else if ( cmd.command == "endrepeat" ) {
			TR.Debug << "CMD:  Repeat: " << repeat_count << std::endl;
			repeat_count -- ;
			total_repeat_count ++ ;
			if ( repeat_count <= 0 ) {}
			else {
				ncmd = repeat_step;
			}
		} else if ( cmd.command == "dump" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }

			for ( core::Size index=0; index < relax_decoys.size(); ++ index ) {
				relax_decoys[index].current_struct->fill_pose( pose );
				pose.dump_pdb( "dump_" + right_string_of( index, 4, '0' ) + "_" + right_string_of( (int) cmd.param1, 4, '0' ) );
			}
		} else if ( cmd.command == "repack" ) {
			for ( auto & relax_decoy : relax_decoys ) {
				if ( !relax_decoy.active ) continue;
				relax_decoy.current_struct->fill_pose( pose );
				if ( input_csts ) pose.constraint_set( input_csts );
				pack_full_repack_->apply( pose );
				core::Real score = (*local_scorefxn)(pose);
				relax_decoy.current_score = score;
				relax_decoy.current_struct->fill_struct( pose );
			}

		} else if ( cmd.command == "min" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }

			for ( auto & relax_decoy : relax_decoys ) {
				if ( !relax_decoy.active ) continue;
				relax_decoy.current_struct->fill_pose( pose );
				if ( input_csts ) pose.constraint_set( input_csts );
				do_minimize( pose, cmd.param1, local_movemap, local_scorefxn  );
				core::Real score = (*local_scorefxn)(pose);
				relax_decoy.current_score = score;
				relax_decoy.current_struct->fill_struct( pose );
			}
		} else if ( cmd.command.substr(0,5) == "scale" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(6));
			local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * cmd.param1 );
		}   else if ( cmd.command.substr(0,6) == "rscale" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
			local_scorefxn->set_weight( scale_param, full_weights[ scale_param ] * ((cmd.param2 - cmd.param1 ) * numeric::random::uniform() + cmd.param1 ));
		}   else if ( cmd.command.substr(0,6) == "switch" ) {
			// no input validation as of now, relax will just die
			if ( cmd.command.substr(7) == "torsion" ) {
				TR << "Using AtomTreeMinimizer"  << std::endl;
				cartesian( false );
			} else if ( cmd.command.substr(7) == "cartesian" ) {
				TR << "Using CartesianMinimizer"  << std::endl;
				cartesian( true );
			}
		}   else if ( cmd.command.substr(0,6) == "weight" ) {
			// no input validation as of now, relax will just die
			scoring::ScoreType scale_param = scoring::score_type_from_name(cmd.command.substr(7));
			local_scorefxn->set_weight( scale_param, cmd.param1 );
			// I'm not too sure if the changing the default weight makes sense
			full_weights[ scale_param ] = cmd.param1;
		}   else if ( cmd.command == "show_weights" ) {
			local_scorefxn->show(TR, pose);
		}   else if ( cmd.command == "ramp_repack_min" ) {
			if ( cmd.nparams < 2 ) { utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }
			local_scorefxn->set_weight( scoring::fa_rep, full_weights[ scoring::fa_rep ] * cmd.param1 );


			for ( auto & relax_decoy : relax_decoys ) {
				try {
					clock_t starttime = clock();
					if ( !relax_decoy.active ) continue;
					relax_decoy.current_struct->fill_pose( pose );
					if ( input_csts ) pose.constraint_set( input_csts );
					if ( total_repeat_count > 1 && repeat_count > 2 ) {
						if ( cmd.param1 < 0.2 ) {
							if ( do_rama_repair ) {
								fix_worst_bad_ramas( pose, ramady_num_rebuild_, 0.0, ramady_cutoff_, ramady_rms_limit_);
							}
						}
					}

					pack_full_repack_->apply( pose );
					do_minimize( pose, cmd.param2, local_movemap, local_scorefxn  );
					relax_decoy.current_score = (*local_scorefxn)(pose);
					relax_decoy.current_struct->fill_struct( pose );

					clock_t endtime = clock();
					TR.Debug << "time:" << endtime - starttime << " Score: " << relax_decoy.current_score << std::endl;
				} catch ( utility::excn::Exception& excn ) {
					std::cerr << "Ramp_repack_min exception: " << std::endl;
					excn.show( std::cerr );
					// just deactivate this pose
					relax_decoy.active = false;
					// and need to "reset scoring" of the pose we reuse
					TR << "Throwing out one structure due to scoring problems!" << std::endl;
					pose.scoring_end(*local_scorefxn);
				}
			}

		} else if ( cmd.command == "batch_shave" ) {
			if ( cmd.nparams < 1 ) { utility_exit_with_message( "More parameters expected after : " + cmd.command  ); }
			core::Real reduce_factor = cmd.param1;
			if ( (reduce_factor <= 0) ||  (reduce_factor >= 1.0)  ) { utility_exit_with_message( "Parameter after : " + cmd.command + " should be > 0 and < 1 " ); }

			TR.Debug << "SHAVE FACTOR: " << reduce_factor << std::endl;

			std::vector < core::Real > energies;
			for ( auto & relax_decoy : relax_decoys ) {
				if ( !relax_decoy.active ) continue;
				TR.Debug << "SHAVE: " << relax_decoy.current_score << std::endl;
				energies.push_back( relax_decoy.current_score );
				TR.Debug << relax_decoy.current_score << std::endl;
			}

			if ( energies.size() < 1 ) {
				TR.Error << "Cannot shave off structures - there are not enough left" << std::endl;
				continue;
			}

			std::sort( energies.begin(), energies.end() );
			auto cutoff_index = core::Size( floor(core::Real(energies.size()) * (reduce_factor)) );
			TR.Debug << "Energies: cutoff index " << cutoff_index << std::endl;
			core::Real cutoff_energy = energies[ cutoff_index ];
			TR.Debug << "Energies: cutoff " << cutoff_energy << std::endl;
			for ( core::Size index=0; index < relax_decoys.size(); ++ index ) {
				if ( !relax_decoys[index].active ) continue;
				TR.Debug << "Shaving candidate: " << index << "  "  << relax_decoys[index].current_score  << "  " << cutoff_energy << std::endl;
				if ( relax_decoys[index].current_score > cutoff_energy ) relax_decoys[index].active = false;
			}


		} else if ( cmd.command == "accept_to_best" ) {  // grab the score and remember the pose if the score is better then ever before.


			std::vector < core::Real > energies;

			for ( core::Size index=0; index < relax_decoys.size(); ++ index ) {
				if ( !relax_decoys[index].active ) continue;
				relax_decoys[index].current_struct->fill_pose( pose );
				if ( input_csts ) pose.constraint_set( input_csts );
				core::Real score = 0;
				try {
					score = (*local_scorefxn)( pose );
				} catch ( utility::excn::Exception& excn ) {
					std::cerr << "Accept_to_best scoring exception: " << std::endl;
					excn.show( std::cerr );
					// just deactivate this pose
					relax_decoys[index].active = false;
					// and need to "reset scoring" of the pose we reuse
					TR << "Throwing out one structure due to scoring problems!" << std::endl;
					pose.scoring_end(*local_scorefxn);
					continue;
				}
				TR.Debug << "Comparison: " << score << " " << relax_decoys[index].best_score << std::endl;

				if ( ( score < relax_decoys[index].best_score) || (relax_decoys[index].accept_count == 0) ) {
					relax_decoys[index].best_score = score;
					relax_decoys[index].best_struct->fill_struct( pose );
#ifdef DEBUG
					core::pose::Pose pose_check;
					relax_decoys[index].best_struct->fill_pose( pose_check );
					if ( input_csts ) pose.constraint_set( input_csts );
					core::Real score_check = (*local_scorefxn)( pose_check );
					TR.Debug << "Sanity: "<< score << " ==  " << score_check << std::endl;
#endif
				}

				if ( get_native_pose() ) {
					if ( core::pose::symmetry::is_symmetric( pose ) && option[ basic::options::OptionKeys::evaluation::symmetric_rmsd ].user() ) {
						core::Real rms = CA_rmsd_symmetric( *get_native_pose() , pose );
						TR << "MRP: " << index << " " << relax_decoys[index].accept_count << "  " << score << "  " << relax_decoys[index].best_score << "  "
							<< rms << "  "
							<< std::endl;

					} else {
						core::Real rms =  native_CA_rmsd( *get_native_pose() , pose );
						TR << "MRP: " << index << " " << relax_decoys[index].accept_count << "  " << score << "  " << relax_decoys[index].best_score << "  "
							<< rms << "  "
							<< std::endl;
					}
				} else {
					TR << "MRP: " << index << " " << relax_decoys[index].accept_count << "  " << score << "  " << relax_decoys[index].best_score << "  "
						<< std::endl;
				}

				relax_decoys[index].curr_score_log.push_back( score );
				relax_decoys[index].accept_count ++ ;
				energies.push_back( relax_decoys[index].best_score );
			}


			// now decide who to keep relaxing and who to abandon!
			if ( energies.size() < 1 ) {
				TR.Debug << "Cannot shave off structures - there are not enough left" << std::endl;
				continue;
			}
			std::sort( energies.begin(), energies.end() );
			core::Real cutoff_energy = energies[ core::Size( floor(core::Real(energies.size()) * (decay_rate)) )];
			for ( auto & relax_decoy : relax_decoys ) {
				if ( !relax_decoy.active ) continue;
				if ( relax_decoy.best_score > cutoff_energy ) relax_decoy.active = false;
			}

		} else if ( cmd.command == "exit" ) {
			utility_exit_with_message( "EXIT INVOKED FROM SEQUENCE RELAX SCRIPT" );
		} else {

			utility_exit_with_message( "Unknown command: " + cmd.command );
		}

		TR.Debug << "CMD: " <<  cmd.command << "  Rep_Wght:"
			<< local_scorefxn->get_weight( scoring::fa_rep )
			<< std::endl;
	}


	input_structs.clear();
	// Finally return all the scores;
	for ( auto & relax_decoy : relax_decoys ) {
		relax_decoy.best_struct->fill_pose( pose );
		if ( input_csts ) pose.constraint_set( input_csts );
		setPoseExtraScore( pose, "giveup", relax_decoy.accept_count  );
		core::Real rms = 0;
		core::Real gdtmm = 0;
		if ( get_native_pose() ) {
			rms =  native_CA_rmsd( *get_native_pose() , pose );
			gdtmm =  native_CA_gdtmm( *get_native_pose() , pose );
			TR.Trace << "BRELAXRMS: " << rms << "  " << gdtmm << std::endl;
		}
		//pose.dump_pdb("best_relax_step_"+string_of( index )+".pdb" );
		core::Real score = (*local_scorefxn)( pose );
		TR << "BRELAX: Rms: "<< rms
			<< " CenScore: "  << relax_decoy.initial_score
			<< " CenRMS: "    << relax_decoy.initial_rms
			<< " FAScore: "   << score
			<< " Check: "     << relax_decoy.best_score
			<< " Acc: "       << relax_decoy.accept_count
			<< " Go: "        << (relax_decoy.active ? "  1" : "  0") <<  std::endl;

		if ( !relax_decoy.active ) continue;

		core::io::silent::SilentFileOptions opts;
		SilentStructOP new_struct =  force_nonideal_?
			SilentStructFactory::get_instance()->get_silent_struct("binary", opts) :
			SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		new_struct->fill_struct( pose );
		new_struct->copy_scores( *(relax_decoy.start_struct) );
		new_struct->energies_from_pose( pose );
		new_struct->add_energy("rms", rms, 1.0 );
		new_struct->add_energy("gdtmm", gdtmm, 1.0 );
		input_structs.push_back( new_struct );
	}


}

void
FastRelax::set_constraint_weight( core::scoring::ScoreFunctionOP local_scorefxn,
	core::scoring::EnergyMap const & full_weights,
	core::Real const weight,
	core::pose::Pose & ) const {
	local_scorefxn->set_weight( scoring::coordinate_constraint,  full_weights[ scoring::coordinate_constraint ] * weight );
}

std::string FastRelax::get_name() const {
	return mover_name();
}

std::string FastRelax::mover_name() {
	return "FastRelax";
}

void FastRelax::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_fast_relax( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Performs the fast relax protocol" )
		.write_complex_type_to_schema( xsd );

}

utility::tag::XMLSchemaComplexTypeGeneratorOP
FastRelax::complex_type_generator_for_fast_relax( utility::tag::XMLSchemaDefinition & xsd){
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function(attlist);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"disable_design", xsct_rosetta_bool,
		"Do not perform design even if a resfile is specified",
		"true");

	rosetta_scripts::attributes_for_parse_task_operations_w_factory(attlist);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"repeats", xs_integer,
		"Same as cmd-line FR. Number of fastrelax repeats to perform",
		"5");

	attlist + XMLSchemaAttribute(
		"name", xs_string,
		"Name of the mover");

	attlist + XMLSchemaAttribute(
		"relaxscript", xs_string,
		"a filename for a relax script, as described in the documentation for the Relax application; "
		"the default relax script is used if not specified");

	attlist + XMLSchemaAttribute(
		"cst_file", xs_string,
		"Add constraints from the constraint file");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"batch", xsct_rosetta_bool,
		"Perform BatchRelax?",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"cartesian", xsct_rosetta_bool,
		"Use cartesian minimization",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"dualspace", xsct_rosetta_bool,
		"Use dualspace minimization",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"ramp_down_constraints", xsct_rosetta_bool,
		"Ramp down constraint weights over the course of the protocol",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"bondangle", xsct_rosetta_bool,
		"Minimize bond angles",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"bondlength", xsct_rosetta_bool,
		"Minimize bond lengths",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_type", xsct_minimizer_type,
		"Minimizer type to use",
		"lbfgs_armijo_nonmonotone");

	attlist + XMLSchemaAttribute(
		"delete_virtual_residues_after_FastRelax", xsct_rosetta_bool,
		"Should virtual residues be deleted when the protocol completes?");

	core::select::movemap::attributes_for_parse_movemap_factory_default_attr_name( attlist,
		"The name of the already-defined MoveMapFactory that will be used to alter the"
		" default behavior of the MoveMap. By default, all backbone, chi, and jump DOFs"
		" are allowed to change. A MoveMapFactory can be used to change which of those"
		" DOFs are actually enabled. Be warned that combining a MoveMapFactory with"
		" a MoveMap can result in unexpected behavior. The MoveMap provided as a"
		" subelement of this element will be generated, and then the DOF modifications"
		" specified in the MoveMapFactory will be applied afterwards. E.g., if the"
		" MoveMap says to allow bb flexibility for res 1-10, and then the MoveMapFactory"
		" says to disable all bb flexibility, the MoveMapFactory modification will happen"
		" second, and this second modification will completely override the MoveMap's"
		" initial instructions. Pay attention to the precedence rules of the MoveMap. If"
		" high-level instructions (e.g. 'apply to all') are performed after low-level"
		" instructions, then the low-level instructions are erased. If low-level instructions"
		" (e.g. 'apply to a particular residue') are performed after a high-level instructions,"
		" the low-level instructions override the high-level instructions. E.g., if there are"
		" 100 residues, then the instruction 'all-all-bb' followed by 'disable-bb-for-residues-10-to-30'"
		" followed by the instruction 'allow-bb-for-bb-dihedral-number-3-for-residue-15' will mean"
		" that residues 1-9 will be allowed backbone flexibility, and residues 31-100, too, and"
		" that only backbone dihedral 3 on residue 15 will have flexibility among the residues from"
		" 10 to 30. However, if a final instruction 'disallow-all-bb' is given, then no backbone"
		" dihedrals will be free, not even dihedral 3 on residue 15." );

	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->add_attributes( attlist )
		.set_subelements_repeatable( subelements, 0, 1)
		.complex_type_naming_func( & moves::complex_type_name_for_mover );
	return ct_gen;
}




std::string FastRelaxCreator::keyname() const {
	return FastRelax::mover_name();
}

protocols::moves::MoverOP
FastRelaxCreator::create_mover() const {
	return protocols::moves::MoverOP( new FastRelax );
}

void FastRelaxCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FastRelax::provide_xml_schema( xsd );
}


} // namespace relax

} // namespace protocols
