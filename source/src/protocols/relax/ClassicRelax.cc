// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_protocols
/// @brief protocols that are specific to relax
/// @details
/// @author Mike Tyka, Monica Berrondo


#include <protocols/relax/ClassicRelax.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.hh>



#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/RampingMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>

// Symmetry
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <protocols/simple_moves/GunnCost.hh>
#include <core/io/raw_data/ScoreMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/simple_moves/WobbleMover.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif


//#include <utility/io/mpistream.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.relax.ClassicRelax" );

using namespace core;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {

std::string
ClassicRelaxCreator::keyname() const
{
	return ClassicRelaxCreator::mover_name();
}

protocols::moves::MoverOP
ClassicRelaxCreator::create_mover() const {
	return protocols::moves::MoverOP( new ClassicRelax );
}

std::string
ClassicRelaxCreator::mover_name()
{
	return "ClassicRelax";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ClassicRelax::ClassicRelax(
	core::scoring::ScoreFunctionOP scorefxn_in
) :
	parent("ClassicRelax", scorefxn_in ),
	checkpoints_("ClassicRelax")
{
	set_default();
	// set these to true - by deafult we're using the default
	// types for these objects. if the user chooses to set their own
	// these will be set to false to indicate that the user has done so
	use_default_pack_full_repack_ = true;
	use_default_pack_rottrial_ = true;
	use_default_mc_ = true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief constructor taking both ScoreFunction and MoveMap
ClassicRelax::ClassicRelax( core::scoring::ScoreFunctionOP scorefxn_in, core::kinematics::MoveMapOP movemap ) :
	parent( "ClassicRelax",scorefxn_in ),
	checkpoints_("ClassicRelax")
{
	use_default_pack_full_repack_ = true;
	use_default_pack_rottrial_ = true;
	use_default_mc_ = true;
	set_movemap( movemap );
	set_default( false  /* use_default_movemap */ );
}
////////////////////////////////////////////////////////////////////////////////////////////////////

ClassicRelax::ClassicRelax() :
	parent( std::string("ClassicRelax") ),
	checkpoints_( std::string("ClassicRelax") )
{
	set_default();

	use_default_pack_full_repack_ = true;
	use_default_pack_rottrial_ = true;
	use_default_mc_ = true;
}

ClassicRelax::ClassicRelax( ClassicRelax const & other ) :
	//utility::pointer::ReferenceCount(),
	parent( other ),
	min_mover_( other.min_mover_ ),
	checkpoints_( other.checkpoints_ ),
	mc_( other.mc_ ),
	use_default_mc_( other.use_default_mc_ ),
	pack_full_repack_( other.pack_full_repack_ ),
	use_default_pack_full_repack_( other.use_default_pack_full_repack_ ),
	pack_rottrial_( other.pack_rottrial_ ),
	use_default_pack_rottrial_( other.use_default_pack_rottrial_ ),
	m_Temperature( other.m_Temperature ),
	nmoves_( other.nmoves_ ),
	energycut( other.energycut ),
	min_type( other.min_type ),
	nb_list( other.nb_list ),
	min_tolerance( other.min_tolerance ),
	moveset_phase1_( other.moveset_phase1_ ),
	moveset_phase2_( other.moveset_phase2_ ),
	moveset_phase3_( other.moveset_phase3_ ),
	lj_ramp_cycles( other.lj_ramp_cycles ),
	lj_ramp_inner_cycles( other.lj_ramp_inner_cycles ),
	start_rep_weight( other.start_rep_weight ),
	end_rep_weight( other.end_rep_weight ),
	st_rep_( other.st_rep_ ),
	st_atr_( other.st_atr_ ),
	st_sol_( other.st_sol_ ),
	stage2_repack_period( other.stage2_repack_period ),
	stage2_cycles( other.stage2_cycles ),
	stage3_cycles( other.stage3_cycles ),
	score_stage2_beginning( other.score_stage2_beginning ),
	score_stage2_quarter( other.score_stage2_quarter ),
	score_stage2_half( other.score_stage2_half ),
	score_stage2_end( other.score_stage2_end ),
	filter_stage2_beginning( other.filter_stage2_beginning ),
	filter_stage2_quarter( other.filter_stage2_quarter ),
	filter_stage2_half( other.filter_stage2_half ),
	filter_stage2_end( other.filter_stage2_end )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP ClassicRelax::clone() const {
	return protocols::moves::MoverOP( new ClassicRelax(*this) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
ClassicRelax::~ClassicRelax()= default;
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details registering of options that are relevant for AbrelaxApplication
void ClassicRelax::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	parent::register_options();

	option.add_relevant( OptionKeys::relax::wobblemoves );
	option.add_relevant( OptionKeys::relax::constrain_relax_to_native_coords );
	option.add_relevant( OptionKeys::relax::constrain_relax_to_start_coords );
	option.add_relevant( OptionKeys::relax::constrain_relax_segments );
	option.add_relevant( OptionKeys::relax::ramp_constraints );
	option.add_relevant( OptionKeys::relax::energycut );
	option.add_relevant( OptionKeys::relax::stage1_ramp_cycles );
	option.add_relevant( OptionKeys::relax::stage1_ramp_inner_cycles );
	option.add_relevant( OptionKeys::relax::stage2_repack_period );
	option.add_relevant( OptionKeys::relax::stage2_cycles );
	option.add_relevant( OptionKeys::relax::min_tolerance );
	option.add_relevant( OptionKeys::relax::stage3_cycles );
	option.add_relevant( OptionKeys::relax::cycle_ratio );
	option.add_relevant( OptionKeys::relax::filter_stage2_beginning );
	option.add_relevant( OptionKeys::relax::filter_stage2_quarter );
	option.add_relevant( OptionKeys::relax::filter_stage2_half );
	option.add_relevant( OptionKeys::relax::filter_stage2_end );
}


void ClassicRelax::set_default( core::scoring::ScoreFunctionOP scorefxn_in ) {
	set_scorefxn( scorefxn_in );
	set_default();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::set_default( bool const use_default_movemap ){
	TR <<  "Setting up default relax setting" << std::endl;
	// minimization
	//min_type = std::string("lbfgs_armijo_nonmonotone");
	min_type = RelaxProtocolBase::min_type();                 // use the base class min_type value
	nb_list = true;


	st_rep_ = core::scoring::fa_rep;
	st_atr_ = core::scoring::fa_atr;
	st_sol_ = core::scoring::fa_sol;

	m_Temperature = 0.8;

	nmoves_ = 5;
	energycut            = basic::options::option[ basic::options::OptionKeys::relax::energycut ];

	lj_ramp_cycles       = basic::options::option[ basic::options::OptionKeys::relax::stage1_ramp_cycles];
	lj_ramp_inner_cycles = basic::options::option[ basic::options::OptionKeys::relax::stage1_ramp_inner_cycles];
	start_rep_weight     = 0.02;
	end_rep_weight       = 1.0;

	// PHASE2 stuff
	stage2_repack_period =    basic::options::option[ basic::options::OptionKeys::relax::stage2_repack_period];
	stage2_cycles        =    basic::options::option[ basic::options::OptionKeys::relax::stage2_cycles ];
	min_tolerance     =    basic::options::option[ basic::options::OptionKeys::relax::min_tolerance ]; //0.00025; //as in stage2 in rosetta++

	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
		stage2_cycles = 1;
		stage3_cycles = 1;
		lj_ramp_cycles = 1;
		min_tolerance     = 0.2;
	}

	// PHASE3 stuff
	stage3_cycles =           basic::options::option[ basic::options::OptionKeys::relax::stage3_cycles ];

	filter_stage2_beginning = basic::options::option[ basic::options::OptionKeys::relax::filter_stage2_beginning ];
	filter_stage2_quarter   = basic::options::option[ basic::options::OptionKeys::relax::filter_stage2_quarter   ];
	filter_stage2_half      = basic::options::option[ basic::options::OptionKeys::relax::filter_stage2_half      ];
	filter_stage2_end       = basic::options::option[ basic::options::OptionKeys::relax::filter_stage2_end       ];

	score_stage2_beginning = 0;
	score_stage2_quarter = 0;
	score_stage2_half = 0;
	score_stage2_end = 0;

	//chu move these two calls after values are set, otherwise it causes bad bug with uninitialized variables.

	if ( use_default_movemap ) {
		set_default_movemap();
	}
	set_default_minimizer();

	set_default_moveset_phase1();
	set_default_moveset_phase2();
	set_default_moveset_phase3();
}

void ClassicRelax::set_tolerance( core::Real new_tolerance ){
	TR.Info << "Setting min tolerance: " << new_tolerance << std::endl;
	min_tolerance = new_tolerance;
	set_default_minimizer();
}


// sets up the default minimizer object with all the options
void ClassicRelax::set_default_minimizer() {
	// options for minimizer
	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
		min_mover_ = protocols::simple_moves::MinMoverOP( new simple_moves::symmetry::SymMinMover( get_movemap(), get_scorefxn(), min_type, min_tolerance, nb_list ) );
	} else {
		min_mover_ = protocols::simple_moves::MinMoverOP( new protocols::simple_moves::MinMover( get_movemap(), get_scorefxn(), min_type, min_tolerance, nb_list ) );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details At stage 1 we're only doing small and shear moves
void ClassicRelax::set_default_moveset_phase1()
{
	// setup the move objects
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( get_movemap(), m_Temperature, nmoves_ ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );

	// setup the move objects
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover( get_movemap(), m_Temperature, nmoves_ ) );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 2.0 );
	shear_mover->angle_max( 'L', 3.0 );

	// create a Random Mover, fill it with individual moves
	moves::RandomMoverOP moveset_phase1_temp( new moves::RandomMover() );
	//moveset_phase1_temp ->add_mover( small_mover );
	//moveset_phase1_temp ->add_mover( shear_mover );

	// then set our internal moveset
	moveset_phase1_ = moveset_phase1_temp;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details At stage 2 we're doing small, shear, wobble and crank moves (the latter two are works in progress)
void ClassicRelax::set_default_moveset_phase2()
{
	// setup the move objects
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( get_movemap(), m_Temperature, nmoves_ ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );

	// setup the move objects
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover( get_movemap(), m_Temperature, nmoves_ ) );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 2.0 );
	shear_mover->angle_max( 'L', 3.0 );

	// create a Random Mover, fill it with individual moves
	moves::RandomMoverOP moveset_phase2_temp( new moves::RandomMover() );
	moveset_phase2_temp ->add_mover( small_mover );
	moveset_phase2_temp ->add_mover( shear_mover );

	// setup the move object
	if ( basic::options::option[ basic::options::OptionKeys::relax::wobblemoves ].user() &&
			basic::options::option[ basic::options::OptionKeys::in::file::frag3 ].user() ) {
		std::string frag3_file  = basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]();
		core::fragment::ConstantLengthFragSetOP fragset3mer( new core::fragment::ConstantLengthFragSet( 3 ) );
		fragset3mer->read_fragment_file( frag3_file );
		protocols::simple_moves::WobbleMoverOP wobble_mover( new protocols::simple_moves::WobbleMover( fragset3mer, get_movemap(), protocols::simple_moves::FragmentCostOP( new protocols::simple_moves::GunnCost ) ) );

		moveset_phase2_temp ->add_mover( wobble_mover );
		moveset_phase2_temp ->add_mover( wobble_mover );
	}

	// then set our internal moveset
	moveset_phase2_ = moveset_phase2_temp;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details At stage 3 we're only doing small and shear moves
void ClassicRelax::set_default_moveset_phase3()
{
	// setup the move objects
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( get_movemap(), m_Temperature, nmoves_ ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );

	// setup the move objects
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover( get_movemap(), m_Temperature, nmoves_ ) );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 2.0 );
	shear_mover->angle_max( 'L', 3.0 );

	// create a Random Mover, fill it with individual moves
	moves::RandomMoverOP moveset_phase3_temp( new moves::RandomMover() );
	moveset_phase3_temp ->add_mover( small_mover );
	moveset_phase3_temp ->add_mover( shear_mover );

	// then set our internal moveset
	moveset_phase3_ = moveset_phase3_temp;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::check_default_mc( core::pose::Pose &pose) {
	if ( use_default_mc_ ) {
		mc_ = moves::MonteCarloOP( new moves::MonteCarlo( pose , *get_scorefxn() , 0.8 ) );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
moves::MonteCarloOP ClassicRelax::get_mc( core::pose::Pose &pose ) {
	check_default_mc( pose );
	return mc_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::set_mc ( moves::MonteCarloOP new_mc_ ){
	mc_ = new_mc_;
	use_default_mc_ = false; // user has set his own
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::check_default_full_repacker( core::pose::Pose & pose, core::kinematics::MoveMap & movemap ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task::operation;

	if ( use_default_pack_full_repack_ ) {

		// Add TF behavior Jadolfbr 5/2/2013

		core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );

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
				//Though, as mentioned in the doc, movemap now overrides chi_move as it should.

				local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
				PreventRepackingOP turn_off_packing( new PreventRepacking() );
				for ( Size pos = 1; pos <= pose.size(); ++pos ) {
					if ( ! movemap.get_chi(pos) ) {
						turn_off_packing->include_residue(pos);
					}
				}
				local_tf->push_back(turn_off_packing);
			}
		}
		//Include current rotamer by default - as before.
		local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));

		if ( limit_aroma_chi2() ) {
			local_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::LimitAromaChi2Operation() ));
		}

		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			pack_full_repack_ = protocols::simple_moves::PackRotamersMoverOP( new simple_moves::symmetry::SymPackRotamersMover( get_scorefxn()) );
		} else {
			pack_full_repack_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( get_scorefxn()) );
		}
		pack_full_repack_->task_factory(local_tf);

		(*get_scorefxn())( pose );
	}


}
////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::set_full_repack( protocols::simple_moves::PackRotamersMoverOP new_pack_full_repack ) {
	pack_full_repack_ = new_pack_full_repack;
	use_default_pack_full_repack_ = false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::check_default_rottrial( core::pose::Pose & pose, core::kinematics::MoveMap & movemap ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task::operation;

	if ( use_default_pack_full_repack_ ) {

		// Add TF behavior Jadolfbr 5/2/2013

		core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );

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
				//Though, as mentioned in the doc, movemap now overrides chi_move as it should.

				local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
				PreventRepackingOP turn_off_packing( new PreventRepacking() );
				for ( Size pos = 1; pos <= pose.size(); ++pos ) {
					if ( ! movemap.get_chi(pos) ) {
						turn_off_packing->include_residue(pos);
					}
				}
				local_tf->push_back(turn_off_packing);
			}
		}
		//Include current rotamer by default - as before.
		local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));

		if ( limit_aroma_chi2() ) {
			local_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::LimitAromaChi2Operation() ));
		}
		(*get_scorefxn())( pose );
		/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies( pose ); // fix this
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			pack_rottrial_ = protocols::simple_moves::RotamerTrialsMoverOP( new simple_moves::symmetry::SymEnergyCutRotamerTrialsMover( get_scorefxn(), local_tf, mc_, energycut ) );
		} else {
			pack_rottrial_ = protocols::simple_moves::RotamerTrialsMoverOP( new protocols::simple_moves::EnergyCutRotamerTrialsMover( get_scorefxn(), local_tf, mc_, energycut ) );
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::set_rottrial ( protocols::simple_moves::RotamerTrialsMoverOP new_pack_rottrial ){
	pack_rottrial_ = new_pack_rottrial;
	use_default_pack_full_repack_ = false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

void ClassicRelax::setPoseExtraScore( pose::Pose &pose ){
	core::pose::setPoseExtraScore( pose, "Filter_Stage2_dEnd", score_stage2_end);
	core::pose::setPoseExtraScore( pose, "Filter_Stage2_cHalf", score_stage2_half);
	core::pose::setPoseExtraScore( pose, "Filter_Stage2_bQuarter",score_stage2_quarter);
	core::pose::setPoseExtraScore( pose, "Filter_Stage2_aBefore", score_stage2_beginning);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassicRelax::apply( core::pose::Pose & pose ){
	using namespace moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	runtime_assert( get_scorefxn() != nullptr );
	(*get_scorefxn())(pose);

	/// Invoke parent local_movemap initialization routines
	core::kinematics::MoveMapOP local_movemap = get_movemap()->clone();
	initialize_movemap( pose, *local_movemap );

	// Make sure we only allow symmetrical degrees of freedom to move
	// if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
	//  core::pose::symmetry::make_symmetric_movemap( pose, *local_movemap );
	// }

	// remember the original pose before refinement kicks in, we'll need it later
	core::pose::Pose prerefine_pose = pose;

	// Set up any internal constraints that need to be set
	set_up_constraints( pose, *local_movemap );

	check_default_mc( pose );
	check_default_full_repacker( pose, *local_movemap);
	check_default_rottrial( pose, *local_movemap);

	// Make sure we only allow symmetrical degrees of freedom to move
	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *local_movemap );
	}

	min_mover_->movemap(local_movemap);

	get_scorefxn()->show(TR.Debug , pose);

	apply_disulfides(pose);

	if ( ( lj_ramp_cycles > 0 ) &&  (!checkpoints_.recover_checkpoint( pose, "stage_1", get_current_tag(),  true, true  )) ) {

		// part 1 ----------------------------------------
		TR.Info  << std::endl << std::endl <<  "===================================================================" << std::endl;
		TR.Info  <<  "   Stage 1                                                         " << std::endl;
		TR.Info  <<  "   Ramping repulsives with " << lj_ramp_cycles << " outer cycles and "
			<< lj_ramp_inner_cycles << " inner cycles" << std::endl;

		// setup the repack cycle
		moves::CycleMoverOP repack_cycle( new moves::CycleMover() );
		repack_cycle->add_mover( pack_full_repack_ );
		for ( int i=1; i<lj_ramp_inner_cycles; ++i ) {
			repack_cycle->add_mover( pack_rottrial_ );
		}

		moves::SequenceMoverOP phase1_cycle( new moves::SequenceMover() );
		//phase1_cycle->add_mover( moveset_phase1_ );
		phase1_cycle->add_mover( repack_cycle );

		moves::JumpOutMoverOP phase1_min( new moves::JumpOutMover( phase1_cycle, min_mover_, get_scorefxn(), 1E6 ) );

		moves::TrialMoverOP phase1_trial( new moves::TrialMover( phase1_min, mc_ ) );

		moves::RampingMoverOP full_cycle_phase1_;
		if ( ramp_down_constraints() ) {
			core::scoring::EnergyMap starting_weights, final_weights;
			starting_weights = final_weights = get_scorefxn()->weights();
			starting_weights[ fa_rep             ] = start_rep_weight * starting_weights[ fa_rep ]; // ramp up repulsion
			final_weights[ coordinate_constraint ] = 0;
			final_weights[ atom_pair_constraint ] = 0;
			final_weights[ angle_constraint ] = 0;
			final_weights[ dihedral_constraint ] = 0;

			full_cycle_phase1_ = moves::RampingMoverOP( new moves::RampingMover(
				phase1_trial, get_scorefxn(),
				starting_weights, final_weights,
				lj_ramp_cycles, lj_ramp_inner_cycles, mc_ ) );
			full_cycle_phase1_->set_func_for_weight( coordinate_constraint, RampingFuncOP( new moves::FastLinearFunc( 0, 0.6 ) ) );
			full_cycle_phase1_->set_func_for_weight( atom_pair_constraint, RampingFuncOP( new moves::FastLinearFunc( 0, 0.6 ) ) );
		} else {
			full_cycle_phase1_ = moves::RampingMoverOP( new moves::RampingMover(
				phase1_trial, get_scorefxn(), scoring::fa_rep,
				lj_ramp_cycles, lj_ramp_inner_cycles, mc_ ) );
			full_cycle_phase1_->start_weight( start_rep_weight * get_scorefxn()->weights()[ scoring::fa_rep ] );
			full_cycle_phase1_->end_weight( end_rep_weight * get_scorefxn()->weights()[ scoring::fa_rep ]  );
		}

		full_cycle_phase1_->apply( pose );

		checkpoints_.checkpoint( pose, "stage_1", get_current_tag(), true );

	}

	mc_->reset_scorefxn( pose, *get_scorefxn() );
	(*get_scorefxn())(pose);
	get_scorefxn()->show(TR.Debug , pose);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		// save time if the pose is symmetric
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
		// need to copy virtuals first
		for ( Size ii = prerefine_pose.size(); ii>=1; --ii ) {
			if ( symm_info->fa_is_independent(ii) ) {
				prerefine_pose.replace_residue( ii, pose.residue( ii ), false);
			}
		}
	} else {
		for ( Size ii = 1; ii <= prerefine_pose.size(); ++ii ) {
			prerefine_pose.replace_residue( ii, pose.residue( ii ), false);
		}
	}
	pose = prerefine_pose;
	(*get_scorefxn())(pose);

	check_default_mc( pose );
	check_default_full_repacker( pose, *local_movemap );
	check_default_rottrial( pose, *local_movemap );


	/// Stage 2 is
	/// Stage 2 is mysterious, just like the comment above.
	/// Stage 2 is <useless> ! Use fastrelax instead.
	if ( stage2_cycles < 0 ) {
		stage2_cycles = pose.size() * 4;
	}


	if ( stage2_cycles > 0 ) {
		if ( !checkpoints_.recover_checkpoint( pose, "stage_2", get_current_tag(), true, true) ) {

			core::Real temp_score  =  (*get_scorefxn())(pose);
			get_scorefxn()->show( TR.Debug , pose );
			score_stage2_beginning = temp_score;
			score_stage2_quarter = temp_score;
			score_stage2_half =  temp_score;
			score_stage2_end = temp_score;
			setPoseExtraScore( pose );
			mc_->reset(pose);

			if ( filter_stage2_beginning <  temp_score ) {
				TR.Info << "Structure failed filter_stage2_beginning: " <<  temp_score << " > " << filter_stage2_beginning << std::endl;
				return;
			}

			// part 2 ----------------------------------------
			TR.Info  << std::endl << "===================================================================" << std::endl;
			TR.Info  <<  "   Stage 2                                                         " << std::endl;
			TR.Info  <<  "   Mainmintrial for " << stage2_cycles << " with a full repack every "
				<< stage2_repack_period << " cycles" << std::endl;

			// setup the repack cycle
			moves::CycleMoverOP repack_cycle( new moves::CycleMover() );
			for ( int i=1; i<stage2_repack_period; ++i ) {
				repack_cycle->add_mover( pack_rottrial_ );
			}
			repack_cycle->add_mover( pack_full_repack_ );

			moves::SequenceMoverOP phase2_cycle( new moves::SequenceMover() );
			phase2_cycle->add_mover( moveset_phase2_ );
			phase2_cycle->add_mover( repack_cycle );

			moves::JumpOutMoverOP phase2_min( new moves::JumpOutMover( phase2_cycle, min_mover_, get_scorefxn(), 25.0 ) );

			moves::TrialMoverOP phase2_trial( new moves::TrialMover( phase2_min, mc_ ) );

			moves::RepeatMoverOP full_cycle_phase2_;
			full_cycle_phase2_ = moves::RepeatMoverOP( new moves::RepeatMover( phase2_trial, int( core::Real(stage2_cycles) * basic::options::option[ basic::options::OptionKeys::relax::cycle_ratio] / 4.0 ) ) );

			full_cycle_phase2_->apply( pose );

			score_stage2_quarter = mc_->lowest_score();
			score_stage2_half = mc_->lowest_score();
			score_stage2_end = mc_->lowest_score();
			if ( filter_stage2_quarter <  mc_->lowest_score() ) {
				TR.Info << "Structure failed filter_stage2_quarter: " <<  mc_->lowest_score() << " > " << filter_stage2_quarter << std::endl;
				mc_->recover_low( pose );
				setPoseExtraScore( pose );
				return;
			}

			full_cycle_phase2_->apply( pose );

			score_stage2_half = mc_->lowest_score();
			score_stage2_end = mc_->lowest_score();
			if ( filter_stage2_half <  mc_->lowest_score() ) {
				TR.Info << "Structure failed filter_stage2_half: " <<  mc_->lowest_score() << " > " << filter_stage2_quarter << std::endl;
				mc_->recover_low( pose );
				setPoseExtraScore( pose );
				return;
			}

			full_cycle_phase2_->apply( pose );
			full_cycle_phase2_->apply( pose );

			score_stage2_end = mc_->lowest_score();
			if ( filter_stage2_end<  mc_->lowest_score() ) {
				TR.Info << "Structure failed score_stage2_end: " <<  mc_->lowest_score() << " > " << filter_stage2_quarter << std::endl;
				mc_->recover_low( pose );
				setPoseExtraScore( pose );
				return;
			}

			checkpoints_.checkpoint( pose, "stage_2", get_current_tag(), true);
		}
	}

	mc_->recover_low( pose );
	(*get_scorefxn())(pose);
	get_scorefxn()->show(TR.Debug, pose);
	setPoseExtraScore( pose );
	mc_->reset(pose);

	// output_debug_structure( pose, "rl_stage2" );

	if ( stage3_cycles < 0 ) {
		stage3_cycles = pose.size();
	}
	if ( stage3_cycles > 0 ) {
		if ( !checkpoints_.recover_checkpoint( pose, "stage_3", get_current_tag(), true, true ) ) {

			// part 3 ----------------------------------------
			TR.Info  << std::endl << "===================================================================" << std::endl;
			TR.Info  <<  "   Stage 3                                                         " << std::endl;
			TR.Info  <<  "   Mainmintrial for " << stage3_cycles << std::endl;
			moves::SequenceMoverOP phase3_cycle( new moves::SequenceMover() );
			phase3_cycle->add_mover( moveset_phase3_ );
			phase3_cycle->add_mover( pack_rottrial_ );

			moves::JumpOutMoverOP phase3_min( new moves::JumpOutMover( phase3_cycle, min_mover_, get_scorefxn(), 15.0 ) );

			moves::TrialMoverOP phase3_trial( new moves::TrialMover( phase3_min, mc_ ) );

			moves::RepeatMoverOP full_cycle_phase3_;
			full_cycle_phase3_ = moves::RepeatMoverOP( new moves::RepeatMover( phase3_trial, int( core::Real(stage3_cycles) * basic::options::option[ basic::options::OptionKeys::relax::cycle_ratio])  ) );

			full_cycle_phase3_->apply( pose );

			checkpoints_.checkpoint( pose, "stage_3", get_current_tag(), true );

		}
	}

	// output_debug_structure( pose, "rl_stage3" );
	mc_->recover_low( pose );

	// add scores to map for output
	std::map < std::string, core::Real > score_map;
	core::io::raw_data::ScoreMap::nonzero_energies( score_map, get_scorefxn(), pose );
	if ( get_native_pose()!=nullptr ) {
		score_map["rms"] = CA_rmsd( *get_native_pose(), pose );
	}

	get_scorefxn()->show( TR.Debug , pose );
	TR << std::endl;

	// cache the score map to the pose
	// why does this obliterate any scores that were already there?
	using namespace basic::datacache;
	pose.data().set(CacheableDataType::SCORE_MAP, DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData(score_map) ));


	(*get_scorefxn())(pose);
}

std::string
ClassicRelax::get_name() const {
	return "ClassicRelax";
}

}
} // namespace protocols
