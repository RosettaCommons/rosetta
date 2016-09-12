// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/AnchoredDesign/AnchorMoversData.cc
/// @brief AnchorMoversData methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/anchored_design/AnchorMoversData.hh>
#include <protocols/anchored_design/Anchor.hh>

// Project Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID.hh> //used for movemap to set omega angle false
#include <core/id/types.hh> //used for movemap to set omega angle false

#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/chemical/ResidueType.hh>

#include <core/conformation/Conformation.hh>

#include <core/id/AtomID.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// C++ Headers
#include <string>

// option key includes
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.AnchoredDesign.AnchorMoversData" );

namespace protocols {
namespace anchored_design {

std::string const EMPTY_STRING("");
core::Real const VDW_WEIGHT(2.0);
core::Real const CHAINBREAK_WEIGHT(10.0);
std::string const INTERFACE_CALC("AnchoredDesign_InterfaceNeighborDefinitionCalculator");
std::string const NEIGHBORHOOD_CALC("AnchoredDesign_NeighborhoodByDistanceCalculator");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////AnchorMoversData//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////mechanicals (ctor, dtor)////////////////////////////////////////////////////////////

/// @brief empty, useless ctor.  You'll need to manually set all the data later.  Why did you make me waste my time writing this ctor?
protocols::anchored_design::AnchorMoversData::AnchorMoversData() :
	utility::pointer::ReferenceCount(),
	anchor_( /* 0 */ ),
	anchor_loop_index_( 0 ),
	movemap_fa_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	movemap_cen_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	//loops_and_fa_mms_( ), //er, not really possible to initialize
	//loops_and_cen_mms_( ),
	//loops_( ),
	fragset_( /* 0 */ ),
	task_factory_( /* 0 */ ),
	late_factory_( /* 0 */ ),
	fullatom_scorefunction_( /* 0 */ ),
	centroid_scorefunction_( /* 0 */ ),
	centroid_scorefunction_min_( /* 0 */ ),
	interface_calc_( INTERFACE_CALC ),
	neighborhood_calc_( NEIGHBORHOOD_CALC ),
	akash_dyepos_( 0 ),
	unbound_mode_( false ),
	anchor_via_constraints_( false ),
	VDW_weight_( VDW_WEIGHT ),
	chainbreak_weight_( CHAINBREAK_WEIGHT ),
	allow_anchor_repack_( false ),
	resfile_1_( EMPTY_STRING ),
	resfile_2_( EMPTY_STRING ),
	loop_file_( EMPTY_STRING ),
	frag3_( EMPTY_STRING ),
	no_frags_( false ),
	anchor_noise_constraints_mode_( false ),
	super_secret_fixed_interface_mode_( false )
{}

/// @details constructor takes an anchor and loop object and sets up internals reasonably; options boolean is optionally optional.  If you use this constructor, you will need to later manually give fragments to the AnchorMoversData object or use one of its fragments functions to determine what type of fragments to use.
protocols::anchored_design::AnchorMoversData::AnchorMoversData(
	protocols::anchored_design::AnchorCOP anchor,
	protocols::loops::Loops const & loops,
	bool const options) :
	utility::pointer::ReferenceCount(),
	anchor_(anchor), //redundant, given that set_loops_and_anchor takes care of it
	anchor_loop_index_( 0 ),
	movemap_fa_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	movemap_cen_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	//loops_and_fa_mms_( ), //er, not really possible to initialize
	//loops_and_cen_mms_( ),
	//loops_( ),
	fragset_( /* 0 */ ),
	task_factory_( /* 0 */ ),
	late_factory_( /* 0 */ ),
	fullatom_scorefunction_( /* 0 */ ),
	centroid_scorefunction_( /* 0 */ ),
	centroid_scorefunction_min_( /* 0 */ ),
	interface_calc_( INTERFACE_CALC ),
	neighborhood_calc_( NEIGHBORHOOD_CALC ),
	akash_dyepos_( 0 ),
	unbound_mode_( false ),
	anchor_via_constraints_( false ),
	VDW_weight_( VDW_WEIGHT ),
	chainbreak_weight_( CHAINBREAK_WEIGHT ),
	allow_anchor_repack_( false ),
	resfile_1_( EMPTY_STRING ),
	resfile_2_( EMPTY_STRING ),
	loop_file_( EMPTY_STRING ),
	frag3_( EMPTY_STRING ),
	no_frags_( false ),
	anchor_noise_constraints_mode_( false ),
	super_secret_fixed_interface_mode_( false )
{
	//TR << "loops/anchor ctor" << std::endl;


	//read commandline options
	if ( options ) read_options();

	//set up loops and anchor internal data
	set_loops_and_anchor(anchor, loops);
	//get other defaults set as necessary
	set_unset_defaults();

}

/// @details use pose and option system to set up loops and create anchor, then call set_loops_and_anchor and set_unset_defaults to set up internal data to reasonable defaults, with or without option system help
protocols::anchored_design::AnchorMoversData::AnchorMoversData( core::pose::Pose const & pose ) :
	utility::pointer::ReferenceCount(),
	anchor_( /* 0 */ ),
	anchor_loop_index_( 0 ),
	movemap_fa_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	movemap_cen_all_( core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ) ),
	//loops_and_fa_mms_( ), //er, not really possible to initialize
	//loops_and_cen_mms_( ),
	//loops_( ),
	fragset_( /* 0 */ ),
	task_factory_( /* 0 */ ),
	late_factory_( /* 0 */ ),
	fullatom_scorefunction_( /* 0 */ ),
	centroid_scorefunction_( /* 0 */ ),
	centroid_scorefunction_min_( /* 0 */ ),
	interface_calc_( INTERFACE_CALC ),
	neighborhood_calc_( NEIGHBORHOOD_CALC ),
	akash_dyepos_( 0 ),
	unbound_mode_( false ),
	anchor_via_constraints_( false ),
	VDW_weight_( VDW_WEIGHT ),
	chainbreak_weight_( CHAINBREAK_WEIGHT ),
	allow_anchor_repack_( false ),
	resfile_1_( EMPTY_STRING ),
	resfile_2_( EMPTY_STRING ),
	loop_file_( EMPTY_STRING ),
	frag3_( EMPTY_STRING ),
	no_frags_( false ),
	anchor_noise_constraints_mode_( false ),
	super_secret_fixed_interface_mode_( false )

{
	//TR << "pose ctor" << std::endl;

	//read commandline options; you don't get a choice with this ctor, because we need the loops file and anchor spec
	read_options();

	//set up loops and anchor internal data

	// read anchor file - automatically reads anchor file name from options system
	anchor_ = protocols::anchored_design::AnchorCOP( protocols::anchored_design::AnchorOP( new protocols::anchored_design::Anchor(pose) ) );

	// read loops file
	loops_ = protocols::loops::Loops( loop_file_ );
	//we do not need to check if the anchor and loop are compatible here; pick_new_cutpoints will take care of it later

	//this will overwrite anchor_ and loops_ with modified versions of themselves, and fill in movemaps, etc
	set_loops_and_anchor(anchor_, loops_);

	//get other defaults set as necessary
	set_unset_defaults();

	//handle fragments - can touch other data, thus needs to occur last
	autogenerate_frags(pose);

}

protocols::anchored_design::AnchorMoversData::~AnchorMoversData()= default;

/// @brief copy ctor
AnchorMoversData::AnchorMoversData( AnchorMoversData const & rhs ) :
	utility::pointer::ReferenceCount(rhs),
	interface_calc_( INTERFACE_CALC ), //const; must be initialized here
	neighborhood_calc_( NEIGHBORHOOD_CALC )//const; must be initialized here
{
	*this = rhs;
}

/// @brief assignment operator
AnchorMoversData & AnchorMoversData::operator=( AnchorMoversData const & rhs ){

	//abort self-assignment
	if ( this == &rhs ) return *this;

	//some data does not need to be copied; handled by set_loops_and_anchor below; required because loops_and_x_mms_ objects cannot be cloned
	//anchor_ = ;
	//anchor_loop_index_ = rhs.anchor_loop_index_;
	//movemap_fa_all_ = rhs.movemap_fa_all()->clone();
	//movemap_cen_all_ = rhs.movemap_cen_all()->clone();
	//no way to copy for these next two objects; handled at end
	//loops_and_fa_mms_ = rhs.loops_and_fa_mms_.clone();
	//loops_and_cen_mms_ = rhs.loops_and_cen_mms_.clone();
	//loops_ = rhs.loops(); //returns by value
	fragset_ = rhs.get_frags()->clone();
	task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*(rhs.get_task_factory())) );
	late_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*(rhs.get_late_factory())) );
	fullatom_scorefunction_ = rhs.get_fullatom_scorefunction()->clone();
	centroid_scorefunction_ = rhs.get_centroid_scorefunction()->clone();
	centroid_scorefunction_min_ = rhs.get_centroid_scorefunction_min()->clone();
	akash_dyepos_ = rhs.get_akash_dyepos();
	unbound_mode_ = rhs.get_unbound_mode();
	anchor_via_constraints_ = rhs.get_anchor_via_constraints();
	VDW_weight_ = rhs.get_VDW_weight();
	chainbreak_weight_ = rhs.get_chainbreak_weight();
	allow_anchor_repack_ = rhs.get_allow_anchor_repack();
	resfile_1_ = rhs.get_resfile_1();
	resfile_2_ = rhs.get_resfile_2();
	loop_file_ = rhs.get_loop_file();
	frag3_ = rhs.get_frag3();
	no_frags_ = rhs.get_no_frags();
	anchor_noise_constraints_mode_ = rhs.get_anchor_noise_constraints_mode();
	super_secret_fixed_interface_mode_ = rhs.get_super_secret_fixed_interface_mode();

	set_loops_and_anchor(rhs.anchor_, rhs.loops());

	return *this;
}

protocols::anchored_design::AnchorMoversDataOP protocols::anchored_design::AnchorMoversData::clone() const {
	return protocols::anchored_design::AnchorMoversDataOP( new AnchorMoversData(*this) );
}

/// @brief randomly reset loop cutpoints.  Useful only when starting structure is well-closed.  Best for MPI-style runs
void protocols::anchored_design::AnchorMoversData::pick_new_cutpoints(bool reset_always){
	TR << "checking/(re)setting loop cutpoints; loops become:" << std::endl;

	core::Size const numloop = num_loops();
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size const loopstart(loop(i).start()), loopend(loop(i).stop()), loopcut(loop(i).cut());
		if ( reset_always || (loopcut < loopstart) || (loopcut > loopend) || ((anchor_->start() <= loopcut) && (anchor_->end() >= loopcut)) ) {
			core::Size cut(pick_new_cutpoint(loopstart, loopend));
			loops_and_fa_mms_[ i ].key1().set_cut(cut);
			loops_and_cen_mms_[ i ].key1().set_cut(cut);
		}

		TR << loop(i) << std::endl;
	}

}

/// @brief randomly reset just one cutpoint; used by pick_new_cutpoints
core::Size protocols::anchored_design::AnchorMoversData::pick_new_cutpoint( core::Size const loopstart, core::Size const loopend ){
	core::Size newcutpoint(0);
	do{
		newcutpoint = (loopstart) + int( numeric::random::rg().uniform()*(loopend-loopstart+1) );
	} while( ((anchor_->start()) <= newcutpoint) && ((anchor_->end()) >= newcutpoint) ); //the cutpoint is in the anchor
	return newcutpoint;
}

///////////////////set functions if defaults aren't acceptable//////////////////////
//set fragments object
void protocols::anchored_design::AnchorMoversData::set_frags( core::fragment::FragSetOP in ) { fragset_ = in; }

//set packertask factory
void protocols::anchored_design::AnchorMoversData::set_task_factory( core::pack::task::TaskFactoryOP in ) { task_factory_ = in; }

//set fullatom scorefunction
void protocols::anchored_design::AnchorMoversData::set_fullatom_scorefunction( core::scoring::ScoreFunctionOP in )
{ fullatom_scorefunction_ = in; }

//set centroid scorefunction
void protocols::anchored_design::AnchorMoversData::set_centroid_scorefunction( core::scoring::ScoreFunctionOP in )
{ centroid_scorefunction_ = in; }

/// @details set up kinematics' loops and anchors; these are combined because loop setup depends on anchor
void protocols::anchored_design::AnchorMoversData::set_loops_and_anchor( protocols::anchored_design::AnchorCOP anchor,
	protocols::loops::Loops loops)
{
	anchor_ = anchor;

	//get loops and movemaps set up - most work is encoded in other functions
	loops.sequential_order();
	input_loops_into_tuples(loops);
	locate_anchor_loop();
	setup_movemaps();
	loops_ = loops; //ordered

	return;
}

/////////////////loops and movemap functions//////////////////////////////////////////////
//accessor for anchored loop
protocols::loops::Loop const & protocols::anchored_design::AnchorMoversData::anchored_loop() const
{ return loops_and_fa_mms_[ anchor_loop_index_ ].key1(); } //ISSUE: WHY NOT USE CEN HERE?  resolve w/ vectors?

//access for movemap that covers all loops; centroid
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_cen_all() const {return movemap_cen_all_;}

//access for movemap that covers all loops; fullatom may allow anchor movement w/constraints
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_fa_all() const {return movemap_fa_all_;}

//accessor for a loop //ISSUE: WHY NOT USE CEN HERE?  resolve w/ vectors?
protocols::loops::Loop const & protocols::anchored_design::AnchorMoversData::loop( core::Size i ) const
{ return loops_and_fa_mms_[ i ].key1(); }

//accessor for omega-variable movemap (most movers); fullatom may allow anchor movement w/constraints
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_fa( core::Size i ) const
{ return loops_and_fa_mms_[ i ].key2(); }

//accessor for omega-fixed movemap (appropriate for CCD movers); fullatom may allow anchor movement w/constraints
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_fa_omegafixed( core::Size i ) const
{ return loops_and_fa_mms_[ i ].key3(); }

//accessor for omega-variable movemap (most movers); centroid
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_cen( core::Size i ) const
{ return loops_and_cen_mms_[ i ].key2(); }

//accessor for omega-fixed movemap (appropriate for CCD movers); centroid
core::kinematics::MoveMapOP protocols::anchored_design::AnchorMoversData::movemap_cen_omegafixed( core::Size i ) const
{ return loops_and_cen_mms_[ i ].key3(); }

//accessor for loops object
protocols::loops::Loops const & protocols::anchored_design::AnchorMoversData::loops() const
{ return loops_; }


/////////////////////////////anchor functions//////////////////////////////////////////////
/// @details access for anchor start
core::Size protocols::anchored_design::AnchorMoversData::anchor_start() const { return anchor_->start(); }

/// @details access for anchor end
core::Size protocols::anchored_design::AnchorMoversData::anchor_end() const { return anchor_->end(); }

////////////////////////scfxn, packertask, fragments accessors//////////////////////////////
//access fragments object
core::fragment::FragSetCOP protocols::anchored_design::AnchorMoversData::get_frags() const
{ return fragset_; }
//access packertask factory
core::pack::task::TaskFactoryCOP protocols::anchored_design::AnchorMoversData::get_task_factory() const
{ return task_factory_; }
//access packertask factory
core::pack::task::TaskFactoryCOP protocols::anchored_design::AnchorMoversData::get_late_factory() const
{ return late_factory_; }
//access fullatom scorefunction
core::scoring::ScoreFunctionOP protocols::anchored_design::AnchorMoversData::get_fullatom_scorefunction() const
{ return fullatom_scorefunction_; }
//access centroid scorefunction
core::scoring::ScoreFunctionOP protocols::anchored_design::AnchorMoversData::get_centroid_scorefunction() const
{ return centroid_scorefunction_; }
//access centroid scorefunction for minimization
core::scoring::ScoreFunctionOP protocols::anchored_design::AnchorMoversData::get_centroid_scorefunction_min() const
{ return centroid_scorefunction_min_; }
/// @details runs the member factory to create a task
core::pack::task::PackerTaskOP protocols::anchored_design::AnchorMoversData::get_task( core::pose::Pose const & pose) const
{ return task_factory_->create_task_and_apply_taskoperations( pose ); }


////////////////////////private functions generate internal data from input and defaults///////////////
/// @details rearranges input loops data structure into the class's internal data structure
void protocols::anchored_design::AnchorMoversData::input_loops_into_tuples(protocols::loops::Loops const & loops )
{
	loops_and_fa_mms_.clear();
	loops_and_cen_mms_.clear();
	for ( auto const & loop : loops ) {
		//instantiate tuple
		loops_and_fa_mms_.push_back( Loop_mm_tuple(
			loop,
			core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ),
			core::kinematics::MoveMapOP( new core::kinematics::MoveMap() )
			) );

		loops_and_cen_mms_.push_back( Loop_mm_tuple(
			loop,
			core::kinematics::MoveMapOP( new core::kinematics::MoveMap() ),
			core::kinematics::MoveMapOP( new core::kinematics::MoveMap() )
			) );
	}
}

/// @details determines which loop contains the anchor (this loop may be treated differently)
void protocols::anchored_design::AnchorMoversData::locate_anchor_loop()
{
	core::Size const numloop = num_loops();
	core::Size const anchorstart(anchor_->start()), anchorend(anchor_->end());
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size const loopstart(loop(i).start()), loopend(loop(i).stop());
		if ( (loopstart < anchorstart) && (loopend > anchorend) ) {
			TR << "anchor start/end " << anchorstart << "/" << anchorend
				<< " fits within loop start/end " << loopstart << "/" << loopend << std::endl;
			anchor_loop_index_ = i;
			return;
		}// end if-anchor-is-in-this-loop
	}//end iterate over all loops

	if ( get_super_secret_fixed_interface_mode() ) {
		anchor_loop_index_ = 0; //I don't know what this will do, but it should cause bounds errors if something untoward happens (good except bad)
		TR << "anchor start/end " << anchorstart << "/" << anchorend
			<< " not within a loop; appropriate for super_secret_fixed_interface_mode" << std::endl;
		return;
	}

	//if we did not return above
	Error() << "Anchor does not reside completely within a loop.  Anchor start and end: " << anchorstart
		<< " " << anchorend << std::endl;
	utility_exit();
}

/// @details sets all movemaps for the class - paired ones and single movemap_all_, also sets loop_or_notloop_
void protocols::anchored_design::AnchorMoversData::setup_movemaps()
{
	movemap_fa_all_->clear();
	movemap_cen_all_->clear();
	core::Size const numloop = num_loops();

	//for each loop, modify that pair's mm and modify movemap_all_
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size loopstart(loop(i).start()), loopend(loop(i).stop());
		movemap_fa(i)->clear();
		movemap_fa_omegafixed(i)->clear();
		movemap_cen(i)->clear();
		movemap_cen_omegafixed(i)->clear();

		//set this loop's mm, and in movemap_all_
		for ( core::Size j = loopstart; j <= loopend; ++j ) {
			set_movemap( movemap_fa( i ), j );
			set_movemap( movemap_fa_omegafixed( i ), j, true );
			set_movemap( movemap_fa_all_, j );

			set_movemap( movemap_cen( i ), j );
			set_movemap( movemap_cen_omegafixed( i ), j, true );
			set_movemap( movemap_cen_all_, j );
		}

	}//for each loop

	if ( !get_super_secret_fixed_interface_mode() ) {
		fix_anchor( movemap_fa( anchor_loop_index_ ), false );
		fix_anchor( movemap_fa_omegafixed( anchor_loop_index_ ), false );
		fix_anchor( movemap_fa_all_, false );

		fix_anchor( movemap_cen( anchor_loop_index_ ), true );
		fix_anchor( movemap_cen_omegafixed( anchor_loop_index_ ), true );
		fix_anchor( movemap_cen_all_, true );
	}

	return;
}

void protocols::anchored_design::AnchorMoversData::set_movemap(core::kinematics::MoveMapOP movemap, core::Size seqpos, bool omega)
{
	//this if statement protects noncanonical dye residues; see AnchoredDesign application documentation.
	if ( seqpos == akash_dyepos_ ) return;
	using namespace core::id;
	movemap->set_bb(seqpos, true); //backbone mobile
	if ( omega ) movemap->set( TorsionID(seqpos, BB, omega_torsion), false ); //fixes omega angle
	movemap->set_chi(seqpos, true); // chi of loop residues
	//movemap->set_jump(false); // call in fix_anchor (more efficient)
}

void protocols::anchored_design::AnchorMoversData::fix_anchor( core::kinematics::MoveMapOP movemap, bool const centroid )
{
	//skip this function in unbound mode or super_secret_fixed_interface_mode_ or anchors_via_constraints mode
	if ( unbound_mode_ || get_super_secret_fixed_interface_mode() ) return;
	if ( anchor_via_constraints_ && !centroid ) { //if we have constraints, and this is NOT for the centroid phase...
		movemap->set_jump(1, true); //magic number: the anchor-target jump is 1; TODO: make this a named variable
		//do NOT re-fix the anchor residues
		return;
	}

	for ( core::Size i = anchor_->start(); i <= anchor_->end(); ++i ) {
		movemap->set_bb(i, false); //fix the backbones
		movemap->set_chi(i, false); // and the chi
	}

	using namespace core::id;
	movemap->set( TorsionID(anchor_->start(), BB, phi_torsion), true);
	movemap->set( TorsionID(anchor_->end(), BB, psi_torsion), true);
	const bool omega(movemap->get(TorsionID(anchor_->end()+1, BB, omega_torsion)));//are omegas free in this mm?
	movemap->set( TorsionID(anchor_->end(), BB, omega_torsion), omega);

	movemap->set_jump(false); //perhaps this call should be in set_movemap()?
}

void protocols::anchored_design::AnchorMoversData::set_unset_defaults()// std::string const & frag_file )
{
	set_unset_scorefunctions();
	set_unset_packertask_factory();
}

void protocols::anchored_design::AnchorMoversData::set_unset_scorefunctions()
{
	using namespace core::scoring;

	//handle centroid scorefunction
	centroid_scorefunction_ = core::scoring::ScoreFunctionOP( new ScoreFunction );
	centroid_scorefunction_->set_weight( env,         1.0 ); //no derivative
	centroid_scorefunction_->set_weight( cbeta,       1.0 ); //no derivative
	centroid_scorefunction_->set_weight( vdw,         VDW_weight_ ); //defines eval_atom_derivative
	centroid_scorefunction_->set_weight( pair,        1.0 ); //no derivative
	centroid_scorefunction_->set_weight( cenpack,     1.0 ); //no derivative
	centroid_scorefunction_->set_weight( rama,        1.0 ); //defines eval_dof_derivative
	centroid_scorefunction_->set_weight( chainbreak,  2.0 ); //defines eval_atom_derivative
	centroid_scorefunction_->set_weight( hbond_lr_bb, 1.0 );
	centroid_scorefunction_->set_weight( hbond_sr_bb, 1.0 );
	//DO NOT add constraints in centroid mode - unless all bb-bb constraints it will break due to csts using blind AtomIDs...
	//core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *centroid_scorefunction_ ); //protected if(option) internally
	TR << "Using default centroid scorefunction\n" << *centroid_scorefunction_ << std::flush;

	//handle centroid minimization scorefunction
	centroid_scorefunction_min_ = centroid_scorefunction_->clone();
	centroid_scorefunction_->set_weight( env,         0. ); //no derivative
	centroid_scorefunction_->set_weight( cbeta,       0. ); //no derivative
	centroid_scorefunction_->set_weight( pair,        0. ); //no derivative
	centroid_scorefunction_->set_weight( cenpack,     0. ); //no derivative


	fullatom_scorefunction_ = get_score_function();
	fullatom_scorefunction_->set_weight( chainbreak, chainbreak_weight_ );
	core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ ); //protected if(option) internally
	TR << "Using default fullatom scorefunction (TALARIS_2013 plus chainbreak @ a lot)\n"
		<< *fullatom_scorefunction_ << std::flush;

	return;
}


void protocols::anchored_design::AnchorMoversData::set_unset_packertask_factory()
{
	using namespace core::pack::task;
	using namespace core::pose::metrics;

	//set up interface-plus-neighbors-positions operation
	core::Size const numloop = num_loops();
	std::set< core::Size > loop_posns;
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size loopstart(loop(i).start()), loopend(loop(i).stop());
		for ( core::Size j = loopstart; j <= loopend; ++j ) {
			loop_posns.insert(j);
		}//for each residue in loop
	}//for each loop

	//std::string const interface_calc("AnchoredDesign_InterfaceNeighborDefinitionCalculator");
	//std::string const neighborhood_calc("AnchoredDesign_NeighborhoodByDistanceCalculator");
	//NOTE that these calculators ASSUME your interface is between chains 1 and 2, because most of the rest of AnchoredDesign carries the same assumption

	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( interface_calc_ ) ) {
		Warning() << "In AnchoredDesign, calculator " << interface_calc_ << " already exists.  "
			<< "Given the two-chain restriction, this is hopefully correct for your purposes" << std::endl;
	} else {
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2) ) ) );
	}

	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( neighborhood_calc_ ) ) {
		Warning() << "In AnchoredDesign, calculator " << neighborhood_calc_ << " already exists.  If you have multiple instances of AnchoredDesign with different loops coexisting in the same program, this is going to cause problems, because you will have the wrong loop definitions for determining what residues to pack" << std::endl;
	} else {
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc_, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) ) );
	}

	//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
	utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
	calcs_and_calcns.push_back(std::make_pair(interface_calc_, "interface_residues"));
	calcs_and_calcns.push_back(std::make_pair(neighborhood_calc_, "neighbors"));

	using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
	operation::TaskOperationOP rbcop( new RestrictByCalculatorsOperation( calcs_and_calcns ) );

	//command line operation
	operation::InitializeFromCommandlineOP initop( new operation::InitializeFromCommandline );

	//operation to protect anchor
	operation::PreventRepackingOP prop( new operation::PreventRepacking );
	if ( !allow_anchor_repack_ && !get_super_secret_fixed_interface_mode() ) {
		TR << "autogenerated TaskFactory will prevent repacking for anchor" << std::endl;
		for ( core::Size i(anchor_->start()); i<= anchor_->end(); ++i ) {
			prop->include_residue(i);
		}
	}

	//'early' Factory
	task_factory_ = core::pack::task::TaskFactoryOP( new TaskFactory() );
	//resfile operation
	if ( resfile_1_ != EMPTY_STRING ) {
		operation::ReadResfileOP rrop1( new operation::ReadResfile );
		rrop1->filename( resfile_1_ );
		task_factory_->push_back( rrop1 );
	}
	task_factory_->push_back( initop );
	task_factory_->push_back( rbcop );
	if ( !allow_anchor_repack_ ) task_factory_->push_back( prop );
	TR << "Using default TaskFactory.  Inits from command line and second resfile, then restricts "
		<< "to the loop/interface area" << std::endl;

	//'late' factory
	if ( (resfile_2_ != EMPTY_STRING) ) {
		late_factory_ = core::pack::task::TaskFactoryOP( new TaskFactory() );
		operation::ReadResfileOP rrop2( new operation::ReadResfile );
		rrop2->filename( resfile_2_ );//second resfile
		late_factory_->push_back( rrop2 );

		late_factory_->push_back( initop );
		late_factory_->push_back( rbcop );
		if ( !allow_anchor_repack_ ) late_factory_->push_back( prop );

		TR << "Using default late TaskFactory.  Inits from command line and second resfile, then restricts "
			<< "to the loop/interface area" << std::endl;
	} else {
		late_factory_ = task_factory_;
		TR << "Not using separate late TaskFactory." << std::endl;
	}
	return;
}

/// @details autogenerate_design_frags will use Andrew Ban's fragment picker to read the Vall and find fragments.  It looks for fragments of loop secondary structure (for now, just 3mers.).  This function is used in a design context and thus uses only secondary structure to define the fragments.
core::fragment::FragSetCOP protocols::anchored_design::AnchorMoversData::autogenerate_design_frags(){

	core::Size const frags_length(3); //magic number: 3mer fragments!!
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( frags_length ) );

	std::string ss_string(frags_length, 'L');
	core::fragment::FragDataOPs list;
	list =  core::fragment::picking_old::vall::pick_fragments_by_ss( ss_string, 4000, false ); //magic number: 4000 fragments

	core::Size const numloop = num_loops();
	//std::set< core::Size > loop_posns;
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size loopstart(loop(i).start()), loopend(loop(i).stop());
		for ( core::Size j = loopstart; j <= loopend - frags_length+1; ++j ) {
			if ( !( (j >anchor_->start() - frags_length) && (j <= anchor_->end()) ) ) {
				TR << "adding frame, start at " << j << " go for " << frags_length << " to " << j+frags_length << std::endl;
				core::fragment::FrameOP frame;
				frame = core::fragment::FrameOP( new core::fragment::Frame( j ) );
				frame->add_fragment( list );
				fragset->add( frame );
			}//not in anchor
		}//for each residue in loop
	}//for each loop

	TR << "recovering memory?" << std::endl;
	// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
	core::fragment::picking_old::FragmentLibraryManager::get_instance()->clear_Vall();

	set_frags(fragset);
	return fragset;
}//autogenerate_design_frags

/// @details autogenerate_constseq_frags will use Andrew Ban's fragment picker to read the Vall and find fragments.  It looks for fragments of loop secondary structure and known, fixed sequence.  The string argument is the whole pose sequence.
core::fragment::FragSetCOP protocols::anchored_design::AnchorMoversData::autogenerate_constseq_frags(std::string const & seq){

	core::Size const frags_length(3); //magic number: 3mer fragments!!
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( frags_length ) );

	std::string ss_string(frags_length, 'L');
	core::fragment::FragDataOPs list;

	core::Size const numloop = num_loops();
	//std::set< core::Size > loop_posns;
	for ( core::Size i = 1; i <= numloop; ++i ) {
		core::Size loopstart(loop(i).start()), loopend(loop(i).stop());
		for ( core::Size j = loopstart; j <= loopend - frags_length+1; ++j ) {
			if ( !( (j >anchor_->start() - frags_length) && (j <= anchor_->end()) ) ) {
				std::string const seqsubstr(seq, j-1, frags_length);
				TR << "adding frame, start at " << j << " go for " << frags_length << " to " << j+frags_length << " seq " << seqsubstr << std::endl;
				list =  core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_string, seqsubstr, 200, false ); //magic number: 200 fragments per position (not duplicated - this will be like robetta server fragments)

				core::fragment::FrameOP frame;
				frame = core::fragment::FrameOP( new core::fragment::Frame( j ) );
				frame->add_fragment( list );
				fragset->add( frame );
			}//not in anchor
		}//for each residue in loop
	}//for each loop

	TR << "recovering memory?" << std::endl;
	// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
	core::fragment::picking_old::FragmentLibraryManager::get_instance()->clear_Vall();

	set_frags(fragset);
	return fragset;
}//autogenerate_design_frags


/// @details autogenerate_frags will determine from this object's state what fragments to use, and generate fragments as needed.  The logic was originally at the executeable level, src/apps/pilot/smlewis/AnchoredDesign.cc:116-136 SVN 40529
void protocols::anchored_design::AnchorMoversData::autogenerate_frags( core::pose::Pose const & pose ){

	core::fragment::ConstantLengthFragSetOP fragset3mer;
	bool const frags_file( !(frag3_ == EMPTY_STRING) );//basic::options::option[ basic::options::OptionKeys::in::file::frag3].user() );
	//bool const no_frags( basic::options::option[ basic::options::OptionKeys::AnchoredDesign::no_frags ].value() );
	if ( frags_file && no_frags_ ) {
		utility_exit_with_message("you've specified a fragments file and requested no_frags - please choose only one");
	} else if ( frags_file ) {
		TR << "reading from fragments file " << frag3_ << std::endl;
		fragset3mer = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 3 ) );
		fragset3mer->read_fragment_file( frag3_ );//basic::options::option[ basic::options::OptionKeys::in::file::frag3].value() );
		set_frags( fragset3mer );
	} else if ( !no_frags_ ) {
		//strange logic: are we designing or not?  To determine, we must run a TaskFactory to get a Task and ask it!
		if ( get_task_factory()->create_task_and_apply_taskoperations(pose)->design_any() ) {
			TR << "creating sequence-generic loop fragments (LLLLL...)" << std::endl;
			autogenerate_design_frags();
		} else {
			TR << "creating constant-sequence fragments for loops" << std::endl;
			autogenerate_constseq_frags(pose.sequence());
		}
	} else { TR << "using no fragments" << std::endl; }
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////option system replacement/////////////////////////////
/// @brief dye position used in dye modeling publication
void protocols::anchored_design::AnchorMoversData::set_akash_dyepos(core::Size const akash_dyepos){ akash_dyepos_ = akash_dyepos;}
/// @brief used for unbound mode
void protocols::anchored_design::AnchorMoversData::set_unbound_mode(bool unbound_mode){ unbound_mode_ = unbound_mode;}
/// @brief used to test anchoring via constraints
void protocols::anchored_design::AnchorMoversData::set_anchor_via_constraints(bool anchor_via_constraints){ anchor_via_constraints_ = anchor_via_constraints;}
/// @brief VDW weight in centroid scorefunction
void protocols::anchored_design::AnchorMoversData::set_VDW_weight(core::Real VDW_weight){ VDW_weight_ = VDW_weight;}
/// @brief chainbreak weight in fullatom scorefunction
void protocols::anchored_design::AnchorMoversData::set_chainbreak_weight(core::Real chainbreak_weight){ chainbreak_weight_ = chainbreak_weight;}
/// @brief allow anchor to repack
void protocols::anchored_design::AnchorMoversData::set_allow_anchor_repack(bool allow_anchor_repack){ allow_anchor_repack_ = allow_anchor_repack;}
/// @brief resfile for design
void protocols::anchored_design::AnchorMoversData::set_resfile_1(std::string const & resfile_1){ resfile_1_ = resfile_1;}
/// @brief later-stage resfile if desired
void protocols::anchored_design::AnchorMoversData::set_resfile_2(std::string const & resfile_2){ resfile_2_ = resfile_2;}
/// @brief loop file
//void protocols::anchored_design::AnchorMoversData::set_loop_file(std::string const & loop_file){ loop_file_ = loop_file;}
/// @brief loop file
void protocols::anchored_design::AnchorMoversData::set_frag3(std::string const & frag3){ frag3_ = frag3;}
/// @brief do not use fragments?
void protocols::anchored_design::AnchorMoversData::set_no_frags(bool const no_frags) { no_frags_= no_frags;}
/// @brief special anchor_noise_constraints_mode
void protocols::anchored_design::AnchorMoversData::set_anchor_noise_constraints_mode(bool const anchor_noise_constraints_mode) { anchor_noise_constraints_mode_= anchor_noise_constraints_mode;}
/// @brief special super_secret_fixed_interface_mode
void protocols::anchored_design::AnchorMoversData::set_super_secret_fixed_interface_mode(bool const super_secret_fixed_interface_mode) { super_secret_fixed_interface_mode_= super_secret_fixed_interface_mode;}

/// @brief dye position used in dye modeling publication
core::Size protocols::anchored_design::AnchorMoversData::get_akash_dyepos() const { return akash_dyepos_;}
/// @brief used for unbound mode
bool protocols::anchored_design::AnchorMoversData::get_unbound_mode() const { return unbound_mode_;}
/// @brief used to test anchoring via constraints
bool protocols::anchored_design::AnchorMoversData::get_anchor_via_constraints() const { return anchor_via_constraints_;}
/// @brief VDW weight in centroid scorefunction
core::Real protocols::anchored_design::AnchorMoversData::get_VDW_weight() const { return VDW_weight_;}
/// @brief chainbreak weight in fullatom scorefunction
core::Real protocols::anchored_design::AnchorMoversData::get_chainbreak_weight() const { return chainbreak_weight_;}
/// @brief allow anchor to repack
bool protocols::anchored_design::AnchorMoversData::get_allow_anchor_repack() const { return allow_anchor_repack_;}
/// @brief resfile for design
std::string const & protocols::anchored_design::AnchorMoversData::get_resfile_1() const { return resfile_1_;}
/// @brief later-stage resfile if desired
std::string const & protocols::anchored_design::AnchorMoversData::get_resfile_2() const { return resfile_2_;}
/// @brief loop file
std::string const & protocols::anchored_design::AnchorMoversData::get_loop_file() const { return loop_file_;}
/// @brief frag3 file
std::string const & protocols::anchored_design::AnchorMoversData::get_frag3() const { return frag3_;}
/// @brief do not use fragments?
bool protocols::anchored_design::AnchorMoversData::get_no_frags() const { return no_frags_;}
/// @brief special anchor_noise_constraints_mode
bool protocols::anchored_design::AnchorMoversData::get_anchor_noise_constraints_mode() const { return anchor_noise_constraints_mode_;}
/// @brief special super_secret_fixed_interface_mode
bool protocols::anchored_design::AnchorMoversData::get_super_secret_fixed_interface_mode() const { return super_secret_fixed_interface_mode_;}

void protocols::anchored_design::AnchorMoversData::read_options() {

	TR << "initializing from options system" << std::endl;

	using namespace basic::options;
	using namespace basic::options::OptionKeys::AnchoredDesign;

	//core::Size akash_dyepos_;
	if ( option[akash::dyepos].user() ) {
		akash_dyepos_ =  option[akash::dyepos].value();
	} else {
		akash_dyepos_ = 0;
	}

	//bool unbound_mode_;
	unbound_mode_ = option[ unbound_mode ].value();

	//bool anchor_via_constraints_;
	anchor_via_constraints_ = option[ testing::anchor_via_constraints ].value();

	//core::Real VDW_weight_;
	VDW_weight_ = option[testing::VDW_weight].value();

	//core::Real chainbreak_weight_;
	chainbreak_weight_ = option[chainbreak_weight].value();

	//bool allow_anchor_repack_;
	allow_anchor_repack_ = option[allow_anchor_repack].value();

	//std::string resfile_1_;
	if ( option[OptionKeys::packing::resfile].user() ) {
		resfile_1_ = option[OptionKeys::packing::resfile].value().at(1);
	} else {
		resfile_1_ = EMPTY_STRING; //empty string flag for don't use resfile
	}

	//std::string resfile_2_;
	if ( option[OptionKeys::packing::resfile].user() && (option[ OptionKeys::packing::resfile ].value().size() > 1) ) {
		resfile_2_ = option[OptionKeys::packing::resfile].value().at(2);
	} else {
		resfile_2_ = EMPTY_STRING; //empty string flag for don't use resfile
	}

	//std::string loop_file_;
	if ( option[OptionKeys::loops::loop_file].user() ) {
		loop_file_ = option[OptionKeys::loops::loop_file].value().at(1);
	}

	//std::string frag3_;
	if ( option[OptionKeys::in::file::frag3].user() ) {
		frag3_ = option[OptionKeys::in::file::frag3].value();
	}

	//bool no_frags_;
	no_frags_ = option[ no_frags ].value();

	//bool anchor_noise_constraints_mode_;
	anchor_noise_constraints_mode_ = option[ testing::anchor_noise_constraints_mode ].value();

	super_secret_fixed_interface_mode_ = option[ testing::super_secret_fixed_interface_mode ].value();

	return;
}

/// @brief get string name for neighborhood_calc_
std::string const & protocols::anchored_design::AnchorMoversData::neighborhood_calc() const { return neighborhood_calc_;}
/// @brief get string name for interface_calc_
std::string const & protocols::anchored_design::AnchorMoversData::interface_calc() const { return interface_calc_;}

/// @details This function automatically generates constraints for an anchor from and for the given pose.  It is intended for use only with single-residue anchors.  The constraints are between the residue CA and the four closest cross-chain CA.  They are parabolic constraints with scores of 0.5 units at 1 Angstrom deviation.
void protocols::anchored_design::AnchorMoversData::anchor_noise_constraints_setup( core::pose::Pose & pose )
{
	//exclusion condition: can't use this mode with multi-residue anchor
	if ( anchor_start() != anchor_end() ) {
		utility_exit_with_message("You can only use one-residue anchors with anchor_noise_constraints_mode");
	}

	//Calculate CA-CA distances between anchor residue and all other atoms
	core::Size const anchor(anchor_start());
	core::Size const CA(pose.residue_type(anchor).atom_index("CA"));
	core::id::AtomID const anchor_ID(core::id::AtomID(CA, anchor));
	core::Vector const & anchor_xyz(pose.xyz(anchor_ID));

	//we'll need to know which is the opposite chain (we are assuming there is only one)
	core::Size const chain1end(pose.conformation().chain_end(1));
	if ( pose.conformation().num_chains() != 2 ) { //2 is the value we expect
		utility_exit_with_message("cannot use anchor_noise_constraints_mode with more than two chains");
	}
	core::Size const oppchain( anchor > chain1end ? 1 : 2);

	//Data object type to hold resid, distance pairs
	typedef std::pair< core::Real, core::id::AtomID> dist_resid;
	std::set< dist_resid > distances;

	//loop to generate distances
	for ( core::Size i(pose.conformation().chain_begin(oppchain)), end(pose.conformation().chain_end(oppchain)); i<=end; ++i ) {
		//if the friend function version of distance worked, I bet I could fit this all on one line of code
		core::id::AtomID const target_atomID(CA, i);
		distances.insert(std::make_pair( anchor_xyz.distance(pose.xyz(target_atomID)), target_atomID) );
	}

	//generate N constraints; add directly to pose
	core::Size const this_many_constraints(4);
	core::Real const sd(sqrt(0.5));
	auto opp_CA(distances.begin());
	for ( core::Size i(1); i<=this_many_constraints; ++i, ++opp_CA ) {
		using namespace core::scoring::constraints;

		core::scoring::func::HarmonicFuncOP anchor_func( new core::scoring::func::HarmonicFunc(opp_CA->first, sd) );

		TR << "anchor_noise_constraints constraint, opp_CA, center, sd: " << opp_CA->second << " " << opp_CA->first << " " << sd << std::endl;

		pose.add_constraint(core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new AtomPairConstraint(anchor_ID, opp_CA->second, anchor_func) ) ));
		//++opp_CA;
	}

	//report constraints chosen
	pose.constraint_set()->show_definition(TR, pose);

	//add constraints to cen and fa scorefunctions
	fullatom_scorefunction_->set_weight( core::scoring::atom_pair_constraint, 1 );
	centroid_scorefunction_->set_weight( core::scoring::atom_pair_constraint, 1 );
	centroid_scorefunction_min_->set_weight( core::scoring::atom_pair_constraint, 1 );

	//free other widgets in AnchorMoversData to allow constraint freedom
	movemap_fa_all_->set_jump(1, true);
	movemap_cen_all_->set_jump(1, true);
	core::Size const numloop = num_loops();
	for ( core::Size i = 1; i <= numloop; ++i ) {
		movemap_fa(i)->set_jump(1, true);
		movemap_fa_omegafixed(i)->set_jump(1, true);
		movemap_cen(i)->set_jump(1, true);
		movemap_cen_omegafixed(i)->set_jump(1, true);
	}

	return;
}

}//namespace anchored_design
}//namespace protocols
