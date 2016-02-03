// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mp_quick_relax.cc
/// @brief   Do a quick relax run for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <core/kinematics/MoveMap.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <protocols/membrane/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/membrane/util.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.mp_quick_relax" );

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
using namespace core::kinematics;
using namespace core::scoring;
using namespace protocols::membrane;

class TestQuickRelaxMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: dih angle_max = 1, nmoves = "nres", movemap = bb and all chi
	TestQuickRelaxMover();

	/// @brief Custom Constructor
	/// @details nmoves is a string because it can either be "nres" or an integer
	TestQuickRelaxMover( Real angle_max, std::string nmoves );

	/// @brief Custom Constructor
	TestQuickRelaxMover( Real angle_max, std::string nmoves, MoveMapOP movemap );

	/// @brief Copy Constructor
	TestQuickRelaxMover( TestQuickRelaxMover const & src );

	/// @brief Assignment Operator
	TestQuickRelaxMover & operator = ( TestQuickRelaxMover const & src );

	/// @brief Destructor
	virtual ~TestQuickRelaxMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	//  /// @brief Create a Clone of this mover
	// virtual protocols::moves::MoverOP clone() const;
	//
	//  /// @brief Create a Fresh Instance of this Mover
	// virtual protocols::moves::MoverOP fresh_instance() const;
	//
	//  /// @brief Pase Rosetta Scripts Options for this Mover
	// void parse_my_tag(
	//       utility::tag::TagCOP tag,
	//       basic::datacache::DataMap &,
	//       protocols::filters::Filters_map const &,
	//       protocols::moves::Movers_map const &,
	//       core::pose::Pose const &
	//       );

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TestQuickRelaxMover)
	virtual std::string get_name() const;

	/// @brief Flip the downstream partner in the membrane
	virtual void apply( Pose & pose );

	/// @brief Run AddMembraneMover before?
	/// @details If you want to keep your anchor point for MEM, then pick no
	void add_membrane_again( bool yesno );

	/// @brief Run MembranePositionFromTopology again?
	/// @details Will change the starting membrane position
	void membrane_from_topology( bool yesno );

	/// @brief Optimize membrane position before relax?
	void optimize_membrane( bool yesno );

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Set default values
	void set_defaults();

	/// @brief Initialize from commandline
	void init_from_cmd();

	/// @brief Initialize from commandline
	utility::vector1< bool > get_repack_residues( Pose & pose, core::SSize center1, core::SSize center2, core::SSize halfrange );

private: // data

	/// @brief Native
	PoseOP native_;

	/// @brief Maximum allowed dihedral angle change for Small and ShearMover
	Real angle_max_;

	/// @brief Number of moves Small and ShearMover can make
	/// @details moves_ is a string and can take 'nres' as well as a number
	///    nmoves_ is the actual number that is taken after conversion
	std::string moves_;
	Size nmoves_;

	/// @brief Movemap for Small and ShearMover
	MoveMapOP movemap_;

	/// @brief Scorefxn
	ScoreFunctionOP sfxn_;

	/// @brief constraint filename
	std::string cst_file_;

	/// @brief constraint weight
	core::Real cst_weight_;

	/// @brief Run AddMembraneMover again?
	/// @details This is stupid: we need this as a workaround for attaching the
	///   MEM at the correct anchor point, neither sliding the jump nor
	///   removing the MEM and re-attaching it works. :(
	bool addmem_;

	/// @brief Run MembranePositionFromTopology?
	bool mem_from_topo_;


	/// @brief Optimize membrane before?
	bool opt_mem_;

	/// @brief Additional round of repack which packs all residues simultaneously
	bool repack_again_;

};

//////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< TestQuickRelaxMover > TestQuickRelaxMoverOP;

//////////////////////////////////////////////////////////////////////

/// @brief Default Constructor
/// @details Defaults: dih angle_max = 1, nmoves = nres, movemap = bb and all chi
TestQuickRelaxMover::TestQuickRelaxMover() :
	protocols::moves::Mover(),
	movemap_( new MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
TestQuickRelaxMover::TestQuickRelaxMover( Real angle_max, std::string nmoves ) :
	protocols::moves::Mover(),
	movemap_( new MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();

	angle_max_ = angle_max;
	moves_ = nmoves;
}

/// @brief Custom constructor
TestQuickRelaxMover::TestQuickRelaxMover( Real angle_max, std::string nmoves, MoveMapOP movemap ) :
	protocols::moves::Mover(),
	movemap_( new MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();

	angle_max_ = angle_max;
	moves_ = nmoves;
	movemap_ = movemap;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TestQuickRelaxMover::TestQuickRelaxMover( TestQuickRelaxMover const & src ) : protocols::moves::Mover( src ),
	native_( src.native_ ),
	angle_max_( src.angle_max_ ),
	moves_( src.moves_ ),
	nmoves_( src.nmoves_ ),
	movemap_( src.movemap_ ),
	sfxn_( src.sfxn_ ),
	cst_file_( src.cst_file_ ),
	cst_weight_( src.cst_weight_ ),
	addmem_( src.addmem_ ),
	mem_from_topo_ ( src.mem_from_topo_ ),
	opt_mem_( src.opt_mem_ )
{}

/// @brief Assignment Operator
TestQuickRelaxMover & TestQuickRelaxMover::operator = ( TestQuickRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new TestQuickRelaxMover( *this ) );
}

/// @brief Destructor
TestQuickRelaxMover::~TestQuickRelaxMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//TestQuickRelaxMover::clone() const {
// return ( protocols::moves::MoverOP( new TestQuickRelaxMover( *this ) ) );
//}
//
// /// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//TestQuickRelaxMover::fresh_instance() const {
// return protocols::moves::MoverOP( new TestQuickRelaxMover() );
//}
//
// /// @brief Pase Rosetta Scripts Options for this Mover
//void
//TestQuickRelaxMover::parse_my_tag(
//        utility::tag::TagCOP /*tag*/,
//        basic::datacache::DataMap &,
//        protocols::filters::Filters_map const &,
//        protocols::moves::Movers_map const &,
//        core::pose::Pose const &
//        ) {
//
//  // TODO: implement this
//
//}
//
// /// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//TestQuickRelaxMoverCreator::create_mover() const {
// return protocols::moves::MoverOP( new TestQuickRelaxMover() );
//}
//
// /// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//TestQuickRelaxMoverCreator::keyname() const {
// return TestQuickRelaxMoverCreator::mover_name();
//}
//
// /// @brief Mover name for Rosetta Scripts
//std::string
//TestQuickRelaxMoverCreator::mover_name() {
// return "TestQuickRelaxMover";
//}
//
/// @brief Get the name of this Mover (TestQuickRelaxMover)
std::string
TestQuickRelaxMover::get_name() const {
	return "TestQuickRelaxMover";
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Run a quick relax on a membrane protein
/// @details Requires AddMembraneMover to be run beforehand
void TestQuickRelaxMover::apply( Pose & pose ) {

	// get time
	std::clock_t start = std::clock();

	using namespace protocols::idealize;
	using namespace protocols::membrane;
	using namespace protocols::simple_moves;
	using namespace protocols::moves;
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	TR << "Running MPQuickRelax protocol..." << std::endl;

	// add membrane to pose
	if ( addmem_ == true ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// get number of residues and number of moves
	Size nres( nres_protein( pose ) );
	if ( moves_ == "nres" ) {
		nmoves_ = nres;
	} else {
		nmoves_ = utility::string2Size( moves_ );
	}

	// find membrane position around the protein
	if ( mem_from_topo_ ) {
		MembranePositionFromTopologyMoverOP mempos( new MembranePositionFromTopologyMover() );
		mempos->anchor_at_res1( false );
		mempos->apply( pose );
	}

	// optimize membrane position using smooth scorefunction
	if ( opt_mem_ ) {
		OptimizeMembranePositionMoverOP optmem( new OptimizeMembranePositionMover() );
		optmem->apply( pose );
	}

	// creating constraint set
	if ( cst_file_.size() > 0 && cst_weight_ > 0 ) {

		ConstraintSetOP cst( ConstraintIO::get_instance()->read_constraints( cst_file_, ConstraintSetOP( new ConstraintSet() ), pose ) );
		pose.constraint_set( cst );

		// set constraints in scorefunction
		sfxn_->set_weight( atom_pair_constraint, cst_weight_ );
		sfxn_->set_weight( angle_constraint, cst_weight_ );
		sfxn_->set_weight( dihedral_constraint, cst_weight_ );
	}

	// starting position: Shake up the protein
	TR << "Choosing starting position: shake up the protein - Waka waka ..." << std::endl;

	// set small and shearmover and run it to mess up a pose
	Real kT = 1.0;
	SmallMoverOP small( new SmallMover( movemap_, kT, nmoves_ ) );
	small->angle_max( angle_max_ );
	small->apply( pose );

	ShearMoverOP shear( new ShearMover( movemap_, kT, nmoves_ ) );
	shear->angle_max( angle_max_ );
	shear->apply( pose );

	// create MC object
	MonteCarloOP mc( new MonteCarlo( pose, *sfxn_, 1.0 ) );

	// initialize AtomTreeMinimizer
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;

	//    fa_rep min_tol cst_wts
	// repeat 5
	// ramp_repack_min 0.02  0.01     1.0
	// ramp_repack_min 0.250 0.01     0.5
	// ramp_repack_min 0.550 0.01     0.0
	// ramp_repack_min 1     0.00001  0.0

	// sfxn_->set_weight( fa_rep, 0.5 );

	// run small and shearmover a second time
	// from then on it's only repacking and minimization
	TR << "SmallMover and ShearMover - Waka waka ..." << std::endl;
	small->apply( pose );
	shear->apply( pose );

	// do this for a certain number of iterations
	// since both the packer and especially minimization takes a while, just do one for now
	Size breakpoint( 0 );

	// score the pose
	Real tot_score = ( *sfxn_ )( pose );
	Real fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
	Real fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];

	// if fa_rep exceeds threshold (number of residues of the protein, as empirically determined)
	// or score is greater than 0, keep continuing
	while ( tot_score > 0 || fa_rep > 1.3 * nres_protein( pose ) ) {

		pose.energies().show_total_headers( TR );
		TR << std::endl;
		pose.energies().show_totals( TR );
		TR << std::endl;
		TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

		++breakpoint;

		if ( breakpoint == 3 ) {
			TR << "Counter == 3, score: " << tot_score << ", fa_rep: " << fa_rep << ", FAIL_DO_NOT_RETRY..." << std::endl;
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return;
		}

		////////////////////////////////////////////////////////////////////////////////

		// PACKING SIDECHAINS FROM RESIDUE CLOSEST TO TM COM OUTWARDS
		// get residue closest to TM COM
		core::SSize res_tm_com( rsd_closest_to_pose_tm_com( pose ) );

		// start 5 residues lower and repack for 10 residues
		core::SSize halfrange = 4;
		core::SSize range = 2 * halfrange;
		core::SSize nres = static_cast< core::SSize >( pose.total_residue() );

		// create packer task - will be re-used
		PackerTaskOP repack = TaskFactory::create_packer_task( pose );
		repack->restrict_to_repacking();

		// pack the center
		TR << "Packing center..." << std::endl;
		utility::vector1< bool > repack_center( get_repack_residues( pose, res_tm_com, res_tm_com, halfrange ) );
		repack->restrict_to_residues( repack_center );
		core::pack::pack_rotamers( pose, *sfxn_, repack );

		// go backwards and forwards by 10 residues and repack and minimize those
		// get the center points first for the lower and upper ranges, then get
		// the repack residues from that
		for ( core::SSize i = res_tm_com, j = res_tm_com;
				i >= 1 || j <= nres;
				i -= range, j += range ) {

			// pack both lower and upper ranges
			TR << "Packing..." << std::endl;
			utility::vector1< bool > repack_res( get_repack_residues( pose, i, j, halfrange ) );
			//   for ( core::Size k = 1; k <= repack_res.size(); ++k ) {
			//    TR << "k: " << k << " bool " << repack_res[ k ] << std::endl;
			//   }
			repack->restrict_to_residues( repack_res );
			core::pack::pack_rotamers( pose, *sfxn_, repack );

		}

		// score pose and print scores
		tot_score = ( *sfxn_ )( pose );
		fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
		fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];
		TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

		// evaluate Boltzmann
		mc->boltzmann( pose );
		pose = mc->lowest_score_pose();

		// packing
		if ( repack_again_ == true ) {

			TR << "Packing all rotamers one more time..." << std::endl;
			PackerTaskOP repack_all = TaskFactory::create_packer_task( pose );
			repack_all->restrict_to_repacking();
			core::pack::pack_rotamers( pose, *sfxn_, repack_all );

			// score pose and print scores
			tot_score = ( *sfxn_ )( pose );
			fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
			fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];
			TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

			// evaluate Boltzmann
			mc->boltzmann( pose );
			pose = mc->lowest_score_pose();

			// score pose and print scores
			tot_score = ( *sfxn_ )( pose );
			fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
			fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];
			TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;
		}

		//  // MINIMIZE OVER WINDOWS
		//  TR << "minimize center" << std::endl;
		//  MoveMapOP mm_cntr( new MoveMap() );
		//  mm_cntr->set_bb( false );
		//  mm_cntr->set_chi( false );
		//  mm_cntr->set_bb_true_range( res_tm_com - halfrange, res_tm_com + halfrange );
		//  mm_cntr->set_chi_true_range( res_tm_com - halfrange, res_tm_com + halfrange );
		//  atm.run( pose, *mm_cntr, *sfxn_, min_opts );
		//
		//  // minimize from center outwards
		//  for ( core::SSize i = res_tm_com, j = res_tm_com;
		//      i >= 1 || j <= nres;
		//      i -= range, j += range ) {
		//
		//   TR << "Minimizing range ..." << std::endl;
		//   MoveMapOP mm_range( new MoveMap() );
		//   mm_range->set_bb( false );
		//   mm_range->set_chi( false );
		//   mm_range->set_bb_true_range( i - halfrange, i + halfrange );
		//   mm_range->set_chi_true_range( i - halfrange, i + halfrange );
		//   mm_range->set_bb_true_range( j - halfrange, j + halfrange );
		//   mm_range->set_chi_true_range( j - halfrange, j + halfrange );
		//   atm.run( pose, *mm_range, *sfxn_, min_opts );
		//
		//  }

		//  // idealize decoy
		//  IdealizeMoverOP idealize( new IdealizeMover() );
		//  idealize->fast( true );
		//  idealize->apply( pose );

		//  // evaluate Boltzmann
		//  mc->boltzmann( pose );
		//  pose = mc->lowest_score_pose();

		// minimize all
		TR << "Minimizing final round..." << std::endl;
		atm.run( pose, *movemap_, *sfxn_, min_opts );

		////////////////////////////////////////////////////////////////////////////////

		// evaluate Boltzmann
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

		// set pose to lowest scoring pose from MC simulation
		pose = mc->lowest_score_pose();

		// score the pose
		tot_score = ( *sfxn_ )( pose );
		fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
		pose.energies().show_total_headers( TR );
		TR << std::endl;
		pose.energies().show_totals( TR );
		TR << std::endl;
		TR << "tot: " << tot_score << "fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;
		TR << "Iteration: " << breakpoint << " score: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

		//  pose.dump_scored_pdb("relax_test_" + utility::to_string( breakpoint ) + ".pdb", *sfxn_ );

	} // number of iterations for search

	// reset the weight in the scorefunction to what it was before (for mpsmooth)
	// sfxn_->set_weight( core::scoring::fa_rep, 0.44);
	tot_score = ( *sfxn_ )( pose );
	fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
	TR << "Final score: " << tot_score << " final fa_rep: " << fa_rep << std::endl;

	// superimpose poses with native
	SuperimposeMoverOP super( new SuperimposeMover( *native_, 1, nres, 1, nres, true ) );
	super->apply( pose );

	// get job for adding rmsd to scorefile
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// calculate and store the rmsd in the score file
	job->add_string_real_pair("rms", core::scoring::bb_rmsd( pose, *native_ ));

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

	// get time
	std::clock_t stop = std::clock();
	core::Real duration = ((double) stop - start )/CLOCKS_PER_SEC;

	// // setting mover status to fail-retry if something is going wrong
	// if ( duration > pose.total_residue() ) {
	//  TR << "Rosetta is thinking for too long, setting Mover status to FAIL_RETRY." << std::endl;
	//  TR << "Try running this app in release mode." << std::endl;
	//  set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	//  return;
	// }

	TR << "Run took " << duration << " seconds" << std::endl;

	// dump pose
	// pose.dump_scored_pdb("test_relax.pdb", *sfxn_ );

}// apply

////////////////////////////////////////////////////////////////////////////////

/// @brief Run AddMembraneMover again?
/// @details If you want to keep your anchor point for MEM, then pick no
void TestQuickRelaxMover::add_membrane_again( bool yesno ) {
	addmem_ = yesno;
}

/// @brief Run MembranePositionFromTopology again?
/// @details Will change the starting membrane position
void TestQuickRelaxMover::membrane_from_topology( bool yesno ) {
	mem_from_topo_ = yesno;
}

/// @brief Optimize membrane position before relax?
void TestQuickRelaxMover::optimize_membrane( bool yesno ) {
	opt_mem_ = yesno;
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void TestQuickRelaxMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::mp::quickrelax::angle_max );
	option.add_relevant( OptionKeys::mp::quickrelax::nmoves );
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::mp::quickrelax::repack_again );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Set default values
void TestQuickRelaxMover::set_defaults() {

	// maximum change in dihedral angles
	angle_max_ = 1.0;

	// number of moves
	moves_ = "nres";

	// Movemap
	movemap_->set_bb( true );
	movemap_->set_chi( true );

	// create scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

	// default constraint weight of 1
	cst_weight_ = 1.0;

	// run AddMembraneMover again?
	addmem_ = true;

	// starting MembranePositionFromTopology?
	mem_from_topo_ = false;

	// optimize membrane before relax?
	opt_mem_ = false;

	// additional round of repack?
	repack_again_ = false;

}// set_defaults

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void TestQuickRelaxMover::init_from_cmd() {

	using namespace basic::options;

	// read native and attach membrane to it
	if ( option[ OptionKeys::in::file::native ].user() ) {
		native_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( *native_ );
	}

	if ( option[ OptionKeys::mp::quickrelax::angle_max ].user() ) {
		angle_max_ = option[ OptionKeys::mp::quickrelax::angle_max ]();
	}

	if ( option[ OptionKeys::mp::quickrelax::nmoves ].user() ) {
		if ( option[ OptionKeys::mp::quickrelax::nmoves ]() == "nres" ) {
			moves_ = utility::to_string( nres_protein( *native_ ) );
		} else {
			moves_ = option[ OptionKeys::mp::quickrelax::nmoves ]();
		}
	}

	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		cst_file_ = option[ OptionKeys::constraints::cst_file ]()[1];
	}

	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
	}

	if ( option[ OptionKeys::mp::quickrelax::repack_again ].user() ) {
		repack_again_ = option[ OptionKeys::mp::quickrelax::repack_again ]();
	}

}// init from cmdline

//////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
utility::vector1< bool > TestQuickRelaxMover::get_repack_residues( Pose & pose, core::SSize center1, core::SSize center2, core::SSize halfrange ){

	core::SSize nres = static_cast< core::SSize >( pose.total_residue() );
	core::SSize m, n;

	// initialize vector with false
	utility::vector1< bool > repack_res( nres, false );

	TR << "center1 " << center1 << " center2 " << center2 << std::endl;

	// go through residues in the lower range
	for ( core::SSize i = center1, j = center1; i >= center1-halfrange || j <= center1+halfrange; --i, ++j ) {

		// reset counters when out of the pose
		i <= 1 || i >= nres-1 ? m = 1 : m = i;
		j >= nres-1 || j <= 1 ? n = nres-1 : n = j;

		//  TR << "i " << i << " j " << j << std::endl;
		//  TR << "m " << m << " n " << n << std::endl;

		repack_res[ m ] = true;
		repack_res[ n ] = true;
	}

	// go through residues in the upper range
	for ( core::SSize i = center2, j = center2; i >= center2-halfrange || j <= center2+halfrange; --i, ++j ) {
		// reset counters when out of the pose
		i <= 1 || i >= nres-1 ? m = 1 : m = i;
		j >= nres-1 || j <= 1 ? n = nres-1 : n = j;

		//  TR << "m " << m << " n " << n << std::endl;

		repack_res[ m ] = true;
		repack_res[ n ] = true;
	}

	return repack_res;

}// get pack residues

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::jd2;
		using namespace protocols::relax::membrane;


		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		MPRangeRelaxMoverOP mprr( new MPRangeRelaxMover() );
		JobDistributor::get_instance()->go(mprr);

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
