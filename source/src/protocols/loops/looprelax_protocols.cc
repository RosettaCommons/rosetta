// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to LoopRebuild
/// @detailed
/// @author Mike Tyka
/// @author James Thompson
/// @author Srivatsan Raman

#include <protocols/loops/looprelax_protocols.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/file/FileName.hh>
#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <numeric/model_quality/rms.hh>

// External library headers
#include <utility/exit.hh>
#include <numeric/random/random.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>


#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>


using basic::T;
using basic::Error;
using basic::Warning;

//using namespace ObjexxFCL;

namespace protocols {

static thread_local basic::Tracer TR( "protocols.looprelax_protocols" );

using namespace core;
using io::pdb::dump_pdb;




using namespace protocols::loops;

	LoopRebuild::LoopRebuild(
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::loops::Loops Loops_in
	) : Mover(),
		scorefxn_( scorefxn ),
		Loops_in_( Loops_in )
	{
		protocols::loops::read_loop_fragments( frag_libs_ );
		Mover::type("LoopRebuild");
		set_default_settings();
	}

LoopRebuild::~LoopRebuild() {}

/// @brief Clone this object
protocols::moves::MoverOP
LoopRebuild::clone() const {
	return protocols::moves::MoverOP( new LoopRebuild(*this) );
}


/// @details  Apply the loop rebuilding protocol to a pose.
///    The loops to be rebuilt and fragments are defined upon object construction.
void LoopRebuild::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//core::Real score_filter_cutoff;
	bool passed_score_filter( false );
	int const n_score_filter_fail_tol( 2 );
	int n_score_filter_fail( 0 );
	int n_loop_fail_tol( 2 );
	int loop_build_round( 1 ); //build full loops this many times
	success = false;

	//setting mc object
	if( !mc_created )
		set_default_mc( pose );

	//initialize fragments
	//		std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	//		initialize_fragments( frag_libs );

	//		found_loops = read_loop_file( ori_loops_begin, ori_loops_end );
	TR.Debug << "Loop file size : " << Loops_in_.size() << std::endl;
	//read_coord_cst(); //include this function later !

	Loops_in_.verify_against( pose );
	protocols::loops::set_secstruct_from_psipred_ss2( pose );
	//score_filter_cutoff = get_score_filter_cutoff();  // set but never used ~Labonte

	TR.Debug << "Loop file size (after filtering) : " << Loops_in_.size() << std::endl;
	//save the initial pose  in case failed scorefilter
	pose::Pose init_pose_obj;
	init_pose_obj = pose;

	bool loop_done( false );
	while( !passed_score_filter && n_score_filter_fail < n_score_filter_fail_tol ){
		n_score_filter_fail++;

		//reset pose
		pose = init_pose_obj;


		int loop_nfail( 0 );
		loop_done = false;
		while( !loop_done && loop_nfail <= n_loop_fail_tol ) {
			loop_nfail++;
			for( int i = 0; i < loop_build_round; ++i ) {
				loop_done = build_random_loops( pose );
			}
		}

		if (!loop_done) TR.Error << "build_random_loops() returned failure." << std::endl;

		passed_score_filter = true; //change this to a call to score_filter function, once it is written
	}

	success = (loop_done && passed_score_filter);
	if (!success) TR.Error << "LoopRebuild::apply() returned false" << std::endl;
}

std::string
LoopRebuild::get_name() const {
	return "LoopRebuild";
}


///////////////////////////////////////////////////////////////////////
core::Real LoopRebuild::get_rmsd_tolerance() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static core::Real rmsd_tol = { 10000.0 };
	static bool init = { false };

	if ( !init ) {
		rmsd_tol = option[ OptionKeys::loops::rmsd_tol ]();
		if( rmsd_tol <= 0 ) rmsd_tol = 10000.0;
		init = true;
	}
	return rmsd_tol;
}


////////////////////////////////////////////////////////////////////////
core::Real LoopRebuild::get_chain_break_tolerance() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static core::Real chain_break_tol = { 0.2 };
	static bool init = { false };

	if ( !init ) {

		chain_break_tol = option[ OptionKeys::loops::chain_break_tol ]();
		init = true;
	}
	return chain_break_tol;
}


///////////////////////////////////////////////////////////////////////
/// @details  Rebuild all the loops in the pose, one at a time, choosing each in random order.
bool LoopRebuild::build_random_loops(	core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	int const nres( pose.total_residue() );

	std::vector< int > free_res; // stores residue numbers in real loops
  for( Loops::const_iterator it=Loops_in_.begin(), it_end=Loops_in_.end(); it != it_end; ++it ) {
		TR.Debug << "Loop res " <<  it->start() << " " <<  it->stop() << std::endl;
		for ( int k = (int)it->start(); k <= (int)it->stop(); ++k )
			free_res.push_back(k);
	}

	std::vector<int> inter_res; // residues in-between loop_file defined loops

	// now go through all the loop regions
	core::pose::Pose stage_pose; // stores a pose at certain stage
	std::vector< int > folded_loops; // loops folded
	int num_desired_loops ( get_desired_loops_exist()? std::min( desired_loops(), Loops_in_.size() ) :  Loops_in_.size()  );
	int cutpoint(0);

	if ( loop_model() ) num_desired_loops = Loops_in_.size(); // build all obligate loops


	for ( int loop_counter = 0; loop_counter < num_desired_loops; ++loop_counter ){

		// save the original pose at this stage
		stage_pose = pose;

		// randomly select one loop from the loop regions
		inter_res.clear();
		int selected_loop, def_loop_begin, def_loop_end, combine_number;
		bool loops_combined( false );
		bool extend_this_loop = option[ OptionKeys::loops::extended ].user();
		bool use_selected_loop =
			select_one_loop( nres, selected_loop, folded_loops,
											 inter_res, def_loop_begin, def_loop_end, cutpoint, extend_this_loop,
											 loops_combined, combine_number, loop_counter );

		if ( !use_selected_loop ) continue; // skip this loop

		// params record the behaviors of individual loop modeling process
		int   nfail( 0 );
		int   n_chain_break_fail( 1 );
		bool  is_chain_break( true );
		bool  rmsd_acceptable( false );
		int   barcst_extend_begin( 0 );
		int   barcst_extend_end( 0 );
		int total_combine( 0 );
		if ( loops_combined ) total_combine = combine_number;
		int backward_combine( 0 );

		// select further random stems which will be constrainted by barcode cst
		if ( get_random_loop_flag() ){
			barcst_extend_begin = static_cast <int> ( numeric::random::uniform() * 2 );
			barcst_extend_end   = static_cast <int> ( numeric::random::uniform() * 2 );
		}

		// store starting fold tree
	 	kinematics::FoldTree f_orig=pose.fold_tree();

		int final_loop_begin=0;
		int final_loop_end=0;
		int time_start = time(NULL);
		float time_per_build = 0.0;
		int nclosurefail = 0;
		int nrmsfail = 0;

		while( is_chain_break || !rmsd_acceptable ) {
			// epic fail
			if( nfail++ > get_allowed_failure_before_stop() ) {
				// restore simple fold tree
				TR.Error << "Failed to rebuild loop " << nfail << " times!  Aborting." << std::endl;
				remove_cutpoint_variants( pose );
				pose.fold_tree( f_orig );

				if( get_abort_on_failed_loop() ) return false;
				else                             	break;
			}

			// refresh pose
			pose = stage_pose;

			// combine next loop if this loop can't be closed after 20 iterations
			if ( get_random_loop_flag() && get_combine_if_fail_exist() ){
				if ( n_chain_break_fail % get_allowed_failure_before_extend() == 0 ){
					if ( core::Size(selected_loop + total_combine + 1) < Loops_in_.size() ){
						def_loop_end = Loops_in_[ selected_loop + total_combine + 1 ].stop();
						extend_this_loop |= Loops_in_[ selected_loop + total_combine + 1 ].is_extended();
						if( !loop_model() ) folded_loops.push_back( selected_loop + total_combine + 1);
						total_combine++; // one more combined loop
					} else if ( selected_loop - backward_combine > 1 ){
						def_loop_begin = Loops_in_[ selected_loop - backward_combine - 1 ].start();
						extend_this_loop |= Loops_in_[ selected_loop - backward_combine - 1 ].is_extended();
						if( !loop_model() ) folded_loops.push_back( selected_loop - backward_combine - 1 );
						backward_combine++;
					}
					else
						break;// return false; // give up if no loop to combine

					if( !loop_model() )
						loop_counter++;

					barcst_extend_begin = static_cast <int> ( numeric::random::uniform() * 2 );
					barcst_extend_end = static_cast <int> ( numeric::random::uniform() * 2 );
				}
			}

			// further extend loop regions, and add barcode constraints on these
			// extended regions
			int loop_begin = def_loop_begin - 1;
			int loop_end = def_loop_end;
			if( loop_begin < 1 ) loop_begin = 1;
			if( loop_end > nres ) loop_end = nres;

			if( get_random_loop_flag() ){
				barcode_extend_stems( pose, barcst_extend_begin, barcst_extend_end,
			 												loop_begin, loop_end, def_loop_begin, def_loop_end,
			 												nres, selected_loop, total_combine,
			 												backward_combine );
			}

			if( option[ OptionKeys::loops::extended_beta ].user()  ){
				Real expfactor =  exp( - option[ OptionKeys::loops::extended_beta ]() * fabs( (double) loop_begin - loop_end ) );
				bool stochastic_extend_this_loop = ( numeric::random::uniform() < expfactor );
				extend_this_loop = extend_this_loop || stochastic_extend_this_loop;
			}

			build_loop_with_ccd_closure( pose, loop_begin, loop_end, cutpoint, extend_this_loop );

			// extend more barcode regions
			extend_barcode_regions_if_chain_break( pose, loop_begin, loop_end,
			                                       n_chain_break_fail, is_chain_break,
			                                       barcst_extend_begin, barcst_extend_end );
			rmsd_acceptable = acceptable_rmsd_change( stage_pose, pose  );

			if( ! rmsd_acceptable ){ nclosurefail++; TR.Info << "WARNING: Loop Built, rms not acceptable - trying again" << std::endl; }
			if(   is_chain_break  ){ nrmsfail++;     TR.Info << "WARNING: Chain_break_remains - tryign again" << std::endl; }

			final_loop_begin = loop_begin;
			final_loop_end   = loop_end;
		}//accepted this loop

		int time_end = time(NULL);
		time_per_build = float(time_end - time_start) / float(nfail);

		using namespace ObjexxFCL::format;
		TR  << "Loopstat: "
		    << "  " << I(3,def_loop_begin)
		    << "  " << I(3,def_loop_end    )
		    << "  " << I(3,final_loop_begin )
		    << "  " << I(3,final_loop_end    )
				<< "  " << I(3,final_loop_end - final_loop_begin )
				<< "  time " << F(5,1,time_per_build )
				<< "  " << I(3,nfail)
				<< "  " << I(3,nclosurefail )
				<< "  " << I(3,nrmsfail)
				<< "  " << (extend_this_loop  ? std::string(" ext ") : std::string(" noext " ))  << std::endl;


		remove_cutpoint_variants( pose ); // remove cutpoint variants
		pose.fold_tree( f_orig );  // restore simple fold tree

	} // all loops folded

	return true;
}


////////////////////////////////////////////////////////////////////////////////
void LoopRebuild::set_looprlx_allow_move_map(
																					 int const & loop_begin,
																					 int const & loop_end,
																					 core::kinematics::MoveMap & mm
																					 ) {
	using namespace core::id;

	mm.set_bb( false );
	mm.set_jump( false );

	for ( int i = loop_begin; i <= loop_end; ++i ) {
		mm.set_bb( i, true );
		mm.set( TorsionID(i,BB,3), false ); //omega is fixed
		mm.set_chi( i, true );
	}
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  Rebuild a loop via fragment insertion + ccd closure + minimization
void LoopRebuild::build_loop_with_ccd_closure(
	core::pose::Pose & pose,
	int const & loop_begin,
	int const & loop_end,
	int & cutpoint,
	bool const & extend_this_loop
) {
	using namespace kinematics;
	using namespace scoring;
	using namespace optimization;
	using namespace protocols::simple_moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace numeric::random;

	int const nres =  pose.total_residue();


	kinematics::MoveMap mm_one_loop;
	set_looprlx_allow_move_map( loop_begin, loop_end, mm_one_loop );
	if ( cutpoint== 0 ){
		Loop myloop( loop_begin+1, loop_end, 0 );
		myloop.choose_cutpoint( pose );
		cutpoint = myloop.cut();
	}
	int cut_orig = cutpoint;
	TR.Info << "Loop and cutpoint: " << loop_begin
													 << "  " << cutpoint
													 << "  " << loop_end
													 << (extend_this_loop ? " extended " : " nonextended" ) << std::endl;

	/////
	// ensure (cutpoint+1) is not a proline (this is a horrible hack -fpd)
	int ntries = 0;
	while (cutpoint != nres &&
	       pose.residue(cutpoint+1).aa() == core::chemical::aa_pro &&
	       ntries < 100) {
		Loop myloop( loop_begin, loop_end );
		myloop.choose_cutpoint( pose );
		cutpoint = myloop.cut();
		ntries++;
	}
	// if not sussessfully chosen cutpoint, just chose whatever is not a
	// profile in the loop !
	ntries = 0;
	while (cutpoint != nres &&
	       pose.residue(cutpoint+1).aa() == core::chemical::aa_pro &&
	       ntries < 100) {
		cutpoint = loop_begin + int( float(loop_end - loop_begin) * uniform() );
		ntries++;
	}
	if (cutpoint != nres && pose.residue(cutpoint+1).aa() == core::chemical::aa_pro ) {
		TR << "  Unable to rebuild loop [" << loop_begin << ", " << loop_end << "] " << std::endl;
		return;
	}
	if (cutpoint != cut_orig)
		TR << "  Changing cutpoint to " << cutpoint << "" << std::endl;


	int const loop_size( loop_end - loop_begin + 1 );
	core::Real const cycle_ratio( std::max( 0.2 , get_looprlx_cycle_ratio() ) );  // 0.5->0.2, fpd
	int cycles2;
	int cycles3;

	if( get_ccd_closure_exist() ) {
		cycles2 = loop_model() ? 3:2 ;
		int base_cycles( std::max( 15, static_cast<int>( 5*loop_size*cycle_ratio )));
		cycles3 = loop_model() ? 2*base_cycles:base_cycles;
	} else {
		cycles2 =  10;
		cycles3 =  std::max( 30, static_cast<int>( 10*loop_size*cycle_ratio));
	}

	TR << "Number of cycles: cycles2 and cycles3 " << cycles2 << " " << cycles3 << std::endl;


	/// ----- Set up repacker -------------------------------------------------------------
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( pose ));
	pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
	utility::vector1< bool > allow_repack( nres, false );
	for ( int i = loop_begin; i <= loop_end; ++i )
		allow_repack[ i ] = true;
	this_packer_task->restrict_to_residues( allow_repack );

	set_single_loop_fold_tree( pose, Loop(loop_begin, loop_end, cutpoint ) );

	bool chainbreak_present =  ( loop_begin != 1 && loop_end != nres );

	// special case ... vrt res at last position
	chainbreak_present &= (loop_end != nres-1 || pose.residue( nres ).aa() != core::chemical::aa_vrt );

	// set cutpoint variant residue for chainbreak score if chanbreak is present
	if( chainbreak_present ){
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint+1 );
	}

	(*scorefxn_)(pose);
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies( pose );
	core::pose::Pose start_pose = pose;


	// either extend or at least idealize the loop (just in case).
	if( extended_loop() || extend_this_loop)   set_extended_torsions( pose, Loop( loop_begin, loop_end ) );
	else                                       idealize_loop(  pose, Loop( loop_begin, loop_end ) );


	// prepare fragment movers
	MoveMapOP movemap( new MoveMap() );
	movemap->set_bb_true_range(loop_begin, loop_end);

	std::vector< FragmentMoverOP > fragmover;
	for ( std::vector< core::fragment::FragSetOP >::const_iterator
				it = frag_libs_.begin(), it_end = frag_libs_.end();
				it != it_end; it++ ) {
		ClassicFragmentMoverOP cfm( new ClassicFragmentMover( *it, movemap ) );
		cfm->set_check_ss( false );
		cfm->enable_end_bias_check( false );
		fragmover.push_back( cfm );
	}

	if( !option[OptionKeys::loops::no_randomize_loop] ){
		// insert random fragment as many times as the loop is long (not quite the exact same as the old code)
		for ( int i = loop_begin; i <= loop_end; ++i ) {
			for ( std::vector< FragmentMoverOP >::const_iterator
						it = fragmover.begin(),it_end = fragmover.end(); it != it_end; it++ ) {
				(*it)->apply( pose );
			}
		}
	}

	// setup the move objects
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( movemap, 1.0, 1 ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );

	// setup the move objects
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover( movemap, 1.0, 1 ) );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 2.0 );
	shear_mover->angle_max( 'L', 3.0 );

	// create a Random Mover, fill it with individual moves
	moves::RandomMoverOP refine_mover( new moves::RandomMover() );
	refine_mover->add_mover( small_mover );
	refine_mover->add_mover( shear_mover );

	// Set up Minimizer object
	AtomTreeMinimizer mzr;
	core::Real const dummy_tol( 0.001 );
	bool const use_nblist( true ), deriv_check( false );
	MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check );


	// Set up MonteCarlo Object
	core::Real const init_temp = 2.0;
	core::Real temperature = init_temp;
 	mc_->reset( pose );
	mc_->set_temperature( temperature );


	// --- Figure out constraints (ConstraintSet)	  ------
	core::scoring::constraints::ConstraintSetOP pose_cst( new core::scoring::constraints::ConstraintSet(	*pose.constraint_set() ) );


	// --- Figure out constraints (ConstraintSet)	  ------
	if(  option[OptionKeys::loops::coord_cst ]()  > 0.0 ){
		core::scoring::constraints::ConstraintSetOP new_csts( new core::scoring::constraints::ConstraintSet );
		for( int ir=loop_begin; ir <= loop_end; ir++ ){
			Real middlefactor = 1.0 - (float(std::min( std::abs(ir - loop_begin), std::abs(loop_end - ir) ))/ float(std::abs(loop_end-loop_begin))) * 2 * 0.8;

			core::scoring::func::FuncOP CA_cst( new core::scoring::func::HarmonicFunc( 0, option[OptionKeys::loops::coord_cst ]() * middlefactor ) );
			core::scoring::constraints::ConstraintOP newcst( new core::scoring::constraints::CoordinateConstraint(
								core::id::AtomID( pose.residue_type(ir).atom_index("CA") , ir ),
								core::id::AtomID( pose.residue_type( 1).atom_index("CA") , 1  ),
								start_pose.residue(ir).xyz( "CA" ),
								CA_cst
							) );
			new_csts->add_constraint( newcst );
		}
		pose.constraint_set( new_csts );
	}

	float final_constraint_weight          = option[ basic::options::OptionKeys::constraints::cst_weight ]();

	//bool 	has_residue_pair_constraints     = pose.constraint_set()->has_residue_pair_constraints();
	//bool	has_intra_residue_constraints    = pose.constraint_set()->has_intra_residue_constraints();
	//bool 	has_non_residue_pair_constraints = pose.constraint_set()->has_non_residue_pair_constraints();
	bool 	has_constraints                  = pose.constraint_set()->has_constraints();

	core::scoring::constraints::ConstraintSetCOP full_cst =  pose.constraint_set(); // pose.constraint_set()

	if( has_constraints ){
		// do we want to disable all constraints that are not part of the loop for speed
		// reasons ? (Constraints are really quite slow!)

		// ConstraintSet
		//filter out constraints that have at least one "anchor" in the loop
		scorefxn_->show(  TR , pose );
		core::scoring::constraints::ConstraintSetOP loop_cst( new core::scoring::constraints::ConstraintSet( *full_cst, loop_begin, loop_end ) );
		pose.constraint_set(  loop_cst );
		scorefxn_->show(  TR , pose );
	}


	// OK! Let's go !
	int   starttime    = time(NULL);
	int   frag_count   = 0;
	scorefxn_->show_line_headers( TR );

	if ( get_ccd_closure_exist() ) {
		TR << "***** CCD CLOSURE *****" << std::endl;
		core::Real const final_temp( 1.0 );
		core::Real const gamma = std::pow( (final_temp/init_temp), (1.0/(cycles2*cycles3)) );

		for( int c2 = 1; c2 <= cycles2; ++c2 ) {
			mc_->recover_low( pose );
			(*scorefxn_)(pose);

			// ramp up constraints
			if( has_constraints ) {
				if( c2 != cycles2 ) {
					scorefxn_->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight*float(c2)/float( cycles2 ) );
				} else {
					scorefxn_->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight * 0.2 );
				}
			}

			scorefxn_->show_line( TR , pose );
			TR << std::endl;
			mc_->score_function( *scorefxn_ );

			for( int c3 = 1; c3 <= cycles3; ++c3 ){
				temperature *= gamma;
				mc_->set_temperature( temperature );
				for ( std::vector< FragmentMoverOP >::const_iterator
							it = fragmover.begin(),it_end = fragmover.end(); it != it_end; it++ ) {
					(*it)->apply( pose );
					if( chainbreak_present ) fast_ccd_close_loops( pose, loop_begin, loop_end, cutpoint, mm_one_loop );
					mzr.run( pose, mm_one_loop, *scorefxn_, options );
					mc_->boltzmann( pose, "ccd_closure" );
					frag_count++;
				}
			}
		}

		// now ensure the loop is properly closed!
		if( false && chainbreak_present ){
			mc_->recover_low( pose );
			TR << "--" << std::endl;
			(*scorefxn_)(pose);
			scorefxn_->show_line( TR , pose );
			TR << std::endl;
			fast_ccd_close_loops( pose, std::max(cutpoint-3,loop_begin),  std::min(cutpoint+3,loop_end), cutpoint, mm_one_loop );
			(*scorefxn_)(pose);
			scorefxn_->show_line( TR , pose );
			TR << std::endl;
			mc_->reset( pose );
		}

	} else {
		TR << "***** DOING CCD MOVES *****" << std::endl;
		float const final_chain_break_weight = 5.0;
		float const delta_weight( final_chain_break_weight/cycles2 );
		scorefxn_->set_weight( chainbreak, 0.0 ); //not evaluating quadratic chainbreak
		mc_->score_function( *scorefxn_ );

		//core::Real const final_temp( 1.0 );
		//core::Real const gamma = std::pow( (final_temp/init_temp), (1.0/(cycles2*cycles3)) );

		for ( int c2 = 1; c2 <= cycles2; ++c2 ) {
			mc_->recover_low( pose );

			// ramp up constraints
			if( has_constraints ) {
				if( c2 != cycles2 ) {
					scorefxn_->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight*float(c2)/float( cycles2 ) );
				} else {
					scorefxn_->set_weight( core::scoring::atom_pair_constraint, final_constraint_weight * 0.2 );
				}
			}
			if ( chainbreak_present ) {
				scorefxn_->set_weight( linear_chainbreak, c2*delta_weight );
			}
			mc_->score_function( *scorefxn_ );

			// score and print an info line
			(*scorefxn_)(pose);
			scorefxn_->show_line( TR , pose );
			TR << std::endl;

			for ( int c3 = 1; c3 <= cycles3; ++c3 ) {
				//temperature *= gamma;
				//mc_->set_temperature( temperature );
				if(( !chainbreak_present || uniform()*cycles2 > c2 ))
				{
					//do fragment moves here

					if( !option[OptionKeys::loops::refine_only ]() ){
						for ( std::vector< FragmentMoverOP >::const_iterator
									it = fragmover.begin(),it_end = fragmover.end(); it != it_end; it++ ) {


							if( ((*it)->fragments()->max_frag_length() == 1 ) && (uniform() < option[OptionKeys::loops::skip_1mers ]() ) ) continue;
							if( ((*it)->fragments()->max_frag_length() == 3 ) && (uniform() < option[OptionKeys::loops::skip_3mers ]() ) ) continue;
							if( ((*it)->fragments()->max_frag_length() == 9 ) && (uniform() < option[OptionKeys::loops::skip_9mers ]() ) ) continue;

							(*it)->apply( pose );
							frag_count++;
						}
					}else{
						refine_mover->apply( pose );
					}

				} else {
					//do ccd_moves here
					if( ! option[OptionKeys::loops::skip_ccd_moves ]() ){
						loop_closure::ccd::CCDLoopClosureMover ccd_mover(
								Loop( loop_begin, loop_end, cutpoint ),
								MoveMapCOP( new MoveMap( mm_one_loop ) ) );
						ccd_mover.max_cycles( 25 );  // Used to be 5 moves, which would result in 25 "tries" in the old code. ~Labonte
						ccd_mover.apply( pose );
					}
				}
				mc_->boltzmann( pose, "ccd_moves" );
			} // cycles3
			mzr.run( pose, mm_one_loop, *scorefxn_, options );
		} //cycles2

		// now ensure the loop is properly closed!
		if( false && chainbreak_present ){
			mc_->recover_low( pose );
			TR << "--" << std::endl;
			(*scorefxn_)(pose);
			scorefxn_->show_line( TR , pose );
			TR << std::endl;
			fast_ccd_close_loops( pose, std::max(cutpoint-3,loop_begin),  std::min(cutpoint+3,loop_end), cutpoint, mm_one_loop );
			(*scorefxn_)(pose);
			scorefxn_->show_line( TR , pose );
			TR << std::endl;
			mc_->reset( pose );
		}
	} // if-else get_ccd_closure_exist

	int looptime = time(NULL) - starttime;
	TR << "FragCount: " << frag_count << std::endl;
	TR << "Looptime " << looptime << std::endl;

	//v	pose = mc.lowest_score_pose();
	pose = mc_->lowest_score_pose();
	scorefxn_->show(  TR , pose );
	TR << "-------------------------" << std::endl;
	mc_->show_counters();

	pose.constraint_set( pose_cst );
}


//////////////////////////////////////////////////////////////////////////////////
/// @details  CCD close the loop [loop_begin,loop_end].
/// Wraps protocols::loops::loop_closure::ccd::CCDLoopClosureMover.apply() using most of its default options.
/// rama scores are not checked, however, and the secondary structure is "fixed" afterward.
/// @remark   This is a misnomer; it actually closes a single loop only. ~Labonte
void LoopRebuild::fast_ccd_close_loops(
		core::pose::Pose & pose,
		int const & loop_begin,
		int const & loop_end,
		int const & cutpoint,
		kinematics::MoveMap & mm )
{
	loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover(
			Loop( loop_begin, loop_end, cutpoint ), kinematics::MoveMapCOP( new kinematics::MoveMap( mm ) ) );
	ccd_loop_closure_mover.check_rama_scores( false );
	ccd_loop_closure_mover.apply( pose );

	// fix secondary structure??
	for (int i=loop_begin; i<=loop_end; ++i) {
		char ss_i = pose.conformation().secstruct( i );
		if ( ss_i != 'L' && ss_i != 'H' && ss_i != 'E')
			pose.set_secstruct( i , 'L' );
	}
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::select_one_loop(
																int nres,
																int & selected_loop,
																std::vector< int > & folded_loops,
																std::vector< int > & inter_res,
																int & loop_begin,
																int & loop_end,
																int & cutpoint,
																bool & extend_this_loop,
																bool & are_loops_combined,
																int & combine_interval,
																int & loop_counter
																) {
	int const num_loops = (int)(Loops_in_.size() );

	// make how many loop combination(s) per structure on average
	core::Real const loop_combine_rate( get_loop_combine_rate() );
	core::Real const num_of_loops_to_combine( numeric::random::uniform()*(num_loops - 1)*loop_combine_rate );

	TR << "LoopsToCombine: " << num_of_loops_to_combine << "  " << num_loops << std::endl;

	selected_loop = 0;
	do{
		// do loops in order if random_loop is not set !
		if ( !get_random_loop_flag() ){
			selected_loop += 1;
		}else{
			selected_loop = (int)( numeric::random::uniform()*num_loops + 1);
		}
	}while( std::find( folded_loops.begin(), folded_loops.end(), selected_loop ) !=
				 folded_loops.end() );
	folded_loops.push_back(selected_loop);

	loop_begin  = Loops_in_[selected_loop].start();
	loop_end = Loops_in_[selected_loop].stop();
	cutpoint = Loops_in_[selected_loop].cut();
	extend_this_loop |= Loops_in_[selected_loop].is_extended();

	if (loop_begin <= 0) loop_begin = 1;
	if (loop_end >= nres) loop_end = nres;

	if ( !get_random_loop_flag() ) return true;

	// LOOPSETS APPLIED OUTSIDE THE PROTOCOL
	//if ( numeric::random::uniform() < get_loop_skip_rate() )
	//			return false;

	combine_interval = 1;
	// combine consecutive loop regions
	if ( selected_loop < num_loops - combine_interval && num_of_loops_to_combine > 0 ) {
		int const longer_loop_size ( Loops_in_[ selected_loop + combine_interval ].stop() -
																 Loops_in_[ selected_loop ].start());
		int loop_limit;
		int const terminal_loop(12);
		int const internal_loop(25);
		if( (int)Loops_in_[ selected_loop + combine_interval ].stop() >= nres ||
				(int)Loops_in_[ selected_loop ].start() <= 1 )
			loop_limit = terminal_loop;
		else  loop_limit = internal_loop;
		// don't merge internal loops longer than 25,
		// or terminal loops longer than 12



		if ( numeric::random::uniform() < num_of_loops_to_combine/(num_loops-1+1e-10 ) ) {

			loop_begin = Loops_in_[selected_loop].start();
			loop_end   = Loops_in_[selected_loop].stop() + combine_interval;
			TR.Info << "Combining loops: " << loop_begin << "  " <<  loop_end <<  std::endl;
			// if loop exceeds maximum length - truncate it randomly on either side
			if( longer_loop_size > loop_limit ){
				if(  numeric::random::uniform() < 0.5 ){
					loop_begin = loop_end - (loop_limit - 1);
				}else{
					loop_end = loop_begin + (loop_limit - 1);
				}
			}

			for ( int ll = 0; ll < combine_interval; ++ll ){
				for ( int k  = (int)Loops_in_[ selected_loop +ll ].stop();
									k <= (int)Loops_in_[ selected_loop +ll + 1 ].start(); ++k )
					inter_res.push_back(k);
				if( !loop_model() ) {
					folded_loops.push_back(selected_loop+ll+1);
					loop_counter++;
				}
			}
			are_loops_combined = true;
		}
	}


	// choose random taking off points
	//int const old_loop_begin( loop_begin );
	//int const old_loop_end( loop_end );
	int rand_extension = 0;  //fpd 1->0
	int rand_limit = 4;
	if ( loop_model() || extended_loop() || extend_this_loop ){
		rand_extension = 0;
		rand_limit = 3;
	}

	do {
		if ( loop_begin > 1 )
			loop_begin = loop_begin - static_cast <int> ( numeric::random::uniform() * rand_limit )
				+ rand_extension;
		if ( loop_end != nres )
			loop_end = loop_end + static_cast <int> ( numeric::random::uniform() * rand_limit )
				- rand_extension;
	//	TR << "Fraying" << loop_begin << "  " << loop_end << std::endl;
	} while (  numeric::random::uniform() < 0.3 );

	if ( loop_end - loop_begin < 0 ){
		loop_end++;
		loop_begin--;
	}
	loop_begin = std::max( 1, loop_begin );
	loop_end   = std::min( nres, loop_end );

	return true;
}


//////////////////////////////////////////////////////////////////////////////////
core::Real LoopRebuild::get_score_filter_cutoff() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static core::Real score_filter_cutoff = { 1.0 };
	static bool init = { false };

	if ( !init ) {
		score_filter_cutoff = option[OptionKeys::loops::score_filter_cutoff ]();
		init = true;
	}
	return score_filter_cutoff;
}

//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::loop_model()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool loop_model = { false };
	static bool init = { false };

	if ( !init ) {
		loop_model = option[ OptionKeys::loops::loop_model ];
		init = true;
	}
	return loop_model;
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::get_ccd_closure_exist() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool ccd_closure_exist = { false };
	static bool init = { false };

	if ( !init ) {
		ccd_closure_exist = option[ OptionKeys::loops::ccd_closure];
		init = true;
	}
	return ccd_closure_exist;
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::get_desired_loops_exist() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool desired_loops_exist = { false };
	static bool init = { false };

	if ( !init ) {
		desired_loops_exist = option[ OptionKeys::loops::loops_subset];
		init = true;
	}
	return desired_loops_exist;
}


//////////////////////////////////////////////////////////////////////////////////
core::Size LoopRebuild::desired_loops() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool init = false;
	static int desired_loop = 1;

	if ( !init ) {
		desired_loop = option[ OptionKeys::loops::num_desired_loops]();
		init = true;
	}
	return desired_loop;
}


//////////////////////////////////////////////////////////////////////////////////
core::Real LoopRebuild::get_loop_combine_rate() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool init = { false };
	static core::Real loop_combine_rate = 0.0;

	if ( !init ) {
		loop_combine_rate = option[ OptionKeys::loops::loop_combine_rate ]();
		init = true;
	}
	return loop_combine_rate;
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::get_combine_if_fail_exist() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool combine_loop = { true };
	static bool init = { false };

	if ( !init ) {
		combine_loop = option[ OptionKeys::loops::no_combine_if_fail];
		init = true;
	}
	return combine_loop;
}


///////////////////////////////////////////////////////////////////////////////
void LoopRebuild::barcode_extend_stems(
	core::pose::Pose & pose,
	int & barcst_extend_begin,
	int & barcst_extend_end,
	int & loop_begin,
	int & loop_end,
	int const & old_loop_begin,
	int const & old_loop_end,
	int const & nres,
	int const & selected_loop,
	int const & total_combine,
	int const & backward_combine
) {
	int const num_loop( static_cast<int>( Loops_in_.size() ) );
	// further extend loop regions, and add barcode constraints on these
	// extended regions
	int barcst_begin ( old_loop_begin );
	int barcst_end   ( old_loop_end   );

	barcst_extend_begin = std::min( 5, barcst_extend_begin );// extend 5 res at most
	barcst_extend_end   = std::min( 5, barcst_extend_end );// extend 5 res at most
	//		loop_update_active_cst_list ( barcst_begin, barcst_end, free_res,
	//														barcst_extend_begin, barcst_extend_end );
	barcst_begin -= barcst_extend_begin;
	barcst_end   += barcst_extend_end;
	while ( barcst_begin < 1    ) barcst_begin++;
	while ( barcst_end   > nres ) barcst_end--;

	// make sure not to extend the stems to loop regions
	if ( selected_loop - backward_combine > 1 )
		while ( barcst_begin <= (int) Loops_in_[ selected_loop - backward_combine -1 ].stop() + 1 )
			barcst_begin++;
	if ( selected_loop + total_combine != num_loop )
		while ( barcst_end >= (int) Loops_in_[ selected_loop + total_combine +1 ].start() - 1 )
			barcst_end--;

	// shorten loops when do looprlx ( loop_modeling)
	if( !loop_model() && shorten_long_terminal_loop() ){
		if( barcst_end == nres && barcst_end - barcst_begin > 10 )
			barcst_begin = barcst_end - 10;
		if( barcst_begin == 1 && barcst_end - barcst_begin > 10 )
			barcst_end = barcst_begin + 10;
	}

	// make sure loop length is greater than 3
	if ( barcst_end - barcst_begin < 2 && barcst_begin != 1 && barcst_end != nres ){
		barcst_begin -= 1;
		barcst_end   += 1;
	}
	//vats hack REMOVE THIS
	loop_begin = barcst_begin;
	loop_end   = barcst_end;
	//	loop_begin = loops_begin.at( 0 );
	//	loop_end = loops_end.at( 0 );

	// show the loop sequence for sanity-checking
	TR.Debug << "barcode_extend_stems():  loop_begin = " << loop_begin <<
							" loop_end =" << loop_end   <<
							" loop_sequence =";
	for ( int i=loop_begin; i<= loop_end; ++i )
		TR.Debug << " " << pose.residue( i ).name3();
	TR.Debug << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::shorten_long_terminal_loop() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool init = { false };
	static bool short_loop = { false };

	if ( !init ) {
		short_loop = option[ OptionKeys::loops::shorten_long_terminal_loop];
		init = true;
	}
	return short_loop;
}



//////////////////////////////////////////////////////////////////////////////////
void LoopRebuild::extend_barcode_regions_if_chain_break(
	core::pose::Pose & pose,
	int const & loop_begin,
	int const & loop_end,
	int & n_chain_break_fail,
	bool & is_chain_break,
	int & barcst_extend_begin, // output
	int & barcst_extend_end
) {
	static int n_small_chain_break_fail = { 1 };
	static bool barcst_flip_flop = { true };
	static bool barcst_small_flip_flop = { true };
	is_chain_break = false;

	// params control the behaviors of all loop modeling processes
	core::Real const chain_break_tol( get_chain_break_tolerance() );
	core::Real chain_break_score=0;

	int const nres = pose.total_residue();

	(*scorefxn_)( pose );
	if ( loop_begin != 1 && loop_end != nres ){
		chain_break_score = std::max( (float)pose.energies().total_energies()[ scoring::chainbreak ],
		                              (float)pose.energies().total_energies()[ scoring::linear_chainbreak ] );
		//			chain_break_score = pose.get_0D_score( pose_ns::CHAINBREAK );
	}

	TR.Debug << "loop_begin and loop_end  " << loop_begin << " " << loop_end << " " << chain_break_score << std::endl;
	scorefxn_->show(  TR , pose );
	TR.Debug << "*****\n*****\n*****\n";
	TR.Debug << pose.fold_tree() << std::endl;;
	TR.Debug << "*****\n*****\n*****\n";

	// extend more barcode regionsi
	TR.Info << "Chainbreak: " << chain_break_score << " Max: " << chain_break_tol << std::endl;
	if ( chain_break_score > chain_break_tol ){
		is_chain_break = true;
		n_chain_break_fail++;

		if( chain_break_score <= 10* chain_break_tol )
			n_small_chain_break_fail++;

		if ( get_random_loop_flag() ){
			if ( chain_break_score > 10* chain_break_tol ) {
				if( barcst_flip_flop ){
					//					if (loop_begin != 1 && loop_end != nres &&
					//							pose.secstruct( loop_begin - 1) != 'H' &&
					//							pose.secstruct( loop_begin - 1) != 'E' )
					barcst_extend_begin++;
					barcst_flip_flop = false;
				} else {
					barcst_extend_end++;
					barcst_flip_flop = true;
				}
			}
			else if ( n_small_chain_break_fail % 3 == 0 ) {
				if( barcst_small_flip_flop ){
					barcst_extend_begin++;
					barcst_small_flip_flop = false;
				} else {
					barcst_extend_end++;
					barcst_small_flip_flop = true;
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::acceptable_rmsd_change(
																			 core::pose::Pose & pose1,
																			 core::pose::Pose & pose2
																			 ) {
	runtime_assert( pose1.total_residue() == pose2.total_residue() );
	int natoms( pose1.total_residue() );
	using ObjexxFCL::FArray2D;
	FArray2D< core::Real > p1a( 3, natoms );
	FArray2D< core::Real > p2a( 3, natoms );
	//I am not sure, but I think '3' is CA. Check with someone if this is correct
	int CA_pos( 3 );
	core::Real rmsd_tol ( get_rmsd_tolerance() );

	int atom_count( 0 );
	for( int i = 1; i <= int( pose1.total_residue() ); ++i ) {
		const numeric::xyzVector< Real > & vec1( pose1.residue( i ).xyz( CA_pos ) );
		const numeric::xyzVector< Real > & vec2( pose2.residue( i ).xyz( CA_pos ) );
		atom_count++;
		for( int k = 0; k < 3; ++k ){
			p1a( k+1, atom_count ) = vec1[ k ];
			p2a( k+1, atom_count ) = vec2[ k ];
		}
	}
	core::Real rmsd( numeric::model_quality::rms_wrapper( natoms, p1a, p2a ) );
	TR.Debug << "RMSd change " << rmsd << std::endl;
	if ( rmsd <= rmsd_tol ) return true;

	return false;
}


///////////////////////////////////////////////////////////////////////
core::Real LoopRebuild::get_looprlx_cycle_ratio() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool init = { false };
	static core::Real cycle_ratio = 1.0;

	if ( !init ) {
		cycle_ratio = option[ OptionKeys::loops::looprlx_cycle_ratio ]();
		init = true;
	}
	return cycle_ratio;
}


//////////////////////////////////////////////////////////////////////////////////
bool LoopRebuild::extended_loop() {
	if ( use_default_extend_loops ) {
		// from now on, use stored value
		extend_loops = basic::options::option[ basic::options::OptionKeys::loops::extended ];
		use_default_extend_loops = false;
	}

	return extend_loops;
}


//////////////////////////////////////////////////////////////////////////////////
void LoopRebuild::set_extended_loop( bool val ) {
	extend_loops = val;
	use_default_extend_loops = false;
}


//////////////////////////////////////////////////////////
///@brief Helper function tokenizes a str
/////////////////////////////////////////////////////////
// don't use use string_util.hh:split or string_split instead
//
void Tokenize(const std::string              &str,
              utility::vector1<std::string>  &tokens,
              const std::string              &delimiters=" ") {
	tokens.clear();
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}


///////////////////////////////////////////////////////////////////////
/// @brief returns mc object
moves::MonteCarloOP LoopRebuild::get_mc( core::pose::Pose & pose )
{

	set_default_mc( pose );
	mc_created = true;
	return mc_;

}

void LoopRebuild::set_default_settings(){
	use_default_extend_loops = true;
	mc_created = false;
	random_loop_flag_ = basic::options::option[ basic::options::OptionKeys::loops::random_loop ];
	allowed_failure_before_extend_ = 15;
	allowed_failure_before_stop_ = random_loop_flag_ ? 20:5;
	abort_on_failed_loop_ = true;
}


///////////////////////////////////////////////////////////////////////
///@brief sets mc object
void LoopRebuild::set_default_mc( core::pose::Pose & pose )
{
	m_Temperature_ = 2.0;
	mc_ = protocols::moves::MonteCarloOP( new moves::MonteCarlo( pose, *scorefxn_, m_Temperature_ ) );
}
///////////////////////////////////////////////////////////////////////
///@brief Full atom loop refinement
void LoopRefine::apply(
	core::pose::Pose & pose
)
{

	TR << "***** Starting full-atom loop refinement protocol  ****" << std::endl;

    protocols::loops::LoopsOP LoopsToRefine( new protocols::loops::Loops() );

    for( Loops::const_iterator it=Loops_in_.begin(), it_end=Loops_in_.end(); it != it_end; ++it ) {
		Loop refine_loop( *it );
		refine_loop.choose_cutpoint( pose );
		LoopsToRefine->add_loop( refine_loop );
	}

 	core::kinematics::FoldTree f;
	protocols::loops::fold_tree_from_loops( pose, *LoopsToRefine, f );
	pose.fold_tree( f );

	//protocols::loops::refine_loops_with_ccd( pose, pose, LoopsToRefine );
	loop_mover::refine::LoopMover_Refine_CCD refine_ccd( LoopsToRefine );
	refine_ccd.set_native_pose( PoseCOP( new pose::Pose ( pose ) ) );
	refine_ccd.apply( pose );

}

std::string
LoopRefine::get_name() const {
	return "LoopRefine";
}


///////////////////////////////////////////////////////////////////////

} // namespace protocols
