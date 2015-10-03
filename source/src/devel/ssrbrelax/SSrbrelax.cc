// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to looprelax
/// @details
/// @author Vatsan Raman
#include <devel/ssrbrelax/SSrbrelax.hh>
#include <devel/ssrbrelax/SSRbClass.hh>
#include <devel/loops/looprelax_protocols.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>

#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <protocols/frags/VallData.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/ccd_closure.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

// External library headers
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>


//C++ headers
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <utility>
#include <set>
#include <list>
#include <cstdlib>

// option key includes

#include <basic/options/keys/SSrbrelax.OptionKeys.gen.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using io::pdb::dump_pdb;

namespace devel {
namespace ssrbrelax {

	using namespace core;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	using io::pdb::dump_pdb;

	typedef numeric::xyzVector < Real > Vector;
// sets up the default mc object
/*
	void LoopRelax::set_default_mc() {
		m_Temperature_ = 0.8;
		mc_ = new protocols::moves::MonteCarlo( *pose_, *scorefxn_ , m_Temperature_ );
	}
	*/
	///////////////////////////////////////////////////////////////////////
	/*	protocols::moves::MonteCarloOP LoopRelax::get_mc() {
		return mc_;
	}

	void LoopRelax::set_default_move_map(){
		movemap_ = new core::kinematics::MoveMap();
		movemap_->set_bb ( true );
		movemap_->set_chi( true );
		}*/
	///////////////////////////////////////////////////////////////////////
	void RbRelax::apply( core::pose::Pose & pose ){

		using namespace protocols::moves;
		using namespace scoring;

		T("protocols.rbrelax") <<  "I don't know how to do anything yet!";
		rbrelax_main( pose );
	}
	///////////////////////////////////////////////////////////////////////
	void RbRelax::apply() {
		apply( *pose_ );
	}
	///////////////////////////////////////////////////////////////////////
	void RbRelax::rbrelax_main(
														core::pose::Pose & pose
														)
	{
		std::cout << "Inside rbrelax_main " << std::endl;
		//The pose will always get read in as a full-atom pose; even if a centroid pose is supplied
		pose::Pose save_input_pose;
		save_input_pose = pose;
		RbSegments rbsegments;

		//converting to centroid atom types
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		//std::string rb_file;
		rbsegments.read_segments_from_file();
		rbsegments.read_param_file();
		protocols::loops::set_secstruct_from_psipred_ss2( pose );


		//Initialzing fragments for loop building
		std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
		initialize_fragments( frag_libs );
		segment_rb_move( pose, rbsegments, frag_libs );

	}

	///////////////////////////////////////////////////////////////////////

	void RbRelax::segment_rb_move(
																core::pose::Pose & pose,
																protocols::ssrbrelax::RbSegments & rbsegments,
																std::map< Size, devel::frags::TorsionFragmentLibraryOP > & frag_libs
																)
	{
		using namespace kinematics;
		kinematics::FoldTree f, f_reset;
		int const nres( pose.total_residue() );
		f_reset.add_edge( 1, nres, Edge::PEPTIDE );
		while( rbsegments.num_rb() > 0 ){
			RbSegments::iterator it ( rbsegments.one_random_segment() ); //pick one random segment
			RbSegments this_segment; //one segment object to make one segment fold tree
			this_segment.add_segment( it );
			rbsegments.delete_segment( it );
			this_segment.one_segment_fold_tree( f, nres );
			pose.fold_tree( f );
			perturb_segment_and_close_loops( pose, rbsegments, this_segment, frag_libs );
			std::cout << "Finished this segment " << std::endl;
			std::cout << "f_reset " << f_reset << std::endl;
			pose.fold_tree( f_reset );

		}
		pose.dump_pdb("all_segments_done.pdb");
		std::cout << "Do nothing " << std::endl;
	}

	///////////////////////////////////////////////////////////////////////

	void RbRelax::perturb_segment_and_close_loops(
																								core::pose::Pose & pose,
																								devel::ssrbrelax::RbSegments & rbsegments,
																								devel::ssrbrelax::RbSegments & this_segment,
																								std::map< Size, protocols::frags::TorsionFragmentLibraryOP > & frag_libs
																								)
	{

		using namespace conformation;

		std::cout << "Inside perturb segment" << std::endl;
		int index( this_segment.index() );
		std::cout << "this_segment.index() " << this_segment.index() << std::endl;
		int n_dof( 6 ); // all six degrees of freedom
		std::vector < int > dof_vec;

		for ( int i = 1; i <= n_dof; ++i ){
			if ( rbsegments.distribution_exists( index, i ) ){
				std::cout << "dof_vec.push_back " << index << " " << i << std::endl;
				dof_vec.push_back( i );
			}
		}
		std::random_shuffle( dof_vec.begin(), dof_vec.end() ); //shuffle the dof vector

		while( dof_vec.size() > 0 ){
			int size( dof_vec.size() );
			int this_dof( dof_vec.at( size - 1 ) ); //get the last element of the shuffled vector
			std::cout << "dof_vec.size before " << dof_vec.size() << std::endl;
			dof_vec.pop_back(); //remove the last element of the vector
			std::cout << "dof_vec.size after " << dof_vec.size() << " " << this_dof <<  std::endl;
			float stddev = rbsegments.get_gaussian_parameters( index, this_dof );
			std::cout << "stddev " << rbsegments.get_gaussian_parameters( index, this_dof ) << " " << stddev << std::endl;
		debug_assert ( stddev > 0.0 );

			std::cout << "Input arguments " <<  this_dof << " " << stddev << std::endl;
			perturb_segment( pose, this_segment, this_dof, stddev );
			pose.dump_pdb("after_segment_perturb.pdb");
		}
		close_both_loops( pose, this_segment, frag_libs );
		pose.dump_pdb("one_segment.pdb");


		//		refine_segment( pose, this_segment );
	}


	///////////////////////////////////////////////////////////////////////

	void RbRelax::perturb_segment(
																core::pose::Pose & pose,
																devel::ssrbrelax::RbSegments & this_segment,
																int const & dof,
																float const & stddev
																)
	{

		using namespace kinematics;
		using namespace scoring;
		using namespace optimization;
		using namespace conformation;

		assert( dof > 0 && dof <= 6 );

		int const cycles1( 2 ); //arbirary now
		int const cycles2( 4 ); //arbitrary as well

		int const flexible_jump_number( 1 ); //Ensure that this is one jump fold tree  !!
		float const step_size( 0.2*stddev ); //This is a fraction of the std deviation
		//just to suppress compiler warning
		std::cout << "step_size " << step_size << std::endl;
		kinematics::MoveMap rb_move_map;
		set_rbrelax_allow_move_map( rb_move_map, flexible_jump_number );

		core::kinematics::Jump flexible_jump( pose.jump( flexible_jump_number ) );

		std::string function_tag("cen_std"),patch_tag("score3");
		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( function_tag, patch_tag ) );
		( *scorefxn )( pose );
		/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
		scorefxn->show( std::cout );

		//Monte Carlo object
		float init_temp( 2.0 ), final_temp( 1.0 );
		protocols::moves::MonteCarlo mc( pose, *scorefxn, init_temp );
		mc.show_scores();

		float const gamma =
			std::pow( (final_temp/init_temp), (1.0f/(cycles1*cycles2)) );

		//Minimizer object
		AtomTreeMinimizer rb_mzr;
		float const dummy_tol( 0.001 );
		bool const use_nblist( true ), deriv_check( false );
		MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check );

		float temperature = init_temp;
		float rb_tmp( 0.0 );


		//Coordinate frame set up for rotation/translation for this_segment
		int const anchor_pos( this_segment.anchor_res() );
		//		if ( this_segment.get_flex_jump_pos( this_segment.seg_start(), this_segment.seg_stop() ) != 0 )
		//			int const flex_jump_pos( this_segment.get_flex_jump_pos( this_segment.seg_start(), this_segment.seg_stop() ) );

		Residue const & rsd_anchor( pose.residue( anchor_pos ) );
		//		Residue const & rsd_flex_jump_pos( pose.residue( flex_jump_pos ) );

		//		kinematics::Stub Stub_fixed_pos( rsd_anchor.xyz("CA"),rsd_anchor.xyz("N"),
		//																		 rsd_anchor.xyz("CA"),rsd_anchor.xyz("C") );

		//		kinematics::Stub Stub_fixed_pos( rsd_anchor.xyz("C"),rsd_anchor.xyz("N"),
		//															 rsd_anchor.xyz("CA"),rsd_anchor.xyz("C") );
		kinematics::Stub test_stub( pose.conformation().upstream_jump_stub( flexible_jump_number ));

		//DoF Nomenclature   1    2    3         4        5         6
		// TRANS along       X    Y    Z   ROT about X  about Y  about Z

		//		numeric::xyzVector < Real > X_axis( this_segment.X_axis( pose ) );
		Vector X_axis( this_segment.alt_X_axis( pose ) );
		Vector Y_axis( this_segment.alt_Y_axis( pose ) );
		Vector Z_axis( this_segment.alt_Z_axis( pose ) );
		Vector Perturbation_axis;

		Vector Center( rsd_anchor.xyz("CA") );//center is set at Calpha of
		//		pose.dump_pdb("before_perturbation.pdb");
		if ( dof == 1 || dof == 4 )
			Perturbation_axis = X_axis;
		else if ( dof == 2 || dof == 5 )
			Perturbation_axis = Y_axis;
		else if ( dof == 3 || dof == 6 )
			Perturbation_axis = Z_axis;
		else {
			Error() << "The degree of freedom of perturbation should lie between 1-6. It is " << dof << "\n";
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		float cumulative_rb_tmp( 0.0 );
		for( int c1 = 1; c1 <= cycles1; ++c1 ) {
			mc.recover_low( pose );
			( *scorefxn )( pose );
			mc.show_scores();

			for( int c2 = 1; c2 <= cycles2; ++c2 ) {
				temperature *= gamma;
				mc.set_temperature( temperature );
				rb_tmp = (-1.0)*rb_tmp + numeric::random::gaussian()*stddev;
				//**********BEWARE STEP SIZE TURNED OFF ***************
				//						rb_tmp = (-1.0)*rb_tmp + numeric::random::gaussian()*stddev*step_size;
				//set up calls translation_along_axis and rotation_by_axis
				//				flexible_jump.set_rb_delta( dof, dir, rb_tmp );
				cumulative_rb_tmp += rb_tmp;
				if ( dof <= 3 )
					flexible_jump.translation_along_axis( test_stub, Perturbation_axis, rb_tmp );
				else if ( dof > 3 )
					flexible_jump.rotation_by_axis( test_stub, Perturbation_axis, Center, rb_tmp );
				std::cout << "dof and rb_tmp and cumulative_rb_tmp " << dof << " " << rb_tmp << " " << cumulative_rb_tmp << std::endl;
				pose.set_jump( flexible_jump_number, flexible_jump );
				pose.dump_pdb("perturb_"+string_of(c1)+"_"+string_of(c2)+".pdb");
				rb_mzr.run( pose, rb_move_map, *scorefxn, options );
				//The final change in angle or distance has to be done after minimization. Stated differently, rb_tmp has to be computed after minimization not after perturbation.
				mc.boltzmann( pose );
				mc.show_scores();

			}

		}
		//		pose.dump_pdb("after_perturbation.pdb");
		//		pose = mc.lowest_score_pose();      ****** BE SURE TO UNCOMMENT THIS *********
		rb_move_map.set_jump( flexible_jump_number, false );
	}

	///////////////////////////////////////////////////////////////////////

		void RbRelax::set_rbrelax_allow_move_map(
																						 core::kinematics::MoveMap & rb_move_map,
																						 int const & flexible_jump
																						 )
		{
			rb_move_map.set_bb( false );
			rb_move_map.set_chi( false );
			rb_move_map.set_jump( false );

			rb_move_map.set_jump( flexible_jump, true );

		}


	///////////////////////////////////////////////////////////////////////

	void RbRelax::initialize_fragments(
																		 std::map< Size, protocols::frags::TorsionFragmentLibraryOP > & frag_libs
																		 )

	{

		//fragment initialization : This has to be changed to the new variable fragment length implementation. This is a good place to start.

		using namespace basic::options;
		utility::vector1<int> frag_sizes( option[ OptionKeys::SSrbrelax::frag_sizes ] );
		FileVectorOption frag_files( option[OptionKeys::SSrbrelax::frag_files] );
		assert( frag_sizes.size() == frag_files.size() );

		std::map< Size, bool > frag_libs_init;
		for ( Size i = 1; i <= frag_sizes.size(); ++i ) {
			Size const frag_size = Size(frag_sizes[i]);
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			std::cout << "Frag libraries debug " << frag_files[i] << " " << frag_size << std::endl;
			frag_libs_init.insert( std::make_pair(frag_size, frag_lib_op->read_file( frag_files[i], frag_size, 3 ) ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
		}

		Size prev_size(10000);
		protocols::frags::TorsionFragmentLibraryOP prev_lib_op(0);
		for ( std::map<Size, bool>::const_reverse_iterator it = frag_libs_init.rbegin(),
						it_end = frag_libs_init.rend(); it != it_end; ++it ) {
			Size const frag_size( it->first );
			bool const frag_lib_init( it->second );
			assert( frag_size < prev_size );
			if ( (!frag_lib_init) && prev_lib_op ) {
				std::cout << "set up " << frag_size << "-mer library from " << prev_size << "-mer library" << std::endl;
				protocols::frags::TorsionFragmentLibraryOP current_lib_op( frag_libs.find(frag_size)->second );
				frag_libs_init[frag_size] = current_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
			}
			prev_size = frag_size;
			prev_lib_op = frag_libs[frag_size];
			std::cout << "frag_libs_init: " << frag_size << " " << frag_libs_init[frag_size] << std::endl;
		}

	}

	///////////////////////////////////////////////////////////////////////

	void RbRelax::close_both_loops(
																 core::pose::Pose & pose,
																 devel::ssrbrelax::RbSegments & this_segment,
																 std::map< Size, protocols::frags::TorsionFragmentLibraryOP > & frag_libs
																 )
	{

		bool nterm_loop_built( false );
		bool cterm_loop_built( false );

		devel::looprelax::LoopRelax thisLoopRelax;

		//now to close loops on both sides with ccd closure
		while ( !nterm_loop_built || !cterm_loop_built ) {

			if ( ( numeric::random::uniform() < 0.5 && !nterm_loop_built )  || cterm_loop_built ) { //start with n or c-term loop closure stochastically
				int loop_begin( this_segment.n_term_loop() );
				int loop_end( this_segment.seg_start() );
				int cutpoint( this_segment.seg_start() );
				std::cout << "nterm loop " << loop_begin << " " << loop_end << " " << cutpoint << std::endl;
				thisLoopRelax.build_loop_with_ccd_closure( pose, loop_begin, loop_end, cutpoint, frag_libs );
				nterm_loop_built = true;
				std::cout << "Built n-term loop" << std::endl;
				pose.dump_pdb("built_nterm_loop.pdb");

			} else if ( !cterm_loop_built || nterm_loop_built ) {
				int loop_begin( this_segment.seg_stop() );
				int loop_end( this_segment.c_term_loop() );
				int cutpoint( this_segment.seg_stop() );
				std::cout << "cterm loop " << loop_begin << " " << loop_end << " " << cutpoint << std::endl;
				thisLoopRelax.build_loop_with_ccd_closure( pose, loop_begin, loop_end, cutpoint, frag_libs );
				cterm_loop_built = true;
				std::cout << "Built c-term loop" << std::endl;
				pose.dump_pdb("built_cterm_loop.pdb");

			}
			std::cout << "n-term loop built " << nterm_loop_built << " c-term loop built " << cterm_loop_built << std::endl;

		}


	}

	///////////////////////////////////////////////////////////////////////
	void RbRelax::refine_segment(
															 core::pose::Pose & pose,
															 devel::ssrbrelax::RbSegments & this_segment
															 )
	{
		using namespace kinematics;
		using namespace scoring;
		using namespace optimization;

		int const nres( pose.total_residue() );
		int const seg_length( this_segment.c_term_loop() - this_segment.n_term_loop() );

		core::kinematics::MoveMapOP refine_move_map ( new kinematics::MoveMap );
		( *refine_move_map ).set_bb( false );
		( *refine_move_map ).set_chi( false );
		( *refine_move_map ).set_jump( false );

		for ( int i = 1; i <= nres; ++i ) {
			if ( i >= this_segment.n_term_loop() && i <= this_segment.c_term_loop() ) {
				( *refine_move_map ).set_bb( i, true );
			}
		}

		//Local looprelax object to access its member functions
		devel::looprelax::LoopRelax localLPObj;
		int const cutpoint ( localLPObj.choose_cutpoint( pose, this_segment.n_term_loop(), this_segment.c_term_loop() ) );

		//set cutpoint variant residue for chainbreak score
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint+1 );

		std::string function_tag("cen_std"), patch_tag("score4L");
		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( function_tag, patch_tag ) );

		float const start_chainbreak_weight( scorefxn->get_weight ( chainbreak ) );
		float const max_chainbreak_weight( 3.0 );

	debug_assert ( max_chainbreak_weight >= start_chainbreak_weight );
		float const chainbreak_wt_stepsize( ( max_chainbreak_weight - start_chainbreak_weight ) / seg_length );

		(*scorefxn)(pose);
		/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
		scorefxn->show( std::cout );

		//Monte Carlo
		float const init_temp( 2.0 ), final_temp( 1.0 );
		float const gamma = std::pow( (final_temp/init_temp), 1.0f/(seg_length) );
		float temperature = init_temp;
		protocols::moves::MonteCarlo mc( pose, *scorefxn, init_temp );
		mc.show_scores();

		//Minimizer
		AtomTreeMinimizer refine_mzr;
		float const dummy_tol( 0.001 );
		bool const use_nblist( true ), deriv_check( false );
		MinimizerOptions options( "linmin", dummy_tol, use_nblist, deriv_check );

		// small/shear move parameters
		Size const nmoves = { 1 };
		std::map< char, Real > angle_max;
		angle_max.insert( std::make_pair('H', 0.0) );
		angle_max.insert( std::make_pair('E', 5.0) );
		angle_max.insert( std::make_pair('L', 6.0) );

		for ( int i = 1; i <= seg_length; ++i ) {
			temperature *= gamma;
			if ( i%2 == 0 ) { //reduce temperature and raise weight every two steps
				mc.set_temperature( temperature );
				float chainbreak_weight( start_chainbreak_weight + float ( chainbreak_wt_stepsize*i ) );
				scorefxn->set_weight( chainbreak, chainbreak_weight );
			}
			protocols::simple_moves::SmallMover small_moves( refine_move_map, temperature, nmoves );
			small_moves.apply( pose );
			protocols::simple_moves::ShearMover shear_moves( refine_move_map, temperature, nmoves );
			shear_moves.apply( pose );
			localLPObj.fast_ccd_close_loops( pose, this_segment.n_term_loop(),
																			 this_segment.c_term_loop(), cutpoint, *refine_move_map );
			refine_mzr.run( pose, *refine_move_map, *scorefxn, options );
			std::string move_type = "small-shear-ccd-min";
			mc.boltzmann( pose, move_type );
			mc.show_scores();
		}

		mc.show_counters();
		pose = mc.lowest_score_pose();

	}

	///////////////////////////////////////////////////////////////////////
} // namespace rbrelax
} // namespace devel
