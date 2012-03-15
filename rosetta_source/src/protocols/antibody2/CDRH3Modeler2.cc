// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c)University of Washington UW TechTransfer, email:license@u.washington.edu.

/// @file antibody2/CDRH3Modeler2.cc
/// @brief models CDR H3 loop using loop modeling
/// @detailed
///// @author Jianqing Xu (xubest@gmail.com)
//

// Rosetta Headers
#include <protocols/antibody2/CDRH3Modeler2.hh>

#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/chemical/VariantType.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/id/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
//#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/util.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>

#include <protocols/antibody2/Ab_util.hh>

#include <protocols/antibody2/Ab_H3_cter_insert_mover.hh>
#include <protocols/antibody2/Ab_Relax_a_CDR_FullAtom.hh>



static numeric::random::RandomGenerator RG(21141980);

static basic::Tracer TR("protocols.antibody2.CDRH3Modeler2");
using namespace core;

namespace protocols {
namespace antibody2 {

CDRH3Modeler2::CDRH3Modeler2() : Mover()
{
	user_defined_ = false;
	init( false, true, true, false, false );
}

CDRH3Modeler2::CDRH3Modeler2(
	bool model_h3,
	bool apply_centroid_mode,
	bool apply_fullatom_mode,
	bool camelid,
	bool benchmark
	) : Mover()
{
	user_defined_ = true;
	init( model_h3, apply_centroid_mode, apply_fullatom_mode, camelid, benchmark );
}

CDRH3Modeler2::CDRH3Modeler2(
	utility::vector1< fragment::FragSetOP > cdr_h3_frags) : Mover()
{
	user_defined_ = false;
	init( false, true, true, false, false );
	cdr_h3_frags_ = cdr_h3_frags;
}

void CDRH3Modeler2::init(
	bool model_h3,
	bool apply_centroid_mode,
	bool apply_fullatom_mode,
	bool camelid,
	bool benchmark
	)
{
	Mover::type( "CDRH3Modeler2" );

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
	set_default();
//	register_options();
//	init_from_options();
	if ( user_defined_ ) {
		model_h3_ = model_h3;
		apply_centroid_mode_ = apply_centroid_mode;
		apply_fullatom_mode_ = apply_fullatom_mode;
		set_camelid( camelid );
		benchmark_ = benchmark;
	}

	lowres_scorefxn_ = scoring::ScoreFunctionFactory::
		create_score_function( "cen_std", "score4L" );
	lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
	// adding constraints
	lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );

	highres_scorefxn_ = scoring::ScoreFunctionFactory::
		create_score_function("standard", "score12" );
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
	// adding constraints
	highres_scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );

	// set up objects based on the boolean values defined above
	setup_objects();
}

    
void CDRH3Modeler2::setup_objects(){
    ab_h3_cter_insert_mover_ = new Ab_H3_cter_insert_mover(antibody_in_);
    //TODO:  
    //JQX: right now, just want the code to work, whether this antibody_in_ is at the right position or 
    // wehther it has been initialized or not? I don't know. Will come back to address this
}
    
    
    
    
    
    
// CDRH3Modeler2 default destructor
CDRH3Modeler2::~CDRH3Modeler2() {}

void CDRH3Modeler2::set_default()
{
	benchmark_ = false;
	model_h3_ = false;
	c_ter_stem_ = 3;
	cen_cst_ = 10.0;
	high_cst_ = 100.0; // if changed here, please change at the end of AntibodyModeler as well
	max_cycle_ = 20;

	apply_centroid_mode_ = false;
	apply_fullatom_mode_ = false;

	current_loop_is_H3_ = true;
	H3_filter_ = true;
	antibody_refine_ = true;
	snug_fit_ = true;
	loops_flag_ = true;
	docking_local_refine_ = true;
	dle_flag_ = true;
	refine_input_loop_ = true;
	is_camelid_ = false;
	// size of loop above which 9mer frags are used
	cutoff_9_ = 16; // default 16
	// size of loop above which 3mer frags are used
	cutoff_3_ = 6; // default 6

	TR << "H3M Finished Setting Defaults" << std::endl;
} // CDRH3Modeler2 set_default

    
    
    
    
    
		void CDRH3Modeler2::set_lowres_score_func(
      scoring::ScoreFunctionOP lowres_scorefxn ) {
			lowres_scorefxn_ = lowres_scorefxn;
		} // set_lowres_score_func

		void CDRH3Modeler2::set_highres_score_func(
      scoring::ScoreFunctionOP highres_scorefxn) {
			highres_scorefxn_ = highres_scorefxn;
		} // set_highres_score_func

    
    
    
        //################################################
        //###########  apply function ####################
        //################################################

		void CDRH3Modeler2::apply( pose::Pose & pose_in )
		{
			if( !model_h3_ )
				return;

			TR << "H3M Applying CDR H3 modeler" << std::endl;

			using namespace core::pose;
			using namespace core::scoring;
			using namespace protocols::moves;

			start_pose_ = pose_in;
			antibody_in_.setup_CDR_loops( pose_in, is_camelid_ );
			setup_packer_task( pose_in, tf_ );
			pose::Pose start_pose = pose_in;

			if( is_camelid_ && !antibody_in_.is_extended() && !antibody_in_.is_kinked() )
				c_ter_stem_ = 0;

			Size framework_loop_begin( antibody_in_.get_CDR_loop("h3")->start() );
			Size frmrk_loop_end_plus_one( antibody_in_.get_CDR_loop("h3")->stop() );
			//Size framework_loop_size = ( frmrk_loop_end_plus_one -
			//														framework_loop_begin ) + 1;
			Size cutpoint = framework_loop_begin + 1;
			loops::Loop cdr_h3( framework_loop_begin, frmrk_loop_end_plus_one,
													cutpoint,	0, true );
			simple_one_loop_fold_tree( pose_in, cdr_h3 );

			// switching to centroid mode
			simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
			simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

			// Building centroid mode loop
			if( apply_centroid_mode_ ) {
                //#############################
				to_centroid.apply( pose_in );
                //#############################
				build_centroid_loop( pose_in );
				if( is_camelid_ )
					loop_centroid_relax( pose_in, antibody_in_.get_CDR_loop("h1")->start(),
															 antibody_in_.get_CDR_loop("h1")->stop() );
				to_full_atom.apply( pose_in );

				utility::vector1<bool> allow_chi_copy( pose_in.total_residue(),
																							 true );
				for( Size ii = antibody_in_.get_CDR_loop("h3")->start();
						 ii <= ( antibody_in_.get_CDR_loop("h3")->stop() ); ii++ )
					allow_chi_copy[ii] = false;
				//recover sidechains from starting structures
				protocols::simple_moves::ReturnSidechainMover recover_sidechains(
					start_pose_, allow_chi_copy );
				recover_sidechains.apply( pose_in );

				// Packer
				protocols::simple_moves::PackRotamersMoverOP packer;
				packer = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
				packer->task_factory(tf_);
				packer->apply( pose_in );
			}

			if( apply_fullatom_mode_ ) {
                //##############################
                Ab_Relax_a_CDR_FullAtom relax_a_cdr_high_res(current_loop_is_H3_, H3_filter_); 
                relax_a_cdr_high_res.pass_start_pose(start_pose_);
                relax_a_cdr_high_res.apply(pose_in);
				//build_fullatom_loop( pose_in );
                //##############################
				if( !benchmark_ ) {
					Size repack_cycles(1);
					if( antibody_refine_ && !snug_fit_ )
						repack_cycles = 3;
					protocols::simple_moves::PackRotamersMoverOP packer;
					packer = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
					packer->task_factory(tf_);
					packer->nloop( repack_cycles );
					packer->apply( pose_in );
				}
			}

			// Minimize CDR H2 loop if this is a camelid
            
			if( is_camelid_ ) {
                //##############################
                Ab_Relax_a_CDR_FullAtom relax_a_cdr_high_res(false, false, is_camelid_); // because of h2
                relax_a_cdr_high_res.apply(pose_in);
                //##############################

                //JQX: remove the duplicated code, camelid H2 will be automatically taken care of
                //     see the code in Ab_Relax_a_CDR_FullAtom file
			}



			TR << "H3M Finished applying CDR H3 modeler" << std::endl;

			return;
		} // CDRH3Modeler2::apply()

    
    
    
std::string
CDRH3Modeler2::get_name() const {
	return "CDRH3Modeler2";
}


		void CDRH3Modeler2::build_centroid_loop( core::pose::Pose & pose ) {
			using namespace core::pose;
			using namespace core::scoring;
			using namespace protocols::moves;

			if( !apply_centroid_mode_ )
				return;

			TR <<  "H3M Modeling Centroid CDR H3 loop" << std::endl;

			Size frmrk_loop_end_plus_one( antibody_in_.get_CDR_loop("h3")->stop() );
			Size framework_loop_size = ( frmrk_loop_end_plus_one - antibody_in_.get_CDR_loop("h3")->start() ) + 1;
			Size cutpoint = antibody_in_.get_CDR_loop("h3")->start() + 1;
			loops::Loop cdr_h3( antibody_in_.get_CDR_loop("h3")->start(), frmrk_loop_end_plus_one,
													cutpoint,	0, true );
			simple_one_loop_fold_tree( pose, cdr_h3 );

			// silly hack to make extended loops to work
			loops::LoopsOP cdr_h3_loop_list = new loops::Loops();
			cdr_h3_loop_list->add_loop( cdr_h3 );
            /* Commented out by BDW with JX's consent
			loops::loop_mover::LoopMoverOP my_loop_move = new loops::loop_mover::LoopMover( cdr_h3_loop_list );
			my_loop_move->set_extended_torsions( pose, cdr_h3 );
			my_loop_move->apply( pose );
            */
			Size unaligned_cdr_loop_begin(0), unaligned_cdr_loop_end(0);
			std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
			core::import_pose::pose_from_pdb( template_pose_, path+"hfr.pdb" );
			std::string template_name = "h3";
			antibody2::Ab_Info hfr_template( template_pose_, template_name );
			unaligned_cdr_loop_begin = hfr_template.current_start;
			unaligned_cdr_loop_end = hfr_template.current_end;

			pose.set_psi( antibody_in_.get_CDR_loop("h3")->start() - 1,
				template_pose_.psi( unaligned_cdr_loop_begin - 1 ) );
			pose.set_omega(antibody_in_.get_CDR_loop("h3")->start() - 1,
				template_pose_.omega( unaligned_cdr_loop_begin - 1 ) );

			Size modified_framework_loop_end = frmrk_loop_end_plus_one - c_ter_stem_;
			loops::Loop trimmed_cdr_h3( antibody_in_.get_CDR_loop("h3")->start(),
				modified_framework_loop_end, cutpoint, 0, true );

			antibody2::Ab_Info starting_antibody;
			starting_antibody = antibody_in_;
			bool closed_cutpoints( false );

			Size cycle ( 1 );
			while( !closed_cutpoints && cycle < max_cycle_) {
				antibody_in_ = starting_antibody;
				if( framework_loop_size > 5 ){
                    ab_h3_cter_insert_mover_->apply(pose);
                }
				scored_frag_close( pose, trimmed_cdr_h3 );
				if( trimmed_cdr_h3.size() > cutoff_9_  ) { // aroop_temp default cutoff_9_
					Size saved_cutoff_9 = cutoff_9_;
					cutoff_9_ = 100; // never going to reach
					scored_frag_close( pose, trimmed_cdr_h3 );
					cutoff_9_ = saved_cutoff_9; // restoring
				}
				closed_cutpoints = cutpoints_separation( pose, antibody_in_ );
				++cycle;
			} // while( ( cut_separation > 1.9 )

			TR <<  "H3M Finished Modeling Centroid CDR H3 loop" << std::endl;

			return;

		} // build_centroid_loop


		void CDRH3Modeler2::set_offset_frags(
			utility::vector1< core::fragment::FragSetOP > & offset_frags ) {
			cdr_h3_frags_ = offset_frags;
		}
    
    
      
    
    

		///////////////////////////////////////////////////////////////////////////
		/// @begin scored_frag_close
		///
		/// @brief builds a loop from fragments file.
		///
		/// @detailed Loop is built by a monte carlo simulation using fragments
		///           from a fragment files. CCD moves are used to close loops
		///           with gaps at cutpoint.H3_check is enforced if H3_filter flag
		///           is set in command line. Loop building results in many files
		///           containing differnt conformations of the same loop in
		///           phi-psi-omega angles. Parallel processing is utilized.
		///
		/// @param[in] weight_map: in this case its a centroid weight
		///            pose_in: loop to be built on this template provided
		///            loop_begin/loop_end: loop termini definition
		///            frag_size: 3-mer or 9-mer
		///            frag_offset:agreement in frag file numbering & pose numberng
		///            cycles1: max cycles to be spent building loops
		///            cycles2: # of fragment swaps for each loop(depends on size)
		///            do_ccd_moves: should ccd moves be used to close gaps
		///
		/// @global_read benchmark_
		///
		/// @global_write
		///
		/// @remarks
		///
		/// @references
		///
		/// @authors Aroop 02/04/2010
		///
		/// @last_modified 02/04/2010
		///////////////////////////////////////////////////////////////////////////
		void CDRH3Modeler2::scored_frag_close (
			pose::Pose & pose_in,
			loops::Loop const trimmed_cdr_h3 ) {
			using namespace fragment;
			using namespace protocols;
			using namespace protocols::simple_moves;
			using namespace protocols::loops;
            using loop_closure::ccd::CcdMover;
            using loop_closure::ccd::CcdMoverOP;
            using loop_closure::ccd::CcdLoopClosureMover;
            using loop_closure::ccd::CcdLoopClosureMoverOP;

			TR <<  "H3M Fragments based centroid CDR H3 loop building" << std::endl;

			if( trimmed_cdr_h3.size() <= 2)
				utility_exit_with_message("Loop too small for modeling");

			// set cutpoint variants for correct chainbreak scoring
			if( !pose_in.residue( trimmed_cdr_h3.cut() ).is_upper_terminus() ) {
				if( !pose_in.residue( trimmed_cdr_h3.cut() ).has_variant_type(chemical::CUTPOINT_LOWER))
					core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, trimmed_cdr_h3.cut() );
				if( !pose_in.residue( trimmed_cdr_h3.cut() + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
					core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, trimmed_cdr_h3.cut() + 1 );
			}


			Size cycles1(10);
			// aroop_temp default 25 * loop size
			Size cycles2(25 * trimmed_cdr_h3.size() );

			// params
			Real const ccd_threshold( 0.1);
			Size h3_attempts(0);
			Real h3_fraction = 0.75; // 75% of loops are required to be H3's
			Real current_h3_prob = RG.uniform();;
			bool H3_found_ever(false);
			bool loop_found(false);
			Size total_cycles(0);
			Size frag_size(0);
			FragSetOP frags_to_use;
			{
				if( trimmed_cdr_h3.size() > cutoff_9_ ) {
					frags_to_use = cdr_h3_frags_[1]->empty_clone();
					frags_to_use = cdr_h3_frags_[1];
					frag_size = 9;
				}
				else {
					frags_to_use = cdr_h3_frags_[2]->empty_clone();
					frags_to_use = cdr_h3_frags_[2];
					frag_size = 3;
				}
			}

			// Storing Fold Tree
			kinematics::FoldTree old_fold_tree = pose_in.fold_tree();
			// New Fold Tree
			simple_one_loop_fold_tree( pose_in, trimmed_cdr_h3 );

			//setting MoveMap
			kinematics::MoveMapOP cdrh3_map;
			cdrh3_map = new kinematics::MoveMap();
			cdrh3_map->clear();
			cdrh3_map->set_chi( true );
			cdrh3_map->set_bb( false );
			for( Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop(); ii++ )
				cdrh3_map->set_bb( ii, true );
			cdrh3_map->set_jump( 1, false );

			// setup monte_carlo
			Real temp( 2.0);
			protocols::moves::MonteCarloOP mc, outer_mc;
			mc = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
			outer_mc = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
			Size buffer( (is_camelid_ && antibody_in_.is_extended()) ? 2 : 0 );
			while( !loop_found && ( total_cycles++ < cycles1) ) {
				// insert random fragments over the whole loop
				for(Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop()
							- ( buffer + (frag_size - 1 ) ); ii++ ) {
					ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
					cfm->set_check_ss( false );
					cfm->enable_end_bias_check( false );
					cfm->define_start_window( ii );
					cfm->apply( pose_in );
				}
				if( total_cycles == 1 )
					mc->reset( pose_in );

				Size local_h3_attempts(0);
				for ( Size c2 = 1; c2 <= cycles2; ++c2 ) {
					// apply a random fragment
					ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
					cfm->set_check_ss( false );
					cfm->enable_end_bias_check( false );
					cfm->apply( pose_in );

					bool H3_found_current(false);
					if( current_loop_is_H3_ && H3_filter_ &&
							( local_h3_attempts++ < (50 * cycles2) ) ) {
						H3_found_current = CDR_H3_filter(pose_in,
                                                         antibody_in_.get_CDR_loop("h3")->start(),
                                                         ( antibody_in_.get_CDR_loop("h3")->stop()-1 - antibody_in_.get_CDR_loop("h3")->start() ) + 1,
                                                         H3_filter_,
                                                         is_camelid_);
						if( !H3_found_ever && !H3_found_current) {
							--c2;
							mc->boltzmann( pose_in );
						}
						else if( !H3_found_ever && H3_found_current ) {
							H3_found_ever = true;
							mc->reset( pose_in );
						}
						else if( H3_found_ever && !H3_found_current ) {
							--c2;
							continue;
						}
						else if( H3_found_ever && H3_found_current )
							mc->boltzmann( pose_in );
					}
					else
						mc->boltzmann( pose_in );

					if ( (c2 > cycles2/2 && RG.uniform() * cycles2 < c2) ||
							 ( trimmed_cdr_h3.size() <= 5) ) {
						// in 2nd half of simulation, start trying to close the loop:
						CcdMoverOP ccd_moves = new CcdMover( trimmed_cdr_h3, cdrh3_map );
						protocols::moves::RepeatMoverOP ccd_cycle;
						if( trimmed_cdr_h3.size() <= 5 ) {
							ccd_cycle = new protocols::moves::RepeatMover(ccd_moves,500*trimmed_cdr_h3.size());
							ccd_cycle->apply( pose_in );
						}
						else {
							ccd_cycle = new protocols::moves::RepeatMover(ccd_moves, 10*trimmed_cdr_h3.size());
							ccd_cycle->apply( pose_in );
						}
						mc->boltzmann( pose_in );
					}
				}

				mc->recover_low( pose_in );
				CcdLoopClosureMoverOP ccd_closure = new CcdLoopClosureMover(
																            trimmed_cdr_h3, cdrh3_map );
				ccd_closure->set_tolerance( ccd_threshold );
				ccd_closure->set_ccd_cycles( 500 );
				ccd_closure->apply( pose_in );

				if( total_cycles == 1 )
					outer_mc->reset( pose_in );

				if ( ccd_closure->forward_deviation() <= ccd_threshold &&
						 ccd_closure->backward_deviation() <= ccd_threshold ) {
					// CDR-H3 filter for antibody mode
					// introduce enough diversity
					outer_mc->boltzmann( pose_in );
					if( current_loop_is_H3_ && H3_filter_ &&
							(current_h3_prob < h3_fraction) && (h3_attempts++<50) )
						if( !CDR_H3_filter(pose_in, 
                                           antibody_in_.get_CDR_loop("h3")->start(),
                                           ( antibody_in_.get_CDR_loop("h3")->stop()-1 - antibody_in_.get_CDR_loop("h3")->start() ) + 1,
                                           H3_filter_,
                                           is_camelid_) 
                           )
                        {
                           continue; 
                        }
					loop_found = true;
				}
				else if( H3_filter_ )
					h3_attempts++;
			}
			outer_mc->recover_low( pose_in );

			// Restoring Fold Tree
			pose_in.fold_tree( old_fold_tree );

			TR <<  "H3M Finished Fragments based centroid CDR H3 loop building"
				 << std::endl;

			return;
		} // scored_frag_close




		///////////////////////////////////////////////////////////////////////////
		/// @begin loop_centroid_relax
		///
		/// @brief actually relaxes the region specified
		///
		/// @detailed This is all done in low resolution. Intention was to give
		///           camelid CDR H1 a larger perturbation.
		///
		/// @param[in] pose, loop begin position, loop end position
		///
		/// @global_read none
		///
		/// @global_write none
		///
		/// @remarks
		///
		/// @references
		///
		/// @authors Aroop 05/07/2010
		///
		/// @last_modified 05/07/2010
		///////////////////////////////////////////////////////////////////////////
		void CDRH3Modeler2::loop_centroid_relax(
			pose::Pose & pose_in,
			Size const loop_begin,
			Size const loop_end )
		{
			using namespace protocols;
			using namespace protocols::simple_moves;
			using namespace protocols::loops;
			using namespace protocols::moves;
			using namespace pack;
			using namespace pack::task;
			using namespace pack::task::operation;
            using loop_closure::ccd::CcdMover;
            using loop_closure::ccd::CcdMoverOP;

			TR << "H3M Centroid Relaxing Loop" << std::endl;

			// storing starting fold tree
			kinematics::FoldTree tree_in( pose_in.fold_tree() );

			//setting MoveMap
			kinematics::MoveMapOP loop_map;
			loop_map = new kinematics::MoveMap();
			loop_map->clear();
			loop_map->set_chi( false );
			loop_map->set_bb( false );
			utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
			for( Size ii = loop_begin; ii <= loop_end; ii++ )
				allow_bb_move[ ii ] = true;
			loop_map->set_bb( allow_bb_move );
			loop_map->set_jump( 1, false );


			Size loop_size = ( loop_end - loop_begin ) + 1;
			Size cutpoint = loop_begin + Size(loop_size/2);

			loops::Loop one_loop( loop_begin, loop_end,	cutpoint,	0, false );
			simple_one_loop_fold_tree( pose_in, one_loop );

			// set cutpoint variants for correct chainbreak scoring
			if( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
				if( !pose_in.residue( cutpoint ).has_variant_type(
						chemical::CUTPOINT_LOWER))
					core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, cutpoint );
				if( !pose_in.residue( cutpoint + 1 ).has_variant_type(
						chemical::CUTPOINT_UPPER ) )
					core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, cutpoint + 1 );
			}



			Real min_tolerance = 0.001;
			if( benchmark_ ) min_tolerance = 1.0;
			std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
			bool nb_list = true;
			MinMoverOP loop_min_mover = new MinMover( loop_map, lowres_scorefxn_, min_type, min_tolerance, nb_list );

			// more params
			Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
			Size inner_cycles( loop_size );
			Size outer_cycles( 1 );
			if( antibody_refine_ || refine_input_loop_ )
				outer_cycles = 5;
			if( antibody_refine_ && snug_fit_ )
				outer_cycles = 2;
			if( benchmark_ ) {
				n_small_moves = 1;
				inner_cycles = 1;
				outer_cycles = 1;
			}

			Real high_move_temp = 2.00;
			// minimize amplitude of moves if correct parameter is set
			BackboneMoverOP small_mover = new SmallMover( loop_map, high_move_temp, n_small_moves );
			BackboneMoverOP shear_mover = new ShearMover( loop_map, high_move_temp, n_small_moves );
			small_mover->angle_max( 'H', 2.0 );
			small_mover->angle_max( 'E', 5.0 );
			small_mover->angle_max( 'L', 6.0 );

			shear_mover->angle_max( 'H', 2.0 );
			shear_mover->angle_max( 'E', 5.0 );
			shear_mover->angle_max( 'L', 6.0 );

			CcdMoverOP ccd_moves = new CcdMover( one_loop, loop_map );
			RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);

			SequenceMoverOP wiggle_cdr_h3( new SequenceMover() );
			wiggle_cdr_h3->add_mover( small_mover );
			wiggle_cdr_h3->add_mover( shear_mover );
			wiggle_cdr_h3->add_mover( ccd_cycle );


			loop_min_mover->apply( pose_in );

			Real const init_temp( 2.0 );
			Real const last_temp( 0.5 );
			Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
			Real temperature = init_temp;

			MonteCarloOP mc;
			mc = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temperature );
			mc->reset( pose_in ); // monte carlo reset

			// outer cycle
			for(Size i = 1; i <= outer_cycles; i++) {
				mc->recover_low( pose_in );

				// inner cycle
				for ( Size j = 1; j <= inner_cycles; j++ ) {
					temperature *= gamma;
					mc->set_temperature( temperature );
					wiggle_cdr_h3->apply( pose_in );
					loop_min_mover->apply( pose_in );

					mc->boltzmann( pose_in );

				} // inner cycles
			} // outer cycles
			mc->recover_low( pose_in );

			// minimize
			if( !benchmark_ )
				loop_min_mover->apply( pose_in );

			// Restoring pose stuff
			pose_in.fold_tree( tree_in ); // Tree

			TR << "H3M Finished Centroid Relaxing Loop" << std::endl;

			return;
		} // loop_centroid_relax





} // namespace antibody2
}  // namespace protocols
