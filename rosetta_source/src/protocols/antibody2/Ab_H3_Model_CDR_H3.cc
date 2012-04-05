// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c)University of Washington UW TechTransfer, email:license@u.washington.edu.

/// @file protocols/antibody2/Ab_H3_Model_CDR_H3.cc
/// @brief models CDR H3 loop using loop modeling
/// @detailed
///// @author Jianqing Xu ( xubest@gmail.com )
//


#include <protocols/antibody2/Ab_H3_Model_CDR_H3.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>

#include <core/id/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>

//Auto Headers

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>

#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/Ab_H3_perturb_ccd_build.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/Ab_H3_cter_insert_mover.hh>
#include <protocols/moves/PyMolMover.hh>



static basic::Tracer TR("protocols.antibody2.Ab_H3_Model_CDR_H3");
using namespace core;

namespace protocols {
namespace antibody2 {

Ab_H3_Model_CDR_H3::Ab_H3_Model_CDR_H3() : Mover()
{

}

Ab_H3_Model_CDR_H3::Ab_H3_Model_CDR_H3(bool camelid, bool benchmark, AntibodyInfoOP antibody_info) : Mover()
{
	user_defined_ = true;
	init(camelid, benchmark, antibody_info );
}



void Ab_H3_Model_CDR_H3::init(bool camelid, bool benchmark, AntibodyInfoOP antibody_info)
{
	Mover::type( "Ab_H3_Model_CDR_H3" );


	set_default();
    
	if ( user_defined_ ) {
		set_camelid( camelid );
		benchmark_ = benchmark;
        ab_info_= antibody_info;
	}

    if( is_camelid_ && !ab_info_->is_extended() && !ab_info_->is_kinked() ){
        c_ter_stem_ = 0;
    }
    
    
	lowres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
	// adding constraints
	lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );

	highres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function("standard", "score12" );
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
	// adding constraints
	highres_scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );

    
	setup_objects();
    
    
}



    
Ab_H3_Model_CDR_H3::~Ab_H3_Model_CDR_H3() {}
    
    
void Ab_H3_Model_CDR_H3::set_default()
{
    do_cter_insert_ = true;
	benchmark_ = false;
	cen_cst_ = 10.0;
	high_cst_ = 100.0; // if changed here, please change at the end of AntibodyModeler as well
	current_loop_is_H3_ = true;
	antibody_refine_ = true;
	snug_fit_ = true;
	loops_flag_ = true;
	dle_flag_ = true;
	refine_input_loop_ = true;
	is_camelid_ = false;
    use_pymol_diy_ =false;
    c_ter_stem_ = 3;
    max_cycle_ = 20;
    
    // size of loop above which 9mer frags are used
	cutoff_9_ = 16; // default 16
    
    // size of loop above which 3mer frags are used
	cutoff_3_ = 6; // default 6
    

    
    //TODO:
    //JQX:
    //if one decides to insert c_terminal first, it means the h3 loop has 3 less residues
    //should one change the creteria of cutoff_3_ and cutoff_9_?
    
    
	TR << "Finished Setting Defaults" << std::endl;
} 

    
    
    
void Ab_H3_Model_CDR_H3::setup_objects(){
    h3_cter_insert_mover_ = new Ab_H3_cter_insert_mover(ab_info_, is_camelid_);
    h3_perturb_ccd_build_ = new Ab_H3_perturb_ccd_build(current_loop_is_H3_,is_camelid_, ab_info_ );        
}
    
    
    
    
    
void Ab_H3_Model_CDR_H3::set_lowres_score_func( scoring::ScoreFunctionOP lowres_scorefxn ) {
    lowres_scorefxn_ = lowres_scorefxn;
} 

void Ab_H3_Model_CDR_H3::set_highres_score_func(scoring::ScoreFunctionOP highres_scorefxn) {
    highres_scorefxn_ = highres_scorefxn;
}

void Ab_H3_Model_CDR_H3::turn_off_H3_filter(){
    h3_perturb_ccd_build_->turn_off_H3_filter();
}    

void Ab_H3_Model_CDR_H3::turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol){        
    use_pymol_diy_ = true;
    pymol_ = pymol;
    h3_cter_insert_mover_->turn_on_and_pass_the_pymol(pymol);
    h3_perturb_ccd_build_->turn_on_and_pass_the_pymol(pymol);
}
    

void Ab_H3_Model_CDR_H3::apply( pose::Pose & pose_in )
{


    TR << "Applying CDR H3 modeler" << std::endl;

    using namespace core::pose;
    using namespace core::scoring;
    using namespace protocols::moves;

    start_pose_ = pose_in;

    pose::Pose start_pose = pose_in;



    Size framework_loop_begin( ab_info_->get_CDR_loop("h3")->start() );
    Size framework_loop_end  ( ab_info_->get_CDR_loop("h3")->stop()  );
    Size cutpoint = framework_loop_begin + 1;
    Size framework_loop_size = (framework_loop_end - framework_loop_begin) + 1;

    loops::Loop cdr_h3( framework_loop_begin, framework_loop_end, cutpoint, 0, true );
    loops::Loop trimmed_cdr_h3(framework_loop_begin, framework_loop_end - c_ter_stem_, cutpoint, 0, true );
    loops::Loop input_loop;
    
    if (do_cter_insert_){
        //JQX: the h3 loop removing the cterminal 3 residues
        input_loop = trimmed_cdr_h3;
    }
    else{
        //JQX: the original h3 loop
        input_loop = cdr_h3;
    }
        

    simple_one_loop_fold_tree( pose_in, cdr_h3 );

    // switching to centroid mode
    simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
    simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

    
    
    // Building centroid mode loop
    to_centroid.apply( pose_in );
//    pymol_->apply(pose_in);

    
    // some initialization before you do h3 loop modeling
//    my_LoopMover my_loop_mover;
//    my_loop_mover.set_extended_torsions( pose_in, cdr_h3 );
        set_extended_torsions( pose_in, cdr_h3 );
       pose_in.dump_pdb("extend_centroid.pdb");
        //JQX:  this function is in loops_main.cc file
        //      firstly, idealize the loop (indealize bonds as well)
        //      phi(-150),  all the residue, except the first one
        //      psi(150),   all the residue, except the last one
        //      omega(180), all the residue, except the first & last one
        //JQX:  in R2: the function is called "insert_init_frag", which is 
        //      in the file "jumping_util.cc". All the phi, psi, omega are 
        //      assigned to all the residues. "L" secondary structure is 
        //      also assinged. The bonds are idealized using 
        //      framework_pose.insert_ideal_bonds(begin-1, end)
    
    
    
    /*  JQX: the following code is probably not ncessary*/
    
    Size unaligned_cdr_loop_begin(0), unaligned_cdr_loop_end(0);
    std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
    core::import_pose::pose_from_pdb( hfr_pose_, path+"hfr.pdb" );
    std::string cdr_name = "h3";
    AntibodyInfoOP hfr_info =  new AntibodyInfo ( hfr_pose_, cdr_name );
    unaligned_cdr_loop_begin = hfr_info->current_start;
    unaligned_cdr_loop_end   = hfr_info->current_end;
    
    if(framework_loop_size > 4){  //JQX: add this if statement to match R2_antibody
        pose_in.set_psi  (framework_loop_begin - 1, hfr_pose_.psi( unaligned_cdr_loop_begin - 1 )   );
        pose_in.set_omega(framework_loop_begin - 1, hfr_pose_.omega( unaligned_cdr_loop_begin - 1 ) );
    }

        pose_in.dump_pdb("after_copying_nter.pdb");

    
    antibody2::AntibodyInfoOP starting_antibody;
    starting_antibody = ab_info_;
    bool closed_cutpoints( false );
    
    h3_perturb_ccd_build_->pass_the_loop(input_loop);
    
    Size cycle ( 1 );
    while( !closed_cutpoints && cycle < max_cycle_) {
        ab_info_ = starting_antibody;
        if (do_cter_insert_){
            if( framework_loop_size > 6 ){ //JQX: replace 5 by 6 to match R2_antibody
                h3_cter_insert_mover_->apply(pose_in);
            }
            else{
                utility_exit_with_message("Loop Size is Less than 6");
            }
        }
        

        
        pose_in.dump_pdb("after_c_insert.pdb");


        h3_perturb_ccd_build_->apply(pose_in);

        
        closed_cutpoints = cutpoints_separation( pose_in, ab_info_ );
        ++cycle;
    } // while( ( cut_separation > 1.9 )
    
    TR <<  "Finished Modeling Centroid CDR H3 loop" << std::endl;
    

    
    
    //#############################  //JQX: this should not be here
    if( is_camelid_ ){
        loop_centroid_relax( pose_in, ab_info_->get_CDR_loop("h1")->start(), ab_info_->get_CDR_loop("h1")->stop() );
    }
    //#############################
    
    
    
    
    to_full_atom.apply( pose_in );

    utility::vector1<bool> allow_chi_copy( pose_in.total_residue(), true );
    for( Size ii = ab_info_->get_CDR_loop("h3")->start();
						 ii <= ( ab_info_->get_CDR_loop("h3")->stop() ); ii++ )
					allow_chi_copy[ii] = false;
    //recover sidechains from starting structures
    protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose_, allow_chi_copy );
    recover_sidechains.apply( pose_in );
    
    

    // Packer
    protocols::simple_moves::PackRotamersMoverOP packer;
    packer = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
    packer->task_factory(tf_);
    packer->apply( pose_in );



    pose_in.dump_pdb("finish_apply_CDR_H3_modeler.pdb");
    TR << "Finished applying CDR H3 modeler" << std::endl;
    

    return;
} // Ab_H3_Model_CDR_H3::apply()

    
    
    
std::string Ab_H3_Model_CDR_H3::get_name() const {
	return "Ab_H3_Model_CDR_H3";
}


    

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
		void Ab_H3_Model_CDR_H3::loop_centroid_relax(
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
				if( !pose_in.residue( cutpoint ).has_variant_type(chemical::CUTPOINT_LOWER))
					core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, cutpoint );
				if( !pose_in.residue( cutpoint + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
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
			if( antibody_refine_ || refine_input_loop_ ){
				outer_cycles = 5;
            }
			if( antibody_refine_ && snug_fit_ ){
				outer_cycles = 2;
            }
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

			TR << "Finished Centroid Relaxing Loop" << std::endl;

			return;
		} // loop_centroid_relax





} // namespace antibody2
}  // namespace protocols



