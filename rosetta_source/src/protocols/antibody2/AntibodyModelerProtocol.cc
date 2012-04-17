// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/AntibodyModelerProtocol.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>


#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>

#include <protocols/jd2/ScoreMap.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyModelerProtocol.hh>
#include <protocols/antibody2/ModelCDRH3.hh>

#include <protocols/antibody2/RefineCDRH3HighRes.hh>
#include <protocols/antibody2/RefineBetaBarrel.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::fmt;




using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.AntibodyModelerProtocol");
using namespace core;

namespace protocols {
namespace antibody2 {

// default constructor
AntibodyModelerProtocol::AntibodyModelerProtocol() : Mover() {
	user_defined_ = false;
	init();
}

// default destructor
AntibodyModelerProtocol::~AntibodyModelerProtocol() {}

//clone
protocols::moves::MoverOP
AntibodyModelerProtocol::clone() const {
	return( new AntibodyModelerProtocol() );
}

    
    
void AntibodyModelerProtocol::init() 
{
	Mover::type( "AntibodyModelerProtocol" );
    
	set_default();
	init_from_options();
	setup_objects();
}

    
    
void AntibodyModelerProtocol::set_default()
{
	TR <<  "Setting Up defaults.........." << std::endl;
    model_h3_  = true;
    extreme_repacking_ =  true;
	snugfit_   = true;
    refine_h3_ = true;
    
	benchmark_ = false;
	camelid_   = false;
	camelid_constraints_ = false;
    cst_weight_ = 0.0;
    cen_cst_ = 10.0;
    high_cst_ = 100.0; // if changed here, please change at the end of AntibodyModeler as well
    H3_filter_ = true;
    cter_insert_ = true;
    use_pymol_diy_ = true;
    LH_repulsive_ramp_ = true;
}

    
    
void AntibodyModelerProtocol::register_options()
{
	using namespace basic::options;

    option.add_relevant( OptionKeys::antibody::model_h3 );
	option.add_relevant( OptionKeys::antibody::snugfit );
    option.add_relevant( OptionKeys::run::benchmark );
	option.add_relevant( OptionKeys::antibody::camelid );
	option.add_relevant( OptionKeys::antibody::camelid_constraints );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::in::file::native );
    //option.add_relevant( OptionKeys::antibody::H3_filter );
    //option.add_relevant( OptionKeys::antibody::cter_insert );
}

    
    
    
    
void AntibodyModelerProtocol::init_from_options() 
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
    
	TR <<  "Start Reading and Setting Options ..." << std::endl;
    
    if ( option[OptionKeys::antibody::model_h3].user() ){
        set_ModelH3(option[OptionKeys::antibody::model_h3]() );
    }
	if ( option[ OptionKeys::antibody::snugfit ].user() ){
        set_SnugFit( option[ OptionKeys::antibody::snugfit ]() );
    }
	if ( option[ OptionKeys::antibody::camelid ].user() ){
        set_camelid( option[ OptionKeys::antibody::camelid ]() );
    }
	if ( option[ OptionKeys::antibody::camelid_constraints ].user() ){
        set_camelid_constraints( option[ OptionKeys::antibody::camelid_constraints ]() );
    }
	if ( option[ OptionKeys::run::benchmark ].user() ){
        set_BenchMark( option[ OptionKeys::run::benchmark ]() );
    }
    if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		set_cst_weight( option[ OptionKeys::constraints::cst_weight ]() );
	}
    //if ( option[ OptionKeys::antibody::H3_filter ].user() ) {
	//	set_H3Filter( option[ OptionKeys::antibody::H3_filter ]() );
	//}
    //if ( option[ OptionKeys::antibody::cter_insert ].user() ) {
	//	set_CterInsert( option[ OptionKeys::antibody::cter_insert ]() );
	//}

	//set native pose if asked for
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, option[ OptionKeys::in::file::native ]() );
		set_native_pose( native_pose );
	}
	else{
		set_native_pose(NULL);
	}
    

    
	if( camelid_ ) {
		snugfit_ = false;
	}
    
    TR <<  "Finish Reading and Setting Options !!!" << std::endl;
}


    
void
AntibodyModelerProtocol::setup_objects() {
    
	sync_objects_with_flags();
    
    // setup all the fold trees
    
    // setup all the move maps
    
    // setup all the task factories
    tf_ = new pack::task::TaskFactory;

    // setup all the scoring functions
    pack_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("standard" );
    dock_scorefxn_highres_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
        dock_scorefxn_highres_->set_weight( core::scoring::chainbreak, 1.0 );
        dock_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
    loop_scorefxn_centroid_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
        loop_scorefxn_centroid_->set_weight( scoring::chainbreak, 10./3. );
        loop_scorefxn_centroid_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
    loop_scorefxn_highres_ = scoring::ScoreFunctionFactory::create_score_function("standard", "score12" );
        loop_scorefxn_highres_->set_weight( scoring::chainbreak, 1.0 );
        loop_scorefxn_highres_->set_weight( scoring::overlap_chainbreak, 10./3. );
        loop_scorefxn_highres_->set_weight( scoring::atom_pair_constraint, high_cst_ );
    
    // miscellaneous
    pymol_ = new protocols::moves::PyMolMover;
    pymol_->keep_history(true);
}
    
void AntibodyModelerProtocol::sync_objects_with_flags() 
{
	using namespace protocols::moves;
	flags_and_objects_are_in_sync_ = true;
	first_apply_with_current_setup_ = true;
}


std::string AntibodyModelerProtocol::get_name() const 
{        
    return "AntibodyModelerProtocol";
}

    
    

void AntibodyModelerProtocol::finalize_setup( pose::Pose & frame_pose ) 
{
	TR<<"AAAAAAAA     cst_weight: "<<cst_weight_<<std::endl;
	if(  cst_weight_ != 0.00  ) {
		simple_moves::ConstraintSetMoverOP cdr_constraint = new simple_moves::ConstraintSetMover();
		cdr_constraint->apply( frame_pose );
	}

	// check for native and input pose
	if ( !get_input_pose() ) {
		pose::PoseOP input_pose = new pose::Pose(frame_pose); 
		set_input_pose( input_pose );   // JQX: pass the input_pose to the mover.input_pose_
	}


	pose::PoseOP native_pose;
	if ( !get_native_pose() ) {
		TR << "Danger Will Robinson! Native is an impostor!" << std::endl;
        TR << "   'native_pose' is just a copy of the 'input_pose'    " << std::endl;
        TR << "    since you didn't sepcifiy the native pdb name"<<std::endl;
		native_pose = new pose::Pose(frame_pose);
	} else {
		native_pose = new pose::Pose( *get_native_pose() );
	}

	pose::set_ss_from_phipsi( *native_pose ); // JQX: this is the secondary structure from the native pose

	set_native_pose( native_pose ); // pass the native pose to the mover.native_pose_

    ab_info_ = new AntibodyInfo(frame_pose,camelid_);
    TR<<*ab_info_<<std::endl;
    
    model_cdrh3_        = new ModelCDRH3( ab_info_, loop_scorefxn_centroid_, loop_scorefxn_highres_);
    refine_beta_barrel_ = new RefineBetaBarrel(ab_info_, dock_scorefxn_highres_, pack_scorefxn_);


    
    
    TR << "Utility: Setting Up Packer Task" << std::endl;
    using namespace pack::task;
    using namespace pack::task::operation;
    tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
    tf_->push_back( new InitializeFromCommandline );
    tf_->push_back( new IncludeCurrent );
    tf_->push_back( new RestrictToRepacking );
    tf_->push_back( new NoRepackDisulfides );
    
    // incorporating Ian's UnboundRotamer operation.
    // note that nothing happens if unboundrot option is inactive!
    pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new pack::rotamer_set::UnboundRotamersOperation();
    unboundrot->initialize_from_command_line();
    operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
    tf_->push_back( unboundrot_operation );
    // adds scoring bonuses for the "unbound" rotamers, if any
    core::pack::dunbrack::load_unboundrot( frame_pose );
        
    TR << "Utility: Done: Setting Up Packer Task" << std::endl;

}



//APPLY
void AntibodyModelerProtocol::apply( pose::Pose & frame_pose ) {

    using namespace chemical;
    using namespace id;
    using namespace scoring;
    using namespace core::scoring::constraints;
    using namespace protocols::moves;


    // the default inital secstruct is all "L" loop!
    start_pose_ = frame_pose;
    


    if ( !flags_and_objects_are_in_sync_ ){ 
       sync_objects_with_flags(); 
    }
    
    if ( first_apply_with_current_setup_ ){ 
        finalize_setup(frame_pose);  
        first_apply_with_current_setup_=false; 
    }



	basic::prof_reset();
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	pose::set_ss_from_phipsi( frame_pose );
    

	// display constraints and return
	if( camelid_constraints_ ) {
		display_constraint_residues( frame_pose );
		return;
	}




    // Step 1: model the cdr h3
    // JQX notes: pay attention to the way it treats the stems when extending the loop
    if(use_pymol_diy_) pymol_->apply(frame_pose);
    if(model_h3_){        
        if(cter_insert_ ==false) { model_cdrh3_->turn_off_cter_insert(); }
        if(H3_filter_   ==false) { model_cdrh3_->turn_off_H3_filter();   }
        model_cdrh3_->set_task_factory(tf_);
        if(use_pymol_diy_) model_cdrh3_->turn_on_and_pass_the_pymol(pymol_);
        model_cdrh3_->apply( frame_pose );
    }
    
    
    // Step 2: packing the CDRs
    if(extreme_repacking_) { relax_cdrs( frame_pose );    }
    if(use_pymol_diy_) pymol_->apply(frame_pose);
    
    
    
	// Step 3: SnugFit: relieve the clashes between L-H
	if ( snugfit_ ) { 
        if (LH_repulsive_ramp_) {refine_beta_barrel_-> turn_off_repulsive_ramp();};
        refine_beta_barrel_->apply(frame_pose);
	}


    
	// Step 4: Full Atom Relax 
    if(refine_h3_){
        //$$$$$$$$$$$$$$$$$$$$$$$$
        cdr_highres_refine_ = new RefineCDRH3HighRes(ab_info_, "h3", loop_scorefxn_highres_); 
        cdr_highres_refine_ -> set_task_factory(tf_);
        cdr_highres_refine_ -> pass_start_pose(start_pose_);
        if(use_pymol_diy_) cdr_highres_refine_ -> turn_on_and_pass_the_pymol(pymol_);
        cdr_highres_refine_ -> apply(frame_pose);
        frame_pose.dump_pdb("finish_h3_refinement.pdb");


        //$$$$$$$$$$$$$$$$$$$$$$$$$
        if( !benchmark_ ) 
        {
            Size repack_cycles(1);
            if( antibody_refine_ && !snugfit_ ){repack_cycles = 3;}
            protocols::simple_moves::PackRotamersMoverOP packer;
            packer = new protocols::simple_moves::PackRotamersMover( loop_scorefxn_highres_ );
            packer->task_factory(tf_);
            packer->nloop( repack_cycles );
            packer->apply( frame_pose );
        }
        if(use_pymol_diy_) pymol_->apply(frame_pose);

        return;
        // Minimize CDR H2 loop if this is a camelid
    
        if( camelid_ ) {
            //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            RefineCDRH3HighRes cdr_highres_refine( ab_info_); // because of h2
            //cdr_highres_refine.turn_off_h3_default();
            cdr_highres_refine.turn_off_h3_filter();
            cdr_highres_refine.apply(frame_pose);
            //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
            //JQX: remove the duplicated code, camelid H2 will be automatically taken care of
            //     see the code in RefineCDRH3HighRes file
        }
    }
    
    

    
    
    // align pose to native pose
    pose::Pose native_pose = *get_native_pose();
    antibody2::AntibodyInfo native_ab( native_pose, camelid_ );
    //		ab_info_.align_to_native( pose, native_ab, native_pose );
    
    
	// Step 5: Store the homolgy models
    
	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( frame_pose, true );

	// Define CDR H3 loop
	Size frag_size   = (ab_info_->get_CDR_loop("h3")->stop()  - ab_info_->get_CDR_loop("h3")->start()) + 3;
	Size cutpoint    =  ab_info_->get_CDR_loop("h3")->start() + int( frag_size / 2 );
	loops::Loop cdr_h3( ab_info_->get_CDR_loop("h3")->start(), ab_info_->get_CDR_loop("h3")->stop(), cutpoint, 0, false );

	// Fold Tree
	antibody2::simple_one_loop_fold_tree( frame_pose, cdr_h3 );

	// Redefining CDR H3 cutpoint variants
	loops::add_single_cutpoint_variant( frame_pose, cdr_h3 );

    // reset the score weight here, to make it match to Aroop's R3 antibody code
    loop_scorefxn_highres_->set_weight( core::scoring::chainbreak, 1.0 );
    loop_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
    loop_scorefxn_highres_->set_weight( core::scoring::atom_pair_constraint, 1.00 );
    
    // add scores to map for outputting constraint score
	( *loop_scorefxn_highres_ )( frame_pose );
	Real constraint_score = frame_pose.energies().total_energies()[ core::scoring::atom_pair_constraint ];

	// removing constraint score
	loop_scorefxn_highres_->set_weight( core::scoring::atom_pair_constraint, 0.00 );
	// add scores to map for output
	( *loop_scorefxn_highres_ )( frame_pose );

	job->add_string_real_pair("AA_H3", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("h3") ));
	job->add_string_real_pair("AB_H2", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("h2") ));
	job->add_string_real_pair("AC_H1", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("h1") ));
	if( !camelid_ ) {
		job->add_string_real_pair("AC_L3", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("l3")));
		job->add_string_real_pair("AD_L2", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("l2") ));
		job->add_string_real_pair("AE_L1", global_loop_rmsd( frame_pose, *get_native_pose(), ab_info_->get_CDR_loop("l1") ));
	}
	job->add_string_real_pair("AF_constraint", constraint_score);

	set_last_move_status( protocols::moves::MS_SUCCESS );   

	basic::prof_show();


}// end apply





    




///////////////////////////////////////////////////////////////////////////
/// @begin relax_cdrs     //JQX: packing, minimization, and mintrial
///
/// @brief relaxes all cdrs simultaneously
///
/// @detailed based on the all_cdrs loop definiton, minimizes only those
///           regions. A standard dfpmin is utilized with score12 and chain
///           -break and chain-overlap set. The allow_bb/chi arrays are
///           changed accordingly but then are reset to their initial
///           states before exiting the routine. Similarly the fold tree
///           and jump movements are restored to their initial states
///
///
/// @authors Aroop 02/15/2010
///
/// @last_modified 02/15/2010
///////////////////////////////////////////////////////////////////////////
void AntibodyModelerProtocol::relax_cdrs( core::pose::Pose & pose )
{
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::moves;
	// Storing initial fold tree
	kinematics::FoldTree const input_tree( pose.fold_tree() );

	// changing to all cdr fold tree
	ab_info_->all_cdr_fold_tree( pose );

	// adding cutpoint variants for chainbreak score computation
	loops::add_cutpoint_variants( pose );
    
    // score functions
    loop_scorefxn_highres_->set_weight( core::scoring::chainbreak, 10. / 3. );
    loop_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10. / 3. );
    ( *loop_scorefxn_highres_ )( pose );
    
	//setting MoveMap
	kinematics::MoveMapOP allcdr_map;
	allcdr_map = new kinematics::MoveMap();
	allcdr_map->clear();
	allcdr_map->set_chi( false );
	allcdr_map->set_bb( false );
	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
    utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	select_loop_residues( pose, ab_info_->all_cdr_loops_, false /*include_neighbors*/, bb_is_flexible );
	allcdr_map->set_bb( bb_is_flexible );
	select_loop_residues( pose, ab_info_->all_cdr_loops_, true /*include_neighbors*/, sc_is_flexible );
	allcdr_map->set_chi( sc_is_flexible );
	for( Size ii = 1; ii <= ab_info_->all_cdr_loops_.num_loop(); ii++ ){
		allcdr_map->set_jump( ii, false );
    }



	Real min_tolerance = 0.1;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
	bool nb_list = true;
    simple_moves::MinMoverOP all_cdr_min_moves = new simple_moves::MinMover( allcdr_map,
                                                                    loop_scorefxn_highres_, min_type, min_tolerance, nb_list );
    all_cdr_min_moves->apply( pose );

    if( !benchmark_ ) {
        simple_moves::PackRotamersMoverOP repack=new simple_moves::PackRotamersMover( loop_scorefxn_highres_ );
        ( *loop_scorefxn_highres_ )( pose );
        tf_->push_back( new RestrictToInterface( sc_is_flexible ) );
        repack->task_factory( tf_ );
        repack->apply( pose );

        simple_moves::RotamerTrialsMinMoverOP rtmin = new simple_moves::RotamerTrialsMinMover( loop_scorefxn_highres_, tf_ );
        rtmin->apply( pose );
    }

    // Restoring pose fold tree
    pose.fold_tree( input_tree );
} // relax_cdrs

    
    
    
    
    







    
    
    
    
    
    
    

void AntibodyModelerProtocol::display_constraint_residues( core::pose::Pose & pose ) 
{

    // Detecting di-sulfide bond

    Size H1_Cys(0), H3_Cys(0);

    if(      pose.residue( pose.pdb_info()->pdb2pose('H',32 ) ).name3() == "CYS" ){
        H1_Cys = pose.pdb_info()->pdb2pose( 'H', 32 );
    }
    else if( pose.residue( pose.pdb_info()->pdb2pose('H',33 ) ).name3() == "CYS" ){
        H1_Cys = pose.pdb_info()->pdb2pose( 'H', 33 );
    }

    for( Size ii = ab_info_->get_CDR_loop("h3")->start(); ii <= ab_info_->get_CDR_loop("h3")->stop(); ii++ ){
        if( pose.residue(ii).name3() == "CYS" ) {
            H3_Cys = ii;
        }
    }

    if( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) ){
        TR << "CONSTRAINTS: "<< "AtomPair CA " << H1_Cys << " CA " << H3_Cys
           << " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;
    }

    // Specifying extended kink

    Size hfr_46(0), h3_closest(0);
    hfr_46 = pose.pdb_info()->pdb2pose( 'H', 46 );
    if( ab_info_->is_extended() ) h3_closest = ab_info_->get_CDR_loop("h3")->stop() - 5;
    if( h3_closest != 0 ) {
        TR << "CONSTRAINTS: " << "AtomPair CA " << hfr_46 << " CA " << h3_closest
           << " BOUNDED 6.5 9.1 0.7 DISTANCE; mean 8.0 sd 0.7" << std::endl;
    }

    return;
} // display_constraint_residues

    
    
    
    
    
    
    
    
    
    
/// @details  Show the complete setup of the antibody modeler protocol
void AntibodyModelerProtocol::show( std::ostream & out ) {
    if ( !flags_and_objects_are_in_sync_ ){
        sync_objects_with_flags();
    }
    out << *this;
}
    
std::ostream & operator<<(std::ostream& out, const AntibodyModelerProtocol & ab_m_2 ){
    using namespace ObjexxFCL::fmt;
        
    // All output will be 80 characters - 80 is a nice number, don't you think?
    std::string line_marker = "///";
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    out << line_marker << A( 47, "Rosetta 3 Antibody Modeler" ) << space( 27 ) << line_marker << std::endl;
    out << line_marker << space( 74 ) << line_marker << std::endl;

    // Display the state of the antibody modeler protocol that will be used
    out << line_marker << "  camelid                : " << ab_m_2.camelid_     << std::endl;
    out << line_marker << "  model_h3               : " << ab_m_2.model_h3_    << std::endl;
    out << line_marker << "     cter_insert         : " << ab_m_2.cter_insert_ << std::endl;
    out << line_marker << "     H3_filter           : " << ab_m_2.H3_filter_   << std::endl;
    out << line_marker << "  extreme_repacking_     : " << ab_m_2.extreme_repacking_   << std::endl;
    out << line_marker << "  snugfit                : " << ab_m_2.snugfit_     << std::endl;
    out << line_marker << "     LH_repulsive_ramp   : " << ab_m_2.LH_repulsive_ramp_ << std::endl;
    out << line_marker << "  refine_h3              : " << ab_m_2.refine_h3_     << std::endl;

    // Close the box I have drawn
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    return out;
}
    


    
    
    
        
    

} // end antibody2
} // end protocols

