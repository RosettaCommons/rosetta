// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/AnchoredGraftMover.cc
/// @brief   Method definitions for AnchoredGraftMover
/// @author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author  Steven Lewis (smlewi@gmail.com)

//Unit Headers
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/GraftMoverBase.hh>
#include <protocols/moves/Mover.hh>

//Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

//Protocol Headers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <protocols/moves/MonteCarlo.hh>

//Option Headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Basic and Utility
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <map>
#include <utility/PyAssert.hh>

static basic::Tracer TR("protocols.grafting.AnchoredGraftMover");
namespace protocols {
namespace grafting {
    using namespace core::scoring;
    using namespace core::pose;
    using namespace core::pack::task;
	using namespace basic::options;
    using namespace protocols::loops;

AnchoredGraftMover::AnchoredGraftMover(const Size start, Size const end)
	:GraftMoverBase(start, end, "AnchoredGraftMover")
{
	//moves::Mover("GraftMover"),
    set_defaults();
}


AnchoredGraftMover::~AnchoredGraftMover() {}

void
AnchoredGraftMover::set_defaults(){

	set_skip_sampling(false);
	set_use_smooth_centroid_settings(false);
	set_mintype(option [OptionKeys::relax::min_type]());

	//Arm types - default is single_loop_single_arm
	set_use_single_loop_double_CCD_arms(false);
	set_use_double_loop_double_CCD_arms(false);
	set_use_double_loop_quad_CCD_arms(false);

	set_test_control_mode(false);
	set_use_default_movemap_from_flexibility(true);
    set_scaffold_flexibility(2, 2);
    set_insert_flexibility(0, 0);
    set_cycles(300);

}

void
AnchoredGraftMover::set_scaffold_flexibility(const Size Nter_scaffold_flexibility, const Size Cter_scaffold_flexibility){
    //Insert flexibility stays at default values (0)
    Nter_scaffold_flexibility_=Nter_scaffold_flexibility;
    Cter_scaffold_flexibility_=Cter_scaffold_flexibility;
}

Size
AnchoredGraftMover::get_nterm_scaffold_flexibility(){
	return Nter_scaffold_flexibility_;
}

Size
AnchoredGraftMover::get_cterm_scaffold_flexibility(){
	return Cter_scaffold_flexibility_;
}

void
AnchoredGraftMover::set_insert_flexibility(Size Nter_insert_flexibility, Size Cter_insert_flexibility)
{
	Nter_insert_flexibility_=Nter_insert_flexibility;
	Cter_insert_flexibility_=Cter_insert_flexibility;
}

void
AnchoredGraftMover::set_cycles(Size cycles){
    cycles_=cycles;
}

void
AnchoredGraftMover::set_default_cen_scorefunction(){
	cen_scorefxn_ = new ScoreFunction();
	cen_scorefxn_->set_weight( chainbreak,  20.00);
	cen_scorefxn_->set_weight( cbeta,       1.0 );
	cen_scorefxn_->set_weight( vdw,         1.0 );
	cen_scorefxn_->set_weight( pair,        1.0 );
	cen_scorefxn_->set_weight( cenpack,     1.0 );
	cen_scorefxn_->set_weight( rama,        5.0 );
	cen_scorefxn_->set_weight( hbond_lr_bb, 1.0 );
	cen_scorefxn_->set_weight( hbond_sr_bb, 1.0 );
	cen_scorefxn_->set_weight( omega,       5.0 );
}

void
AnchoredGraftMover::set_cen_scorefunction(ScoreFunctionOP score){
	cen_scorefxn_ = score->clone();
	if (cen_scorefxn_->get_weight(chainbreak) == 0.0){
		cen_scorefxn_->set_weight(chainbreak, 20.00);
	}
}

void
AnchoredGraftMover::set_mintype(std::string mintype){
	mintype_ = mintype;
}

void
AnchoredGraftMover::set_skip_sampling(bool skip_sampling){
	skip_sampling_ = skip_sampling;
}

void
AnchoredGraftMover::set_use_smooth_centroid_settings(bool use_smooth){
	if (use_smooth){
		std::string def_score = option [OptionKeys::relax::centroid::weights]();
		cen_scorefxn_ = ScoreFunctionFactory::create_score_function(def_score);
		cen_scorefxn_->set_weight(chainbreak, 20.00);
	} else {
		set_default_cen_scorefunction();
	}
}

void
AnchoredGraftMover::set_use_single_loop_double_CCD_arms(bool single_loop_double_arm){
	single_loop_double_arm_ = single_loop_double_arm;
}

void
AnchoredGraftMover::set_use_double_loop_double_CCD_arms(bool double_loop_double_arm){
	double_loop_double_arm_=double_loop_double_arm;
}

void
AnchoredGraftMover::set_use_double_loop_quad_CCD_arms(bool double_loop_quad_arm){
	double_loop_quad_arm_=double_loop_quad_arm;
}

core::Size
AnchoredGraftMover::get_Cter_loop_end(){
	return Cter_loop_end_;
}

void
AnchoredGraftMover::set_test_control_mode(bool test_control_mode){
	test_control_mode_=test_control_mode;
}

void
AnchoredGraftMover::set_use_default_movemap_from_flexibility(bool def){
    use_default_movemap_=def;
}

void
AnchoredGraftMover::set_movemaps(const MoveMapOP scaffold_mm, const MoveMapOP insert_mm){
    
    scaffold_movemap_ = scaffold_mm;
    insert_movemap_ = insert_mm;
    use_default_movemap_=false;

}

void
AnchoredGraftMover::apply(Pose & pose){
	//No Local copy since we will be messing with scaffold_pose anyway.

    setup_movemap_and_regions(pose);
    
	//Run the insertion.
	Pose combined = insert_piece(pose);
	core::kinematics::FoldTree original_ft = combined.fold_tree();
	//Setup for the remodeling
	core::Size const insert_start(start_+1); //this will be the first residue of the insert
	core::Size const insert_end(start_+insertion_length_); //this will be the last residue of the insert


	
	///Add variants, create the loops and set the foldtree that will be used for CCD.
	using core::pose::add_variant_type_to_pose_residue;
	
    
	Loop Nter_loop;
	Loop Cter_loop;
    LoopsOP loop_set = new Loops();
    std::map< Loop, loop_closure::ccd::CcdLoopClosureMoverOP > loop_set_map; //Would not work without owning pointer.

    //Assert that none or only one boolean for the close type is set.

    if (double_loop_double_arm_){

		Nter_loop = Loop(Nter_loop_start_, Nter_loop_end_+1, Nter_loop_end_);//(LEFT LOOP)
		Cter_loop = Loop(Cter_loop_start_-1, Cter_loop_end_, Cter_loop_start_-1);//(RIGHT LOOP)
		core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_LOWER, Nter_loop.cut() );
		core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_UPPER, Nter_loop.cut()+1 );


        loop_set->add_loop(Nter_loop);
        loop_set->add_loop(Cter_loop);
		FoldTreeFromLoops ft_loop = FoldTreeFromLoops();
		ft_loop.loops(loop_set);
		ft_loop.apply(combined);

        loop_closure::ccd::CcdLoopClosureMover close_for_Nter_loop(Nter_loop, movemap_);
        loop_closure::ccd::CcdLoopClosureMover close_for_Cter_loop(Cter_loop, movemap_);
        loop_set_map[Nter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Nter_loop, movemap_);
        loop_set_map[Cter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Cter_loop, movemap_);
        
	}

	else if (double_loop_quad_arm_){
		Nter_loop = Loop(Nter_loop_start_, Nter_loop_end_+1, start_);//(LEFT LOOP)
		Cter_loop = Loop(Cter_loop_start_-1, Cter_loop_end_, end_-1);//(RIGHT LOOP)
		core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_LOWER, Nter_loop.cut() );
		core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_UPPER, Nter_loop.cut()+1 );


        loop_set->add_loop(Nter_loop);
        loop_set->add_loop(Cter_loop);
		FoldTreeFromLoops ft_loop = FoldTreeFromLoops();
		ft_loop.loops(loop_set);
		ft_loop.apply(combined);

        loop_closure::ccd::CcdLoopClosureMover close_for_Nter_loop(Nter_loop, movemap_);
        loop_closure::ccd::CcdLoopClosureMover close_for_Cter_loop(Cter_loop, movemap_);
        loop_set_map[Nter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Nter_loop, movemap_);
        loop_set_map[Cter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Cter_loop, movemap_);
	}
	else if (single_loop_double_arm_){
	    
		setup_single_loop_double_arm_remodeling_foldtree(combined, Nter_loop_start_, Cter_loop_end_);
		Cter_loop = Loop(Nter_loop_start_, Cter_loop_end_, insert_end);
        loop_set->add_loop(Cter_loop);
        
        loop_closure::ccd::CcdLoopClosureMover close_for_Cter_loop(Cter_loop, movemap_);
        loop_set_map[Cter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Cter_loop, movemap_);
        
	}
	else{   //single_loop_single_arm

		setup_single_loop_single_arm_remodeling_foldtree(combined, Nter_loop_start_, Cter_loop_end_);
		Cter_loop = Loop(Nter_loop_start_, Cter_loop_end_, Cter_loop_end_-1);
        loop_set->add_loop(Cter_loop);
        
        loop_closure::ccd::CcdLoopClosureMover close_for_Cter_loop(Cter_loop, movemap_);
        loop_set_map[Cter_loop]=new loop_closure::ccd::CcdLoopClosureMover(Cter_loop, movemap_);
        
	}
	core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_LOWER, Cter_loop.cut() );
	core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_UPPER, Cter_loop.cut()+1 );






	///////////////////////////////////////Idealize////////////////////////////////////////////////////////
	//this code also resets conformation variables: omegas to 180, newly made connections phi or psi to reasonable
	//edges of insert will be somewhat mobile inside minimization (small and CCD moves will ignore it)
	using namespace core::id;

	//Idealize the All residues in Movemap since they will all help close the loop(s).
    //combined.dump_pdb("before_idealize.pdb");

	//Order of idealization matters here.
	//Nter regions
	core::Size i(Nter_loop_start_);
	for(; i<=start_; ++i) {
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		if (movemap_->get_bb(i)){
			combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
			combined.set_omega(i, 180);
			TR << "mobile " << i << std::endl;
		}
	}
	TR << "ideal " << insert_start << std::endl;

	//Set individual torsions ON in the movemap for the start and end of the insert
	movemap_->set( TorsionID(insert_start, BB, phi_torsion), true);
	combined.set_phi(insert_start, -60);

	//insert regions
	i=start_+1;
	for(; i<=insert_end; ++i) {
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		if (movemap_->get_bb(i)){
			combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
			combined.set_omega(i, 180);
			TR << "mobile " << i << std::endl;
		}
	}

	movemap_->set( TorsionID(insert_end, BB, psi_torsion), true);
	combined.conformation().insert_ideal_geometry_at_polymer_bond(insert_end);

	TR << "ideal " << insert_end << std::endl;
	combined.set_omega(insert_end, 180);//Note, upper Nter idealize loop takes care of omega for insert_start-1
	combined.set_psi(insert_end, -40);

	//Cter regions
	i=insert_end+1;
	for(; i<=Cter_loop_end_; ++i) {
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		if (movemap_->get_bb(i)){
			combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
			combined.set_omega(i, 180);
			TR << "mobile " << i << std::endl;
		}
	}



	//centroidize the pose before we do stuff to it - sidechains are expensive and unnecessary
	protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap(core::chemical::CENTROID);
	protocols::simple_moves::ReturnSidechainMover return_sidechains( combined );
	typeset_swap.apply( combined );

	//combined.dump_pdb("combined_preclose_cen.pdb");

	//[Note- Consider having the ability to set a small mover.]
	protocols::simple_moves::SmallMover small(movemap_, 10, 200); //huge moves for sampling
	small.angle_max( 'H', 180.0 );
	small.angle_max( 'E', 180.0 );
	small.angle_max( 'L', 180.0 );


	protocols::simple_moves::MinMover min_mover(movemap_, cen_scorefxn_, mintype_, 0.01, true /*use_nblist*/ );
    

	/////////////////////////TESTING/////////////////////////////////////////////////////////////

	if (test_control_mode_){perturb_backbone_for_test(combined, movemap_);}

	/////////////////////////Monte Carlo////////////////////////////////////////////////////////
	using protocols::moves::MonteCarlo;
	using protocols::moves::MonteCarloOP;
	MonteCarlo mc(combined, (*cen_scorefxn_), 0.8);

	/////////////////////////Protocol//////////////////////////////////////////////////////////
	TR << "start " << ((*cen_scorefxn_))(combined) << std::endl;

	for( core::Size i(1); i<=cycles_; ++i){

		if (!skip_sampling_){small.apply(combined);}

        for (protocols::loops::Loops::const_iterator it=loop_set->begin(), it_end=loop_set->end(); it!=it_end; ++it){

        	loop_set_map[*it]->apply(combined);

			combined.conformation().insert_ideal_geometry_at_polymer_bond(it->cut());

            min_mover.apply(combined);

			combined.conformation().insert_ideal_geometry_at_polymer_bond(it->cut());
        }
        
		if(mc.boltzmann(combined)) TR << i << " " << ((*cen_scorefxn_))(combined) << std::endl;
          
	}
        
	mc.recover_low(combined);
	TR << "finish " << ((*cen_scorefxn_))(combined) << std::endl;
	//combined.conformation().insert_ideal_geometry_at_polymer_bond(Cter_loop.cut());

	return_sidechains.apply( combined );




	//Remove cutpoints that were required for CCD.
	for (protocols::loops::Loops::const_iterator it=loop_set->begin(), it_end=loop_set->end(); it!=it_end; ++it){
		core::pose::remove_variant_type_from_pose_residue(combined, core::chemical::CUTPOINT_LOWER, it->cut());
		core::pose::remove_variant_type_from_pose_residue(combined, core::chemical::CUTPOINT_UPPER, it->cut()+1);
	}

	//Give back foldtree from pose_into_pose.
	combined.fold_tree(original_ft);
	TR << "Complete"<<std::endl;
	pose = combined;

}

void
AnchoredGraftMover::repack_connection_and_residues_in_movemap(Pose & pose, ScoreFunctionOP fa_scorefxn){

	PackerTaskOP task = TaskFactory::create_packer_task( pose );

	task->restrict_to_repacking();
	task->temporarily_fix_everything();

	task->temporarily_set_pack_residue(start_, true);
    task->temporarily_set_pack_residue(start_+1, true);
    
    task->temporarily_set_pack_residue(end_, true);
    task->temporarily_set_pack_residue(end_-1, true);
    
    for (Size i=1; i<=pose.total_residue(); ++i){
        if (movemap_->get_chi(i)){
            task->temporarily_set_pack_residue(i, true);
        }
    }
    
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover(fa_scorefxn, task);
	packer->apply(pose);
}

void
AnchoredGraftMover::repack_connection_and_residues_in_movemap_and_piece(Pose & pose, ScoreFunctionOP fa_scorefxn){

	PackerTaskOP task = TaskFactory::create_packer_task( pose );
	task->restrict_to_repacking();
	task->temporarily_fix_everything();
	task->temporarily_set_pack_residue(start_, true);
    task->temporarily_set_pack_residue(start_+1, true);
    
    task->temporarily_set_pack_residue(end_, true);
    task->temporarily_set_pack_residue(end_-1, true);
    
    for (Size i=1; i<=pose.total_residue(); ++i){
        if (movemap_->get_chi(i)){
            task->temporarily_set_pack_residue(i, true);
        }
    }
    for (Size i=start_+2; i<=end_-2; ++i){
        task->temporarily_set_pack_residue(i, true);
    }
    
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover(fa_scorefxn, task);
	packer->apply(pose);
}


}
}
