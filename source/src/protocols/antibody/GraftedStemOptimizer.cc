// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/GraftedStemOptimizer.cc
/// @brief Optimize the CDR Grafted stem
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody/GraftedStemOptimizer.hh>
#include <basic/Tracer.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/util.hh> //idealize
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/moves/PyMolMover.hh>



static basic::Tracer TRG("protocols.antibody.GraftedStemOptimizer");

namespace protocols {
namespace antibody {
using namespace core;



GraftedStemOptimizer::GraftedStemOptimizer( CDRNameEnum const & cdr_name,
        AntibodyInfoOP antibody_info) : Mover( "GraftedStemOptimizer" ) {
	cdr_name_  = cdr_name;
	ab_info_   = antibody_info;

	init();
}

void
GraftedStemOptimizer::init() {
	scorefxn_       = NULL;
	mc_             = NULL;
	optimize_stems_ = NULL;
	benchmark_ = false;

	stem_size_        = 4;
	deep_optimization_ = true;
	cdr_loop_    = new loops::Loop( ab_info_->get_CDR_loop(cdr_name_) );
}


void
GraftedStemOptimizer::set_stem_size(Size const & setting) {
	stem_size_=setting;
	if( stem_size_ % 2 != 0 ) {
		TRG<<"The stem_size_ must be dividable by 2"<<std::endl;
		exit(-1);
	}
}

GraftedStemOptimizer::~GraftedStemOptimizer() {}


void
GraftedStemOptimizer::setup_protocol(pose::Pose & pose) {

	TRG<<"          Setting Up the Optimizer Mover ..... "<<std::endl;

	/// scoring function
	if(!scorefxn_) {
		scorefxn_=scoring::get_score_function();
		scorefxn_->set_weight( scoring::chainbreak, 30./3. );
		scorefxn_->set_weight( scoring::overlap_chainbreak, 30./3. );
		scorefxn_->set_weight( scoring::dslf_ss_dst, 3.0);
		scorefxn_->set_weight( scoring::dslf_cs_ang, 3.0);
		scorefxn_->set_weight( scoring::dslf_ss_dih, 3.0);
		scorefxn_->set_weight( scoring::dslf_ca_dih, 3.0);

	}

	/// PyMol_Mover
	//moves::PyMolMoverOP pymol = new moves::PyMolMover();
	//pymol->keep_history(true);

	using namespace protocols::simple_moves;

	/// Small_Mover and Shear_Mover
	Real high_move_temp = 2.00;
	Size n_small_moves = 5 ;
	BackboneMoverOP small_mover = new SmallMover( get_stem_movemap(pose,"NC"), high_move_temp, n_small_moves );
	BackboneMoverOP shear_mover = new ShearMover( get_stem_movemap(pose,"NC"), high_move_temp, n_small_moves );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	/// Ccd_Loop_Closure_Mover
	loops::Loop Nter_stem(cdr_loop_->start()-stem_size_, cdr_loop_->start()+stem_size_-1, cdr_loop_->start()-1);
	loops::Loop Cter_stem(cdr_loop_->stop()-stem_size_+1, cdr_loop_->stop()+stem_size_, cdr_loop_->stop() );
	using namespace loops::loop_closure::ccd;
	CcdLoopClosureMoverOP close_Nter_stem=new CcdLoopClosureMover(Nter_stem, get_stem_movemap(pose, "N"));
	close_Nter_stem->set_tolerance(0.001);
	CcdLoopClosureMoverOP close_Cter_stem=new CcdLoopClosureMover(Cter_stem, get_stem_movemap(pose, "C"));
	close_Cter_stem->set_tolerance(0.001);

	/// Rotamer_Trial_Mover
	RotamerTrialsMoverOP rotamer_trial_mover = new RotamerTrialsMover( scorefxn_, get_stem_taskfactory(pose, "NC") );

	/// Pack_Rotamers_Mover
	PackRotamersMoverOP repack = new PackRotamersMover( scorefxn_ );
	repack->task_factory( get_stem_taskfactory(pose, "NC") );

	/// Min_Mover
	MinMoverOP min_mover = new MinMover( get_stem_movemap(pose, "NC"), scorefxn_, "dfpmin_armijo_nonmonotone", 0.001, true );

	/// Sequence_Mover
	optimize_stems_ = new moves::SequenceMover();

	if(deep_optimization_) {
		optimize_stems_ -> add_mover(small_mover);
		optimize_stems_ -> add_mover(shear_mover);
		optimize_stems_ -> add_mover(close_Nter_stem);
		optimize_stems_ -> add_mover(close_Cter_stem);
	}
	optimize_stems_ -> add_mover(rotamer_trial_mover);
	optimize_stems_ -> add_mover(repack);
	optimize_stems_ -> add_mover(min_mover);

	TRG<<"          Finished Setting Up the Optimizer Mover!! "<<std::endl;

}


void
GraftedStemOptimizer::apply( pose::Pose & pose ) {

	TRG<<"Optimizing the stem of the Grafted CDR   "<< ab_info_->get_CDR_name(cdr_name_) <<" ............"<<std::endl;

	setup_protocol(pose);

	pose.fold_tree(*get_N_C_stems_foldtree(pose));

	TRG<<"Score Before Stem Idealization: "<<(*scorefxn_)(pose)<<std::endl;

	/*
	if (deep_optimization_){
	    for(Size i=cdr_loop_->start()-stem_size_; i<=cdr_loop_->start()+stem_size_-1; ++i){
	        conformation::idealize_position( i, pose.conformation() );
	        pose.set_omega(i, 179.8);
	    }
	    for(Size i=cdr_loop_->stop()-stem_size_+1; i<=cdr_loop_->stop()+stem_size_; ++i){
	        conformation::idealize_position( i, pose.conformation() );
	        pose.set_omega(i, 179.8);
	    }
	    // REFERENCE: protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh and .cc
	    // see the   OMEGA_MEAN_(179.8)
	}
	 */


	// adding cutpoint variants for chainbreak score computation
	loops::remove_cutpoint_variants( pose, true ); //remove first
	loops::add_cutpoint_variants( pose );

	TRG<<"Score After Idealization, but Before Stem Optimization: "<<(*scorefxn_)(pose)<<std::endl;




	Size inner_cycles( stem_size_ * 3 );
	Size outer_cycles( 5 );
	if(benchmark_) {
		inner_cycles=1;
		outer_cycles=1;
	}
	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	mc_ = new moves::MonteCarlo( pose, *scorefxn_, temperature );
	mc_->reset( pose ); // monte carlo reset
	moves::TrialMoverOP optimize_trial = new moves::TrialMover(optimize_stems_, mc_);


	for(Size i = 1; i <= outer_cycles; i++) {
		mc_->recover_low( pose );
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc_->set_temperature( temperature );
			optimize_trial->apply( pose );
			//TRG<<pose.fold_tree()<<std::endl;

		} // inner cycles
	} // outer cycles
	mc_->recover_low( pose );


	TRG<<"Score After Stem Optimization: "<<(*scorefxn_)(pose)<<std::endl;


	TRG<<"Stem of CDR "<<ab_info_->get_CDR_name(cdr_name_) <<" has finished!!!"<<std::endl;


}





///
/// 	      ####LLLLLLLLLLLLLLLLLLL####              L1-3.pdb, H1-3.pdb
/// ..ffffffff@@@@                   @@@@ffffffffff..  framework.pdb
///                       ||
///                       || grafting
///                       \/
///   ffffffff@@##LLLLLLLLLLLLLLLLLLL##@@ffffffffff
///
///
///           _____ set up fold tree  _____
///          |     |                 |     |
///          |  |  |                 |  |  |
///   ffffffff@@|##LLLLLLLLLLLLLLLLLLL##|@@ffffffffff
///
///

void
GraftedStemOptimizer::set_scorefxn(core::scoring::ScoreFunctionOP setting) {
	scorefxn_ = setting;
}

kinematics::FoldTreeOP
GraftedStemOptimizer::get_N_C_stems_foldtree( pose::Pose const & pose ) const  {
	using namespace core::kinematics;

	FoldTreeOP ft = new FoldTree();
	ft ->clear();

	const Size  jumppoint1=cdr_loop_->start()-stem_size_-1;
	const Size  cutpoint1=jumppoint1+(stem_size_/2);
	const Size  jumppoint2=cdr_loop_->start();

	const Size  jumppoint3=cdr_loop_->stop();
	const Size  cutpoint2=cdr_loop_->stop()+(stem_size_/2);
	const Size  jumppoint4=cdr_loop_->stop()+stem_size_+1;

	ft->add_edge( 1, jumppoint1, Edge::PEPTIDE );
	ft->add_edge( jumppoint1, jumppoint2, 1 );
	ft->add_edge( jumppoint1, cutpoint1, Edge::PEPTIDE );
	ft->add_edge( jumppoint2, cutpoint1 + 1, Edge::PEPTIDE );

	ft->add_edge( jumppoint2, jumppoint3, Edge::PEPTIDE );

	ft->add_edge( jumppoint3, jumppoint4, 2 );
	ft->add_edge( jumppoint3, cutpoint2, Edge::PEPTIDE );
	ft->add_edge( jumppoint4, cutpoint2 + 1, Edge::PEPTIDE );
	ft->add_edge( jumppoint4, pose.total_residue(), Edge::PEPTIDE );


	//TRG<<"##################################################"<<std::endl;
	//TRG<<ab_info_->get_CDR_name(cdr_name_) <<std::endl;
	//TRG<<"start = "<<cdr_loop_->start() <<std::endl;
	//TRG<<"stop  = "<<cdr_loop_->stop()  <<std::endl;
	//TRG<<*ft<<std::endl;
	//TRG<<"##################################################"<<std::endl;

	return ft;

}

kinematics::FoldTreeOP
GraftedStemOptimizer::get_Nstem_foldtree( pose::Pose const & pose ) const  {
	using namespace core::kinematics;

	FoldTreeOP ft = new FoldTree();
	ft ->clear();

	const Size  jumppoint1=cdr_loop_->start()-stem_size_-1;
	const Size  cutpoint1=jumppoint1+(stem_size_/2);
	const Size  jumppoint2=cdr_loop_->start();

	ft->add_edge( 1, jumppoint1, Edge::PEPTIDE );
	ft->add_edge( jumppoint1, jumppoint2, 1 );
	ft->add_edge( jumppoint1, cutpoint1, Edge::PEPTIDE );
	ft->add_edge( jumppoint2, cutpoint1 + 1, Edge::PEPTIDE );
	ft->add_edge( jumppoint2, pose.total_residue(), Edge::PEPTIDE );

	return ft;
}


kinematics::FoldTreeOP
GraftedStemOptimizer::get_Cstem_foldtree( pose::Pose const & pose ) const  {
	using namespace core::kinematics;

	FoldTreeOP ft = new FoldTree();
	ft ->clear();

	const Size  jumppoint3=cdr_loop_->stop();
	const Size  cutpoint2=cdr_loop_->stop()+(stem_size_/2);
	const Size  jumppoint4=cdr_loop_->stop()+stem_size_+1;

	ft->add_edge( 1, jumppoint3, Edge::PEPTIDE );
	ft->add_edge( jumppoint3, jumppoint4, 1 );
	ft->add_edge( jumppoint3, cutpoint2, Edge::PEPTIDE );
	ft->add_edge( jumppoint4, cutpoint2 + 1, Edge::PEPTIDE );
	ft->add_edge( jumppoint4, pose.total_residue(), Edge::PEPTIDE );

	return ft;
}




kinematics::MoveMapOP
GraftedStemOptimizer::get_stem_movemap( pose::Pose const & pose, std::string const & type, bool const & include_nb_sc ) const  {
	using namespace core::kinematics;

	MoveMapOP mm = new MoveMap();
	mm->clear();

	mm->set_chi( false );
	mm->set_bb( false );

	utility::vector1< bool> bb_is_flexible( pose.total_residue(), false );
	utility::vector1< bool> sc_is_flexible( pose.total_residue(), false );

	if ( (type == "N") || (type == "NC")  ) {
		for(Size i=cdr_loop_->start()-stem_size_; i<=cdr_loop_->start()-1; ++i) {
			bb_is_flexible[i]=true;
			sc_is_flexible[i]=true;
		}
	}
	if ( (type == "C") || (type == "NC")  ) {
		for(Size i=cdr_loop_->stop()+1; i<=cdr_loop_->stop()+stem_size_; ++i) {
			bb_is_flexible[i]=true;
			sc_is_flexible[i]=true;
		}
	}

	if(include_nb_sc) {
		loops::get_tenA_neighbor_residues(pose, sc_is_flexible);
	}
	// a function in the loops namespace to calculate the neighbors.
	// for the "true" values in "sc_is_flexible", the function will find
	// neighbors of these "true" residues within 10A (cb-cb distance).
	// "CYD" residues will be automatically turned off

	mm->set_bb( bb_is_flexible );
	mm->set_chi( sc_is_flexible );

	mm->set_jump( 1, false );
	if(type == "NC") {
		mm->set_jump( 2, false );
	}

	//TRG<<*mm<<std::endl;

	return mm;
}



pack::task::TaskFactoryOP
GraftedStemOptimizer::get_stem_taskfactory( pose::Pose & pose, std::string const & type, bool const & include_nb_sc ) const  {
	pack::task::TaskFactoryOP tf = new pack::task::TaskFactory();
	tf->clear();

	tf= setup_packer_task(pose);

	vector1< bool> sc_is_packable( pose.total_residue(), false );

	if ( (type == "N") || (type == "NC")  ) {
		for(Size i=cdr_loop_->start()-stem_size_; i<=cdr_loop_->start()-1; ++i) {
			sc_is_packable[i]=true;
		}
	}
	if ( (type == "C") || (type == "NC")  ) {
		for(Size i=cdr_loop_->stop()+1; i<=cdr_loop_->stop()+stem_size_; ++i) {
			sc_is_packable[i]=true;
		}
	}

	if(include_nb_sc) {
		loops::get_tenA_neighbor_residues(pose, sc_is_packable);
	}

	using namespace protocols::toolbox::task_operations;
	tf->push_back( new RestrictToInterface(sc_is_packable) );

	//tf->create_task_and_apply_taskoperations(pose);
	//pack::task::PackerTaskOP my_task(tf->create_task_and_apply_taskoperations(pose));
	//TRG<<*my_task<<std::endl;

	return tf;
}



std::string
GraftedStemOptimizer::get_name() const {
	return "GraftedStemOptimizer";
}

// copy ctor
GraftedStemOptimizer::GraftedStemOptimizer( GraftedStemOptimizer const & rhs ) : Mover(rhs) {
	initForEqualOperatorAndCopyConstructor(*this, rhs);
}


///@brief assignment operator
GraftedStemOptimizer & GraftedStemOptimizer::operator=( GraftedStemOptimizer const & rhs ) {
	//abort self-assignment
	if (this == &rhs) return *this;
	Mover::operator=(rhs);
	initForEqualOperatorAndCopyConstructor(*this, rhs);
	return *this;
}

void
GraftedStemOptimizer::initForEqualOperatorAndCopyConstructor(GraftedStemOptimizer & lhs, GraftedStemOptimizer const & rhs) {
	lhs.stem_size_    =  rhs.stem_size_;
	lhs.cdr_name_     =  rhs.cdr_name_;
	lhs.ab_info_      =  rhs.ab_info_;
	lhs.benchmark_    =  rhs.benchmark_;
	lhs.cdr_loop_     =  rhs.cdr_loop_;
	lhs.scorefxn_     =  rhs.scorefxn_;
}





}  // namespace antibody
}  // namespace protocols
