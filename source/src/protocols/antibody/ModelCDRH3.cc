// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/ModelCDRH3.cc
/// @brief models CDR H3 loop using loop modeling
/// @details
///// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/antibody/ModelCDRH3.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <utility/exit.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/H3PerturbCCD.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/H3CterInsert.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/util.hh>


static basic::Tracer TR( "protocols.antibody.ModelCDRH3" );

using namespace core;

namespace protocols {
namespace antibody {

ModelCDRH3::ModelCDRH3() : Mover() {}


ModelCDRH3::~ModelCDRH3() = default;

ModelCDRH3::ModelCDRH3( AntibodyInfoOP antibody_info) : Mover() {
	user_defined_ = false;
	ab_info_ = antibody_info;

	init();
}


ModelCDRH3::ModelCDRH3( AntibodyInfoOP antibody_info,
	core::scoring::ScoreFunctionCOP lowres_scorefxn) : Mover() {
	user_defined_ = true;
	ab_info_ = antibody_info;
	lowres_scorefxn_  = lowres_scorefxn->clone();

	init();
}


void ModelCDRH3::init( ) {
	Mover::type( "ModelCDRH3" );

	set_default();

	//TODO:
	//JQX: need to deal with this
	if ( ab_info_->is_camelid() && ab_info_->get_H3_kink_type()!=Kinked && ab_info_->get_H3_kink_type()!=Extended ) {
		c_ter_stem_ = 0;
	}

	h3_cter_insert_mover_ = H3CterInsertOP( new H3CterInsert(ab_info_) );
	h3_perturb_ccd_build_ = H3PerturbCCDOP( new H3PerturbCCD(ab_info_, lowres_scorefxn_) );
}


void ModelCDRH3::set_default() {
	benchmark_          = false;
	do_cter_insert_     = true;
	loops_flag_         = true;
	dle_flag_           = true;
	bad_nter_           = true;

	remodel_            = "legacy_perturb_ccd";

	c_ter_stem_ = 3;
	max_cycle_ = 20;

	cen_cst_ = 10.0;

	extend_h3_ = true;
	idealize_h3_stems_ = false;

	if ( !user_defined_ ) {
		lowres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
		lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
		lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
	}
}


void ModelCDRH3::set_lowres_score_func(scoring::ScoreFunctionCOP lowres_scorefxn ) {
	lowres_scorefxn_ = lowres_scorefxn->clone();
}


void ModelCDRH3::set_task_factory(pack::task::TaskFactoryCOP tf) {
	tf_ = pack::task::TaskFactoryOP( new pack::task::TaskFactory(*tf) );
}


void ModelCDRH3::turn_off_H3_filter() {
	h3_perturb_ccd_build_->turn_off_H3_filter();
}


void ModelCDRH3::apply( pose::Pose & pose_in ) {


	TR << "Applying CDR H3 modeler" << std::endl;

	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::moves;

	/// FIXME: JQX: very redudent here, just get one loops object
	Size framework_loop_begin( ab_info_->get_CDR_loop(h3).start() - 2 ); // use the AHo loop defs without mucking up the code everywhere.
	Size framework_loop_end  ( ab_info_->get_CDR_loop(h3).stop()  );
	Size cutpoint = ab_info_->get_CDR_loop(h3).cut() ; // keep the cutpoint unchanged
	Size framework_loop_size = (framework_loop_end - framework_loop_begin) + 1;

	loops::Loop cdr_h3( framework_loop_begin, framework_loop_end, cutpoint, 0, true );
	loops::Loop trimmed_cdr_h3(framework_loop_begin, framework_loop_end - c_ter_stem_, cutpoint, 0, true );
	loops::Loop input_loop;

	if ( framework_loop_size <= 6 ) {
		do_cter_insert_ = false;
		TR<<"loop_size <= 6, AUTOMATICALLY TURNING OFF THE C_TERMINAL INSERT"<<std::endl;
	}

	if ( do_cter_insert_ ) {
		//JQX: the h3 loop removing the cterminal 3 residues
		input_loop = trimmed_cdr_h3;
	} else {
		//JQX: the original h3 loop
		input_loop = cdr_h3;
	}

	simple_one_loop_fold_tree( pose_in, cdr_h3 );
	//    TR<<"*******************************************"<<std::endl;
	//    TR<<pose_in.fold_tree()<<std::endl;
	//    TR<<"*******************************************"<<std::endl;

	if ( idealize_h3_stems_ ) {
		set_single_loop_fold_tree( pose_in, cdr_h3 );

		TR<<"Idealizing the H3 stems ........ "<<std::endl;
		core::conformation::idealize_position(framework_loop_begin-1, pose_in.conformation());
		pose_in.set_omega( framework_loop_begin - 1, 179.6 );
		core::conformation::idealize_position(framework_loop_end+1, pose_in.conformation());
		pose_in.set_omega( framework_loop_end + 1, 179.6 );
		// 179.6 is from Daisuke's literature search, it is on Graylab wiki for idealization benchmark

		simple_one_loop_fold_tree( pose_in, cdr_h3 );
	}

	if ( extend_h3_ ) {
		TR<<"Extend the H3 loop ..........."<<std::endl;
		set_extended_torsions( pose_in, cdr_h3 );
	}

	if ( remodel_=="legacy_perturb_ccd" ) {
		//pose_in.dump_pdb("extended_idealized_centroid.pdb");
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


		h3_perturb_ccd_build_->pass_the_loop(input_loop);

	} else {
		/// FIXME: JQX very redudent loops defitions
		// use H3 to define a loops object
		loops::LoopsOP pass_loops( new loops::Loops() );
		pass_loops->add_loop( input_loop   );
		pass_loops->set_extended(true); // the IndepdentLoopMover will extend the loop

		// create a LoopMover type based on the string remode_
		remodel_mover_ = utility::pointer::static_pointer_cast< loops::loop_mover::IndependentLoopMover >
			( loops::LoopMoverFactory::get_instance()->create_loop_mover(remodel_, pass_loops) ) ;
		if ( !remodel_mover_ ) {
			utility_exit_with_message( "Error: no remodel mover defined!" );
		}

		// deal with the fragment files if the remodel_ type is not KIC
		if ( remodel_ != "perturb_kic"  ) {
			utility::vector1< core::fragment::FragSetOP > frag_libs;
			loops::read_loop_fragments( frag_libs );
			runtime_assert( frag_libs.size() > 0 );
			for ( Size i = 1; i <= frag_libs.size(); ++i ) {
				remodel_mover_->add_fragments( frag_libs[i]) ;
			}
		}

		// if you have native structure to compare, do this
		if ( get_native_pose() ) remodel_mover_->set_native_pose(get_native_pose()) ;

		// scoring function, their default scoring function is the same as specified here, but put it anyway
		remodel_mover_->set_scorefxn( lowres_scorefxn_ );

	}


	/*  JQX: the following code is probably not ncessary*/
	if ( bad_nter_ ) {
		//Size unaligned_cdr_loop_begin(0)/*, unaligned_cdr_loop_end(0)*/;
		std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
		core::import_pose::pose_from_file( hfr_pose_, path+"hfr.pdb" , core::import_pose::PDB_file);
		Size const unaligned_cdr_loop_begin( hfr_pose_.pdb_info()->pdb2pose('H', 95) );
		//unaligned_cdr_loop_end   = hfr_pose_.pdb_info()->pdb2pose('H', 103);
		//unaligned_cdr_loop_end -= 1 ;

		if ( framework_loop_size > 4 ) { //JQX: add this if statement to match R2_antibody
			pose_in.set_psi  (framework_loop_begin - 1, hfr_pose_.psi( unaligned_cdr_loop_begin - 1 )   );
			pose_in.set_omega(framework_loop_begin - 1, hfr_pose_.omega( unaligned_cdr_loop_begin - 1 ) );
		}
		//pose_in.dump_pdb("after_copying_nter.pdb");
	}

	antibody::AntibodyInfoOP starting_antibody;
	starting_antibody = antibody::AntibodyInfoOP( new AntibodyInfo(*ab_info_) );
	bool closed_cutpoints( false );


	Size cycle ( 1 );
	while ( !closed_cutpoints && cycle < max_cycle_ ) {
		ab_info_ = starting_antibody;
		if ( do_cter_insert_ ) {
			h3_cter_insert_mover_->apply(pose_in);
		}
		//pose_in.dump_pdb("after_c_insert.pdb");

		if ( remodel_=="legacy_perturb_ccd" ) {
			h3_perturb_ccd_build_->apply(pose_in);
		} else {
			remodel_mover_->apply(pose_in);
		}


		closed_cutpoints = cutpoints_separation( pose_in, ab_info_ );
		++cycle;
	} // while( ( cut_separation > 1.9 )

	TR <<  "Finished Modeling Centroid CDR H3 loop" << std::endl;

	TR << "Finished applying CDR H3 modeler" << std::endl;


	return;
} // ModelCDRH3::apply()


std::string ModelCDRH3::get_name() const {
	return "ModelCDRH3";
}


//basic::Tracer & my_LoopMover::tr() const{
//    return TR;
//}


} // namespace antibody
} // namespace protocols


