// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/AntibodyModelerProtocol.cc
/// @brief Build a homology model of an antibody
/// @details
/// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/jobdist/JobDistributors.hh> // SJF Keep first for mpi and mac boinc build (fails if not first for some odd reason) -dek
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyModelerProtocol.hh>
#include <protocols/antibody/CDRsMinPackMin.hh>
#include <protocols/antibody/ModelCDRH3.hh>
#include <protocols/antibody/RefineBetaBarrel.hh>
#include <protocols/antibody/RefineOneCDRLoop.hh>
#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/util.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

using namespace ObjexxFCL::format;


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.antibody.AntibodyModelerProtocol" );
using namespace core;

namespace protocols {
namespace antibody {

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
	return( protocols::moves::MoverOP( new AntibodyModelerProtocol() ) );
}


void AntibodyModelerProtocol::init() {
	Mover::type( "AntibodyModelerProtocol" );

	set_default();
	init_from_options();
	setup_objects();
}


void AntibodyModelerProtocol::set_default() {
	TR <<  "Setting Up defaults.........." << std::endl;
	model_h3_              = true;
	h3_filter_         = true;
	cter_insert_       = true;
	snugfit_               = true;

	LH_repulsive_ramp_ = true;
	refine_h3_             = true;
	flank_residue_min_ = true;
	middle_pack_min_       = true;

	bad_nter_  = true;
	extend_h3_before_modeling_ = true;
	idealize_h3_stems_before_modeling_ = false;

	benchmark_ = false;
	camelid_   = false;
	camelid_constraints_ = false;
	use_csts_ = false;
	constrain_vlvh_qq_ = false;
	constrain_cter_ = false;
	packonly_after_graft_=false;

	cst_weight_ = 0.0;
	cen_cst_ = 10.0;
	high_cst_ = 100.0; // if changed here, please change at the end of AntibodyModeler as well
	flank_residue_size_ = 2;
	h3_filter_tolerance_ = 20;

	sc_min_ = false;
	rt_min_ = false;

	h3_perturb_type_ = "legacy_perturb_ccd"; // legacy_perturb_ccd, kic, ccd
	h3_refine_type_  = "legacy_refine_ccd"; // legacy_refine, kic, ccd

	cdr_constraint_ = NULL;
}


void AntibodyModelerProtocol::register_options() {
	using namespace basic::options;

	option.add_relevant( OptionKeys::antibody::model_h3 );
	option.add_relevant( OptionKeys::antibody::snugfit );
	option.add_relevant( OptionKeys::antibody::camelid );
	option.add_relevant( OptionKeys::antibody::camelid_constraints );
	option.add_relevant( OptionKeys::run::benchmark );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::antibody::constrain_vlvh_qq );
	option.add_relevant( OptionKeys::antibody::constrain_cter );
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::antibody::refine_h3 );
	option.add_relevant( OptionKeys::antibody::h3_filter );
	option.add_relevant( OptionKeys::antibody::cter_insert );
	option.add_relevant( OptionKeys::antibody::sc_min);
	option.add_relevant( OptionKeys::antibody::rt_min);
	option.add_relevant( OptionKeys::antibody::flank_residue_min);
	//option.add_relevant( OptionKeys::antibody::flank_residue_size);
	option.add_relevant( OptionKeys::antibody::remodel);
	option.add_relevant( OptionKeys::antibody::refine);
	//option.add_relevant( OptionKeys::antibody::middle_pack_min);
	option.add_relevant( OptionKeys::antibody::h3_filter_tolerance);
	option.add_relevant( OptionKeys::antibody::bad_nter);
	option.add_relevant( OptionKeys::antibody::extend_h3_before_modeling);
	option.add_relevant( OptionKeys::antibody::idealize_h3_stems_before_modeling);
	option.add_relevant( OptionKeys::antibody::packonly_after_graft);
}


void AntibodyModelerProtocol::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR <<  "Start Reading and Setting Options ..." << std::endl;

	if ( option[OptionKeys::antibody::model_h3].user() ) {
		set_ModelH3( option[OptionKeys::antibody::model_h3]() );
	}
	if ( option[ OptionKeys::antibody::snugfit ].user() ) {
		set_SnugFit( option[ OptionKeys::antibody::snugfit ]() );
	}
	if ( option[ OptionKeys::antibody::refine_h3 ].user() ) {
		set_refine_h3( option[ OptionKeys::antibody::refine_h3 ]()  );
	}
	if ( option[ OptionKeys::antibody::cter_insert ].user() ) {
		set_CterInsert( option[ OptionKeys::antibody::cter_insert ]() );
	}
	if ( option[ OptionKeys::antibody::h3_filter ].user() ) {
		set_H3Filter ( option[ OptionKeys::antibody::h3_filter ]() );
	}
	if ( option[ OptionKeys::antibody::h3_filter_tolerance ].user() ) {
		set_H3Filter_Tolerance( option[ OptionKeys::antibody::h3_filter_tolerance ]()  );
	}
	if ( option[ OptionKeys::antibody::flank_residue_min ].user() ) {
		set_flank_residue_min ( option[ OptionKeys::antibody::flank_residue_min ]() );
	}
	if ( option[ OptionKeys::run::benchmark ].user() ) {
		set_BenchMark( option[ OptionKeys::run::benchmark ]() );
	}
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		set_cst_weight( option[ OptionKeys::constraints::cst_weight ]() );
	}
	if ( option[ OptionKeys::constraints::cst_file].user() ) {
		set_use_constraints( option[ OptionKeys::constraints::cst_file ].user() ); // .user here??
	}
	if ( option[ OptionKeys::antibody::constrain_cter ].user() ) {
		set_constrain_cter( option[ OptionKeys::antibody::constrain_cter ]() );
	}
	if ( option[ OptionKeys::antibody::constrain_vlvh_qq ].user() ) {
		set_constrain_vlvh_qq( option[ OptionKeys::antibody::constrain_vlvh_qq ]() );
	}
	if ( option[ OptionKeys::antibody::packonly_after_graft ].user()  ) {
		set_packonly_after_graft( option[ OptionKeys::antibody::packonly_after_graft ]()  );
	}
	if ( option[ OptionKeys::antibody::sc_min ].user() ) {
		set_sc_min( option[ OptionKeys::antibody::sc_min ]() );
	}
	if ( option[ OptionKeys::antibody::rt_min ].user() ) {
		set_rt_min( option[ OptionKeys::antibody::rt_min ]() );
	}
	if ( option[ OptionKeys::antibody::remodel ].user() ) {
		set_perturb_type( option[ OptionKeys::antibody::remodel ]() );
	}
	if ( option[ OptionKeys::antibody::refine ].user() ) {
		set_refine_type( option[ OptionKeys::antibody::refine ]() );
	}
	if ( option[ OptionKeys::antibody::bad_nter].user()  ){
		set_bad_nter(option[ OptionKeys::antibody::bad_nter]() );
	}
	if ( option[ OptionKeys::antibody::extend_h3_before_modeling].user()  ){
		set_extend_h3_before_modeling(option[ OptionKeys::antibody::extend_h3_before_modeling]() );
  }
	if ( option[ OptionKeys::antibody::idealize_h3_stems_before_modeling].user()  ){
		set_idealize_h3_stems_before_modeling(option[ OptionKeys::antibody::idealize_h3_stems_before_modeling]() );
	}
	
    //if ( option[ OptionKeys::antibody::middle_pack_min].user() ){
    //  set_middle_pack_min( option[ OptionKeys::loops::refine ] )
    //}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//if( option[ out::levels ].active() ) T("Levels:") << option[ out::levels ]() <<  std::endl;
	option[ out::levels ]("protocols.simple_moves.ConstraintSetMover:debug");  // echo constraints to log
	//if( option[ out::levels ].active() ) T("Levels:") << option[ out::levels ]() <<  std::endl;

	//set native pose if asked for
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_pdb( *native_pose, option[ OptionKeys::in::file::native ]() );
		set_native_pose( native_pose );
	} else {
		set_native_pose(NULL);
	}

	TR <<  "Finish Reading and Setting Options !!!" << std::endl;
}


void
AntibodyModelerProtocol::setup_objects() {

	sync_objects_with_flags();

	if(use_csts_) {
		if ( cst_weight_ == 0.0 ) cst_weight_ = 1.0;
	}

	// setup all the scoring functions
	pack_scorefxn_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	dock_scorefxn_highres_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	dock_scorefxn_highres_->set_weight( core::scoring::chainbreak, 1.0 );
	dock_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
	if(constrain_vlvh_qq_) {
		dock_scorefxn_highres_->set_weight( scoring::atom_pair_constraint, cst_weight_ );
	}
	loop_scorefxn_centroid_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	loop_scorefxn_centroid_->set_weight( scoring::chainbreak, 10./3. );
	//loop_scorefxn_centroid_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
	loop_scorefxn_highres_ = scoring::get_score_function();
	loop_scorefxn_highres_->set_weight( scoring::chainbreak, 1.0 );
	loop_scorefxn_highres_->set_weight( scoring::overlap_chainbreak, 10./3. );
	if(constrain_cter_) {
		loop_scorefxn_highres_->set_weight(scoring::dihedral_constraint, cst_weight_);
	}
}

void AntibodyModelerProtocol::sync_objects_with_flags() {
	using namespace protocols::moves;
	flags_and_objects_are_in_sync_ = true;
}


std::string AntibodyModelerProtocol::get_name() const {
	return "AntibodyModelerProtocol";
}


void AntibodyModelerProtocol::finalize_setup( pose::Pose & pose ) {
	// check for native and input pose
	if ( !get_input_pose() ) {
		pose::PoseOP input_pose( new pose::Pose(pose) );
		set_input_pose( input_pose );   // JQX: pass the input_pose to the mover.input_pose_
	}

	pose::PoseOP native_pose;
	if ( !get_native_pose() ) {
		TR << "Danger Will Robinson! Native is an impostor!" << std::endl;
		TR << "   'native_pose' is just a copy of the 'input_pose'    " << std::endl;
		TR << "    since you didn't sepcifiy the native pdb name"<<std::endl;
		native_pose = pose::PoseOP( new pose::Pose(pose) );
	} else {
		native_pose = pose::PoseOP( new pose::Pose( *get_native_pose() ) );
	}

	pose::set_ss_from_phipsi( *native_pose ); // JQX: this is the secondary structure from the native pose

	set_native_pose( native_pose ); // pass the native pose to the mover.native_pose_

	ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose) );

	pose.fold_tree( * ab_info_->get_FoldTree_AllCDRs(pose) ) ;
	TR<<*ab_info_<<std::endl;

	//AntibodyInfoOP native_ab_info = new AntibodyInfo(*native_pose,camelid_);

	//core::pack::task::PackerTaskOP my_task2(tf_->create_task_and_apply_taskoperations(pose));
	//TR<<*my_task2<<std::endl; exit(-1);
}


//APPLY
void AntibodyModelerProtocol::apply( pose::Pose & pose ) {

	using namespace chemical;
	using namespace id;
	using namespace scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;

	// the default inital secstruct is all "L" loop!
	pose::Pose start_pose_ = pose;

	if ( !flags_and_objects_are_in_sync_ ) {
		sync_objects_with_flags();
	}

	finalize_setup(pose);

	basic::prof_reset();
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	pose::set_ss_from_phipsi( pose );

	// display constraints and return
	if( camelid_constraints_ ) {
		display_constraint_residues( pose );
		return;
	}

	// Step 1: model the cdr h3 in centroid mode
	// JQX notes: pay attention to the way it treats the stems when extending the loop
	if(model_h3_) {

		// switching to centroid mode
		simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
		simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

		// Building centroid mode loop
		pose::Pose start_pose = pose;
		to_centroid.apply( pose );

		// call ConstraintSetMover
		TR<<"Centroid cst_weight: "<<cst_weight_<<std::endl;
		if(  cst_weight_ != 0.0  ) {
			cdr_constraint_ = protocols::simple_moves::ConstraintSetMoverOP( new simple_moves::ConstraintSetMover() );
			cdr_constraint_->apply( pose );
		}

		ModelCDRH3OP model_cdrh3( new ModelCDRH3( ab_info_, loop_scorefxn_centroid_) );
		model_cdrh3->set_perturb_type(h3_perturb_type_); //legacy_perturb_ccd, ccd, kic
		if(cter_insert_ ==false) {
			model_cdrh3->turn_off_cter_insert();
		}
		if(h3_filter_   ==false) {
			model_cdrh3->turn_off_H3_filter();
		}
		model_cdrh3->set_bad_nter(bad_nter_);
		model_cdrh3->set_extend_h3(extend_h3_before_modeling_ );
		model_cdrh3->set_idealize_h3_stems(idealize_h3_stems_before_modeling_);
		model_cdrh3->apply( pose );
		//pose.dump_pdb("1st_finish_model_h3.pdb");

		// back to fullatom
		to_full_atom.apply( pose );

		utility::vector1<bool> allow_chi_copy( pose.total_residue(), true );
		/// FIXME: JQX very redudent loops defition
		for( Size ii=ab_info_->get_CDR_loop(h3).start(); ii<=ab_info_->get_CDR_loop(h3).stop(); ii++ ) {
			allow_chi_copy[ii] = false;
		}
		//recover sidechains from starting structures except H3
		protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose, allow_chi_copy );
		recover_sidechains.apply( pose );
	}

	// call ConstraintSetMover
	TR << "Full-atom cst_weight: " << cst_weight_ << std::endl;
	if(  cst_weight_ != 0.0  ) {
		cdr_constraint_ = protocols::simple_moves::ConstraintSetMoverOP( new simple_moves::ConstraintSetMover() );
		cdr_constraint_->apply( pose );
	}

	//if(middle_pack_min_){
	CDRsMinPackMinOP cdrs_min_pack_min( new CDRsMinPackMin(ab_info_) );
	if(sc_min_) cdrs_min_pack_min->set_sc_min(true);
	if(rt_min_) cdrs_min_pack_min->set_rt_min(true);
	cdrs_min_pack_min -> set_turnoff_minimization(packonly_after_graft_);
	cdrs_min_pack_min -> apply(pose);
	//}

	// Step 2: SnugFit: relieve the clashes between L-H
	if ( snugfit_ ) {
		RefineBetaBarrelOP refine_beta_barrel( new RefineBetaBarrel(ab_info_, dock_scorefxn_highres_, pack_scorefxn_) );
		// it has default movemap, tf, and fold_tree
		if (!LH_repulsive_ramp_) {
			refine_beta_barrel-> turn_off_repulsive_ramp();
		}
		if(sc_min_)              {
			refine_beta_barrel->set_sc_min(true);
		}
		if(rt_min_)              {
			refine_beta_barrel->set_rt_min(true);
		}
		refine_beta_barrel->apply(pose);
		//pose.dump_pdb("2nd_finish_snugfit.pdb");
	}

	// Step 3: Full Atom Relax
	if(refine_h3_) {
		RefineOneCDRLoopOP cdr_highres_refine_( new RefineOneCDRLoop(ab_info_, h3_refine_type_, loop_scorefxn_highres_) );
		cdr_highres_refine_ -> set_refine_mode(h3_refine_type_);
		cdr_highres_refine_ -> set_h3_filter(h3_filter_);
		cdr_highres_refine_ -> set_num_filter_tries(h3_filter_tolerance_);
		cdr_highres_refine_ -> set_flank_relax(flank_residue_min_);
		if (flank_residue_min_) cdr_highres_refine_->set_flank_size(flank_residue_size_);
		cdr_highres_refine_ -> pass_start_pose(start_pose_);
		cdr_highres_refine_ -> apply(pose);
		//pose.dump_pdb("3rd_finish_h3_refine.pdb");
	}

	//FoldTree, MoveMap, TaskFactory and Variants will be taken care of inside
	cdrs_min_pack_min -> apply(pose);
	//pose.dump_pdb("4th_final_min_pack_min.pdb");

	// Step 4: Store the homolgy models
	pose.fold_tree( * ab_info_->get_FoldTree_AllCDRs(pose) ) ;

	// Redefining CDR H3 cutpoint variants
	loops::remove_cutpoint_variants( pose, true );
	loops::add_cutpoint_variants( pose );

	// Final score (with constraints) before jd2 output the results
	if(constrain_vlvh_qq_) {
		loop_scorefxn_highres_->set_weight( scoring::atom_pair_constraint, cst_weight_ );
	}
	( *loop_scorefxn_highres_ )( pose );

	// Finish
	echo_metrics_to_jd2(pose,job);
	set_last_move_status( protocols::moves::MS_SUCCESS );
	basic::prof_show();
	TR<<"Antibody Modeling Protocol Finished!!!!"<<std::endl<<std::endl<<std::endl;

}// end apply


void AntibodyModelerProtocol::echo_metrics_to_jd2(core::pose::Pose & pose, protocols::jd2::JobOP job) {

	// align pose to native pose
	pose::Pose native_pose = *get_native_pose();
	antibody::AntibodyInfoOP native_ab_info( new AntibodyInfo(native_pose) );
	align_to_native( pose, native_pose, ab_info_, native_ab_info );

	// the specific constraint terms for output in the log file
	Real atom_pair_constraint_score = pose.energies().total_energies()[ core::scoring::atom_pair_constraint ];
	Real dihedral_constraint_score = pose.energies().total_energies()[ core::scoring::dihedral_constraint ];
	Real total_score = pose.energies().total_energies()[ core::scoring::total_score ];
	Real unconstrained_score = total_score - atom_pair_constraint_score - dihedral_constraint_score;

	TR << " 		atom_pair_constraint_score = " << atom_pair_constraint_score << std::endl;
	TR << "      dihedral_constraint_score = " << dihedral_constraint_score << std::endl;
	TR << "                    total_score = " << total_score << std::endl;
	TR << "            unconstrained_score = " << unconstrained_score << std::endl;

	job->add_string_real_pair( "unconstr_score", unconstrained_score );

	align_to_native( pose, native_pose, ab_info_, native_ab_info, "H" );
	job->add_string_real_pair("H3_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h3) ));
	job->add_string_real_pair("H2_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h2) ));
	job->add_string_real_pair("H1_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h1) ));
	//pose.dump_pdb("aligned_H.pdb");
	if( camelid_ == false ) {
		align_to_native( pose, native_pose, ab_info_, native_ab_info, "L" );
		job->add_string_real_pair("L3_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l3) ));
		job->add_string_real_pair("L2_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l2) ));
		job->add_string_real_pair("L1_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l1) ));
		//pose.dump_pdb("aligned_L.pdb");
	}
	
	//job->add_string_real_pair("AP_constraint", atom_pair_constraint_score);
	job->add_string_real_pair("VL_VH_distance", vl_vh_orientation_coords( pose, *ab_info_ )[1]);
	job->add_string_real_pair("VL_VH_opening_angle", vl_vh_orientation_coords( pose, *ab_info_ )[2]);
	job->add_string_real_pair("VL_VH_opposite_opening_angle", vl_vh_orientation_coords( pose, *ab_info_ )[3]);
	job->add_string_real_pair("VL_VH_packing_angle", vl_vh_orientation_coords( pose, *ab_info_ )[4]);

	//kink metrics
	job->add_string_real_pair( "kink_RD_HB", kink_RD_Hbond( pose, *ab_info_ ));
	job->add_string_real_pair( "kink_bb_HB", kink_bb_Hbond( pose, *ab_info_ ));
	job->add_string_real_pair( "kink_Trp_HB", kink_Trp_Hbond( pose, *ab_info_ ));

	std::pair<core::Real,core::Real> q = kink_dihedral( pose, *ab_info_);
	job->add_string_real_pair( "kink_q", q.first );
	job->add_string_real_pair( "kink_qbase", q.second );

	std::pair<ParatopeMetric< core::Real >, ParatopeMetric<core::Real> > sasa = paratope_sasa( pose, *ab_info_);
	ParatopeMetric< core::SSize> p_charge = paratope_charge( pose, *ab_info_ );
	job->add_string_real_pair( "CDR_SASA", sasa.first.paratope );
	job->add_string_real_pair( "CDR_SASA_HP", sasa.second.paratope);
	job->add_string_real_pair( "CDR_charge", Real(p_charge.paratope));
}


void AntibodyModelerProtocol::display_constraint_residues( core::pose::Pose & pose ) {

	// Detecting di-sulfide bond

	Size H1_Cys(0), H3_Cys(0);

	if(      pose.residue( pose.pdb_info()->pdb2pose('H',32 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( 'H', 32 );
	} else if( pose.residue( pose.pdb_info()->pdb2pose('H',33 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( 'H', 33 );
	}

	for( Size ii = ab_info_->get_CDR_loop(h3).start(); ii <= ab_info_->get_CDR_loop(h3).stop(); ii++ ) {
		if( pose.residue(ii).name3() == "CYS" ) {
			H3_Cys = ii;
		}
	}

	if( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) ) {
		TR << "CONSTRAINTS: "<< "AtomPair CA " << H1_Cys << " CA " << H3_Cys
		   << " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;
	}

	// Specifying extended kink

	Size hfr_46(0), h3_closest(0);
	hfr_46 = pose.pdb_info()->pdb2pose( 'H', 46 );
	if( ab_info_->get_H3_kink_type() == Extended ) h3_closest = ab_info_->get_CDR_loop(h3).stop() - 5;
	if( h3_closest != 0 ) {
		TR << "CONSTRAINTS: " << "AtomPair CA " << hfr_46 << " CA " << h3_closest
		   << " BOUNDED 6.5 9.1 0.7 DISTANCE; mean 8.0 sd 0.7" << std::endl;
	}

	return;
} // display_constraint_residues


/// @details  Show the complete setup of the antibody modeler protocol
void AntibodyModelerProtocol::show( std::ostream & out ) const {
	/*if ( !flags_and_objects_are_in_sync_ ){
	 sync_objects_with_flags();
	 }*/  // show() should be const
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const AntibodyModelerProtocol & ab_m ) {
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << std::endl;
	out << line_marker << A( 47, "Rosetta 3 Antibody Modeler" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;

	// Display the state of the antibody modeler protocol that will be used
	out << line_marker << "  camelid              : " << ab_m.camelid_                   << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  model_h3  :  "          << ab_m.model_h3_            << std::endl;
	out << line_marker << "         h3_perturb_type     = '"<< ab_m.h3_perturb_type_<<"'"<< std::endl;
	out << line_marker << "         cter_insert         = " << ab_m.cter_insert_         << std::endl;
	out << line_marker << "         h3_filter           = " << ab_m.h3_filter_           << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  snugfit   :  "          << ab_m.snugfit_             << std::endl;
	out << line_marker << "         LH_repulsive_ramp   = " << ab_m.LH_repulsive_ramp_   << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  refine_h3 :  "          << ab_m.refine_h3_           << std::endl;
	out << line_marker << "         h3_refine_type      = '"<< ab_m.h3_refine_type_<<"'" << std::endl;
	out << line_marker << "         h3_filter           = " << ab_m.h3_filter_           << std::endl;
	out << line_marker << "         h3_filter_tolerance = " << ab_m.h3_filter_tolerance_ << std::endl;
	out << line_marker << std::endl;
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}


} // end antibody
} // end protocols
