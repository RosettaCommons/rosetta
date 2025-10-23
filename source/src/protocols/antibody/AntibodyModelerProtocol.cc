// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyModelerProtocol.cc
/// @brief Build a homology model of an antibody
/// @details
/// @author Jianqing Xu ( xubest@gmail.com )


#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/chemical/ResidueType.hh>

#include <ObjexxFCL/format.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyModelerProtocol.hh>
#include <protocols/antibody/AntibodyNumberingConverterMover.hh>
#include <protocols/antibody/CDRsMinPackMin.hh>
#include <protocols/antibody/RefineBetaBarrel.hh>
#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/snugdock/SnugDock.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/jd2/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeler/LoopModeler.hh>

#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue

using namespace ObjexxFCL::format;


static basic::Tracer TR( "protocols.antibody.AntibodyModelerProtocol" );
using namespace core;

namespace protocols {
namespace antibody {

// default constructor
AntibodyModelerProtocol::AntibodyModelerProtocol() : Mover() {
	user_defined_ = false;
	init();
}

// default destructor
AntibodyModelerProtocol::~AntibodyModelerProtocol() = default;

//clone
protocols::moves::MoverOP
AntibodyModelerProtocol::clone() const {
	return( utility::pointer::make_shared< AntibodyModelerProtocol >() );
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
	snugfit_               = true;

	LH_repulsive_ramp_ = true;
	refine_h3_             = true;
	middle_pack_min_       = true;

	extend_h3_before_modeling_ = true;

	benchmark_ = false;
	constrain_vlvh_qq_ = true;
	h3_loop_csts_lr_ = true;
	h3_loop_csts_hr_ = true;
	auto_h3_constraint_ = true;
	packonly_after_graft_=false;

	cst_weight_ = 0.0;
	cen_cst_ = 10.0;
	high_cst_ = 100.0; // if changed here, please change at the end of AntibodyModeler as well

	sc_min_ = false;
	rt_min_ = false;

	cdr_constraint_ = nullptr;
}


void AntibodyModelerProtocol::register_options() {
	using namespace basic::options;

	option.add_relevant( OptionKeys::antibody::model_h3 );
	option.add_relevant( OptionKeys::antibody::snugfit );
	option.add_relevant( OptionKeys::antibody::camelid );
	option.add_relevant( OptionKeys::run::benchmark );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::antibody::constrain_vlvh_qq );
	option.add_relevant( OptionKeys::antibody::h3_loop_csts_lr ); // previously constrain_cter
	option.add_relevant( OptionKeys::antibody::h3_loop_csts_hr ); // previously all_atom_mode_kink_constraint
	option.add_relevant( OptionKeys::antibody::auto_generate_h3_kink_constraint ); // previously auto_generate_kink_constraint
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::antibody::refine_h3 );
	option.add_relevant( OptionKeys::antibody::sc_min);
	option.add_relevant( OptionKeys::antibody::rt_min);
	//option.add_relevant( OptionKeys::antibody::middle_pack_min);
	option.add_relevant( OptionKeys::antibody::extend_h3_before_modeling);
	option.add_relevant( OptionKeys::antibody::packonly_after_graft);
	option.add_relevant( OptionKeys::antibody::run_snugdock);
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
	if ( option[ OptionKeys::run::benchmark ].user() ) {
		set_BenchMark( option[ OptionKeys::run::benchmark ]() );
	}
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		set_cst_weight( option[ OptionKeys::constraints::cst_weight ]() );
	}
	if ( option[ OptionKeys::antibody::h3_loop_csts_lr ].user() ) {
		set_h3_loop_csts_lr( option[ OptionKeys::antibody::h3_loop_csts_lr ]() );
	}
	if ( option[ OptionKeys::antibody::h3_loop_csts_hr ].user() ) {
		set_h3_loop_csts_hr( option[ OptionKeys::antibody::h3_loop_csts_hr ]() );
	}
	if ( option[ OptionKeys::antibody::auto_generate_h3_kink_constraint ].user() ) {
		set_auto_h3_constraint( option[ OptionKeys::antibody::auto_generate_h3_kink_constraint ]() );
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
	if ( option[ OptionKeys::antibody::extend_h3_before_modeling].user()  ) {
		set_extend_h3_before_modeling(option[ OptionKeys::antibody::extend_h3_before_modeling]() );
	}

	run_snugdock_ = option[ OptionKeys::antibody::run_snugdock]();

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
		core::import_pose::pose_from_file( *native_pose, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		set_native_pose( native_pose );
	} else {
		set_native_pose(nullptr);
	}

	TR <<  "Finish Reading and Setting Options !!!" << std::endl;
}


void
AntibodyModelerProtocol::setup_objects() {

	sync_objects_with_flags();

	// setup all the scoring functions
	pack_scorefxn_ = core::scoring::get_score_function(); //JAB - changing this to normal.

	dock_scorefxn_highres_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	dock_scorefxn_highres_->set_weight( core::scoring::chainbreak, 1.0 );
	dock_scorefxn_highres_->set_weight( core::scoring::overlap_chainbreak, 10./3. );

	loop_scorefxn_centroid_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	loop_scorefxn_centroid_->set_weight( scoring::chainbreak, 10./3. );

	loop_scorefxn_highres_ = scoring::get_score_function();
	loop_scorefxn_highres_->set_weight( scoring::chainbreak, 1.0 );
	loop_scorefxn_highres_->set_weight( scoring::overlap_chainbreak, 10./3. );

	// set cst weight, if not already done
	if ( cst_weight_ == 0.0 ) cst_weight_ = 1.0;

	// This adds a new score term (called atom_pair_constraint with weight 1) to the energy function.
	TR << "Scorefunction before: " << *dock_scorefxn_highres_ << std::endl;
	if ( constrain_vlvh_qq_ ) {
		dock_scorefxn_highres_->set_weight( scoring::atom_pair_constraint, cst_weight_ );
	}
	TR << "Scorefunction after: " << *dock_scorefxn_highres_ << std::endl;

	if ( h3_loop_csts_lr_ ) {
		// Always enable constraints in low-resolution mode (i.e. when the sampling is aggressive enough for it to matter)
		loop_scorefxn_centroid_->set_weight( scoring::dihedral_constraint, cst_weight_ );
		loop_scorefxn_centroid_->set_weight( scoring::angle_constraint, cst_weight_ );

		if ( h3_loop_csts_hr_  ) {
			loop_scorefxn_highres_->set_weight( scoring::dihedral_constraint, cst_weight_ );
			loop_scorefxn_highres_->set_weight( scoring::angle_constraint, cst_weight_ );
		}
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
		native_pose = utility::pointer::make_shared< pose::Pose >(pose);
	} else {
		native_pose = utility::pointer::make_shared< pose::Pose >( *get_native_pose() );
	}

	pose::set_ss_from_phipsi( *native_pose ); // JQX: this is the secondary structure from the native pose

	set_native_pose( native_pose ); // pass the native pose to the mover.native_pose_

	ab_info_ = utility::pointer::make_shared< AntibodyInfo >(pose);

	if ( ab_info_->is_camelid() == true ) {
		// no vhvl constrain for camelids
		constrain_vlvh_qq_ = false;
	}

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
	using namespace protocols::analysis;

	// the default inital secstruct is all "L" loop!
	pose::Pose start_pose_ = pose;

	if ( !flags_and_objects_are_in_sync_ ) {
		sync_objects_with_flags();
	}

	finalize_setup(pose);

	basic::prof_reset();
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	pose::set_ss_from_phipsi( pose );

	// Step 1: model the cdr h3 in centroid mode
	// JQX notes: pay attention to the way it treats the stems when extending the loop
	if ( model_h3_ ) {

		// switching to centroid mode
		simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
		simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

		// Building centroid mode loop
		pose::Pose start_pose = pose;
		to_centroid.apply( pose );

		// call ConstraintSetMover
		TR<<"Centroid cst_weight: "<<cst_weight_<<std::endl;
		if (  cst_weight_ != 0.0  ) {
			if ( ! auto_h3_constraint_ ) {
				cdr_constraint_ = utility::pointer::make_shared< constraint_movers::ConstraintSetMover >();
				cdr_constraint_->apply( pose );
			} else {
				// Create constraints on-the-fly here.

				// All of this stuff, and the work that is being done in finalize_setup() for that matter, only needs to happen
				// once per input. I don't trust the 'fresh_instance' and related settings in the other antibody movers well
				// enough to actually rely on that, so I'm going to do this every apply for now.
				antibody::kink_constrain_antibody_H3( pose, ab_info_ );

			}
		}

		// ModelCDRH3-independent approach to loop modeling
		protocols::loop_modeler::LoopModelerOP h3_loop_mover ( new protocols::loop_modeler::LoopModeler() );

		// if fragments are given on the command-line, run fKIC, else just NGK (default)
		if ( basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user() ) {
			h3_loop_mover->setup_kic_with_fragments_config();
		}
		// Fyi, this can also be setup for loophash ...

		// only run centroid
		h3_loop_mover->disable_fullatom_stage();

		// If testing, run reduced cycles
		if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
			h3_loop_mover->centroid_stage()->mark_as_test_run();
		}

		loops::Loop cdr_h3_loop = get_cdr_h3_loop();

		// setting extended to true ensures LoopBuilder will idealize
		// residue bond lengths and angles before closure
		if ( extend_h3_before_modeling_ ) {
			cdr_h3_loop.set_extended( true );
		} else {
			// do not extend loop before modeling
			// do not rebuild loop before modeling
			h3_loop_mover->disable_build_stage();
		}

		// pass loops and score function for both stages
		h3_loop_mover->set_loop( cdr_h3_loop );
		h3_loop_mover->set_cen_scorefxn( loop_scorefxn_centroid_ );

		h3_loop_mover->apply( pose );

		// back to fullatom
		to_full_atom.apply( pose );

		//recover sidechains from starting structures except H3
		utility::vector1<bool> allow_chi_copy( pose.size(), true );
		for ( core::Size ii=cdr_h3_loop.start(); ii<=cdr_h3_loop.stop(); ++ii ) {
			allow_chi_copy[ii] = false;
		}
		protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose, allow_chi_copy );
		recover_sidechains.apply( pose );
	}

	// QQ constraint
	if ( constrain_vlvh_qq_ ) {
		core::Size qq_light_residue ( ab_info_->qq_light_residue(pose) );
		core::Size qq_heavy_residue ( ab_info_->qq_heavy_residue(pose) );

		// Check if both residues are glutamines
		TR << "Light chain residue in question: " << pose.residue_type(qq_light_residue).name() << std::endl;
		TR << "Heavy chain residue in question: " << pose.residue_type(qq_heavy_residue).name() << std::endl;

		std::string gln_type = "GLN";

		if ( pose.residue_type(qq_light_residue).name() == gln_type ) {
			TR << "Got a glutamine in the light chain!" << std::endl;
			if ( pose.residue_type(qq_heavy_residue).name() == gln_type ) {
				TR << "Got a glutamine in the heavy chain as well!" << std::endl;

				// If both are glutamines, make QQ constraint
				antibody::qq_constrain_antibody(pose, qq_heavy_residue, qq_light_residue);
			}
		}
	}


	// Repack loop in full atom
	// call ConstraintSetMover
	TR << "Full-atom cst_weight: " << cst_weight_ << std::endl;
	if (  cst_weight_ != 0.0 && ! auto_h3_constraint_ ) {
		cdr_constraint_ = utility::pointer::make_shared< constraint_movers::ConstraintSetMover >();
		cdr_constraint_->apply( pose );
	}

	CDRsMinPackMinOP cdrs_min_pack_min( new CDRsMinPackMin(ab_info_) );
	if ( sc_min_ ) cdrs_min_pack_min->set_sc_min(true);
	if ( rt_min_ ) cdrs_min_pack_min->set_rt_min(true);
	cdrs_min_pack_min -> set_turnoff_minimization(packonly_after_graft_);
	cdrs_min_pack_min -> apply(pose);

	// Step 2: SnugFit: relieve the clashes between L-H
	//JAB - Why are we refining the LH chain in a class called RefineBetaBarrel?????
	// Seriously.  And we use pre-talaris 2013 to pack the sidechains in it.  Just FYI....

	if ( snugfit_ && ! ab_info_->is_camelid() ) {
		RefineBetaBarrelOP refine_beta_barrel( new RefineBetaBarrel(ab_info_, dock_scorefxn_highres_, pack_scorefxn_) );
		// it has default movemap, tf, and fold_tree
		if ( !LH_repulsive_ramp_ ) {
			refine_beta_barrel-> turn_off_repulsive_ramp();
		}
		if ( sc_min_ )              {
			refine_beta_barrel->set_sc_min(true);
		}
		if ( rt_min_ )              {
			refine_beta_barrel->set_rt_min(true);
		}
		refine_beta_barrel->apply(pose);
		//pose.dump_pdb("2nd_finish_snugfit.pdb");
	}

	// Step 3: Full Atom Relax
	if ( refine_h3_ ) {
		loops::Loop cdr_h3_loop = get_cdr_h3_loop();

		protocols::loop_modeler::LoopModelerOP refine_kic ( new protocols::loop_modeler::LoopModeler() );
		// if fragments are given on the command-line, run fKIC, else just NGK (default)
		if ( basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user() ) {
			refine_kic->setup_kic_with_fragments_config();
		}
		// do not build
		refine_kic->disable_build_stage();
		// pose is fullatom, disable centroid
		refine_kic->disable_centroid_stage();
		// pass cdrh3 loop, defined above
		refine_kic->set_loop( cdr_h3_loop );
		// pass scorefunction from above
		refine_kic->set_fa_scorefxn( loop_scorefxn_highres_ );
		// If testing, run reduced cycles
		if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
			refine_kic->fullatom_stage()->mark_as_test_run();
		}
		// run
		refine_kic->apply( pose );

	}

	//FoldTree, MoveMap, TaskFactory and Variants will be taken care of inside
	cdrs_min_pack_min -> apply(pose);
	//pose.dump_pdb("4th_final_min_pack_min.pdb");

	// Step 4: Store the homolgy models
	pose.fold_tree( * ab_info_->get_FoldTree_AllCDRs(pose) ) ;

	// Redefining CDR H3 cutpoint variants
	// remove this if FT gymnastics are no longer a thing (they shouldn't be)
	loops::remove_cutpoint_variants( pose, true );
	loops::add_cutpoint_variants( pose );

	// Final score (with constraints) before jd2 output the results
	if ( constrain_vlvh_qq_ ) {
		loop_scorefxn_highres_->set_weight( scoring::atom_pair_constraint, cst_weight_ );
	}
	( *loop_scorefxn_highres_ )( pose );

	// Finish
	echo_metrics_to_output(pose);
	set_last_move_status( protocols::moves::MS_SUCCESS );


	if ( basic::options::option[ basic::options::OptionKeys::antibody::output_ab_scheme].user() ) {
		AntibodyNumberingConverterMover converter = AntibodyNumberingConverterMover();
		converter.apply(pose);
	}

	basic::prof_show();
	TR<<"Antibody Modeling Protocol Finished!!!!"<<std::endl<<std::endl<<std::endl;

	if ( ab_info_->antigen_present() ) {

		if ( run_snugdock_ ) {
			core::Real start = pack_scorefxn_->score(pose);
			TR << "Running Snugdock" << std::endl;
			snugdock::SnugDockOP snug( new snugdock::SnugDock() );

			snug->set_scorefxn(pack_scorefxn_);
			snug->set_antibody_info(ab_info_);//Updated info for movemaps, foldtrees, etc.
			snug->apply(pose);

			core::Real post_snugdock = pack_scorefxn_->score(pose);
			TR <<"start:            " << start <<std::endl;
			TR <<"postSD:           " << post_snugdock << std::endl;

		}

		TR << "Running Interface AnalyzerMover" << std::endl;
		core::pose::DockingPartners ab_dock_chains;
		ab_dock_chains.partner1.push_back("A");
		ab_dock_chains.partner2 = ab_info_->get_antibody_chains();
		InterfaceAnalyzerMover analyzer = InterfaceAnalyzerMover(design::get_dock_chains_from_ab_dock_chains(ab_info_, ab_dock_chains), false /* tracer */, pack_scorefxn_, false /* compute_packstat */ , false /* pack_input */,  true /* pack_separated */) ;

		analyzer.apply(pose); //Adds PoseExtraScore_Float to be output into scorefile.
	}

}// end apply


loops::Loop AntibodyModelerProtocol::get_cdr_h3_loop() {
	// get loop from antibody info
	loops::Loop cdr_h3_loop( ab_info_->get_CDR_loop(h3) );
	// Chothia is the default numbering and CDR definition convention
	// under the Chothia definition of the CDR H3, structural divergence occurs
	// at residue 93, not 95 (which is the defined loop start), so we alter this here.
	cdr_h3_loop = loops::Loop( cdr_h3_loop.start() - 2,
		cdr_h3_loop.stop(),
		cdr_h3_loop.cut(),
		0,
		true );
	return cdr_h3_loop;
}

void AntibodyModelerProtocol::echo_metrics_to_output(core::pose::Pose & pose) {

	// align pose to native pose
	pose::Pose native_pose = *get_native_pose();
	antibody::AntibodyInfoOP native_ab_info( new AntibodyInfo(native_pose) );
	align_to_native( pose, native_pose, ab_info_, native_ab_info );

	// the specific constraint terms for output in the log file
	Real atom_pair_constraint_score = pose.energies().total_energies()[ core::scoring::atom_pair_constraint ];
	Real dihedral_constraint_score = pose.energies().total_energies()[ core::scoring::dihedral_constraint ];
	Real angle_constraint_score = pose.energies().total_energies()[ core::scoring::angle_constraint ];
	Real total_score = pose.energies().total_energies()[ core::scoring::total_score ];
	Real unconstrained_score = total_score - atom_pair_constraint_score - dihedral_constraint_score - angle_constraint_score;

	TR << " \t\tatom_pair_constraint_score = " << atom_pair_constraint_score << std::endl;
	TR << "      dihedral_constraint_score = " << dihedral_constraint_score << std::endl;
	TR << "         angle_constraint_score = " << angle_constraint_score << std::endl;
	TR << "                    total_score = " << total_score << std::endl;
	TR << "            unconstrained_score = " << unconstrained_score << std::endl;

	protocols::jd2::add_string_real_pair_to_current_job( "unconstr_score", unconstrained_score );

	align_to_native( pose, native_pose, ab_info_, native_ab_info, "H" );
	protocols::jd2::add_string_real_pair_to_current_job("H3_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h3) ));
	protocols::jd2::add_string_real_pair_to_current_job("H2_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h2) ));
	protocols::jd2::add_string_real_pair_to_current_job("H1_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(h1) ));
	//pose.dump_pdb("aligned_H.pdb");
	if ( ab_info_->is_camelid() == false ) {
		align_to_native( pose, native_pose, ab_info_, native_ab_info, "L" );
		protocols::jd2::add_string_real_pair_to_current_job("L3_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l3) ));
		protocols::jd2::add_string_real_pair_to_current_job("L2_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l2) ));
		protocols::jd2::add_string_real_pair_to_current_job("L1_RMS", global_loop_rmsd( pose, *get_native_pose(), ab_info_->get_CDR_in_loopsop(l1) ));
		//pose.dump_pdb("aligned_L.pdb");
	}

	//protocols::jd2::add_string_real_pair_to_current_job("AP_constraint", atom_pair_constraint_score);
	protocols::jd2::add_string_real_pair_to_current_job("VL_VH_distance", vl_vh_orientation_coords( pose, *ab_info_ )[1]);
	protocols::jd2::add_string_real_pair_to_current_job("VL_VH_opening_angle", vl_vh_orientation_coords( pose, *ab_info_ )[2]);
	protocols::jd2::add_string_real_pair_to_current_job("VL_VH_opposite_opening_angle", vl_vh_orientation_coords( pose, *ab_info_ )[3]);
	protocols::jd2::add_string_real_pair_to_current_job("VL_VH_packing_angle", vl_vh_orientation_coords( pose, *ab_info_ )[4]);

	//kink metrics
	protocols::jd2::add_string_real_pair_to_current_job( "kink_RD_HB", kink_RD_Hbond( pose, *ab_info_ ));
	protocols::jd2::add_string_real_pair_to_current_job( "kink_bb_HB", kink_bb_Hbond( pose, *ab_info_ ));
	protocols::jd2::add_string_real_pair_to_current_job( "kink_Trp_HB", kink_Trp_Hbond( pose, *ab_info_ ));

	std::pair<core::Real,core::Real> q = kink_dihedral( pose, *ab_info_);
	protocols::jd2::add_string_real_pair_to_current_job( "kink_q", q.first );
	protocols::jd2::add_string_real_pair_to_current_job( "kink_qbase", q.second );

	std::pair<ParatopeMetric< core::Real >, ParatopeMetric<core::Real> > sasa = paratope_sasa( pose, *ab_info_);
	ParatopeMetric< core::SSize> p_charge = paratope_charge( pose, *ab_info_ );
	protocols::jd2::add_string_real_pair_to_current_job( "CDR_SASA", sasa.first.paratope );
	protocols::jd2::add_string_real_pair_to_current_job( "CDR_SASA_HP", sasa.second.paratope);
	protocols::jd2::add_string_real_pair_to_current_job( "CDR_charge", Real(p_charge.paratope));
}


void AntibodyModelerProtocol::display_constraint_residues( core::pose::Pose & pose ) {

	// Detecting di-sulfide bond

	core::Size H1_Cys(0), H3_Cys(0);

	if (      pose.residue( pose.pdb_info()->pdb2pose("H",32 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( "H", 32 );
	} else if ( pose.residue( pose.pdb_info()->pdb2pose("H",33 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( "H", 33 );
	}

	for ( core::Size ii = ab_info_->get_CDR_loop(h3).start(); ii <= ab_info_->get_CDR_loop(h3).stop(); ii++ ) {
		if ( pose.residue(ii).name3() == "CYS" ) {
			H3_Cys = ii;
		}
	}

	if ( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) ) {
		TR << "CONSTRAINTS: "<< "AtomPair CA " << H1_Cys << " CA " << H3_Cys
			<< " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;
	}

	// Specifying extended kink

	core::Size hfr_46(0), h3_closest(0);
	hfr_46 = pose.pdb_info()->pdb2pose( "H", 46 );
	if ( ab_info_->get_H3_kink_type() == Extended ) h3_closest = ab_info_->get_CDR_loop(h3).stop() - 5;
	if ( h3_closest != 0 ) {
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
	out << line_marker << std::endl;
	out << line_marker << " ******  model_h3  :  "          << ab_m.model_h3_            << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  snugfit   :  "          << ab_m.snugfit_             << std::endl;
	out << line_marker << "         LH_repulsive_ramp   = " << ab_m.LH_repulsive_ramp_   << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  loop_constraints_lr :  "          << ab_m.h3_loop_csts_lr_           << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  loop_constraints_hr :  "          << ab_m.h3_loop_csts_hr_           << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  h3_kink_constraint :  "          << ab_m.auto_h3_constraint_         << std::endl;
	out << line_marker << std::endl;
	out << line_marker << " ******  vlvh_qq_constraint :  "          << ab_m.constrain_vlvh_qq_         << std::endl;
	out << line_marker << std::endl;
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}


} // end antibody
} // end protocols
