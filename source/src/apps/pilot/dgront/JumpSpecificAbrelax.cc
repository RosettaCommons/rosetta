// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Application level code for relax-type protocols
/// @detailed
///
/// use JumpSpecificAbrelax in the following way:
///
/// AbrelaxAppliaction::register_options();
/// devel::init
/// AbrelaxAppliaction my_app;
/// my_app.add_evaluation( new MySpecialEvaluator );
/// ...
/// my_app.run();
///
/// ---------------------------------------
/// control flow:
/// the run method calls
/// setup()
/// and then either fold() or rerun() (depending on option -rerun )
///
/// each decoy is evaluated by process_decoy() and results are written to the score-file (if specified) or
/// to the silent_output file ( if applicable ).
/// the score file is a silent-file without structural data ( just SCORE lines )
///
/// rerun(): run thru structures in in:file:silent and call process_decoy for each
/// fold(): produce structures with an Abinitio-type protocol and call process_decoy
///
/// options specific to JumpSpecificAbrelax can be found by using -help at the command-line.
/// if you add new options please follow the scheme in the static method register options
///
/// the behaviour of JumpSpecificAbrelax is controlled by comman-line-options. Refer to -help (usage) and the code
///
/// information that is not always present is stored as xxxOP, and the NULL-pointer is interpreted that the respective
/// behaviour is not present. (i.e., native_pose_ is either pointing to the native pose (-native) or to NULL.
/// when you use such pointers ask if they are non-NULL.
///
///
///
/// @author Oliver Lange

// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

// Unit Headers
#include <protocols/abinitio/AbrelaxApplication.hh>

// dgront headers
#include "JumpSpecificAbrelax.hh"
#include "SpecificJumpSetup.cc"
#include "LibraryJumpSetup.cc"


// Package Headers
#include <core/kinematics/util.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/JumpingFoldConstraints.hh>

#include <protocols/abinitio/KinematicTaskControl.hh>
#include <protocols/abinitio/LoopJumpFoldCst.hh>
#include <protocols/abinitio/DoubleLayerKinematicAbinitio.hh>

#include <protocols/abinitio/Templates.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/StrandConstraints.hh>
#include <protocols/abinitio/FragmentMover.hh>

#include <protocols/Protocol.hh>
#include <protocols/relax_protocols.hh>

#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/ResiduePairJumpSetup.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <protocols/jumping/StrandPairing.hh>
#include <protocols/jumping/util.hh>

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/MetricValue.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/conformation/util.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <core/sequence/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#include <protocols/evaluation/PCA.hh>
#include <protocols/evaluation/PoseMetricEvaluator.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/evaluation/EvaluationFactory.hh>

#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/ShortLoopClosure.hh>
#include <protocols/loops/LoopClosure.hh>
#include <protocols/loops/LoopClass.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/util.hh>

//#include <protocols/loops/looprelax_protocols.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/RGFilter.hh>
#include <protocols/simple_filters/COFilter.hh>
#include <protocols/simple_filters/SheetFilter.hh>



//#include <protocols/simple_moves/MinMover.hh>

//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>
// C++ headers
#include <cstdlib>
#include <string>
#include <vector>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>




static basic::Tracer tr("protocols.abinitio.JumpSpecificAbrelax");
static numeric::random::RandomGenerator RG(423465);  // <- Magic number, do not change it!

////////////////////////////////////////////////////////////////////////////////////////////////////
//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, rerun )
OPT_KEY( Boolean, jdist_rerun )
OPT_KEY( Boolean, steal )
OPT_KEY( Boolean, dump_frags )
OPT_KEY( Boolean, start_native )
OPT_KEY( File, pca )
OPT_KEY( File, sf )
OPT_KEY( Boolean, viol )
OPT_KEY( Integer, viol_level )
OPT_KEY( String, viol_type )
OPT_KEY( StringVector, tag_selector )
OPT_KEY( Real, perturb )
OPT_KEY( IntegerVector, rmsd_residues )
OPT_KEY( Boolean, bGDT )
OPT_KEY( Boolean, steal_3mers )
OPT_KEY( Boolean, steal_9mers )

OPT_KEY( Boolean, no_prof_info_in_silentout )

OPT_1GRP_KEY( File,          jumps, pairing_file )
OPT_1GRP_KEY( IntegerVector, jumps, sheets )
OPT_1GRP_KEY( File,          jumps, fix_jumps )
OPT_1GRP_KEY( File,          jumps, jump_lib )
OPT_1GRP_KEY( Boolean,       jumps, fix_chainbreak )
OPT_1GRP_KEY( Boolean,       jumps, evaluate )
OPT_1GRP_KEY( File,          jumps, extra_frags_for_ss )
OPT_1GRP_KEY( File,          jumps, loop_definition_from_file)
OPT_1GRP_KEY( Boolean,       jumps, no_chainbreak_in_relax )
OPT_1GRP_KEY( File,          jumps, residue_pair_jump_file )
OPT_1GRP_KEY( File,          jumps, topology_file )

OPT_1GRP_KEY( IntegerVector,          dgront, jump_cheat )
OPT_1GRP_KEY( IntegerVector,          dgront, jump_from_library )

OPT_1GRP_KEY( Boolean,       loop, close_loops )
OPT_1GRP_KEY( Boolean,       loop, debug_loop_closure )
OPT_1GRP_KEY( Real,          loop, short_frag_cycles )
OPT_1GRP_KEY( Real,          loop, scored_frag_cycles )

OPT_1GRP_KEY( File,          loops, extended_loops )

OPT_1GRP_KEY( Real,       loopfcst, coord_cst_weight )
OPT_1GRP_KEY( Boolean,    loopfcst, coord_cst_all_atom )
OPT_1GRP_KEY( Boolean,    loopfcst, use_general_protocol )
OPT_1GRP_KEY( File,       loopfcst, coord_cst_weight_array )
OPT_1GRP_KEY( File,       loopfcst, dump_coord_cst_weight_array )

OPT_1GRP_KEY( Boolean,       abinitio, relax_with_jumps )
OPT_1GRP_KEY( Boolean,       abinitio, process_store )
OPT_1GRP_KEY( IntegerVector, abinitio, fix_residues_to_native )
OPT_1GRP_KEY( Boolean,        abinitio, return_full_atom )
OPT_1GRP_KEY( Boolean,         abinitio, detect_disulfide_before_relax )

OPT_1GRP_KEY( File,          constraints, forest_file )
OPT_1GRP_KEY( Boolean,       constraints, compute_total_dist_cst )
OPT_1GRP_KEY( IntegerVector,       constraints, evaluate_max_seq_sep )
OPT_1GRP_KEY( Boolean,       constraints, no_linearize_bounded )
OPT_1GRP_KEY( File,          constraints, dump_cst_set )
OPT_1GRP_KEY( Boolean,       constraints, no_cst_in_relax )
OPT_1GRP_KEY( Integer,       constraints, cull_with_native )
OPT_1GRP_KEY( Boolean,        constraints, enable_HA )

OPT_1GRP_KEY( File,         templates, config )
OPT_1GRP_KEY( Boolean,      templates, pairings )
OPT_1GRP_KEY( Integer,      templates, min_nr_large_frags )
OPT_1GRP_KEY( Integer,      templates, min_nr_small_frags )
OPT_1GRP_KEY( Integer,      templates, nr_large_copies )
OPT_1GRP_KEY( Integer,      templates, nr_small_copies )
OPT_1GRP_KEY( Boolean,      templates, vary_frag_size )
OPT_1GRP_KEY( Boolean,      templates, fix_aligned_residues )
OPT_1GRP_KEY( Integer,      templates, fix_margin )
OPT_1GRP_KEY( File,         templates, fix_frag_file )
OPT_1GRP_KEY( Boolean,      templates, no_pick_fragments )
OPT_1GRP_KEY( Boolean,      templates, pick_multiple_sizes )
OPT_1GRP_KEY( Boolean,      templates, strand_constraint )
OPT_1GRP_KEY( Integer,      frags, nr_large_copies )
OPT_1GRP_KEY( Boolean,      frags, annotate )

OPT_1GRP_KEY( Integer,         process, proc_id )
OPT_1GRP_KEY( Integer,         process, nproc )
OPT_1GRP_KEY( Boolean,         process, condor )

////////////////////////////////////////////////////////////////////////////////////////////////////
///@details registering of options that are relevant for JumpSpecificAbrelax
void protocols::abinitio::JumpSpecificAbrelax::register_options() {

	NEW_OPT( dgront::jump_cheat, "read specific jump point from a space-separated values",0);
	NEW_OPT( dgront::jump_from_library, "read specific jump point from a space-separated values",0);

	/// options from the main Option file that are relevant in this context ( and should appear in -help output )
	OPT( in::file::native );
	OPT( in::file::silent ); // input silent file
	OPT( in::file::frag3 );
	OPT( in::file::frag9 );
	OPT( in::file::fasta );
	OPT( in::file::native_exclude_res ); // list for residues to exclude

	OPT( out::file::silent );
	OPT( out::nstruct );

	NEW_OPT( process::proc_id, "give process number... Jobdistributor will only work on proc_id mod nproc part of work ", 0 );
	NEW_OPT( process::nproc, "number of process... needed if proc_id is specified",0);
	NEW_OPT( process::condor,"if condor say yes -- proc_id counting starts at 0",false);


	OPT( abinitio::fastrelax );
	OPT( abinitio::relax );
	OPT( abinitio::multifastrelax );
	NEW_OPT( abinitio::relax_with_jumps, "switch to allow relax even if loops are not closed ", false );
	OPT( abinitio::use_filters );
	NEW_OPT( abinitio::detect_disulfide_before_relax, "run detect_disulfides() before relax", false);
	OPT( abinitio::debug );
	OPT( abinitio::number_3mer_frags );
	OPT( abinitio::number_9mer_frags );
	NEW_OPT( abinitio::process_store, "run process_decoy on each structure in the structure store", false );
	NEW_OPT( abinitio::fix_residues_to_native, "these residues torsions are copied from native and fixed",0);
	NEW_OPT( abinitio::return_full_atom, "return a full-atom structure even if no relax is done", false );
// new options that haven't been defined in global option file

	NEW_OPT( loopfcst::use_general_protocol, "use the new machinery around classes KinematicXXX", false );
	NEW_OPT( loopfcst::coord_cst_weight, "use coord constraints for template", 0.0 );
	NEW_OPT( loopfcst::coord_cst_all_atom, "use coord constraints on all atoms and not just CA", false );
	NEW_OPT( loopfcst::coord_cst_weight_array, "use these weights (per seqpos) for coord cst in rigid regions", "");
	NEW_OPT( loopfcst::dump_coord_cst_weight_array, "dump these weights (per seqpos) for coord cst in rigid regions", "");
	NEW_OPT( rerun, "go through intput structures and evaluate ( pca, rmsd, cst-energy )", false );
	NEW_OPT( jdist_rerun, "go through intput structures and evaluate ( pca, rmsd, cst-energy )", false );

	// use fragments from native structure
	NEW_OPT( steal, "use fragments of native (starting) structure", false );
	NEW_OPT( steal_3mers, "if stealing: use 3mers from native", true );
	NEW_OPT( steal_9mers, "if stealing: use 9mers from native", true );
	NEW_OPT( dump_frags, "for control purposes... dump fragments", false );

	// starting conditions
	NEW_OPT( start_native, "start from extended structure (instead of native)", false );
	NEW_OPT( perturb, "add some perturbation (gaussian) to phi/psi of native", 0.0);

	// evaluation
	NEW_OPT( pca, "compute PCA projections", "");
	NEW_OPT( rmsd_residues, "give start and end residue for rmsd calcul.", -1 );
	NEW_OPT( bGDT, "compute gdtmmm", true );
	NEW_OPT( sf, "filename for score output", "score.fsc" );
	NEW_OPT( no_prof_info_in_silentout, "add column <seconds> to silent output", false );
	// jumping
	NEW_OPT( jumps::fix_jumps, "read jump_file", "" );
	NEW_OPT( jumps::jump_lib, "read jump_library_file for automatic jumps", "" );
	NEW_OPT( jumps::fix_chainbreak, "minimize to fix ccd in re-runs", false );
	NEW_OPT( jumps::pairing_file, "file with pairings", "" );
	NEW_OPT( jumps::sheets, "sheet topology--> replaces -sheet1 -sheet2 ... -sheetN", 1);
	NEW_OPT( jumps::evaluate, "evaluate N-CA-C gemoetry for all jumps in the fold-tree",false);
	NEW_OPT( jumps::extra_frags_for_ss, "use ss-def from this fragset","");
	NEW_OPT( jumps::loop_definition_from_file, "use ss-def from this file","");
	NEW_OPT( jumps::no_chainbreak_in_relax, "don't penalize chainbreak in relax",false );
	NEW_OPT( jumps::residue_pair_jump_file, "a file to define residue pair jump","");
	NEW_OPT( jumps::topology_file, "read a file with topology info ( PairingStats )", "");
	//loop closure
	NEW_OPT( loop::close_loops, "close loops", false );
	NEW_OPT( loop::short_frag_cycles, "cycle-ratio for short_frag_cycles ( loop_size<10 ) after jumping-abinitio",0.2 );
	NEW_OPT( loop::scored_frag_cycles, "cycle-ratio for scored_frag_cycles ( loop_size<10 ) after jumping-abinitio",0.1 );
	NEW_OPT( loop::debug_loop_closure, "dump structures before and after loop closing", false );

	NEW_OPT( loops::extended_loops, "for general protocol: extend structure within these loop-definitions", "" );
	// constraints
	OPT( constraints::cst_file );
	NEW_OPT( constraints::forest_file, "file with constraintforest", "" );
	NEW_OPT( constraints::compute_total_dist_cst, "only relevant for debug: atom_pair_constraints during abinito depends on seq_sep, this computes also the score without regarding seq_sep", false );
	NEW_OPT( constraints::no_linearize_bounded, "don't switch to linearized in BOUNDED func", false );
	NEW_OPT( constraints::dump_cst_set, "dump the cstset_ to file ", "" );
	NEW_OPT( constraints::no_cst_in_relax,"remove constraints for relax", false );
	NEW_OPT( constraints::evaluate_max_seq_sep, "evaluate constraints to this seq-sep [vector]", 0);
	NEW_OPT( constraints::cull_with_native, "if option is set all constraints that violate the native structure with more than X are thrown out! ", 1 );
	NEW_OPT( constraints::enable_HA, "enable constraints to the HA-atom in centroid mode", false);

	NEW_OPT( viol, "show violations", false );
	NEW_OPT( viol_level, "how much detail for violation output", 1 );
	NEW_OPT( viol_type, "work only on these types of constraints", "");
	NEW_OPT( tag_selector, "work only on these tag(s)","");

	// homologs
	NEW_OPT( templates::config,"read a list of templates and alignments","templates.dat");
	NEW_OPT( templates::pairings, "use pairings from templates",false);

	//large default number means all frags are used	if this option is not specified
	NEW_OPT( templates::min_nr_large_frags, "how many large fragments should be present", 100000 );
	NEW_OPT( templates::min_nr_small_frags, "how many small fragments should be present", 100000 );

	NEW_OPT( templates::nr_large_copies, "make N copies of each picked template fragment -- a hacky way to weight them", 4 );
	NEW_OPT( templates::nr_small_copies, "make N copies of each picked template fragment -- a hacky way to weight them", 20 );
	NEW_OPT( templates::vary_frag_size, "pick fragments as long as aligned regions", false );
	NEW_OPT( templates::fix_aligned_residues, "pick only from template fragments and then keep these residues fixed", false );
	NEW_OPT( templates::fix_margin, "keep n residues at edges of fixed fragments moveable", 1 );
	NEW_OPT( templates::fix_frag_file," fragments from this file are picked once in beginning and then kept fixed", "");
	NEW_OPT( templates::no_pick_fragments, "no further fragment picking from templates", false);
	NEW_OPT( templates::pick_multiple_sizes,"pick 9mers, 18mers and 27mers",false);
	NEW_OPT( templates::strand_constraint,"use the template-based strand-constraints",false);
	//	NEW_OPT( homologs::steal,"steal fragments from homologs", true);

	NEW_OPT( frags::nr_large_copies, "make N copies for each standard 9mer (or so) fragment",1 );
	NEW_OPT( frags::annotate, "read the annotation from the rosetta++ fragment file", false );


	// generalized protocol:
	OPT( loops::loop_file );


	// this adds all relevant options from the protocols
	// ClassicAbinitio, FoldConstraints, JumpingFoldConstraints
	KinematicAbinitio::register_options();
	Templates::register_options();
	// here we should have
	// ClassicRelax::register_options();
	// FastRelax::register_options(); etc.
}
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace abinitio {

using core::Size;
using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;
using namespace jumping;
using namespace evaluation;
using namespace basic::options;
//using namespace basic::options::OptionKeys;

// little helper classes for evaluation of generated decoys:
////////////////////////////////////////////////////////////////////////////////////////////////////
// evaluates the PCA
class PcaEvaluator : public PoseEvaluator {
public:
	PcaEvaluator ( PCA_OP pca ) : pca_( pca ) {}
	virtual void apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss ) const;
private:
	PCA_OP pca_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// evaluates the violations of atom-pairconstraints always with the full weight and full sequence separation
using namespace scoring::constraints;
class ShowViolation : public PoseEvaluator {
public:
	ShowViolation( ) : constraints_( NULL ) {}
	virtual void apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss ) const;
private:
	mutable ConstraintSetOP constraints_;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
class ComputeTotalDistCst : public PoseEvaluator {
public:
	ComputeTotalDistCst( ) : constraints_( NULL ) {};
	virtual void apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss ) const;
private:
	mutable ConstraintSetOP constraints_;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail c'stor - nothing special
JumpSpecificAbrelax::JumpSpecificAbrelax() :
	silent_score_file_( NULL ),
	native_pose_( NULL ),
	pca_( NULL ),
	bRelax_ ( false ),
	cstset_( NULL ),
	jump_def_ ( NULL ),
	templates_( NULL ),
	fragset_large_( NULL ),
	fragset_small_( NULL ),
	evaluator_ ( new MetaPoseEvaluator ),
	abrelax_checkpoints_( "Abrelax" )
{}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail add a PoseEvaluator derived instance for decoy-processing
void JumpSpecificAbrelax::add_evaluation( PoseEvaluatorOP eval ) {
	evaluator_->add_evaluation( eval );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@details setup of Application data that is used for both, fold() and run()
/// this is mainly stuff for scoring and evaluation ( process_decoys(), evaluator_ )
void JumpSpecificAbrelax::setup() {
	using namespace basic::options::OptionKeys;

	if ( option[ constraints::no_linearize_bounded ] ) {
		tr.Info << "use fully harmonic potential for BOUNDED " << std::endl;
		ConstraintIO::get_func_factory().add_type("BOUNDED", new BoundFunc(0,0,0,1000,"dummy") );
	}

	silent_score_file_ = new io::silent::SilentFileData;
	silent_score_file_-> set_filename( std::string( option[ sf ]()  ) );

	// read native pose
	if ( option[ in::file::native ].user() ) {
		native_pose_ = new pose::Pose;
		core::import_pose::pose_from_pdb( *native_pose_, option[ in::file::native ]() );

		pose::set_ss_from_phipsi( *native_pose_ );

#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( *native_pose_ );
#endif

		core::util::switch_to_residue_type_set( *native_pose_, chemical::CENTROID ); //so that in do_rerun the native pose is the same as the other poses
		//		for ( Size i = 1; i<=native_pose_->total_residue(); i++ ) {
		//			std::cout << native_pose_->phi(i) << ' ' << native_pose_->psi(i) << std::endl;
		//		}
	}

	// specify sequence -- from fasta file or native_pose
	if ( option[ in::file::fasta ].user() ) {
		sequence_ = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1].sequence();
		tr.Info << "read fasta sequence: " << sequence_.size() << " residues\n"  << sequence_ << std::endl;
	} else if ( native_pose_ ) {
		sequence_ = native_pose_->sequence();
	} else if ( !option[ rerun ]() ) { // if we rerun we don't need sequence or native or anything...
		utility_exit_with_message(
			"Error: can't read sequence! Use -in::file::fasta sequence.fasta or -in::file::native native.pdb!"
		);
	}

	std::string native_tag = "";

	// run with homolog info? -- needed for setup_fragments, and rerun keep it upfront
	setup_templates();

	// set rmsd native
	if ( native_pose_ ) {
		if ( option[ in::file::native_exclude_res ].user() ) {
			add_evaluation( new SelectRmsdEvaluator(
						native_pose_,
						invert_exclude_residues( native_pose_->total_residue(), option[ in::file::native_exclude_res ]()),
						native_tag)
			);
			if ( option[ bGDT ]() ) {
				add_evaluation( new SelectGdtEvaluator(
						native_pose_,
						invert_exclude_residues( native_pose_->total_residue(), option[ in::file::native_exclude_res ]()),
						native_tag)
				);
			}
		} else if ( option[ rmsd_residues ].user() ){
			core::Size start = option[ rmsd_residues ]()[ 1 ];
			Size end = option[ rmsd_residues ]()[ 2 ];
			add_evaluation( new RmsdEvaluator( native_pose_, start, end,  native_tag, option[ bGDT ]() ) );
		} else {
			add_evaluation( new RmsdEvaluator( native_pose_, native_tag, option[ bGDT ]() ) );
		}
	} // if ( native_pose_ )

	//add command-line evaluator stuff
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);

	core::pose::metrics::PoseMetricCalculatorOP
		clash_calculator = new protocols::toolbox::pose_metric_calculators::ClashCountCalculator( 2.0 );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "clashes", clash_calculator );
	add_evaluation( new simple_filters::PoseMetricEvaluator<core::Size>( "clashes", "total" ) );
	add_evaluation( new simple_filters::PoseMetricEvaluator<core::Size>( "clashes", "bb" ) );

	if ( option[ viol ]() ) add_evaluation( new ShowViolation );
	if ( option[ constraints::compute_total_dist_cst ] ) add_evaluation( new ComputeTotalDistCst );
	// read PCA info
	if ( option[ pca ].user() ) {
		pca_ = new PCA;
		pca_->read_eigvec_file( option[ pca ](), *native_pose_, 2 );
		// if ( tr.Trace.visible() ) pca_->show( std::cout );
		add_evaluation( new PcaEvaluator( pca_ ) );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// ///@detail minimizes pose with chainbreak score switched on
// /// this reduces the score immensely if decoys are read from silent file --> somehow the
// /// silent-file reader fucks up the chainbreaks... ( not noticable rmsd-wise )
// /// it is probably due to random placement of the variant atoms -> short minimization and they sit in the right place
// /// only small improvement when called after fold()
// void JumpSpecificAbrelax::fix_chainbreaks( pose::Pose &pose ) {
// 	using namespace kinematics;
// 	using namespace basic::options::OptionKeys;

// 	scorefxn_->set_weight( scoring::linear_chainbreak, 1.0 );
// 	FoldTree const &f ( pose.fold_tree() );

// 	// make a MoveMap
// 	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
// 	// disallow bb moves
// 	movemap->set_bb( false );

// 	for ( int i = 1; i<=f.num_cutpoint(); i++ ) {
// 		Size cut = f.cutpoint( i );
// 		movemap->set_bb( cut, true );
// 		movemap->set_bb( cut+1, true );
// 		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cut );
// 		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cut+1 );
// 	}

// 	protocols::simple_moves::MinMoverOP min_move_ = new protocols::simple_moves::MinMover;
// 	min_move_->movemap( movemap );
// 	min_move_->min_type( "dfpmin" );
// 	//get currently used score_function...
// 	min_move_->score_function( scorefxn_ );
// 	min_move_->apply( pose );
// }
////////////////////////////////////////////////////////////////////////////////////////////////////
bool JumpSpecificAbrelax::close_loops( pose::Pose &pose, scoring::ScoreFunctionOP scorefxn, std::string const& tag /*for checkpoints*/ ) {
	if ( !fragset_small_ ) {
		setup_fragments();
	}

	using namespace basic::options::OptionKeys;
	if ( option[ loop::debug_loop_closure ]() ) pose.dump_pdb(tag+"_pre_closure.pdb");

	bool success( false );
	if ( !(success=abrelax_checkpoints_.recover_checkpoint( pose, tag, "close_loops")) ) {

		// make a MoveMap ... could be coming from somewhere else, though
		kinematics::MoveMapOP movemap = new kinematics::MoveMap;
		movemap->set_bb( true );

		// a weird bug occurs if we make a copy of a copy of the pose and set the new fold-tree
		// the behaviour is very different from setting the same fold-tree into the copy of the pose
		loops::SlidingWindowLoopClosureOP closure_protocol =
			new loops::SlidingWindowLoopClosure( fragset_small_, scorefxn, movemap );
		// set options here if you like
		// closure_protocol-> ... and write the setters/accessors, too, if you have to
		closure_protocol->scored_frag_cycle_ratio( option[ loop::scored_frag_cycles ]() );
		closure_protocol->short_frag_cycle_ratio( option[ loop::short_frag_cycles ]() );

		success = jumping::close_chainbreaks( closure_protocol, pose, abrelax_checkpoints_ , tag );
		//core::Real

		if ( option[ loop::debug_loop_closure ]() ) pose.dump_pdb(tag+"_post_closure.pdb");

		if ( !success ) tr.Warning << "WARNING: no success in close_loops()" << std::endl;
		if ( success ) abrelax_checkpoints_.checkpoint( pose, tag, "close_loops");
	}
	return success;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail run all evaluations on the decoy from this function
/// if you want these evaluations also available during internal stages of the protocols -->
///     put them into a PoseEvaluator and use add_evaluation in setup()
/// otherwise you can also use "manual" code right here in process_decoy --> this will only appear in final
/// silent_out and silent_score - files.
void JumpSpecificAbrelax::process_decoy(
	 pose::Pose &pose,
	 scoring::ScoreFunction const& scorefxn,
	 std::string tag,
	 io::silent::ProteinSilentStruct &pss ) const
{
	using namespace basic::options::OptionKeys;
	//pose.dump_pdb("test.pdb");

	// would like to put the following two also in an PoseEvaluator
	// ScoreEvaluator
	// StructureDumper
	scorefxn( pose );
	// if ( tr.Info.visible() ) scorefxn_->show( std::cout, pose );

	pss.fill_struct( pose, tag );

	// run PoseEvaluators
	evaluator_->apply( pose, tag, pss );

	if ( option[ jumps::evaluate ]() ) {
		if ( !native_pose_ ) utility_exit_with_message(" to evaluate jumps you need to specify a native structure ");
		evaluation::MetaPoseEvaluator eval_jumps;
		native_pose_->fold_tree( pose.fold_tree() );
		for ( Size nj = 1; nj<= pose.num_jump(); ++nj ) {
			eval_jumps.add_evaluation( new simple_filters::JumpEvaluator( *native_pose_, nj) );
		}
		eval_jumps.apply( pose, tag, pss );
	}

} // process_decoy
////////////////////////////////////////////////////////////////////////////////////////////////////
void JumpSpecificAbrelax::initialize_constraint_forest( pose::Pose & pose ) {
	using namespace basic::options::OptionKeys;
	if ( option[ constraints::forest_file ].user() ) {
		if( ! constraint_forest_ ) {
			tr.Info << "read ConstraintForest... : " << std::endl;
			constraint_forest_ = ConstraintIO::get_instance()->
				read_constraint_forest( option[ constraints::forest_file ](), pose );
		}
		constraint_forest_->generate_random_sample();
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail  read constraints file (once) and constraints_set to the pose (each call)
void JumpSpecificAbrelax::add_constraints( pose::Pose & pose ) {
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	bool bFirst( !cstset_ );
	if ( option[ constraints::enable_HA ] ) {
		for ( Size i = 1; i<= pose.total_residue(); i++ ) {
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CENTROID_HA, i );
		}
	}

	if ( bFirst ) {
		if ( option[ constraints::cst_file ].user() ) {
				// reads and sets constraints
				tr.Info << "read constraints... : " << std::endl;
				cstset_ = ConstraintIO::get_instance()->read_constraints( core::scoring::constraints::get_cst_file_option(), new ConstraintSet, pose );
	}
	}
	if ( constraint_forest_ ) {
		if ( !cstset_ ) cstset_ = new ConstraintSet;
		cstset_->add_constraints( constraint_forest_->get_constraints() );
	}
	if ( bFirst && templates_ ) {
		if ( !cstset_ ) cstset_ = new ConstraintSet;
		templates_->add_target_constraints( cstset_, pose );
		if ( option[ templates::strand_constraint ] ) {
			ConstraintCOPs my_strand_cst;
			if ( templates_ ) {
				my_strand_cst = StrandConstraints( templates_->strand_pairing_stats() ).build_constraints( pose );
			} else if ( option[ jumps::topology_file ].user() ) {
				utility::io::izstream is( option[ jumps::topology_file ] );
				if ( is.good() ) {
					PairingStatisticsOP ps = new PairingStatistics;
					is >> *ps;
					tr.Info << *ps << std::endl;
					my_strand_cst = StrandConstraints( *ps ).build_constraints( pose );
				} else {
					utility_exit_with_message(" did not find topology_file: " + std::string( option[ jumps::topology_file ]() ) );
				}
			} else {
					utility_exit_with_message(" strand_constraint nees a topology info: either via templates or -topology_file ");
			}
			cstset_->add_constraints( my_strand_cst );
			add_evaluation( new evaluation::ConstraintEvaluator( "strand", my_strand_cst ) );

			if ( native_pose_ ) {//just a temporary hack to test the StrandConstraint
				pose::Pose test_pose = *native_pose_;
				test_pose.add_constraints( my_strand_cst );

				if ( option[ constraints::dump_cst_set ].user() ) {
					tr.Info << "dump strand constraints to file..." << std::endl;
					utility::io::ozstream dump_cst( "STRAND_CST_DUMP" );
					test_pose.constraint_set()->show_definition( dump_cst, test_pose );
				}

				scoring::ScoreFunction cst_score;
				cst_score.set_weight( scoring::atom_pair_constraint, 1.0 );
				cst_score( test_pose );
				tr.Info << " native pose yields this score for the StrandConstraints: " << cst_score( test_pose ) << std::endl;
				cst_score.show( tr, test_pose );
				test_pose.constraint_set()->show_violations( tr , test_pose, 120);
			}
		}
	}

	if ( option[ constraints::cull_with_native ].user() && native_pose_ ) {
		tr.Warning << "************************************************************************************\n"
						<< "*********************  CULL CONSTRAINTS WITH NATIVE STRUCTURE *********************\n"
						<< "************************************************************************************\n" << std::endl;
		ConstraintCOPs filtered;
		core::scoring::constraints::cull_violators( cstset_->get_all_constraints(),
			filtered, *native_pose_, option[ constraints::cull_with_native ]() );
		cstset_=new ConstraintSet;
		cstset_->add_constraints( filtered );
	}

	pose.constraint_set( cstset_ );

	if ( option[ constraints::dump_cst_set ].user() ) {
		tr.Info << "dump constraints to file..." << std::endl;
		utility::io::ozstream dump_cst( option[ constraints::dump_cst_set ]() );
		cstset_->show_definition( dump_cst, pose );
	}

	if ( option[ constraints::evaluate_max_seq_sep ].user() ) {
		Size const neval ( option[ constraints::evaluate_max_seq_sep ]().size() );
		for ( Size i = 1; i<= neval; i++ ) {
			Size const seq_sep( option[ constraints::evaluate_max_seq_sep ]()[ i ] );
			add_evaluation( new evaluation::ConstraintEvaluator( "seq_sep_"+utility::to_string( seq_sep) , *cstset_, 1, seq_sep ) );
		}
	}

} // add_constraints( pose::Pose & pose )
////////////////////////////////////////////////////////////////////////////////////////////////////
class Stage1Sampler : public ClassicAbinitio {
public:
	Stage1Sampler(
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
	) : ClassicAbinitio( fragset_large, fragset_large, movemap ) {};

	Stage1Sampler( protocols::abinitio::FragmentMoverOP brute_move_large	)
	: ClassicAbinitio( brute_move_large, brute_move_large, brute_move_large, 1 /*dummy*/ ) {};

	virtual void apply( core::pose::Pose &pose );
};

void Stage1Sampler::apply( core::pose::Pose &pose ) {
	prepare_stage1( pose );
	do_stage1_cycles( pose );
	//return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void JumpSpecificAbrelax::insert_template_frags( core::pose::Pose &pose, kinematics::MoveMapOP movemap, std::string tag ) const {
	using namespace basic::options::OptionKeys;
	if ( option[ templates::fix_frag_file ].user() ) {
		FrameList fix_frames;
		fragment::FragmentIO().read( option[ templates::fix_frag_file ](), fix_frames );
		Size const frame_id ( static_cast< int >( RG.uniform() * fix_frames.size() ) + 1 );
		FrameOP frame( fix_frames[ frame_id ] );
		Size const frag_id ( static_cast< int >( RG.uniform() * frame->nr_frags() ) + 1 );
		frame->apply( frag_id, pose );

		std::ofstream out( "big_frags.log", std::ios_base::out | std::ios_base::app );
		out << tag << " " << RJ(10,frame->start()) << RJ( 10, frame->stop() ) << RJ( 10, frag_id ) << std::endl;

		movemap->set_bb( true );
		Size const npadding( option[ OptionKeys::templates::fix_margin ] );
		for ( Size pos = frame->start() + npadding; pos<=frame->end() - npadding; ++pos ) {
			movemap->set_bb( pos, false );
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail loop over structures in silent-input file
/// small trick is used to also have native structure in the set of analysis:
/// it is added to the collection of silent_file-structures manually
/// TODO we need to do something about difference between fullatom and centroid input!
void JumpSpecificAbrelax::do_rerun() {
	using namespace core;
	using namespace io::silent;
	using namespace pose;
	using namespace basic::options::OptionKeys;

	core::io::silent::SilentFileDataOP outsfd( NULL );
	if ( option[ out::file::silent ].user() ) {
		outsfd = new	core::io::silent::SilentFileData();
	}

	scoring::ScoreFunctionOP scorefxn( NULL );
	if ( option[ in::file::silent ].user() ) {
		//read silent file for input
		SilentFileData sfd;
		sfd.read_file( *(option [ in::file::silent ]().begin()) );

		// run thru all structures
		Size ct ( 0 );
		for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
			Pose pose;
			std::string tag = it->decoy_tag();
			if ( option[ tag_selector ].user() == 0 || std::find( option[ tag_selector ]().begin(), option[ tag_selector ]().end(), tag ) != option[ tag_selector ]().end() ) {
				if ( option[ in::file::fullatom ].user() ) {
						it->fill_pose( pose,
							option[ in::file::fullatom ] ?
							*(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD )) :
							*(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) ));
					}	else {
						it->fill_pose( pose );
					}
					//      if ( tag == "NATIVE" ) pose = *native_pose_; // replace structure with NATIVE so that we don't suffer from non-idealized stuff
				add_constraints( pose );
				scorefxn = generate_scorefxn( pose.is_fullatom() );
				//screen output
				if ( sfd.size() < 10 )	{
					tr.Info << tag << " " << std::endl;
				} else {
					if ( (ct++ % 50) == 0 ) {
						std::cout << ".";
						std::cout.flush();
					}
				}

				// set score terms
				//		if ( pose.fold_tree().num_cutpoint() != 0 ) {
					//  activate the CUTPOINTS by setting the pose variants -- is done in SilentFile now...
					//					JumpSample jumps( pose.fold_tree() );
					//					jumps.add_chainbreaks( pose );
					// and set score term
					scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
					scorefxn->set_weight( scoring::overlap_chainbreak, 1.0 );
					//}


				if ( option[ loop::close_loops ] ) {
					add_evaluation( new RmsdEvaluator( new pose::Pose( pose ), std::string("closure"), option[ bGDT ]() ) );
					close_loops( pose, scorefxn, tag );
				}

				basic::MetricValue<core::Size> mr;
				pose.metric("clashes","total",mr);
				tr.Info << "Total clashes " << mr.value() << std::endl;

				bool passes_filters = check_filters( pose );
				if( !passes_filters ) {
					tag = "F_"+tag.substr(2);
				}

				ProteinSilentStruct ss;
				process_decoy( pose, *scorefxn, tag, ss );
				// write this to score-file if applicable
				if ( outsfd ) outsfd->add_structure( new ProteinSilentStruct( ss ) );

				//remove closure-rmsd
				if ( option[ loop::close_loops ] ) {
					evaluator_->pop_back();
				}

			}

		}
	}

	// add native structure to the list
	if ( native_pose_ ) {
		ProteinSilentStruct ss;
		add_constraints( *native_pose_ );
		scorefxn = generate_scorefxn( false /*full_atom*/ );
		scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
		scorefxn->set_weight( scoring::overlap_chainbreak, 1.0 );

		if ( option[ loop::close_loops ] ) { //otherwise the column (needed for non-native decoys) doesn't show up in score-file
			add_evaluation( new RmsdEvaluator( new pose::Pose( *native_pose_ ), std::string("closure"), option[ bGDT ]() ) );
		}

		process_decoy( *native_pose_, *scorefxn,  "NATIVE", ss );
		// write this to score-file if applicable
		if ( silent_score_file_ ) {
			silent_score_file_ -> write_silent_struct( ss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
		}

		if ( outsfd ) outsfd->add_structure( new ProteinSilentStruct( ss ) );

		if ( option[ loop::close_loops ] ) evaluator_->pop_back();
	}

	if ( silent_score_file_ ) {
		outsfd->write_all( silent_score_file_->filename(), true /* bWriteScoresOnly */ );
	}

	if ( outsfd ) outsfd->write_all( option[ out::file::silent ]() );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail loop over structures in silent-input file
/// small trick is used to also have native structure in the set of analysis:
/// it is added to the collection of silent_file-structures manually
/// TODO we need to do something about difference between fullatom and centroid input!
void JumpSpecificAbrelax::do_distributed_rerun() {
	using namespace core;
	using namespace io::silent;
	using namespace pose;
	using namespace basic::options::OptionKeys;

	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;

	//read silent file for input
	bRelax_ = option[ OptionKeys::abinitio::relax ]() ||
		option[ OptionKeys::abinitio::fastrelax ]();

	// setup profiling
	evaluation::TimeEvaluatorOP run_time( NULL );
	if ( !option[ no_prof_info_in_silentout ] ) {
		add_evaluation( run_time = new evaluation::TimeEvaluator ); //just don't use this in integration tests!
	}

	if ( option[  OptionKeys::loops::loop_file ].user() ) {
		std::string filename( protocols::loops::get_loop_file_name() );
		loops_in_.read_loop_file( filename );  // <== TODO: select these using density score
	}


	// get input tags
	SilentFileData sfd;
	typedef utility::vector1< std::string > TagList;
	TagList input_tags;
	input_tags = sfd.read_tags_fast( *(option [ in::file::silent ]().begin())  );

	// read silent data
	sfd.read_file( *(option [ in::file::silent ]().begin()) );

	// determine nstruct
	int const nstruct = std::max( 1, option [ out::nstruct ]() );

	// create jobs
	typedef utility::vector1< BasicJobOP > JobList;
	JobList input_jobs;
	for ( TagList::const_iterator it = input_tags.begin(), eit = input_tags.end(); it!=eit; ++it ) {
		BasicJobOP job = new BasicJob( *it, "rerun", nstruct);
		input_jobs.push_back( job );
	}

	// setup JobDistributor
	PlainSilentFileJobDistributor< BasicJobOP > jobdist( input_jobs );
	int const procid ( option[ process::proc_id ] + ( option[ process::condor ] ? 1 : 0 ) );
	if ( procid > option[ process::nproc ] ) {
		utility_exit_with_message("procid to large " + ObjexxFCL::string_of( procid ) + " run only " + ObjexxFCL::string_of( option[ process::nproc ] ) + " processes");
	}
	jobdist.set_proc_id( procid, option[ process::nproc ] );
	jobdist.startup();

	// production loop
	bool bEndrun = false;
	BasicJobOP curr_job;
	int curr_nstruct;
	while ( jobdist.next_job(curr_job, curr_nstruct) && !bEndrun ) {
		if ( run_time ) run_time->reset(); //reset clock of TimeEvaluator
		tr.Info << "Starting " << jobdist.get_current_output_tag() << " ..." << std::endl;
		tr.Info << "read " << curr_job->input_tag() << "..." << std::endl;
		Pose pose;
		//		ResidueTypeSetOP residue_types;
		//		if ( option[ in:file:fullatom ] ) {
		//		}
		sfd.get_structure( curr_job->input_tag() ).fill_pose( pose );/*,
					*chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )
					);*/
		std::string tag = jobdist.get_current_output_tag();

		if( option[ OptionKeys::loops::random_grow_loops_by ].user() ){
			protocols::loops::LoopMover loopmover( loops_in_ );
			loopmover.grow_all_loops( pose, option[ OptionKeys::loops::random_grow_loops_by ]() );
			loops_in_ = loopmover.loops();
			tr.Info << "Enlarged loops: " << std::endl;
			tr.Info << loops_in_ << std::endl;

		};

		add_constraints( pose ); // needs to come before scorefxn setup to know if constraints are present
		scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false /*fullatom*/ ) );
		scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true /*fullatom*/ ) );

		centroid_scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
		centroid_scorefxn->set_weight( scoring::overlap_chainbreak, 1.0 );

		if ( loops_in_.size() ) {
			utility::vector1< core::Real > vecs;
			loops::Loops rigid( loops_in_.invert( pose.total_residue() ) );
			loops::fix_with_coord_cst( rigid, pose, option[ loopfcst::coord_cst_all_atom ], vecs );
		}

		if ( option[ jumps::no_chainbreak_in_relax ] ) {
			fullatom_scorefxn->set_weight( scoring::linear_chainbreak, 0.0 );
			fullatom_scorefxn->set_weight( scoring::overlap_chainbreak, 0.0 );
		}
		// set score function for processing/relaxing stage

		tr.Info << tag << " " << std::endl;

		bool loop_closure_failed( false );
		if ( option[ loop::close_loops ] ) {
			add_evaluation( new RmsdEvaluator( new pose::Pose( pose ), std::string("closure"), option[ bGDT ]() ) );
			loop_closure_failed = !close_loops( pose, centroid_scorefxn, tag );
		}

		bool passes_filters = check_filters( pose );
		// run relax if applicable
		// don't relax if we failed filters or loop_closing, or if option[ relax_with_jumps ] is true
		bool bCanRelax = passes_filters && ( !loop_closure_failed || option[ OptionKeys::abinitio::relax_with_jumps ]() );
		if ( bRelax_ ) {
			if ( !pose.is_fullatom() ) core::util::switch_to_residue_type_set( pose, chemical::FA_STANDARD );
			if ( bCanRelax ) {
				relax( pose, fullatom_scorefxn, jobdist.get_current_output_tag() );
			} else { //cannot relax
				//need proper atom set to score with full-atom
				//				if ( !pose.is_fullatom() ) core::util::switch_to_residue_type_set( pose, chemical::FA_STANDARD );
			}
		} // if ( bRelax_ )

		// process decoy if this hasn't happened yet
			// analyze result
		std::string output_tag = jobdist.get_current_output_tag();
		if ( !passes_filters  && loop_closure_failed ) {
			output_tag = "X_"+output_tag.substr(2);
		} else if( loop_closure_failed ) {
			output_tag = "C_"+output_tag.substr(2);
		} else if( !passes_filters ) {
			output_tag = "F_"+output_tag.substr(2);
		}

		ProteinSilentStruct pss;
		process_decoy( pose, pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, output_tag, pss );
		// write this to score-file if applicable
		if ( silent_score_file_ ) {
			silent_score_file_ -> write_silent_struct( pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
		}

		// write to silent file
		jobdist.dump_silent( curr_nstruct, pss );
		//remove closure-rmsd
		if ( option[ loop::close_loops ] ) {
			evaluator_->pop_back();
		}

	} // end of production loop
	jobdist.shutdown();
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail called by setup_fold() if option[ start_native ] is active
/// the routine defines a fragment of the length of the structure
/// steals the fragment from the native and applies it to the decoy
/// native needs to be idealized!
void JumpSpecificAbrelax::copy_native_structure( core::pose::Pose & extended_pose ) const {
	// requires that the sequences match at the beginning (1..nmatch_res) -- > use sequence alignment later
	tr.Info << " *** use native structure as starting template -- NEEDS TO BE IDEALIZED !!! *** \n";
	copy_structure( extended_pose, *native_pose_ );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void JumpSpecificAbrelax::copy_structure( core::pose::Pose & extended_pose, core::pose::Pose & desired_pose ) const {
	// requires that the sequences match at the beginning (1..nmatch_res) -- > use sequence alignment later
	tr.Info << " *** use native structure as starting template -- NEEDS TO BE IDEALIZED !!! *** \n";
	// determine length of segment to copy from native
	Size seg_len = std::min(extended_pose.total_residue(), desired_pose.total_residue() );
	// chu workaround when folding with ligand/metal
	Size protein_len = 0;
	for ( Size i = 1; i <= seg_len; ++i ) {
		if( extended_pose.residue(i).is_protein() && desired_pose.residue(i).is_protein() ) {
			protein_len ++;
		}
	}
	seg_len = protein_len;
	fragment::Frame long_frame(1, seg_len);

	//create apropriate length FragData object
	FragData frag( new BBTorsionSRFD, seg_len );
	//	for ( Size pos = 1; pos<= seg_len; pos++ ) {
	//		frag.add_residue( new BBTorsionSRFD( nbb, native_pose_->secstruct( pos ), 'X' ) );
	//	}

	// get torsion angles from native pose
	frag.steal( desired_pose, long_frame );

	// apply native torsions to extended structue
	frag.apply( extended_pose, long_frame );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail called by setup_fold(): setup the decoy pose with correct target sequence and extended structure
///
void JumpSpecificAbrelax::generate_extended_pose( core::pose::Pose &extended_pose, std::string const& sequence ) const {

	core::pose::make_pose_from_sequence(
		 extended_pose,
		 sequence,
		 *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
	 );

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose.total_residue(); pos++ ) {
		if ( ! extended_pose.residue(pos).is_protein() ) continue;
		extended_pose.set_phi( pos, -150 );
		extended_pose.set_psi( pos, 150);
		extended_pose.set_omega( pos, 180 );
	}
	tr.Debug 	<< "CHECK PHI/PSI/OMEGA pos 8: " << extended_pose.phi( 8 ) <<  " "
						<< extended_pose.psi( 8 ) << " " << extended_pose.omega( 8 )
						<< std::endl;

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( extended_pose );
#endif

}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail called by setup_fold(): read fragment libraries, I strongly suggest to use different options than A and B
/// if option[ steal ] fragments from the native structure are added to the set.
/// native structure needs to be idealized for this!
void JumpSpecificAbrelax::setup_fragments() {// FragSetOP& fragsetA, FragSetOP& fragsetB ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string frag_large_file, frag_small_file;
	if (option[ in::file::fragA ].user()) {
		frag_large_file  = option[ in::file::fragA ]();
	} else {
		frag_large_file  = option[ in::file::frag9 ]();
	}

	if (option[ in::file::fragB ].user()) {
		frag_small_file  = option[ in::file::fragB ]();
	} else {
		frag_small_file  = option[ in::file::frag3 ]();
	}


//   ConstantLengthFragSetOP _fragset_large = new ConstantLengthFragSet;
//   ConstantLengthFragSetOP _fragset_small = new ConstantLengthFragSet;
//   _fragset_large->read_fragment_file( frag_large_file, option[ OptionKeys::abinitio::number_9mer_frags ], option[ frags::nr_large_copies ], option[ frags::annotate ] );
//   _fragset_small->read_fragment_file( frag_small_file, option[ OptionKeys::abinitio::number_3mer_frags ] );

// 	//copy pointers to member_variables
// 	fragset_large_ = _fragset_large;
// 	fragset_small_ = _fragset_small;

//	fragset_large_ = FragmentIO().read( frag_large_file );
	fragset_large_ = FragmentIO(
		option[ OptionKeys::abinitio::number_9mer_frags ](),
		option[ OptionKeys::frags::nr_large_copies ](),
		option[ OptionKeys::frags::annotate ]()
	).read( frag_large_file );

	fragset_small_ = FragmentIO(
		option[ OptionKeys::abinitio::number_3mer_frags ],
		1, //nr_copies
		option[ OptionKeys::frags::annotate ]
	).read( frag_small_file );

	if ( templates_ && !option[ templates::no_pick_fragments ]() ) {
		if ( option[ templates::vary_frag_size ] ) {
			fragset_templates_ = new OrderedFragSet;
			templates_->pick_large_frags( *fragset_templates_, new BBTorsionSRFD, option[ templates::nr_large_copies ] );
			tr.Info << " merge template frags with standard library " << std::endl;
			fragset_large_ = merge_frags(
																	 *fragset_templates_,
																	 *fragset_large_,
																	 option[ templates::min_nr_large_frags ],
																	 true /* random selection of fill frags */
			);
		} else { // use old-way of picking:
			//pick torsion fragments fragset_large
			tr.Info << "pick large fragments as 9mers " << std::endl;
			if ( option[ templates::min_nr_large_frags ].user() ) {
				Size const min_nr_frags( option[ templates::min_nr_large_frags ] );
				Size const nr_large_copies( option[ templates::nr_large_copies ] );
				fragset_large_ = templates_->pick_frags(
								fragset_large_,
								new FragData( new BBTorsionSRFD, fragset_large_->max_frag_length() ),
								min_nr_frags,
								nr_large_copies );
			} else {
				Size nr = templates_->pick_frags( *fragset_large_, new FragData( new BBTorsionSRFD, fragset_large_->max_frag_length() ) );
				tr.Info << nr << " " << fragset_large_->max_frag_length() << "mer fragments picked from homolog structures" << std::endl;
			}
			if ( option[ templates::pick_multiple_sizes ] ) {
				Size nr = templates_->pick_frags(
								*fragset_large_,
								new FragData( new BBTorsionSRFD, 18 )
				);
				tr.Info << nr << " 18mer fragments picked from homolog structures" << std::endl;
				nr = templates_->pick_frags(
								*fragset_large_,
								new FragData( new BBTorsionSRFD, 24 )
				);
				tr.Info << nr << " 27mer fragments picked from homolog structures" << std::endl;
			}
		} // !vary_frag_size

		if ( option[ templates::min_nr_small_frags ].user() ) {
			Size const min_nr_frags( option[ templates::min_nr_small_frags ] );
			Size const nr_small_copies( option[ templates::nr_small_copies ] );
			fragset_small_ = templates_->pick_frags(
								fragset_small_,
								 new FragData( new BBTorsionSRFD, fragset_small_->max_frag_length() ),
								 min_nr_frags,
								 nr_small_copies );
		} else {
			//pick torsion fragments fragset_small
			Size nr2 = templates_->pick_frags( *fragset_small_, new FragData( new BBTorsionSRFD, fragset_small_->max_frag_length() ) );
			tr.Info << nr2 << " " << fragset_small_->max_frag_length() << "mer fragments picked from homolog structures" << std::endl;
		}
	} // templates && !templates:no_pick_fragments

	if ( option[ steal ]() && native_pose_ && ( option[ steal_3mers ]() || option[ steal_9mers ]() )) {
		tr.Info << " stealing fragments from native pose: ATTENTION: native pose has to be IDEALIZED!!! " << std::endl;
		//		utility_exit_with_message(" stealing fragments from pose: currently not supported! ask Oliver " );
		if ( option[ steal_9mers ]() ) steal_frag_set_from_pose( *native_pose_, *fragset_large_,
			new FragData( new BBTorsionSRFD, fragset_large_->max_frag_length() ) );
		if ( option[ steal_3mers ]() ) steal_frag_set_from_pose( *native_pose_, *fragset_small_,
			new FragData( new BBTorsionSRFD, fragset_small_->max_frag_length() ) );
	} else if ( option[ steal ]() && !native_pose_ && !templates_ ) {
		tr.Warning << "cannot steal fragments without native pose or homologue structures " << std::endl;
	}

	if ( option[ dump_frags ]() ) { //diagnosis
		utility::io::ozstream dump_frag_small( "fragset_small.dump" );
		for ( FrameIterator it=fragset_small_->begin(), eit=fragset_small_->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_small );
		}
		utility::io::ozstream dump_frag_large( "fragset_large.dump" );
		for ( FrameIterator it=fragset_large_->begin(), eit=fragset_large_->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_large );
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail called by setup_fold(). Read template definitions
void JumpSpecificAbrelax::setup_templates() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool const bTemplates( option[ templates::config ].user() );
	if ( !bTemplates ) { // jump-out if not used
		if ( option[ templates::pairings ].user() )
			tr.Warning << "option templates:pairings ignored... specify templates:config!" << std::endl;
		//...
		return;
	}

	if ( native_pose_ ) tr.Info << "native strand pairings " << jumping::StrandPairingSet( *native_pose_ );
	templates_ = new Templates( option[ templates::config ], native_pose_ );
	templates_->target_sequence() = sequence_; // a hack until class SequenceMapping works better
	// want to pick fragments from templates... make sure they are not initialized yet
	assert( !fragset_large_ );

	if( !templates_->is_good() ){
		utility_exit_with_message("ERRORS occured during template setup. check BAD_SEQUENCES file!");
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail called by setup_fold(). Read jump definitions / barcodes (not yet) etc.
/// if jump_def_ points to an object we will use JumpFoldConstraint-protocol in fold()
void JumpSpecificAbrelax::setup_jumps(	pose::Pose const& extended_pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// setup jumps

	bool bDoubleDef = false;

	jump_def_ = NULL;

	ss_def_ = new jumping::SecondaryStructure( *fragset_small_, false /*no JustUseCentralResidue */ );
	if ( option [ jumps::extra_frags_for_ss ].user() ) {
		FragSetOP ss_frags = FragmentIO().read( option[ jumps::extra_frags_for_ss ]() );
		ss_def_ = new jumping::SecondaryStructure( *ss_frags, false );
	}
	if ( option[ jumps::loop_definition_from_file ].user() ) {
		ss_def_ = new jumping::SecondaryStructure();
		ss_def_->read_from_file( option[ jumps::loop_definition_from_file ]() );
	}

	if ( option[ jumps::fix_jumps ].user() ) {
	  JumpSetup *ptr = new JumpSetup( extended_pose.total_residue() );
	  ptr->read_file( option[ jumps::fix_jumps ]() );
		// initialize jumping
		jumping::JumpSample current_jumps( ptr->create_jump_sample() );
		if ( native_pose_ && !current_jumps.has_orientation_and_pleating() ) {
			tr.Warning << "abinitio:CHEAT JumpingFoldConstraints takes orienation and pleating from native structure !!!" << std::endl;
			current_jumps.steal_orientation_and_pleating( *native_pose_ );
			ptr->set_jump_sample( current_jumps );
		}
	  jump_def_ = ptr;
	}
	if ( option[ dgront::jump_from_library ].user() ) {
		utility::vector1< int> const& jump_point( option[ dgront::jump_from_library ]() );
		std::cout << "abinitio:dgront:jump_from_library Jumping parameters are: [";
		for(core::Size i=1;i<=jump_point.size();i++)
			std::cout<<jump_point[i]<<" ";
		  std::cout<<"]"<<std::endl;
		  if(jump_point.size()==4) {
		    BaseJumpSetup *ptr = new LibraryJumpSetup( extended_pose.total_residue(), ss_def_, native_pose_,jump_point[1], jump_point[2],jump_point[3], jump_point[4]);
		    std::cout << "abinitio:dgront:jump_from_library Jumping position is: " << jump_point[1]<< jump_point[2]<<jump_point[3]<< jump_point[4]<<std::endl;
		    jump_def_ = ptr;
		  }
		  else {
			  int nJumps = (int)(jump_point.size() / 4);
			  int* ir = new int[nJumps];
			  int* jr = new int[nJumps];
			  int* o = new int[nJumps];
			  int* p = new int[nJumps];
			  for(int i=0;i<nJumps;i++) {
				  ir[i] = jump_point[i*4+1];
				  jr[i] = jump_point[i*4+2];
				  o[i] = jump_point[i*4+3];
				  p[i] = jump_point[i*4+4];
			  }
			  BaseJumpSetup *ptr = new LibraryJumpSetup( extended_pose.total_residue(), ss_def_, native_pose_,nJumps,ir, jr,o, p);
			  std::cout << "abinitio:dgront:jump_from_library Jumping with " << nJumps<< " jumps"<<std::endl;
			  jump_def_ = ptr;
			  delete ir;
			  delete jr;
			  delete o;
			  delete p;
		  }
	}

	if ( option[ jumps::jump_lib ].user() ) {
		bDoubleDef = jump_def_ != 0;
		JumpSelector *ptr = new JumpSelector( native_pose_->secstruct() );
		ptr->read_file( option[ jumps::jump_lib ] );
		jump_def_ = ptr;
	}
	if ( option[ jumps::sheets ].user() ) {
		bDoubleDef = jump_def_ != 0;

		// get sheet-topology
		jumping::SheetBuilder::SheetTopology sheets = option[ jumps::sheets ]();

		// get secondary structure info
		assert( fragset_small_ );

		ss_def_->show( tr.Trace );
		// get pairings file
		PairingsList pairings;
		if ( option[ jumps::pairing_file ].user() )
			read_pairing_list( option[ jumps::pairing_file ](), pairings );
		else if ( option[ constraints::forest_file ].user() )
			pairings = constraint_forest_->get_pairings();
			// perhaps get_all_pairings()?

		// done: instantiate sheet-builder
		jump_def_ = new SheetBuilder( ss_def_, pairings, sheets );
	}
	if ( option[ jumps::topology_file ].user() ) {
		utility::io::izstream is( option[ jumps::topology_file ] );
		if ( !is.good() ) {
			utility_exit_with_message(" did not find topology_file: " + std::string( option[ jumps::topology_file ]() ) );
		}
		PairingStatisticsOP ps = new PairingStatistics;
		is >> *ps;
		tr.Info << *ps << std::endl;
		jumping::PairingList helix_pairings; //empty for now
		jump_def_ = new TemplateJumpSetup( NULL, ss_def_, ps, helix_pairings );
	}
	if ( option[ constraints::forest_file ].user() ) {
		if( jump_def_ == 0 && constraint_forest_->get_pairings().size() > 0 ) { //constraint_forest pairings but no -sheet option
			jump_def_ = new JumpsFromConstraintForest(extended_pose.total_residue(),
																							 constraint_forest_,
																							 ss_def_->loop_fraction());
//			jump_def_	= new JumpsFromAllPairings(
//				extended_pose.total_residue(),
//				constraint_forest_->get_pairings(),
//				ss_def_->loop_fraction()
//			);
		}
	}
	if ( option[ templates::pairings ] ) {
		bDoubleDef = false;
		if ( option[ jumps::fix_jumps ].user() ) {
			tr.Info << "use fixed jumps but take jump-geometries from template! " << std::endl;
			jump_def_ = new FixTemplateJumpSetup( *templates_->create_jump_def( ss_def_ ), jump_def_ );
		} else {
			bDoubleDef = jump_def_ != 0;
			jump_def_ = templates_->create_jump_def( ss_def_ );
		}
		//		templates_->read_pairings( option[ templates::pairings ] );

// 		assert( native_pose_ );
// 		SecondaryStructureOP ss_contr = new jumping::SecondaryStructure( *native_pose_ );
// 		utility::io::ozstream dump("ss_def__from_frags");
// 		for ( Size i = 1; i<=ss_def_->total_residue(); i++ ) {
// 			dump << i << " " << ss_def_->loop_fraction()(i) << " " << ss_contr->loop_fraction()(i) << std::endl;
// 		}

		// no argument means ss_def_ is taken from homologues --> loop_fraction() controls cutpoint placement


		if ( option[ constraints::forest_file ].user() ) utility_exit_with_message("can't use constraint-forest pairings with template pairings yet");
	}

	if ( option[ jumps::residue_pair_jump_file ].user() ) {
		bDoubleDef = jump_def_ != 0;
		ResiduePairJumpSetup * ptr = new ResiduePairJumpSetup( extended_pose.total_residue() );
		ptr->read_file( option[ jumps::residue_pair_jump_file ]() );
		jump_def_ = ptr;
	}

	if ( bDoubleDef ) {
		utility_exit_with_message("you can only define one jump mode: choose one of -fix_jumps / -jump_lib / -sheets / -pairings");
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail setup_fold() all initialization that is necessary to run abinitio production loop in fold()
/// read fragments, make pose from sequence, get constraints, set jumps, movemap .etc
/// the two parameters are OUTPUT:
///    extended_pose to run A) with ( might actually contain native starting structure (option!) )
///    prot_ptr: an initialized instance of ClassicAbinitio, FoldConstraints or JumpingFoldConstraints depending on
///   user settings.
void JumpSpecificAbrelax::setup_fold( 	pose::Pose& extended_pose, ProtocolOP& prot_ptr ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// ==========================================================================
	///  --------- fold()-specific setup --------------------------
	// ==========================================================================

	// ---------------------------------------------------------------------------------------------------------------
	// initialize pose
	generate_extended_pose( extended_pose, sequence_ );
	if ( option[ start_native ]() ) {
		copy_native_structure( extended_pose );
	} else if ( option[ in::file::s ].user() ) {
		core::pose::PoseOP tmp_pose( new core::pose::Pose );
		std::string fn = option[ in::file::s ](1);
		core::import_pose::pose_from_pdb( *tmp_pose, fn );
		copy_structure( extended_pose, *tmp_pose );
	}
	// ---------------------------------------------------------------------------------------------------------------
	bRelax_ = option[ OptionKeys::abinitio::relax ]() ||
		option[ OptionKeys::abinitio::fastrelax ]() ||
		option[ OptionKeys::abinitio::multifastrelax ]();

	// FRAGMENTS:
	// stores in FragSetOP fragset_large_, fragset_small_;
	setup_fragments();

	// read in constraint_forest, generate random sample
	initialize_constraint_forest( extended_pose );

	// add constraints if available
	add_constraints( extended_pose );

	// setup jumping... evtl. needs fragset to determine sheet/loop-fractions..
	setup_jumps( extended_pose );

	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	// allow bb moves
	movemap->set_bb( true );
	if ( option[ OptionKeys::abinitio::fix_residues_to_native ].user() ) {
		utility::vector1< int> const& fix_start_ends( option[ OptionKeys::abinitio::fix_residues_to_native ]() );
		for ( Size i=1; i + 1 <= fix_start_ends.size(); i+=2 ) {
			Size const start( fix_start_ends[ i ]);
			Size const end( fix_start_ends[ i + 1 ]);
			if ( !(end >= start) ) utility_exit_with_message("end < start in abinitio:fix_residues_to_native");
			fragment::Frame long_frame(start, end-start+1 );
			//create apropriate length FragData object
			FragData frag( new BBTorsionSRFD, end-start+1 );

			// get torsion angles from native pose
			frag.steal( *native_pose_, long_frame );

			// apply native torsions to extended structue
			frag.apply( extended_pose, long_frame );

			for ( Size pos = start; pos <= end; pos++ ) {
				movemap->set_bb( pos, false);
			}
		}
	}

	if ( option[ OptionKeys::loopfcst::use_general_protocol ] ) {
		// parse loops file Loops
		//utility::vector1< protocols::Loop > loops_in, loops;
		if ( option[  OptionKeys::loops::loop_file ].user() ) {
			std::string filename( protocols::loops::get_loop_file_name() );
			loops_in_.read_loop_file( filename );  // <== TODO: select these using density score
		}
		core::chemical::ResidueTypeSetCAP rsd_set;
		// if full-atom load starting structure as full-atom to recover sidechains later
		if ( option[ OptionKeys::loops::input_pdb ].user() ) {
			if ( option[ in::file::fullatom ]() ) {
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
				core::import_pose::pose_from_pdb( extended_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name() );
				if ( !extended_pose.is_fullatom() ) utility_exit_with_message(" this full-atom pose should be a full-atom pose, no? ");
			} else {
				// centroid starting structure for loop-modeling
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
				core::import_pose::pose_from_pdb( extended_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name() );
				if ( extended_pose.is_fullatom() ) utility_exit_with_message(" this centroid pose should not be a full-atom pose, no? ");
			}
			if( option[ OptionKeys::loops::random_grow_loops_by ].user() ){
				protocols::loops::LoopMover loopmover( loops_in_ );
				loopmover.grow_all_loops( extended_pose, option[ OptionKeys::loops::random_grow_loops_by ]() );
				loops_in_ = loopmover.loops();
				tr.Info << "Enlarged loops: " << std::endl;
				tr.Info << loops_in_ << std::endl;
			};
		}
		add_constraints( extended_pose );
		core::kinematics::simple_visualize_fold_tree( extended_pose.fold_tree(), tr.Debug );
		KinematicAbinitioOP sampler = new KinematicAbinitio( fragset_small_, fragset_large_, movemap /*this movemap will be ignored*/ );
		ResolutionSwitcher res_switch(
																	extended_pose,
																	extended_pose.is_fullatom(),
																	sampler->start_from_centroid(),
																	sampler->return_centroid()
		);
		if ( native_pose_ ) sampler->set_native_pose( native_pose_ );
		sampler->set_show_viol_level( option[ viol_level ] );
		sampler->init( res_switch.start_pose() );

		if ( option[ loop::close_loops ]() ) {
			loops::SlidingWindowLoopClosureOP closure_protocol = new loops::SlidingWindowLoopClosure;
			// set options here if you like
			// closure_protocol-> ... and write the setters/accessors, too, if you have to
			closure_protocol->scored_frag_cycle_ratio( option[ loop::scored_frag_cycles ]() );
			closure_protocol->short_frag_cycle_ratio( option[ loop::short_frag_cycles ]() );
			sampler->closure_protocol( closure_protocol );
		}

		LoopJumpFoldCstOP controller;
		if ( option[  OptionKeys::loops::extended_loops ].user() ) {
			loops::Loops extended_loops_in;
			std::string filename( option[ OptionKeys::loops::extended_loops ]().name() );
			extended_loops_in.read_loop_file( filename );  // <== TODO: select these using density score

			KinematicAbinitioOP stage1_sampler = new KinematicAbinitio( *sampler );
			sampler->bSkipStage1_ = true;
	//	if ( sampling_protocol_ ) sampling_protocol_->init(  res_switch.start_pose() ); this overwrites all settings with defaults...

			stage1_sampler->init( res_switch.start_pose() ); //sets default options
			stage1_sampler->bSkipStage3_ = true;
			stage1_sampler->bSkipStage4_ = true;
			stage1_sampler->closure_protocol( NULL );
			loops::Loops rigid_core( loops_in_.invert( extended_pose.total_residue() ) );
			controller = new DoubleLayerKinematicAbinitio(
						 jump_def_,
						 extended_loops_in,
						 rigid_core,
						 sampler,
						 stage1_sampler,
						 ss_def_,
						 option[ loopfcst::coord_cst_weight ],
						 option[ loopfcst::coord_cst_all_atom ]
			);
		} else {
			controller = new LoopJumpFoldCst(
							 jump_def_,
							 loops_in_,
							 sampler,
							 ss_def_,
							 option[ loopfcst::coord_cst_weight ],
							 option[ loopfcst::coord_cst_all_atom ]
			);
		}

		controller->set_input_pose_is_fa( extended_pose.is_fullatom() );
		prot_ptr = controller;
		if ( evaluator_->size() ) sampler->set_evaluation( evaluator_ );

		if ( option[ loopfcst::coord_cst_weight_array ].user() ) {
			utility::io::izstream file( option[ loopfcst::coord_cst_weight_array ]() );
			utility::vector1< core::Real > weights;
			read_vector( file, weights );
			controller->set_coord_cst_weight_array( weights );
		}

		if ( option[ loopfcst::dump_coord_cst_weight_array ].user() ) {
			controller->set_dump_weights_file( option[ loopfcst::dump_coord_cst_weight_array ] );
		}


	} else {
		// setup abinitio protocol: one of either, ClassicAbinitio/ FoldConstraints/ JumpingFoldConstraints
		if ( jump_def_ ) {
			tr.Info << "run JumpingFoldConstraints....." << std::endl;
			// it doesn't matter if we have no constraints the extra FoldConstraints part in the Jumping protocl
			// won't do anything
			JumpingFoldConstraintsWrapper* pp;
			pp = new JumpingFoldConstraintsWrapper( fragset_small_, fragset_large_, movemap, jump_def_ );
			if ( native_pose_ ) pp->set_native_pose( native_pose_ ); //to steal native jumps
			pp->set_show_viol_level( option[ viol_level ] );
			prot_ptr = pp;
		}	else {
			if ( extended_pose.constraint_set()->has_residue_pair_constraints() ) {
				// We have constraints: run xxxFoldConstraints
				tr.Info << "run FoldConstraints....." << std::endl;
				FoldConstraints* pp;
				pp = new FoldConstraints( fragset_small_, fragset_large_, movemap );
				pp->set_show_viol_level( option[ viol_level ] );
				prot_ptr = pp;
			} else {
				/// no constraints ---> ClassicAbinitio
				tr.Info << "run ClassicAbinitio....." << std::endl;
				prot_ptr = new ClassicAbinitio( fragset_small_, fragset_large_, movemap );
			}
		}
	}

	Protocol& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer

	/// initialize protocol
	abinitio_protocol.init( extended_pose );
	abinitio_protocol.return_centroid( !(bRelax_ || option[ OptionKeys::abinitio::return_full_atom ]) );
	if ( evaluator_->size() ) abinitio_protocol.set_evaluation( evaluator_ );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
bool JumpSpecificAbrelax::check_filters( core::pose::Pose & pose ) {
	using namespace protocols::filters;

	// return true if we're not supposed to use filters
	if ( !basic::options::option[ basic::options::OptionKeys::abinitio::use_filters ]() ) return true;
	if ( option[ basic::options::OptionKeys::filters::disable_all_filters ]() ) return true; //makes a lot of sense IMHO

	// apply RG, contact-order and sheet filters
	protocols::simple_filters::RGFilter    rg_filter;
	protocols::simple_filters::COFilter    co_filter;
	protocols::simple_filters::SheetFilter sh_filter;

	if ( !option[ basic::options::OptionKeys::filters::disable_rg_filter ]() && !rg_filter.apply( pose) ) return false;
	if ( !option[ basic::options::OptionKeys::filters::disable_co_filter ]() && !co_filter.apply( pose) ) return false;
	if ( !option[ basic::options::OptionKeys::filters::disable_sheet_filter ]() && !sh_filter.apply( pose) ) return false;
	tr.Info << " passed all filters " << std::endl;
	return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
scoring::ScoreFunctionOP JumpSpecificAbrelax::generate_scorefxn( bool fullatom ) {
	scoring::ScoreFunctionOP scorefxn( NULL );
	if ( fullatom ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	}
	if ( jump_def_  && !option[ OptionKeys::jumps::no_chainbreak_in_relax ] ) {
		scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
		scorefxn->set_weight( scoring::overlap_chainbreak, 1.0 );
	}
	if ( cstset_ && cstset_->has_residue_pair_constraints()  ) {
		scorefxn->set_weight( scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn->set_weight( scoring::angle_constraint, option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn->set_weight( scoring::dihedral_constraint, option[ OptionKeys::constraints::cst_weight ]() );
	}
	if ( option[ OptionKeys::loopfcst::coord_cst_weight ].user() ) {
		scorefxn->set_weight( scoring::coordinate_constraint, option[ OptionKeys::loopfcst::coord_cst_weight ]);
	}
	return scorefxn;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail everything happens in fold()!
/// setup of stuff that is not needed for rerun()
///    read fragments
///    [ optional ] steal fragments ( take fragments from native pose )
///
void JumpSpecificAbrelax::fold() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;

	// setup JobDistributor stuff
	utility::vector1< BasicJobOP > input_jobs;

	int const nstruct = std::max( 1, option [ out::nstruct ]() );
	BasicJobOP job = new BasicJob("" /*no input tag*/, "abinitio_relax", nstruct);
	input_jobs.push_back( job );
	PlainSilentFileJobDistributor< BasicJobOP > jobdist( input_jobs );
	BasicJobOP curr_job;
	int curr_nstruct;
	jobdist.startup();
	bool bEndrun = false;


	// setup pose and abinitio
	ProtocolOP prot_ptr;
	pose::Pose init_pose;
	setup_fold( init_pose, prot_ptr ); //init_pose may or may not be fullatom (depending on flag option[ in:file:fullatom] )
	Protocol& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer

	// setup scorefunctions
	// this is called in setup_fold: add_constraints( extended_pose ); //such that scorefxn setup knows about constraints...
	scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false /*fullatom*/ ) );
	scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true /*fullatom*/ ) );
	abinitio_protocol.set_fullatom_scorefxn( fullatom_scorefxn );
	abinitio_protocol.set_centroid_scorefxn( centroid_scorefxn );

	evaluation::TimeEvaluatorOP run_time( NULL );
	if ( !option[ no_prof_info_in_silentout ]() ) {
		add_evaluation( run_time = new evaluation::TimeEvaluator ); //just don't use this in integration tests!
		abinitio_protocol.set_evaluation( evaluator_ );
	}

	// production loop
	while ( jobdist.next_job(curr_job, curr_nstruct) && !bEndrun ) {
		time_t pdb_start_time = time(NULL);
		if ( run_time ) run_time->reset(); //reset clock of TimeEvaluator

		// retrieve starting pose
		pose::Pose fold_pose ( init_pose );

		//need to save RG states such that choices for constraints and fold-tree are the same.
		// in this case a successful recover doesn't mean that we want to skip folding..
		// rather there will be a checkpoint later on containing the structure to go on: e.g., checkpoint in stage3
		abrelax_checkpoints_.recover_checkpoint( fold_pose, jobdist.get_current_output_tag(), "abrelax_rg_state");
		abrelax_checkpoints_.checkpoint( fold_pose, jobdist.get_current_output_tag(), "abrelax_rg_state");

		// perturb phi/psi randomly -- should be different each run
		if ( option[ perturb ].user() ) {
			Real sig = option[ perturb ];
			for ( Size pos = 1; pos <= fold_pose.total_residue(); pos++ ) {
				fold_pose.set_phi( pos, fold_pose.phi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_psi( pos, fold_pose.psi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_omega( pos, fold_pose.omega( pos ) );
			}
		}

		// run abinitio
		abinitio_protocol.set_current_tag( jobdist.get_current_output_tag() );
		tr.Debug << "fold_pose is " << (fold_pose.is_fullatom() ? " fullatom " : " centroid " ) << "before protocol run "<<std::endl;
		abinitio_protocol.apply( fold_pose );
		bool success = true; // this code is depraciated -- !!!
		std::string output_tag = abinitio_protocol.get_current_tag();
		// if ( tr.Info.visible() ) fold_pose.constraint_set()->show_violations( std::cout, fold_pose, option[ viol_level ] );

		// close loops if wished so
		bool loop_closure_failed( false ); //!fold_pose.fold_tree().num_cutpoint() );
		if ( success && !option[ OptionKeys::loopfcst::use_general_protocol ] &&  option[ loop::close_loops ]() ) {
			loop_closure_failed = !close_loops( fold_pose, centroid_scorefxn, jobdist.get_current_output_tag() );
		}

		if ( bRelax_ && option [ OptionKeys::abinitio::debug ] ) {
			io::silent::SilentFileData outsfd;
			std::string silent_file = option[ basic::options::OptionKeys::out::file::silent ]() + "_" + "before_relax";

			io::silent::ProteinSilentStruct pss;
			process_decoy( fold_pose, fold_pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, jobdist.get_current_output_tag(), pss );
			outsfd.write_silent_struct( pss, silent_file );
		}

		// run relax if applicable. Use filters to decide which structures to relax and which not to relax!
		bool passes_filters = check_filters( fold_pose );
		bool bProcessDecoy( true );
		// don't relax if we failed filters or loop_closing, or if option[ relax_with_jumps ] is true
		bool bCanRelax = success && passes_filters && ( !loop_closure_failed || option[ OptionKeys::abinitio::relax_with_jumps ]() );
		if ( bRelax_ ) {
				if ( !option[ OptionKeys::loopfcst::use_general_protocol ] ) {
					//					utility_exit_with_message("Cannot proceed (about to crash), because pose is fullatom and loopcose expects centroid. Use -use_general_protocol if you want to postrelax foldcst models.");
				//This part of code i think is obsolete and will soon disappear.");
					if ( !fold_pose.is_fullatom() ) core::util::switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD );
				}
				if ( bCanRelax ) {
					if ( option[ basic::options::OptionKeys::abinitio::multifastrelax ]() ) {
						bEndrun = multi_fast_relax( abinitio_protocol, fullatom_scorefxn, jobdist, curr_nstruct, curr_job );
						bProcessDecoy = false;
						if ( bEndrun ) break;
					} else {
						fold_pose.constraint_set( NULL );
						//if ( abinitio_protocol.return_centroid() ) core::util::switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD );
						if ( !fold_pose.is_fullatom() ) core::util::switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD );
						relax( fold_pose, fullatom_scorefxn, jobdist.get_current_output_tag() );
					}
				} else { //cannot relax
					//need proper atom set to score with full-atom
					if ( !fold_pose.is_fullatom() ) {
						core::util::switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD );
					}
					(*fullatom_scorefxn)( fold_pose );
					if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {	//FastRelax adds another two columns, grr
						relax::FastRelax::setPoseExtraScore( fold_pose );
					} else {
						relax::ClassicRelax().setPoseExtraScore( fold_pose ); // so we the same number of columns
					}
				}
		} // if ( bRelax_ )

		// process decoy if this hasn't happened yet
		if ( bProcessDecoy ) {
			// analyze result
			io::silent::SilentFileData outsfd;

			//make sure that number of columns does not change -- ever
			outsfd.strict_column_mode( true );

			io::silent::ProteinSilentStruct pss;


			// abinitio produces n_stored structures -- the last one is the same as the final structure n_stored.
			Size n_stored( abinitio_protocol.structure_store().size() );
			std::string new_output_tag ( output_tag );
			if ( option[ OptionKeys::abinitio::process_store ] ) new_output_tag = "S" + string_of( n_stored ) + "_" + output_tag.substr(2);

			if ( !passes_filters  && loop_closure_failed ) {
				new_output_tag = "X_"+new_output_tag.substr(2);
			} else if( loop_closure_failed ) {
				new_output_tag = "C_"+new_output_tag.substr(2);
			} else if( !passes_filters ) {
				new_output_tag = "F_"+new_output_tag.substr(2);
			}

			// write to silent file
			process_decoy( fold_pose, fold_pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, new_output_tag, pss );
			if ( option[ OptionKeys::abinitio::process_store ] ) pss.add_energy ( "store", n_stored );
			// write this to score-file if applicable
			if ( silent_score_file_ ) {
				silent_score_file_ -> write_silent_struct( pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
			}

			if ( option[ OptionKeys::out::pdb ] ) fold_pose.dump_pdb( std::string(option[ OptionKeys::out::path::path ]())  + "/" + new_output_tag + ".pdb");
			outsfd.add_structure( pss );

			if ( option[ OptionKeys::abinitio::process_store ] ) {
				tr.Info <<" storing intermediate stage3 and stage4 structures: these are currently not relaxed "<< std::endl;
				Size ct( 1 );
				for ( Protocol::StructureStore::iterator it = abinitio_protocol.structure_store().begin(),
								eit = abinitio_protocol.structure_store().end(); ct < n_stored && it!=eit; ++it,++ct ) { //13th structure is identical to fold_pose...
					io::silent::ProteinSilentStruct pss;
					std::string extra_tag( "S" );
					std::string current_tag( extra_tag + string_of( ct ) + "_" + output_tag.substr(2) );

					bool loop_closure_failed( false );
					if ( option[ loop::close_loops ]() ) {
						if ( loop_closure_failed=!close_loops( *it, centroid_scorefxn,  extra_tag + string_of( ct ) + "_" + output_tag ) ) {
							extra_tag = "C";
						}
					}
					bool passes_filters = check_filters( *it );
					if ( !passes_filters  && loop_closure_failed ) {
						current_tag = "X_"+current_tag.substr(2);
					} else if( loop_closure_failed ) {
						current_tag = "C_"+current_tag.substr(2);
					} else if( !passes_filters ) {
						current_tag = "F_"+current_tag.substr(2);
					}
					bool bCanRelax = passes_filters && ( !loop_closure_failed || option[ OptionKeys::abinitio::relax_with_jumps ]() );
					if ( bRelax_ ) {
						if ( bCanRelax ) {
							it->constraint_set( NULL );
							relax( *it, fullatom_scorefxn, current_tag );
						} else { //cannot relax
							//need proper atom set to score with full-atom
							if ( !it->is_fullatom() ) core::util::switch_to_residue_type_set( *it, chemical::FA_STANDARD );
							(*fullatom_scorefxn)( *it );
							if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {	//FastRelax adds another two columns, grr
								relax::FastRelax::setPoseExtraScore(  *it );
							} else {
								relax::ClassicRelax().setPoseExtraScore( *it ); // so we the same number of columns
							}
						}
					}
					it->constraint_set( cstset_ ); // make-sure that all constraints are used for scoring... there might be still the MaxSeqSep... in t

					process_decoy( *it, *centroid_scorefxn,  current_tag, pss );
					pss.add_energy ( "store", ct );
					// write this to score-file if applicable
					if ( silent_score_file_ ) {
						silent_score_file_ -> write_silent_struct( pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
					}
					outsfd.add_structure( pss );
				} // for loop for structure-store
			}
			jobdist.dump_silent( outsfd ); // does the same thing as: outsfd.write_all( filename ); //cool bulk-writing makes clusters happy
		} //bProcessDecoy ( false if multi-fast_relax )
		// clean up
		time_t pdb_end_time = time(NULL);
		tr.Info << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (pdb_end_time - pdb_start_time) << " seconds." << std::endl;
		abinitio_protocol.get_checkpoints().clear_checkpoints();
		abinitio_protocol.structure_store().clear();
		abrelax_checkpoints_.clear_checkpoints();
	} // end of production loop
	jobdist.shutdown();
}
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::pair < core::pose::Pose, core::Real > PoseWithScore;
bool sort_PoseWithScore( const PoseWithScore& left, const PoseWithScore& right )
{
	return left.second < right.second;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail do fast relax on multiple structures that have been visited during abinitio-protocol.
/// MIKE: please give more documentation to this
bool JumpSpecificAbrelax::multi_fast_relax(
	 Protocol& abinitio_protocol,
	 core::scoring::ScoreFunctionOP scorefxn,
	 jobdist::PlainSilentFileJobDistributor< jobdist::BasicJobOP > jobdist,
	 int& curr_nstruct,
	 jobdist::BasicJobOP& curr_job
)
{
	using namespace protocols::jobdist;
	using namespace basic::options::OptionKeys;
	std::vector < PoseWithScore > candidates;
	for(platform::Size i=0; i < abinitio_protocol.structure_store().size(); i ++ ){
		pose::Pose cpose;

		cpose = abinitio_protocol.structure_store()[i];
		cpose.dump_pdb("staged_" + right_string_of(i,3,'0') + ".pdb" );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( cpose );
#endif

		core::util::switch_to_residue_type_set( cpose, chemical::FA_STANDARD );

		FastRelax fast_relax( scorefxn );
		// for making checkpoint tag unique
		fast_relax.set_current_tag( string_of(curr_nstruct) + "_" + right_string_of(i,3,'0') );
		fast_relax.apply( cpose );

		core::pose::setPoseExtraScore( cpose, "extranumber", float(i) );

		Real scorevalue = (*scorefxn)(cpose);

		candidates.push_back( PoseWithScore( cpose, scorevalue ) );

	}

	std::sort( candidates.begin(), candidates.end(), sort_PoseWithScore );

	bool endrun = false;
	for ( platform::Size i=0; i < std::min((size_t)candidates.size(), (size_t)300); i++ ) {
		if ( candidates[i].second > 0 ) continue;
		pose::Pose cpose;
		cpose = candidates[i].first;

		// analyse result
		io::silent::ProteinSilentStruct pss;
		process_decoy( cpose, *scorefxn, jobdist.get_current_output_tag(), pss );
		// write this to score-file if applicable
		if ( silent_score_file_ ) {
			silent_score_file_ -> write_silent_struct( pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
		}

		// write to silent file
		jobdist.dump_silent( curr_nstruct, pss );

		if ( !jobdist.next_job(curr_job, curr_nstruct) ){
		endrun = true;
		break;
		}
	}
	abinitio_protocol.structure_store().clear();
	return endrun;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail full-atom relax of decoys
/// uses either ClassicRelax or FastRelax protocols.
void JumpSpecificAbrelax::relax( pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag ) {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	//bail out if no relax is done
	if ( !option[ OptionKeys::abinitio::relax ]() && !option[ OptionKeys::abinitio::fastrelax ]() ) return;
	// run relax if applicable

	// remove constraints if option is set
	if ( option[ constraints::no_cst_in_relax ] ) {
		pose.constraint_set( NULL );
	}

	add_fa_constraints_from_cmdline( pose, *scorefxn );

	if ( option[ OptionKeys::abinitio::detect_disulfide_before_relax ] ) {
		pose.conformation().detect_disulfides();
	}

	if ( option[ OptionKeys::abinitio::relax ]() ) {
		if ( !abrelax_checkpoints_.recover_checkpoint( pose, tag, "abrelax_relax") ) {
			relax::ClassicRelax relax_protocol;
			relax_protocol.set_default( scorefxn );
			relax_protocol.set_current_tag( tag );
			relax_protocol.apply( pose );
			relax_protocol.get_checkpoints().clear_checkpoints();
			abrelax_checkpoints_.checkpoint( pose, tag, "abrelax_relax" ); //since relax_protocol throws away its checkpoints right here
		}
	}

	if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {
		FastRelax fast_relax( scorefxn );
		fast_relax.set_current_tag( tag );
		fast_relax.apply( pose );
	}
} // JumpSpecificAbrelax::relax
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail after setup() run either fold() or rerun()
void JumpSpecificAbrelax::run() {
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	setup();


	if ( option [ rerun ] ) {
		do_rerun();
		return;
	}

	if ( option [ jdist_rerun ] ) {
		do_distributed_rerun();
		return;
	}


	fold();
	return;
}

//============= Implementation of PoseEvaluators ======================================
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail
void ShowViolation::apply( pose::Pose &pose, std::string, io::silent::SilentStruct& ) const {
	using namespace basic::options::OptionKeys;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		//don't remove its not active ( on BOINC ) if you don't say -viol !!!
		pose.constraint_set()->show_violations(  std::cout, pose, option[ viol_level ] );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///@detail
void ComputeTotalDistCst::apply( pose::Pose &pose, std::string, io::silent::SilentStruct &pss ) const {
	using namespace basic::options::OptionKeys;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		if ( !constraints_ ) {
			constraints_ = new ConstraintSet( *pose.constraint_set() ); //get rid of MAX_SEQ_SEP
		}
		if ( option[ constraints::compute_total_dist_cst ]() ) {
			pose::Pose my_pose ( pose ); //copy of pose, don't want to change any score terms.
			my_pose.constraint_set( constraints_ );
			scoring::ScoreFunction scfxn;
			scfxn.set_weight( scoring::atom_pair_constraint, option[ constraints::cst_weight ]() );
			scfxn( my_pose );
			pss.add_energy("total_dist_cst", my_pose.energies().total_energies()[ scoring::atom_pair_constraint ] );
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void
PcaEvaluator::apply( pose::Pose& pose, std::string , io::silent::SilentStruct &pss ) const {
	PCA::ProjectionVector proj;
	pca_->eval( pose, proj );
	// tr.Info << "PCAEvaluator:  " << proj[1] << " " << proj[2] << std::endl;
	pss.add_energy ( "pca1", proj[1] );
	pss.add_energy ( "pca2", proj[2] );
}
////////////////////////////////////////////////////////////////////////////////////////////////////

} //abinitio
} //protocols
