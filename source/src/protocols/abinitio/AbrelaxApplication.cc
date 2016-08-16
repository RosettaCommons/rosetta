// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Application level code for relax-type protocols
/// @details
///
/// use AbrelaxApplication in the following way:
///
/// AbrelaxApplication::register_options();
/// core::init::init
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
/// options specific to AbrelaxApplication can be found by using -help at the command-line.
/// if you add new options please follow the scheme in the static method register options
///
/// the behaviour of AbrelaxApplication is controlled by comman-line-options. Refer to -help (usage) and the code
///
/// information that is not always present is stored as xxxOP, and the NULL-pointer is interpreted that the respective
/// behaviour is not present. (i.e., native_pose_ is either pointing to the native pose (-native) or to NULL.
/// when you use such pointers ask if they are non-NULL.
///
///
///
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// keep these headers first for compilation with Visual Studio C++
#include <utility/io/izstream.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

// Unit Headers
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/util.hh>
#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loops_main.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

// Package Headers
#include <core/kinematics/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/MembraneAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/KinematicTaskControl.hh>
#include <protocols/abinitio/LoopJumpFoldCst.hh>
#include <protocols/abinitio/DoubleLayerKinematicAbinitio.hh>
#include <protocols/abinitio/Templates.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/StrandConstraints.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/abinitio/Protocol.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/ResiduePairJumpSetup.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <protocols/jumping/util.hh>
#include <protocols/jumping/MembraneJump.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PCA.hh>
#include <protocols/moves/PyMolMover.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/scoring/SS_Killhairpins_Info.hh>
#include <core/scoring/methods/ContactOrderEnergy.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/util.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
#include <protocols/simple_filters/PoseMetricEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>
#include <protocols/loops/loop_closure/ccd/FASelectSlidingWindowLoopClosure.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Exceptions.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/RGFilter.hh>
#include <protocols/simple_filters/COFilter.hh>
#include <protocols/simple_filters/SheetFilter.hh>
#include <protocols/simple_filters/PDDFScoreFilter.hh>
#include <protocols/simple_filters/SAXSScoreFilter.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/RepulsiveOnlyMover.hh>
#include <protocols/forge/methods/util.hh>

//numeric headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>
#include <utility/exit.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameList.hh>
#include <core/id/SequenceMapping.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/abinitio/KinematicAbinitio.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.AbrelaxApplication" );

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details registering of options that are relevant for AbrelaxApplication
void protocols::abinitio::AbrelaxApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	register_common_options();
	// this adds all relevant options from the protocols
	// ClassicAbinitio, FoldConstraints, JumpingFoldConstraints
	KinematicAbinitio::register_options();
	Templates::register_options();
	loops::loop_closure::ccd::WidthFirstSlidingWindowLoopClosure::register_options();
	loops::loop_closure::ccd::FASelectSlidingWindowLoopClosure::register_options();

	// here we should have
	// ClassicRelax::register_options();
	// FastRelax::register_options(); etc.
#ifdef BOINC
	std::cerr << "Registered extra options." << std::endl; std::cerr.flush();
#endif
}

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
	core::Size size() const { return 2; };
	std::string name( core::Size ) const { return "pca1"; };
private:
	PCA_OP pca_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// evaluates the violations of atom-pairconstraints always with the full weight and full sequence separation
using namespace core::scoring::constraints; // has to be core, now that protocols::scoring is visible
class ShowViolation : public PoseEvaluator {
public:
	ShowViolation( ) : constraints_( /* NULL */ ) {}
	virtual void apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss ) const;
	core::Size size() const { return 1; };
	std::string name( core::Size ) const { return "viol"; };
private:
	mutable ConstraintSetOP constraints_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
class ComputeTotalDistCst : public PoseEvaluator {
public:
	ComputeTotalDistCst( ) : constraints_( /* NULL */ ) {};
	virtual void apply( pose::Pose& pose, std::string tag, io::silent::SilentStruct &pss ) const;
	core::Size size() const { return 1; };
	std::string name( core::Size ) const { return "total"; };
private:
	mutable ConstraintSetOP constraints_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail c'stor - nothing special
AbrelaxApplication::AbrelaxApplication() :
	silent_score_file_( /* NULL */ ),
	native_pose_( /* NULL */ ),
	pca_( /* NULL */ ),
	bRelax_ ( false ),
	cstset_( /* NULL */ ),
	jump_def_ ( /* NULL */ ),
	templates_( /* NULL */ ),
	fragset_large_( /* NULL */ ),
	fragset_small_top25_( /* NULL */ ),
	fragset_small_( /* NULL */ ),
	evaluator_ ( evaluation::MetaPoseEvaluatorOP( new MetaPoseEvaluator ) ),
	abrelax_checkpoints_( "Abrelax" )
{}

AbrelaxApplication::~AbrelaxApplication() {}

/// @details Shallow copy to mimic the pre 9/8/09 compiler-generated version of this
/// method.  If you add new
AbrelaxApplication::AbrelaxApplication( AbrelaxApplication const & src ) :
	silent_score_file_( src.silent_score_file_ ),
	native_pose_( src.native_pose_ ),
	loops_in_( src.loops_in_ ),
	pca_( src.pca_ ),
	bRelax_( src.bRelax_ ),
	sequence_( src.sequence_ ),
	cstset_( src.cstset_ ),
	membrane_jumps_( src.membrane_jumps_ ),
	jump_def_ ( src.jump_def_ ),
	ss_def_( src.ss_def_ ),
	templates_( src.templates_ ),
	fragset_large_( src.fragset_large_ ),
	fragset_small_top25_( src.fragset_small_top25_ ),
	fragset_small_( src.fragset_small_ ),
	fragset_templates_( src.fragset_templates_ ),
	evaluator_( src.evaluator_ ),
	abrelax_checkpoints_( src.abrelax_checkpoints_ )
{}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail add a PoseEvaluator derived instance for decoy-processing
void AbrelaxApplication::add_evaluation( evaluation::PoseEvaluatorOP eval ) {
	evaluator_->add_evaluation( eval );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details setup of Application data that is used for both, fold() and run()
/// this is mainly stuff for scoring and evaluation ( process_decoys(), evaluator_ )
void AbrelaxApplication::setup() {
	using namespace basic::options::OptionKeys;
	using core::scoring::func::FuncOP;

	if ( option[ constraints::no_linearize_bounded ] ) {
		tr.Info << "use fully harmonic potential for BOUNDED " << std::endl;
		ConstraintIO::get_func_factory().add_type("BOUNDED", FuncOP( new BoundFunc(0,0,0,1000,"dummy") ) );
	}

	if ( option[ constraints::named ] ) {
		tr.Info << "use named constraints in AtomPairConstraint to avoid problems with cutpoint-variants " << std::endl;
		/// WARNING WARNING WARNING. THREAD UNSAFE. DO NOT USE SINGLETONS THIS WAY.
		core::scoring::constraints::ConstraintFactory::get_instance()->replace_creator(
			ConstraintCreatorCOP( ConstraintCreatorOP( new constraints_additional::NamedAtomPairConstraintCreator ) ) );
	}

	silent_score_file_ = core::io::silent::SilentFileDataOP( new io::silent::SilentFileData );
	silent_score_file_-> set_filename( std::string( option[ out::sf ]()  ) );

	// read native pose
	if ( option[ in::file::native ].user() ) {
		native_pose_ = core::pose::PoseOP( new pose::Pose );
		core::import_pose::pose_from_file( *native_pose_, option[ in::file::native ]() , core::import_pose::PDB_file);

		pose::set_ss_from_phipsi( *native_pose_ );

#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( *native_pose_ );
#endif

		// allow sloppy matches here, because sometimes the Centroid residue set doesn't have all the residue variants
		// that the fullatom set has.
		core::util::switch_to_residue_type_set( *native_pose_, chemical::CENTROID, true ); //so that in do_rerun the native pose is the same as the other poses
	}

	// specify sequence -- from fasta file or native_pose
	if ( option[ in::file::fasta ].user() ) {
		sequence_ = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();
		tr.Info << "read fasta sequence: " << sequence_.size() << " residues\n"  << sequence_ << std::endl;
	} else if ( native_pose_ ) {
		sequence_ = native_pose_->sequence();
		tr.Info << "take sequence from native : " << sequence_ << std::endl;
	} else if ( !option[ OptionKeys::abinitio::rerun ]() && !option[ OptionKeys::abinitio::jdist_rerun ]() ) { // if we rerun we don't need sequence or native or anything...
		utility_exit_with_message(
			"Error: can't read sequence! Use -in::file::fasta sequence.fasta or -in::file::native native.pdb!"
		);
	}

	// run with homolog info? -- needed for setup_fragments, and rerun keep it upfront
	setup_templates();

	//add command-line evaluator stuff
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator_);

	core::pose::metrics::PoseMetricCalculatorOP
		clash_calculator( new protocols::toolbox::pose_metric_calculators::ClashCountCalculator( 2.0 ) );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "clashes", clash_calculator );
	add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::PoseMetricEvaluator<core::Size>( "clashes", "total" ) ) );
	add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::PoseMetricEvaluator<core::Size>( "clashes", "bb" ) ) );

	if ( option[ constraints::viol ]() ) add_evaluation( evaluation::PoseEvaluatorOP( new ShowViolation ) );
	if ( option[ constraints::compute_total_dist_cst ] ) add_evaluation( evaluation::PoseEvaluatorOP( new ComputeTotalDistCst ) );
	// read PCA info
	if ( option[ OptionKeys::in::file::pca ].user() ) {
		pca_ = evaluation::PCA_OP( new PCA );
		pca_->read_eigvec_file( option[ OptionKeys::in::file::pca ](), *native_pose_, 2 );
		add_evaluation( evaluation::PoseEvaluatorOP( new PcaEvaluator( pca_ ) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool AbrelaxApplication::close_loops( pose::Pose &pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag /*for checkpoints*/ ) {
	if ( !fragset_small_ ) {
		setup_fragments();
	}
	runtime_assert( !pose.is_fullatom() );
	(*scorefxn)( pose ); // just to check that we can do this --

	using namespace basic::options::OptionKeys;
	if ( option[ OptionKeys::loops::debug_loop_closure ]() ) pose.dump_pdb(tag+"_pre_closure.pdb");

	bool success( false );

	// Did we already close the loops successfully ?
	if ( (success=abrelax_checkpoints_.recover_checkpoint( pose, tag, "close_loops", false /*fullatom*/, true /*foldtree*/)) )  {
		abrelax_checkpoints_.debug( tag, "close_loops", (*scorefxn)(pose), (core::Real) true );
		return true;

	} else {
		// No ? Have we already tried but failed ?
		if ( abrelax_checkpoints_.recover_checkpoint( pose, tag, "close_loops_failure", false /*fullatom*/, true /*foldtree*/) )  {
			abrelax_checkpoints_.debug( tag, "close_loops", (*scorefxn)(pose), (core::Real) false );
			return false;
		}

		// Oh - we've not tried at all yet - let's go then!
		// make a MoveMap ... could be coming from somewhere else, though
		kinematics::MoveMapOP movemap( new kinematics::MoveMap );
		movemap->set_bb( true );

		// a weird bug occurs if we make a copy of a copy of the pose and set the new fold-tree
		// the behaviour is very different from setting the same fold-tree into the copy of the pose
		loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol( new loops::loop_closure::ccd::SlidingWindowLoopClosure( fragset_small_, scorefxn, movemap ) );

		if ( option[ OptionKeys::loops::alternative_closure_protocol ]() ) {
			closure_protocol = loops::loop_closure::ccd::SlidingWindowLoopClosureOP( new loops::loop_closure::ccd::WidthFirstSlidingWindowLoopClosure( fragset_small_, scorefxn, movemap ) );
		}

		// set options here if you like
		closure_protocol->set_native_pose( native_pose_ );
		closure_protocol->scored_frag_cycle_ratio( option[ OptionKeys::loops::scored_frag_cycles ]() );
		closure_protocol->short_frag_cycle_ratio( option[ OptionKeys::loops::short_frag_cycles ]() );

		//for debugging:
		closure_protocol->set_evaluation( evaluator_ );
		if ( option[ OptionKeys::abinitio::debug ] ) {
			closure_protocol->keep_fragments();
		}

		bool  bIdeal( !option[ OptionKeys::loops::non_ideal_loop_closing ]() );
		closure_protocol->set_bIdealLoopClosing( bIdeal );

		ProtocolOP debug_output( new Protocol );
		debug_output->set_evaluation( evaluator_ );

		success = true;

		try {
			jumping::close_chainbreaks( closure_protocol, pose, abrelax_checkpoints_, tag, kinematics::FoldTree() );
			if ( option[ OptionKeys::loops::debug_loop_closure ]() ) pose.dump_pdb(tag+"_post_closure.pdb");
		} catch ( loops::EXCN_Loop_not_closed& excn ) {
			success = false;
		}

		if ( success && option[ OptionKeys::loops::idealize_after_loop_close ]() ) {
			protocols::idealize::IdealizeMover idealizer;
			idealizer.fast( false );
			pose.constraint_set( NULL );
			idealizer.apply( pose );
			bIdeal = true;
		}

		if ( !bIdeal ) option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary");
		// to know this we'd have to catch the Exception EXCN_Loop_not_closed
		if ( success ) abrelax_checkpoints_.checkpoint( pose, tag, "close_loops", true /*foldtree*/ );
		else           abrelax_checkpoints_.checkpoint( pose, tag, "close_loops_failure", true /*foldtree*/ );

		abrelax_checkpoints_.debug( tag, "close_loops", (*scorefxn)(pose), (core::Real) success );
	}
	return success;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail run all evaluations on the decoy from this function
/// if you want these evaluations also available during internal stages of the protocols -->
///     put them into a PoseEvaluator and use add_evaluation in setup()
/// otherwise you can also use "manual" code right here in process_decoy --> this will only appear in final
/// silent_out and silent_score - files.
void AbrelaxApplication::process_decoy(
	pose::Pose &pose,
	core::scoring::ScoreFunction const& scorefxn,
	std::string tag,
	io::silent::SilentStruct &pss ) const
{
	using namespace basic::options::OptionKeys;
	// would like to put the following two also in an PoseEvaluator
	// ScoreEvaluator
	// StructureDumper
	if ( option[ OptionKeys::abinitio::clear_pose_cache ]() ) {
		tr.Debug << "\n******************************************************** \n"
			<< "              CLEAR POSE CACHE                           \n"
			<< "***********************************************************" << std::endl;
		pose.data().clear();
	}

	scorefxn( pose );
	pss.fill_struct( pose, tag );
	// run PoseEvaluators
	evaluator_->apply( pose, tag, pss );
	if ( option[ jumps::evaluate ]() ) {
		if ( !native_pose_ ) utility_exit_with_message(" to evaluate jumps you need to specify a native structure ");
		evaluation::MetaPoseEvaluator eval_jumps;
		native_pose_->fold_tree( pose.fold_tree() );
		for ( Size nj = 1; nj<= pose.num_jump(); ++nj ) {
			eval_jumps.add_evaluation( PoseEvaluatorOP( new simple_filters::JumpEvaluator( *native_pose_, nj) ) );
		}
		eval_jumps.apply( pose, tag, pss );
	}

} // process_decoy
////////////////////////////////////////////////////////////////////////////////////////////////////
//mjo commenting out 'pose' because it is unused and causes a warning
void AbrelaxApplication::initialize_constraint_forest( pose::Pose & /*pose*/ ) {
	using namespace basic::options::OptionKeys;
	if ( option[ constraints::forest_file ].user() ) {
		tr.Info << "read ConstraintForest... : " << std::endl;
		utility_exit_with_message( "ConstraintForest needs to be revived!" );
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail  read constraints file (once) and constraints_set to the pose (each call)
void AbrelaxApplication::add_constraints( pose::Pose & pose ) {
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;

	bool bFirst( !cstset_ );
	if ( bFirst ) {
		if ( option[ constraints::cst_file ].user() ) {
			// reads and sets constraints
			cstset_ = ConstraintIO::get_instance()->read_constraints(core::scoring::constraints::get_cst_file_option(), ConstraintSetOP( new ConstraintSet ), pose );
		}
	}

	if ( bFirst && templates_ ) {
		if ( !cstset_ ) cstset_ = core::scoring::constraints::ConstraintSetOP( new ConstraintSet );
		templates_->add_target_constraints( cstset_, pose );
		if ( option[ templates::strand_constraint ] ) {
			ConstraintCOPs my_strand_cst;
			if ( templates_ ) {
				my_strand_cst = StrandConstraints( templates_->strand_pairing_stats() ).build_constraints( pose );
			} else if ( option[ jumps::topology_file ].user() ) {
				utility::io::izstream is( option[ jumps::topology_file ] );
				if ( is.good() ) {
					PairingStatisticsOP ps( new PairingStatistics );
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
			add_evaluation( evaluation::PoseEvaluatorOP( new constraints_additional::ConstraintEvaluator( "strand", my_strand_cst ) ) );

			if ( native_pose_ ) { //just a temporary hack to test the StrandConstraint
				pose::Pose test_pose = *native_pose_;
				test_pose.add_constraints( my_strand_cst );

				if ( option[ constraints::dump_cst_set ].user() ) {
					tr.Info << "dump strand constraints to file..." << std::endl;
					utility::io::ozstream dump_cst( "STRAND_CST_DUMP" );
					test_pose.constraint_set()->show_definition( dump_cst, test_pose );
				}

				core::scoring::ScoreFunction cst_score;
				cst_score.set_weight( core::scoring::atom_pair_constraint, 1.0 );
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
		cstset_ = core::scoring::constraints::ConstraintSetOP( new ConstraintSet );
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
			add_evaluation( evaluation::PoseEvaluatorOP( new constraints_additional::ConstraintEvaluator( "seq_sep_"+utility::to_string( seq_sep) , *cstset_, 1, seq_sep ) ) );
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

	Stage1Sampler( protocols::simple_moves::FragmentMoverOP brute_move_large )
	: ClassicAbinitio( brute_move_large, brute_move_large, brute_move_large, 1 /*dummy*/ ) {};

	virtual void apply( core::pose::Pose &pose );
};

void Stage1Sampler::apply( core::pose::Pose &pose ) {
	prepare_stage1( pose );
	do_stage1_cycles( pose );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void AbrelaxApplication::insert_template_frags( core::pose::Pose &pose, kinematics::MoveMapOP movemap, std::string tag ) const {
	using namespace basic::options::OptionKeys;
	if ( option[ templates::fix_frag_file ].user() ) {
		FrameList fix_frames;
		fragment::FragmentIO().read_data( option[ templates::fix_frag_file ](), fix_frames );
		Size const frame_id ( static_cast< int >( numeric::random::rg().uniform() * fix_frames.size() ) + 1 );
		FrameOP frame( fix_frames[ frame_id ] );
		Size const frag_id ( static_cast< int >( numeric::random::rg().uniform() * frame->nr_frags() ) + 1 );
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
/// @detail loop over structures in silent-input file
/// small trick is used to also have native structure in the set of analysis:
/// it is added to the collection of silent_file-structures manually
/// TODO we need to do something about difference between fullatom and centroid input!
void AbrelaxApplication::do_rerun() {
	using namespace core;
	using namespace io::silent;
	using namespace pose;
	using namespace basic::options::OptionKeys;

	core::io::silent::SilentFileDataOP outsfd( NULL );
	if ( option[ out::file::silent ].user() ) {
		outsfd = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData() );
	}

	core::scoring::ScoreFunctionOP scorefxn( NULL );
	if ( option[ in::file::silent ].user() ) {
		//read silent file for input
		SilentFileData sfd;
		sfd.read_file( *(option [ in::file::silent ]().begin()) );

		// run thru all structures
		Size ct ( 0 );
		for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
			Pose pose;
			std::string tag = it->decoy_tag();
			if ( option[ in::file::tags ].user() == 0 || std::find( option[ in::file::tags ]().begin(), option[ in::file::tags ]().end(), tag ) != option[ in::file::tags ]().end() ) {
				if ( option[ in::file::fullatom ].user() ) {
					it->fill_pose( pose,
						option[ in::file::fullatom ] ?
						*(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD )) :
						*(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) ));
				} else {
					it->fill_pose( pose );
				}

				add_constraints( pose );
				scorefxn = generate_scorefxn( pose.is_fullatom() );
				//screen output
				if ( sfd.size() < 10 || option[ constraints::viol ]() ) {
					tr.Info << tag << " " << std::endl;
				} else {
					if ( (ct++ % 50) == 0 ) {
						std::cout << ".";
						std::cout.flush();
					}
				}

				// set score terms
				scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
				scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );

				if ( option[ OptionKeys::abinitio::close_loops ] ) {
					add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::RmsdEvaluator( pose::PoseOP( new pose::Pose( pose ) ), std::string("closure"), option[ OptionKeys::abinitio::bGDT ]() ) ) );
					close_loops( pose, scorefxn, tag );
				}

				basic::MetricValue<core::Size> mr;
				pose.metric("clashes","total",mr);
				tr.Info << "Total clashes " << mr.value() << std::endl;

				bool passes_filters = check_filters( pose );
				if ( !passes_filters ) {
					tag = "F_"+tag.substr(2);
				}

				SilentStructOP ss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
				process_decoy( pose, *scorefxn, tag, *ss );
				// write this to score-file if applicable
				if ( outsfd ) outsfd->add_structure( ss );

				//remove closure-rmsd
				if ( option[ OptionKeys::abinitio::close_loops ] ) {
					evaluator_->pop_back();
				}
			}
		}
	}

	// add native structure to the list
	if ( native_pose_ ) {
		SilentStructOP ss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		add_constraints( *native_pose_ );
		scorefxn = generate_scorefxn( false /*full_atom*/ );
		scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
		scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );

		if ( option[ OptionKeys::abinitio::close_loops ] ) { //otherwise the column (needed for non-native decoys) doesn't show up in score-file
			add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::RmsdEvaluator( pose::PoseOP( new pose::Pose( *native_pose_ ) ), std::string("closure"), option[ OptionKeys::abinitio::bGDT ]() ) ) );
		}

		process_decoy( *native_pose_, *scorefxn,  "NATIVE", *ss );
		// write this to score-file if applicable
		if ( silent_score_file_ ) {
			silent_score_file_ -> write_silent_struct( *ss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
		}

		if ( outsfd ) outsfd->add_structure( ss );
		if ( option[ OptionKeys::abinitio::close_loops ] ) evaluator_->pop_back();
	}

	if ( silent_score_file_ && outsfd ) {
		outsfd->write_all( silent_score_file_->filename(), true /* bWriteScoresOnly */ );
	}

	if ( outsfd ) outsfd->write_all( option[ out::file::silent ]() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail loop over structures in silent-input file
/// small trick is used to also have native structure in the set of analysis:
/// it is added to the collection of silent_file-structures manually
/// TODO we need to do something about difference between fullatom and centroid input!
void AbrelaxApplication::do_distributed_rerun() {
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
	if ( !option[ OptionKeys::run::no_prof_info_in_silentout ] ) {
		add_evaluation( run_time = evaluation::TimeEvaluatorOP( new evaluation::TimeEvaluator ) ); //just don't use this in integration tests!
	}

	loops_in_ = protocols::loops::Loops( true );

	// get input tags
	SilentFileData sfd;
	typedef utility::vector1< std::string > TagList;
	TagList input_tags;
	//WRONG -- these tags are the ones in file --- while after "read_file" tags might be renamed use .tags() after reading
	// input_tags = sfd.read_tags_fast( *(option [ in::file::silent ]().begin())  );

	// read silent data
	sfd.read_file( *(option [ in::file::silent ]().begin()) );
	input_tags = sfd.tags();
	// determine nstruct
	int const nstruct = std::max( 1, option [ out::nstruct ]() );

	// create jobs
	typedef utility::vector1< BasicJobOP > JobList;
	JobList input_jobs;
	for ( TagList::const_iterator it = input_tags.begin(), eit = input_tags.end(); it!=eit; ++it ) {
		BasicJobOP job( new BasicJob( *it, "rerun", nstruct) );
		input_jobs.push_back( job );
	}

	// setup JobDistributor
	PlainSilentFileJobDistributor jobdist( input_jobs );
	if ( option[ run::proc_id ].user() ) {
		int const procid ( option[ run::proc_id ] + ( option[ run::condor ] ? 1 : 0 ) );
		if ( procid > option[ run::nproc ] ) {
			utility_exit_with_message("procid to large " + ObjexxFCL::string_of( procid ) + " run only " + ObjexxFCL::string_of( option[ run::nproc ] ) + " processes");
		}
		jobdist.set_proc_id( procid, option[ run::nproc ] );
	}
	jobdist.startup(); //this will overwrite proc_id settings with mpi_rank if MPI is present.

	// production loop
	bool bEndrun = false;
	BasicJobOP curr_job;
	int curr_nstruct;
	while ( jobdist.next_job(curr_job, curr_nstruct) && !bEndrun ) {
		if ( run_time ) run_time->reset(); //reset clock of TimeEvaluator
		tr.Info << "Starting " << jobdist.get_current_output_tag() << " ..." << std::endl;
		tr.Info << "read " << curr_job->input_tag() << "..." << std::endl;
		Pose pose;

		sfd.get_structure( curr_job->input_tag() ).fill_pose( pose );
		set_ss_from_phipsi( pose );

		//mjo TODO: verify that the disulfides are correct coming out of fill_pose() and then delete this code
		// Fix disulfides if a file is given
		if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
			utility::vector1< std::pair<Size, Size> > disulfides;
			core::io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
			ds_file.disulfides( disulfides, pose);
			pose.conformation().fix_disulfides( disulfides );
		}

		loops_in_.verify_against( pose );

		if ( bRelax_ ) {  //always add all f@@#$ columns so we never-ever have a column mismatch....
			tr.Info << "relax is active... add stupid extra score terms " << std::endl;
			relax::ClassicRelax().setPoseExtraScore( pose );
		}

		std::string tag = jobdist.get_current_output_tag();
		if ( option[ OptionKeys::loops::random_grow_loops_by ].user() ) {
			loops_in_.grow_all_loops( pose, option[ OptionKeys::loops::random_grow_loops_by ]() );
			tr.Info << "Enlarged loops: " << std::endl;
			tr.Info << loops_in_ << std::endl;
		};

		add_constraints( pose ); // needs to come before scorefxn setup to know if constraints are present
		core::scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false /*fullatom*/ ) );
		core::scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true /*fullatom*/ ) );

		centroid_scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
		centroid_scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );

		if ( loops_in_.size() ) {
			utility::vector1< core::Real > vecs;
			loops::Loops rigid( loops_in_.invert( pose.total_residue() ) );
			loops::fix_with_coord_cst( rigid, pose, option[ loopfcst::coord_cst_all_atom ], vecs );
		}

		if ( option[ jumps::no_chainbreak_in_relax ] ) {
			fullatom_scorefxn->set_weight( core::scoring::linear_chainbreak, 0.0 );
			fullatom_scorefxn->set_weight( core::scoring::overlap_chainbreak, 0.0 );
		}
		// set score function for processing/relaxing stage

		tr.Info << tag << " " << std::endl;
		if ( option [ OptionKeys::abinitio::debug ] ) {
			//this functionality is needed for the iterative protocol to have a restart structure with the same tag as the final structure
			io::silent::SilentFileData outsfd;
			std::string silent_file = option[ basic::options::OptionKeys::out::file::silent ]() + "_" + "before_loops";

			io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
			process_decoy( pose, pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, jobdist.get_current_output_tag(), *pss );
			outsfd.write_silent_struct( *pss, silent_file );
		}

		bool loop_closure_failed( false );
		if ( option[ OptionKeys::abinitio::close_loops ] ) {
			add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::RmsdEvaluator( pose::PoseOP( new pose::Pose( pose ) ), std::string("closure"), option[ OptionKeys::abinitio::bGDT ]() ) ) );
			loop_closure_failed = !close_loops( pose, centroid_scorefxn, tag );
		}

		bool passes_filters = check_filters( pose );
		// run relax if applicable
		// don't relax if we failed filters or loop_closing, or if option[ relax_with_jumps ] is true
		bool bCanRelax = passes_filters && ( !loop_closure_failed || option[ OptionKeys::abinitio::relax_with_jumps ]() );
		if ( bRelax_ ) {
			if ( !pose.is_fullatom() ) {
				Pose const centroid_pose ( pose );
				core::util::switch_to_residue_type_set( pose, chemical::FA_STANDARD );
				pose.constraint_set( pose.constraint_set()->remapped_clone( centroid_pose, pose ) );
			}

			if ( bCanRelax ) {
				tr.Info << "relax is active... add stupid extra score terms " << std::endl;
				relax::ClassicRelax().setPoseExtraScore( pose );
				relax( pose, fullatom_scorefxn, jobdist.get_current_output_tag() );
			} else { //cannot relax
				//need proper atom set to score with full-atom
				(*fullatom_scorefxn)( pose );
				if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {}
				else {
					relax::ClassicRelax().setPoseExtraScore( pose ); // ClassicRelax adds four columns
				}
			}
		} // if ( bRelax_ )

		// process decoy if this hasn't happened yet
		// analyze result
		std::string output_tag = jobdist.get_current_output_tag();
		if ( !passes_filters  && loop_closure_failed ) {
			output_tag = "X_"+output_tag.substr(2);
		} else if ( loop_closure_failed ) {
			output_tag = "C_"+output_tag.substr(2);
		} else if ( !passes_filters ) {
			output_tag = "F_"+output_tag.substr(2);
		}

		SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		process_decoy( pose, pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, output_tag, *pss );
		// write this to score-file if applicable
		if ( silent_score_file_ ) {
			silent_score_file_ -> write_silent_struct( *pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
		}

		// write to silent file
		jobdist.dump_silent( curr_nstruct, *pss );

		//remove closure-rmsd
		if ( option[ OptionKeys::abinitio::close_loops ] ) {
			evaluator_->pop_back();
		}

	} // end of production loop
	jobdist.shutdown();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail called by setup_fold() if option[ start_native ] is active
/// the routine defines a fragment of the length of the structure
/// steals the fragment from the native and applies it to the decoy
/// native needs to be idealized!
void AbrelaxApplication::copy_native_structure( core::pose::Pose & extended_pose ) const {
	// requires that the sequences match at the beginning (1..nmatch_res) -- > use sequence alignment later
	tr.Info << " *** use native structure as starting template -- NEEDS TO BE IDEALIZED !!! *** \n";
	copy_structure( extended_pose, *native_pose_ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void AbrelaxApplication::copy_structure( core::pose::Pose & extended_pose, core::pose::Pose & desired_pose ) const {
	// requires that the sequences match at the beginning (1..nmatch_res) -- > use sequence alignment later
	tr.Info << " *** use native structure as starting template -- NEEDS TO BE IDEALIZED !!! *** \n";
	// determine length of segment to copy from native
	Size seg_len = std::min(extended_pose.total_residue(), desired_pose.total_residue() );
	// chu workaround when folding with ligand/metal
	Size protein_len = 0;
	for ( Size i = 1; i <= seg_len; ++i ) {
		if ( extended_pose.residue(i).is_protein() && desired_pose.residue(i).is_protein() ) {
			protein_len ++;
		}
	}
	seg_len = protein_len;
	fragment::Frame long_frame(1, seg_len);

	//create apropriate length FragData object
	FragData frag( SingleResidueFragDataOP( new BBTorsionSRFD ), seg_len );

	// get torsion angles from native pose
	frag.steal( desired_pose, long_frame );

	// apply native torsions to extended structue
	frag.apply( extended_pose, long_frame );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail called by setup_fold(): setup the decoy pose with correct target sequence and extended structure
///
void AbrelaxApplication::generate_extended_pose( core::pose::Pose &extended_pose, std::string const& sequence ) const {
	core::pose::make_pose_from_sequence(
		extended_pose,
		sequence,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
	);

	// Fix disulfides if a file is given
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		utility::vector1< std::pair<Size, Size> > disulfides;
		core::io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		ds_file.disulfides( disulfides, extended_pose);
		extended_pose.conformation().fix_disulfides( disulfides );
	}

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose.total_residue(); pos++ ) {
		if ( ! extended_pose.residue(pos).is_protein() ) continue;
		extended_pose.set_phi( pos, -150 );
		extended_pose.set_psi( pos, 150);
		extended_pose.set_omega( pos, 180 );
	}

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( extended_pose );
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail called by setup_fold(): read fragment libraries, I strongly suggest to use different options than A and B
/// if option[ steal ] fragments from the native structure are added to the set.
/// native structure needs to be idealized for this!
void AbrelaxApplication::setup_fragments() {// FragSetOP& fragsetA, FragSetOP& fragsetB ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string frag_large_file, frag_small_file;
	if ( option[ in::file::fragA ].user() ) {
		frag_large_file  = option[ in::file::fragA ]();
	} else {
		frag_large_file  = option[ in::file::frag9 ]();
	}

	if ( option[ in::file::fragB ].user() ) {
		frag_small_file  = option[ in::file::fragB ]();
	} else {
		frag_small_file  = option[ in::file::frag3 ]();
	}

	fragset_large_ = FragmentIO(
		option[ OptionKeys::abinitio::number_9mer_frags ](),
		option[ OptionKeys::frags::nr_large_copies ](),
		option[ OptionKeys::frags::annotate ]()
		).read_data( frag_large_file );

	if ( option[OptionKeys::abinitio::membrane] ) { //bw stupied way of get top25 frags for 3mer.
		fragset_small_top25_ = FragmentIO(
			option[ OptionKeys::abinitio::number_9mer_frags ],
			1, //nr_copies
			option[ OptionKeys::frags::annotate ]
			).read_data( frag_small_file );
	}
	fragset_small_ = FragmentIO(
		option[ OptionKeys::abinitio::number_3mer_frags ],
		1, //nr_copies
		option[ OptionKeys::frags::annotate ]
		).read_data( frag_small_file );

	if ( templates_ && !option[ templates::no_pick_fragments ]() ) {
		if ( option[ templates::vary_frag_size ] ) {
			fragset_templates_ = core::fragment::FragSetOP( new OrderedFragSet );
			templates_->pick_large_frags( *fragset_templates_, core::fragment::SingleResidueFragDataOP( new BBTorsionSRFD ), option[ templates::nr_large_copies ] );
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
					core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_large_->max_frag_length() ) ) ),
					min_nr_frags,
					nr_large_copies );
			} else {
				Size nr = templates_->pick_frags( *fragset_large_, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_large_->max_frag_length() ) ) ) );
				tr.Info << nr << " " << fragset_large_->max_frag_length() << "mer fragments picked from homolog structures" << std::endl;
			}
			if ( option[ templates::pick_multiple_sizes ] ) {
				Size nr = templates_->pick_frags(
					*fragset_large_,
					core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), 18 ) ) )
				);
				tr.Info << nr << " 18mer fragments picked from homolog structures" << std::endl;
				nr = templates_->pick_frags(
					*fragset_large_,
					core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP(  new BBTorsionSRFD ), 24 ) ) )
				);
				tr.Info << nr << " 27mer fragments picked from homolog structures" << std::endl;
			}
		} // !vary_frag_size

		if ( option[ templates::min_nr_small_frags ].user() ) {
			Size const min_nr_frags( option[ templates::min_nr_small_frags ] );
			Size const nr_small_copies( option[ templates::nr_small_copies ] );
			fragset_small_ = templates_->pick_frags(
				fragset_small_,
				core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_small_->max_frag_length() ) ) ),
				min_nr_frags,
				nr_small_copies );
		} else {
			//pick torsion fragments fragset_small
			Size nr2 = templates_->pick_frags( *fragset_small_, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_small_->max_frag_length() ) ) ) );
			tr.Info << nr2 << " " << fragset_small_->max_frag_length() << "mer fragments picked from homolog structures" << std::endl;
		}
	} // templates && !templates:no_pick_fragments

	if ( native_pose_ && ( option[ OptionKeys::abinitio::steal_3mers ]() || option[ OptionKeys::abinitio::steal_9mers ]() ) ) {
		tr.Info << " stealing fragments from native pose: ATTENTION: native pose has to be IDEALIZED!!! " << std::endl;
		if ( option[ OptionKeys::abinitio::steal_9mers ]() ) {
			steal_frag_set_from_pose( *native_pose_, *fragset_large_,
				core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_large_->max_frag_length() ) ) ) );
		}
		if ( option[ OptionKeys::abinitio::steal_3mers ]() ) {
			steal_frag_set_from_pose( *native_pose_, *fragset_small_,
				core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset_small_->max_frag_length() ) ) ) );
		}
	} else if ( ( option[ OptionKeys::abinitio::steal_3mers ]() || option[ OptionKeys::abinitio::steal_9mers ]() ) && !native_pose_ && !templates_ ) {
		tr.Warning << "cannot steal fragments without native pose or homologue structures " << std::endl;
	}

	if ( option[ OptionKeys::abinitio::dump_frags ]() ) { //diagnosis
		utility::io::ozstream dump_frag_small( "fragset_small.dump" );
		for ( ConstFrameIterator it=fragset_small_->begin(), eit=fragset_small_->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_small );
		}
		utility::io::ozstream dump_frag_large( "fragset_large.dump" );
		for ( ConstFrameIterator it=fragset_large_->begin(), eit=fragset_large_->end(); it!=eit; ++it ) {
			(*it)->show( dump_frag_large );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail called by setup_fold(). Read template definitions
void AbrelaxApplication::setup_templates() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool const bTemplates( option[ templates::config ].user() );
	if ( !bTemplates ) { // jump-out if not used
		if ( option[ templates::pairings ].user() ) {
			tr.Warning << "option templates:pairings ignored... specify templates:config!" << std::endl;
		}
		return;
	}

	if ( native_pose_ ) tr.Info << "native strand pairings " << core::scoring::dssp::StrandPairingSet( *native_pose_ );
	templates_ = TemplatesOP( new Templates( option[ templates::config ], native_pose_ ) );
	templates_->target_sequence() = sequence_; // a hack until class SequenceMapping works better
	// want to pick fragments from templates... make sure they are not initialized yet
	runtime_assert( !fragset_large_ );

	if ( !templates_->is_good() ) {
		utility_exit_with_message("ERRORS occured during template setup. check BAD_SEQUENCES file!");
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail called by setup_fold(). Read jump definitions / barcodes (not yet) etc.
/// if jump_def_ points to an object we will use JumpFoldConstraint-protocol in fold()
void AbrelaxApplication::setup_jumps( pose::Pose const& extended_pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// setup jumps
	bool bDoubleDef = false;
	jump_def_ = NULL;
	ss_def_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure( *fragset_small_, false /*no JustUseCentralResidue */ ) );

	if ( option [ jumps::extra_frags_for_ss ].user() ) {
		FragSetOP ss_frags = FragmentIO().read_data( option[ jumps::extra_frags_for_ss ]() );
		ss_def_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure( *ss_frags, false ) );
	}
	if ( option[ jumps::loop_definition_from_file ].user() ) {
		ss_def_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure() );
		ss_def_->read_from_file( option[ jumps::loop_definition_from_file ]() );
	}

	if ( option[ jumps::fix_jumps ].user() ) {
		JumpSetupOP ptr( new JumpSetup( extended_pose.total_residue() ) );
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
	if ( option[ jumps::jump_lib ].user() ) {
		bDoubleDef = jump_def_ != 0;
		JumpSelectorOP ptr( new JumpSelector( native_pose_->secstruct() ) );
		ptr->read_file( option[ jumps::jump_lib ] );
		jump_def_ = ptr;
	}
	if ( option[ jumps::sheets ].user() || option[ jumps::random_sheets ].user() ) {
		bDoubleDef = jump_def_ != 0;

		// get secondary structure info
		runtime_assert( fragset_small_ != 0 );

		ss_def_->show( tr.Trace );
		// get pairings file
		core::scoring::dssp::PairingsList pairings;
		if ( option[ jumps::pairing_file ].user() ) {
			read_pairing_list( option[ jumps::pairing_file ](), pairings );
		}

		// get sheet-topology
		jumping::SheetBuilder::SheetTopology sheets;
		if ( option[ jumps::sheets ].user() ) {
			sheets = option[ jumps::sheets ]();
			// done: instantiate sheet-builder
			jump_def_ = jumping::BaseJumpSetupOP( new SheetBuilder( ss_def_, pairings, sheets ) );
		} else {
			sheets = option[ jumps::random_sheets ]();
			jump_def_ = jumping::BaseJumpSetupOP( new RandomSheetBuilder( ss_def_, pairings, sheets ) );
		}
	}

	if ( option[ jumps::topology_file ].user() ) {
		utility::io::izstream is( option[ jumps::topology_file ] );
		if ( !is.good() ) {
			utility_exit_with_message(" did not find topology_file: " + std::string( option[ jumps::topology_file ]() ) );
		}
		PairingStatisticsOP ps( new PairingStatistics );
		is >> *ps;
		if ( is.fail() ) {
			utility_exit_with_message(" error reading file: " + std::string( option[ jumps::topology_file ]() ) );
		}
		tr.Info << *ps << std::endl;
		core::scoring::dssp::PairingList helix_pairings; //empty for now
		jump_def_ = jumping::BaseJumpSetupOP( new TemplateJumpSetup( NULL, ss_def_, ps, helix_pairings ) );
	}

	if ( option[ templates::pairings ] ) {
		bDoubleDef = false;
		if ( option[ jumps::fix_jumps ].user() ) {
			tr.Info << "use fixed jumps but take jump-geometries from template! " << std::endl;
			jump_def_ = jumping::BaseJumpSetupOP( new FixTemplateJumpSetup( *templates_->create_jump_def( ss_def_ ), jump_def_ ) );
		} else {
			bDoubleDef = jump_def_ != 0;
			jump_def_ = templates_->create_jump_def( ss_def_ );
		}

		if ( option[ constraints::forest_file ].user() ) utility_exit_with_message("can't use constraint-forest pairings with template pairings yet");
	}

	if ( option[ jumps::residue_pair_jump_file ].user() ) {
		bDoubleDef = jump_def_ != 0;
		ResiduePairJumpSetupOP ptr( new ResiduePairJumpSetup( extended_pose.total_residue() ) );
		ptr->read_file( option[ jumps::residue_pair_jump_file ]() );
		jump_def_ = ptr;
	}

	if ( bDoubleDef ) {
		utility_exit_with_message("you can only define one jump mode: choose one of -fix_jumps / -jump_lib / -sheets / -pairings");
	}

	if ( jump_def_ ) {
		//yields columns named nrjumps
		add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::JumpNrEvaluator ) );
	}
}

void AbrelaxApplication::setup_membrane_topology( pose::Pose & pose, std::string spanfile ) const {
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	core::scoring::MembraneTopologyOP topologyOP( new core::scoring::MembraneTopology );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
	core::scoring::MembraneTopology & topology=*( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
	topology.initialize(spanfile);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail setup_fold() all initialization that is necessary to run abinitio production loop in fold()
/// read fragments, make pose from sequence, get constraints, set jumps, movemap .etc
/// the two parameters are OUTPUT:
///    extended_pose to run A) with ( might actually contain native starting structure (option!) )
///    prot_ptr: an initialized instance of ClassicAbinitio, FoldConstraints or JumpingFoldConstraints depending on
///   user settings.
void AbrelaxApplication::setup_fold( pose::Pose& extended_pose, ProtocolOP& prot_ptr ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// ==========================================================================
	///  --------- fold()-specific setup --------------------------
	// ==========================================================================

	// ---------------------------------------------------------------------------------------------------------------
	// initialize pose
	generate_extended_pose( extended_pose, sequence_ );

	// apply cyclic peptide constraints if the option is selected
	if ( option[ OptionKeys::abinitio::cyclic_peptide ]() ) {
		protocols::relax::cyclize_pose( extended_pose );
	}

	// apply a mover which calculates only repulsive energy on designate residues
	protocols::simple_moves::RepulsiveOnlyMover replonly;
	replonly.apply( extended_pose );


	if ( option[ OptionKeys::abinitio::start_native ]() ) {
		copy_native_structure( extended_pose );
	} else if ( option[ in::file::s ].user() ) {
		core::pose::PoseOP tmp_pose( new core::pose::Pose );
		std::string fn = option[ in::file::s ](1);
		core::import_pose::pose_from_file( *tmp_pose, fn , core::import_pose::PDB_file);
		copy_structure( extended_pose, *tmp_pose );
	}

	// Fix disulfides if a file is given
	if ( option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		io::raw_data::DisulfideFile ds_file( option[ OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, extended_pose);
		extended_pose.conformation().fix_disulfides( disulfides );
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

	if ( option[OptionKeys::abinitio::membrane].user() ) {
		if ( option[in::file::spanfile].user() ) {
			std::string spanfile(option[OptionKeys::in::file::spanfile]());
			std::cout << "Reading spanfile " << spanfile <<"\n";
			setup_membrane_topology(extended_pose, spanfile);
			//read in membrane jumps if available
			membrane_jumps_ = jumping::MembraneJumpOP( new MembraneJump );
			std::cout << "1.TEMPLATE SIZE: " << membrane_jumps_->template_size() << "\n";
			std::cout << "1.PAIRING SIZE:  " << membrane_jumps_->pairings_size() << "\n";
			if ( option[ jumps::pairing_file ].user() && option[ jumps::jump_lib].user() ) {
				membrane_jumps_->init(option[ jumps::jump_lib ], option[ jumps::pairing_file ]); //template_file,pairings_file)
			}
			std::cout << "2.TEMPLATE SIZE: " << membrane_jumps_->template_size() << "\n";
			std::cout << "2.PAIRING SIZE:  " << membrane_jumps_->pairings_size() << "\n";

		}
	} else {
		// setup jumping... evtl. needs fragset to determine sheet/loop-fractions..
		setup_jumps( extended_pose );
	}

	// make a MoveMap
	kinematics::MoveMapOP movemap( new kinematics::MoveMap );
	movemap->set_bb( true );

	if ( option[ OptionKeys::abinitio::fix_residues_to_native ].user() ) {
		utility::vector1< int> const& fix_start_ends( option[ OptionKeys::abinitio::fix_residues_to_native ]() );
		for ( Size i=1; i + 1 <= fix_start_ends.size(); i+=2 ) {
			Size const start( fix_start_ends[ i ]);
			Size const end( fix_start_ends[ i + 1 ]);
			if ( !(end >= start) ) utility_exit_with_message("end < start in abinitio:fix_residues_to_native");
			fragment::Frame long_frame(start, end-start+1 );
			//create apropriate length FragData object
			FragData frag( SingleResidueFragDataOP( new BBTorsionSRFD ), end-start+1 );

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
		if ( option[  OptionKeys::loops::loop_file ].user() ) {
			loops_in_ = protocols::loops::Loops( true );
		}
		// if full-atom load starting structure as full-atom to recover sidechains later
		if ( option[ OptionKeys::loops::input_pdb ].user() ) {
			if ( option[ in::file::fullatom ]() ) {
				core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
				core::import_pose::pose_from_file( extended_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name() , core::import_pose::PDB_file);
				if ( !extended_pose.is_fullatom() ) utility_exit_with_message(" this full-atom pose should be a full-atom pose, no? ");
			} else {
				// centroid starting structure for loop-modeling
				core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
				core::import_pose::pose_from_file( extended_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name() , core::import_pose::PDB_file);
				if ( extended_pose.is_fullatom() ) utility_exit_with_message(" this centroid pose should not be a full-atom pose, no? ");
			}
			loops_in_.verify_against( extended_pose );
			if ( option[ OptionKeys::loops::random_grow_loops_by ].user() ) {
				loops_in_.grow_all_loops( extended_pose, option[ OptionKeys::loops::random_grow_loops_by ]() );
				tr.Info << "Enlarged loops: " << std::endl;
				tr.Info << loops_in_ << std::endl;
			};
		}
		loops_in_.verify_against( extended_pose );
		add_constraints( extended_pose );
		core::kinematics::simple_visualize_fold_tree( extended_pose.fold_tree(), tr.Debug );

		KinematicAbinitioOP sampler( new KinematicAbinitio( fragset_small_, fragset_large_, movemap /*this movemap will be ignored*/ ) );
		ResolutionSwitcher res_switch(
			extended_pose,
			extended_pose.is_fullatom(),
			sampler->start_from_centroid(),
			sampler->return_centroid()
		);
		if ( native_pose_ ) sampler->set_native_pose( native_pose_ );
		sampler->set_show_viol_level( option[ constraints::viol_level ] );
		sampler->init( res_switch.start_pose() );

		if ( option[ OptionKeys::abinitio::close_loops ]() ) {
			loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol( new loops::loop_closure::ccd::SlidingWindowLoopClosure );

			if ( option[ OptionKeys::loops::alternative_closure_protocol ]() ) {
				closure_protocol = loops::loop_closure::ccd::SlidingWindowLoopClosureOP( new loops::loop_closure::ccd::WidthFirstSlidingWindowLoopClosure );
			}

			if ( option[ OptionKeys::loops::fa_closure_protocol ]() ) {
				loops::loop_closure::ccd::FASelectSlidingWindowLoopClosureOP prot;
				closure_protocol = prot = loops::loop_closure::ccd::FASelectSlidingWindowLoopClosureOP( new loops::loop_closure::ccd::FASelectSlidingWindowLoopClosure );
				//spaeter kann man hier einen ResolutionSwitcher uebergeben.
				runtime_assert( extended_pose.is_fullatom() );
				prot->set_fullatom_pose( extended_pose );
			}

			// set options here if you like
			// closure_protocol-> ... and write the setters/accessors, too, if you have to
			closure_protocol->scored_frag_cycle_ratio( option[ OptionKeys::loops::scored_frag_cycles ]() );
			closure_protocol->short_frag_cycle_ratio( option[ OptionKeys::loops::short_frag_cycles ]() );
			closure_protocol->set_native_pose( native_pose_ );
			bool const bIdeal( !option[ OptionKeys::loops::non_ideal_loop_closing ]() );
			closure_protocol->set_bIdealLoopClosing( bIdeal );
			closure_protocol->set_chainbreak_max( option[ OptionKeys::loops::chainbreak_max_accept ]() );
			if ( !bIdeal ) option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary");
			if ( option[ OptionKeys::abinitio::debug ] ) {
				closure_protocol->keep_fragments();
			}

			sampler->closure_protocol( closure_protocol );
		}

		LoopJumpFoldCstOP controller;
		if ( option[  OptionKeys::loops::extended_loop_file ].user() ) {
			std::string filename( option[ OptionKeys::loops::extended_loop_file ]().name() );
			loops::Loops extended_loops_in( filename ); // <== TODO: select these using density score
			extended_loops_in.verify_against( extended_pose );
			KinematicAbinitioOP stage1_sampler( new KinematicAbinitio( *sampler ) );
			sampler->bSkipStage1_ = true;

			stage1_sampler->init( res_switch.start_pose() ); //sets default options
			stage1_sampler->bSkipStage3_ = true;
			stage1_sampler->bSkipStage4_ = true;
			stage1_sampler->closure_protocol( NULL );
			loops::Loops rigid_core( loops_in_.invert( extended_pose.total_residue() ) );
			controller = LoopJumpFoldCstOP( new DoubleLayerKinematicAbinitio(
				jump_def_,
				extended_loops_in,
				rigid_core,
				sampler,
				stage1_sampler,
				ss_def_,
				option[ loopfcst::coord_cst_weight ],
				option[ loopfcst::coord_cst_all_atom ]
				) );
		} else {
			controller = LoopJumpFoldCstOP( new LoopJumpFoldCst(
				jump_def_,
				loops_in_,
				sampler,
				ss_def_,
				option[ loopfcst::coord_cst_weight ],
				option[ loopfcst::coord_cst_all_atom ]
				) );
		}

		controller->set_input_pose_is_fa( extended_pose.is_fullatom() );
		prot_ptr = controller;
		if ( evaluator_->size() ) sampler->set_evaluation( evaluator_ );

		if ( option[ loopfcst::coord_cst_weight_array ].user() ) {
			utility::io::izstream file( option[ loopfcst::coord_cst_weight_array ]() );
			if ( !file.good() ) {
				utility_exit_with_message("ERROR:: Unable to open coord_cst_weight_array file: ");
			}
			utility::vector1< core::Real > weights;
			read_vector( file, weights );
			controller->set_coord_cst_weight_array( weights );
		}

		if ( option[ loopfcst::dump_coord_cst_weight_array ].user() ) {
			controller->set_dump_weights_file( option[ loopfcst::dump_coord_cst_weight_array ] );
		}


	} else { // not -use_general_protocol
		// setup abinitio protocol: one of either, MembraneAbinitio / ClassicAbinitio/ FoldConstraints/ JumpingFoldConstraints
		if ( option[  basic::options::OptionKeys::abinitio::membrane ]() ) {
			tr.Info << "run MembraneAbinitio.... " << std::endl;
			prot_ptr = ProtocolOP( new MembraneAbinitio( fragset_small_, fragset_small_top25_,fragset_large_, movemap) );
			tr.Info << "After new MembraneAbinitio.... " << std::endl;
		} else {
			if ( jump_def_ ) {
				tr.Info << "run JumpingFoldConstraints....." << std::endl;
				// it doesn't matter if we have no constraints the extra FoldConstraints part in the Jumping protocl
				// won't do anything
				JumpingFoldConstraintsWrapperOP pp;
				pp = JumpingFoldConstraintsWrapperOP( new JumpingFoldConstraintsWrapper( fragset_small_, fragset_large_, movemap, jump_def_ ) );
				if ( native_pose_ ) pp->set_native_pose( native_pose_ ); //to steal native jumps
				pp->set_show_viol_level( option[ constraints::viol_level ] );
				prot_ptr = pp;
			} else {
				if ( extended_pose.constraint_set()->has_residue_pair_constraints() ) {
					// We have constraints: run xxxFoldConstraints
					tr.Info << "run FoldConstraints....." << std::endl;
					FoldConstraintsOP pp;
					pp = FoldConstraintsOP( new FoldConstraints( fragset_small_, fragset_large_, movemap ) );
					pp->set_show_viol_level( option[ constraints::viol_level ] );
					prot_ptr = pp;
				} else {
					/// no constraints ---> ClassicAbinitio
					tr.Info << "run ClassicAbinitio....." << std::endl;
					prot_ptr = ProtocolOP( new ClassicAbinitio( fragset_small_, fragset_large_, movemap ) );
				}
			}
		}
	}

	Protocol& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer

	/// initialize protocol
	abinitio_protocol.init( extended_pose );
	if ( option[ OptionKeys::loopfcst::use_general_protocol ] ) {
		abinitio_protocol.return_centroid( !(bRelax_ || option[ OptionKeys::abinitio::return_full_atom ]) );
	} else {
		abinitio_protocol.return_centroid( true );
	}
	if ( evaluator_->size() ) abinitio_protocol.set_evaluation( evaluator_ );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
bool AbrelaxApplication::check_filters( core::pose::Pose & pose ) {
	using namespace protocols::filters;

	// return true if we're not supposed to use filters
	if ( !basic::options::option[ basic::options::OptionKeys::abinitio::use_filters ]() ) return true;
	if ( option[ basic::options::OptionKeys::filters::disable_all_filters ]() ) return true; //makes a lot of sense IMHO

	// apply numeric::random::rg(), contact-order and sheet filters
	protocols::simple_filters::RGFilter    rg_filter;
	protocols::simple_filters::COFilter    co_filter;
	protocols::simple_filters::SheetFilter sh_filter;

	if ( !option[ basic::options::OptionKeys::filters::disable_rg_filter ]() && !rg_filter.apply( pose) ) return false;
	if ( !option[ basic::options::OptionKeys::filters::disable_co_filter ]() && !co_filter.apply( pose) ) return false;
	if ( !option[ basic::options::OptionKeys::filters::disable_sheet_filter ]() && !sh_filter.apply( pose) ) return false;

	if ( ( option[basic::options::OptionKeys::filters::set_pddf_filter ].user() ) &&
			( option[basic::options::OptionKeys::score::saxs::ref_pddf ].user() ) ) {

		protocols::simple_filters::PDDFScoreFilterOP pddf_filter( new protocols::simple_filters::PDDFScoreFilter() );
		bool flag = pddf_filter->apply(pose);
		core::pose::setPoseExtraScore( pose, "pddf_score", pddf_filter->recent_score());
		if ( ! flag ) return false; // We need this flag because filter's score must be set before this if statement
	}

	if ( ( option[basic::options::OptionKeys::filters::set_saxs_filter ].user() ) &&
			( option[basic::options::OptionKeys::score::saxs::ref_spectrum ].user() ) ) {

		protocols::simple_filters::SAXSScoreFilterOP saxs_filter( new protocols::simple_filters::SAXSScoreFilter() );
		bool flag = saxs_filter->apply(pose);
		core::pose::setPoseExtraScore( pose, "saxs_score", saxs_filter->recent_score());
		if ( ! flag ) return false; // We need this flag because filter's score must be set before this if statement
	}

	tr.Info << " passed all filters " << std::endl;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP AbrelaxApplication::generate_scorefxn( bool fullatom ) {
	core::scoring::ScoreFunctionOP scorefxn( NULL );
	if ( fullatom ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		if ( option[  basic::options::OptionKeys::abinitio::membrane ]() ) {
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score_membrane" );
		} else if ( option[ OptionKeys::abinitio::stage4_patch ].user() ) {
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
		} else {
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		}
	}
	if ( jump_def_  && !option[ OptionKeys::jumps::no_chainbreak_in_relax ] ) {
		scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
		scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );
	}
	if ( cstset_ && cstset_->has_residue_pair_constraints()  ) {
		scorefxn->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn->set_weight( core::scoring::angle_constraint, option[ OptionKeys::constraints::cst_weight ]() );
		scorefxn->set_weight( core::scoring::dihedral_constraint, option[ OptionKeys::constraints::cst_weight ]() );
	}
	if ( option[ OptionKeys::loopfcst::coord_cst_weight ].user() ) {
		scorefxn->set_weight( core::scoring::coordinate_constraint, option[ OptionKeys::loopfcst::coord_cst_weight ]);
	}
	return scorefxn;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail everything happens in fold()!
/// setup of stuff that is not needed for rerun()
///    read fragments
///    [ optional ] steal fragments ( take fragments from native pose )
///
void AbrelaxApplication::fold( core::pose::Pose &init_pose, ProtocolOP prot_ptr ) {
	Protocol& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;

	// determine nstruct
	int const nstruct = std::max( 1, option [ out::nstruct ]() );

	// setup JobDistributor stuff
	utility::vector1< BasicJobOP > input_jobs;

	if ( option[ in::file::tags ].user() ) {
		// get input tags
		typedef utility::vector1< std::string > TagList;
		TagList input_tags;
		utility::io::izstream tag_file( option[ in::file::tagfile ]() );

		std::copy( std::istream_iterator< std::string >( tag_file ), std::istream_iterator< std::string >(),
			std::back_inserter( input_tags ) );

		// create jobs
		for ( TagList::const_iterator it = input_tags.begin(), eit = input_tags.end(); it!=eit; ++it ) {
			BasicJobOP job( new BasicJob( *it, "resample", nstruct) );
			tr.Debug << "create resample job" << *it << std::endl;
			input_jobs.push_back( job );
		}
	} else { //default behaviour
		BasicJobOP job( new BasicJob("" /*no input tag*/, "abinitio_relax", nstruct) );
		input_jobs.push_back( job );
	}

	PlainSilentFileJobDistributor jobdist( input_jobs );
	BasicJobOP curr_job;
	int curr_nstruct;
	jobdist.startup();
	bool bEndrun = false;

	// setup scorefunctions
	// this is called in setup_fold: add_constraints( extended_pose ); //such that scorefxn setup knows about constraints...
	core::scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false /*fullatom*/ ) );
	core::scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true /*fullatom*/ ) );
	abinitio_protocol.set_fullatom_scorefxn( fullatom_scorefxn );
	abinitio_protocol.set_centroid_scorefxn( centroid_scorefxn );

	evaluation::TimeEvaluatorOP run_time( NULL );
	if ( !option[ OptionKeys::run::no_prof_info_in_silentout ]() ) {
		add_evaluation( run_time = evaluation::TimeEvaluatorOP( new evaluation::TimeEvaluator ) ); //just don't use this in integration tests!
		abinitio_protocol.set_evaluation( evaluator_ );
	}

	// production loop
	while ( jobdist.next_job(curr_job, curr_nstruct) && !bEndrun ) {
		time_t pdb_start_time = time(NULL);
		if ( run_time && !option[ OptionKeys::abinitio::no_write_failures ]() ) {  //if we omit decoys we want to count from write-event to write-
			run_time->reset(); //reset clock of TimeEvaluator
		}

#ifdef BOINC
		std::cerr << "Starting work on structure: " << curr_job->output_tag(curr_nstruct) << std::endl;
#endif

		// retrieve starting pose
		pose::Pose fold_pose ( init_pose );
		// Can we add the PyMOL mover here?
		if ( option[OptionKeys::run::show_simulation_in_pymol].user()
				&& option[OptionKeys::run::show_simulation_in_pymol].value() > 0.0 ) {
			protocols::moves::AddPyMolObserver(fold_pose,
				option[OptionKeys::run::keep_pymol_simulation_history](),
				option[OptionKeys::run::show_simulation_in_pymol].value());
		}


		//membrane jumping set up the proper foldtree
		if ( membrane_jumps_ && membrane_jumps_->defined() ) {
			Size njumps = option[jumps::njumps]();
			membrane_jumps_->setup_fold_tree(fold_pose,njumps);
		}

		// kill hairpins
		//using core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO;
		if ( option[ OptionKeys::abinitio::kill_hairpins ].user() ) {
			using namespace basic::datacache;
			fold_pose.data().set( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO, DataCache_CacheableData::DataOP( new core::scoring::SS_Killhairpins_Info ) );
			runtime_assert( fold_pose.data().has( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO ) );
			core::scoring::SS_Killhairpins_Info & hairpins=*( utility::pointer::static_pointer_cast< core::scoring::SS_Killhairpins_Info > ( fold_pose.data().get_ptr( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO ) ));
			hairpins.setup_killhairpins();
		}

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( fold_pose );
#endif

		if ( bRelax_ ) {  //always add all f@@#$ columns so we never-ever have a column mismatch....
			relax::ClassicRelax().setPoseExtraScore( fold_pose );
		}

		// need to save numeric::random::rg() states such that choices for constraints and fold-tree are the same.
		if ( ! abrelax_checkpoints_.recover_checkpoint( fold_pose, jobdist.get_current_output_tag(), "rg_state") ) {
			abrelax_checkpoints_.checkpoint( fold_pose, jobdist.get_current_output_tag(), "rg_state");
		}
		abrelax_checkpoints_.debug( jobdist.get_current_output_tag(), "rg_state", numeric::random::rg().uniform() );

		// perturb phi/psi randomly -- should be different each run
		if ( option[ OptionKeys::abinitio::perturb ].user() ) {
			Real sig = option[ OptionKeys::abinitio::perturb ];
			for ( Size pos = 1; pos <= fold_pose.total_residue(); pos++ ) {
				fold_pose.set_phi( pos, fold_pose.phi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_psi( pos, fold_pose.psi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_omega( pos, fold_pose.omega( pos ) );
			}
		}

		// run abinitio

		// set the TAG
		abinitio_protocol.set_current_tag( jobdist.get_current_output_tag() );

		// and the JOB (that's different)
		abinitio_protocol.set_current_job( curr_job );
		tr.Debug << "fold_pose is " << (fold_pose.is_fullatom() ? " fullatom " : " centroid " ) << "before protocol run "<<std::endl;

		std::string output_tag = abinitio_protocol.get_current_tag();
		abinitio_protocol.apply( fold_pose );

		bool loop_closure_failed( false ); //!fold_pose.fold_tree().num_cutpoint() );

		if ( option[ OptionKeys::abinitio::close_loops ]()
				&& abinitio_protocol.get_last_move_status() == moves::MS_SUCCESS
				&& !option[ OptionKeys::loopfcst::use_general_protocol ] ) {
			tr.Info << "OLD PATHWAY: close loops" << std::endl;
			loop_closure_failed = !close_loops( fold_pose, centroid_scorefxn, jobdist.get_current_output_tag() );
		}

		if ( bRelax_ && option [ OptionKeys::abinitio::debug ] ) {
			io::silent::SilentFileData outsfd;
			std::string silent_file = option[ basic::options::OptionKeys::out::file::silent ]() + "_" + "before_relax";

			io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
			process_decoy( fold_pose, fold_pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, jobdist.get_current_output_tag(), *pss );
			outsfd.write_silent_struct( *pss, silent_file );
		}

		// run relax if applicable. Use filters to decide which structures to relax and which not to relax!
		bool passes_filters = check_filters( fold_pose );
		if ( option[OptionKeys::abinitio::membrane] ) {
			passes_filters=true;
		}
		bool bProcessDecoy( true );
		// run relax if applicable
		// don't relax if we failed filters or loop_closing, or if option[ relax_with_jumps ] is true
		bool bCanRelax = abinitio_protocol.get_last_move_status() == moves::MS_SUCCESS;
		if ( option[ basic::options::OptionKeys::abinitio::relax_failures ]() ) bCanRelax = true;
		bCanRelax = bCanRelax && passes_filters
			&& ( !loop_closure_failed || option[ OptionKeys::abinitio::relax_with_jumps ]() );

		if ( bRelax_ ) {
			//fpd  detect disulfides in centroid (or else they will repact as non-disulf and never get detected)
			if ( option[ OptionKeys::abinitio::detect_disulfide_before_relax ] ) {
				fold_pose.conformation().detect_disulfides();
			}

			if ( !fold_pose.is_fullatom() ) {
				pose::Pose const centroid_pose ( fold_pose );
				ResolutionSwitcher res_switch( centroid_pose, false, true, true );

				if ( option[ OptionKeys::constraints::cst_fa_file ].user() )   res_switch.set_map_cst_from_centroid_to_fa( false ); //will override any attempt to map centroid constraints to full-atom constraints -- use user-defined file instead!

				res_switch.apply( fold_pose );

				if ( option[ OptionKeys::constraints::cst_fa_file ].user() )  {
					ConstraintSetOP cstset_ = ConstraintIO::get_instance()->read_constraints( get_cst_fa_file_option(), ConstraintSetOP( new ConstraintSet ), fold_pose );
					fold_pose.constraint_set( cstset_ );
				}
			}

			if ( option[ basic::options::OptionKeys::abinitio::close_loops_by_idealizing ]() ) {
				// record cutpoints
				protocols::loops::LoopsOP cloops( new protocols::loops::Loops() );
				for ( Size ncut = 1; ncut <= (Size) fold_pose.fold_tree().num_cutpoint(); ncut++ ) {
					Size cutpoint = fold_pose.fold_tree().cutpoint( ncut );
					protocols::loops::Loop newloop (
						std::max( (int) 1, int(cutpoint - 5) ),
						std::min( (int) fold_pose.total_residue(), int(cutpoint + 5) ),
						0
					);

					if ( cloops->size() >= 2 ) {
						if ( newloop.start() <= ( *cloops )[cloops->size()-1].stop() ) newloop.set_start( ( *cloops )[cloops->size()-1].stop() +2 );
					}
					newloop.choose_cutpoint( fold_pose );
					cloops->add_loop( newloop );
				}

				cloops->auto_choose_cutpoints( fold_pose );

				// forget about foldtree & cuts
				core::kinematics::FoldTree f_new;
				f_new.simple_tree( fold_pose.total_residue() );
				fold_pose.fold_tree( f_new );

				//idealize
				protocols::idealize::IdealizeMover idealizer;
				idealizer.fast( false );
				idealizer.apply( fold_pose );

				relax( fold_pose, fullatom_scorefxn, jobdist.get_current_output_tag() );

				if ( option[ basic::options::OptionKeys::abinitio::optimize_cutpoints_using_kic ]() ) {
					protocols::loops::fold_tree_from_loops( fold_pose, *cloops, f_new, true /* include terminal cutpoints */);
					fold_pose.fold_tree( f_new );
					core::scoring::ScoreFunctionOP refine_scorefxn = fullatom_scorefxn->clone();
					protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine_kic( cloops, refine_scorefxn );
					refine_kic.apply( fold_pose );

					// Return fold tree to norml state
					f_new.simple_tree( fold_pose.total_residue() );
					fold_pose.fold_tree( f_new );
				}
			}

			if ( bCanRelax ) {
				if ( option[ basic::options::OptionKeys::abinitio::multifastrelax ]() ) {
					bEndrun = multi_fast_relax( abinitio_protocol, fullatom_scorefxn, jobdist, curr_nstruct, curr_job );
					bProcessDecoy = false;
					if ( bEndrun ) break;
				} else {
					relax( fold_pose, fullatom_scorefxn, jobdist.get_current_output_tag() );
				}
			} else { //cannot relax
				(*fullatom_scorefxn)( fold_pose );
				if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {}
				else {
					relax::ClassicRelax().setPoseExtraScore( fold_pose ); // ClassicRelax adds four columns
				}
			}
		} // if ( bRelax_ )

		//Add contact order to score file as an extra column
		core::scoring::methods::ContactOrderEnergy co_energy;
		Real contact_order = co_energy.calculate_contact_order( fold_pose );
		core::pose::setPoseExtraScore( fold_pose, "co", contact_order );

		// process decoy if this hasn't happened yet
		if ( bProcessDecoy ) {
			if ( option[ run::checkpoint ]() ) {
				core::pose::setPoseExtraScore( fold_pose, "ichkpnt",
					abinitio_protocol.get_checkpoints().get_checkpoint_recoveries() +
					abrelax_checkpoints_.get_checkpoint_recoveries()     );
			}

			// analyze result
			io::silent::SilentFileData outsfd;

			//make sure that number of columns does not change -- ever
			outsfd.strict_column_mode( true );

			io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();

			// abinitio produces n_stored structures -- the last one is the same as the final structure n_stored.
			//   Size n_stored( abinitio_protocol.structure_store().size() );
			std::string new_output_tag ( output_tag );
			if ( !passes_filters  && loop_closure_failed ) {
				new_output_tag = "X_"+new_output_tag.substr(2);
			} else if ( loop_closure_failed ) {
				new_output_tag = "C_"+new_output_tag.substr(2);
			} else if ( !passes_filters ) {
				new_output_tag = "F_"+new_output_tag.substr(2);
			} else if ( abinitio_protocol.get_last_move_status() != moves::MS_SUCCESS ) {
				new_output_tag = "P_"+new_output_tag.substr(2);
			}

			if ( !option[ OptionKeys::abinitio::no_write_failures ]()
					|| ( passes_filters && !loop_closure_failed && abinitio_protocol.get_last_move_status() == moves::MS_SUCCESS ) ) {
				// write to silent file
				//process_decoy( fold_pose, fold_pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn, new_output_tag, *pss );
				if ( fold_pose.is_fullatom() ) {
					process_decoy( fold_pose,  *fullatom_scorefxn, new_output_tag, *pss );
				} else if ( fold_pose.is_centroid() ) {
					process_decoy( fold_pose,  *centroid_scorefxn, new_output_tag, *pss );
				} else {
					core::scoring::ScoreFunctionOP cenrot_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_relax");
					process_decoy( fold_pose, *cenrot_scorefxn, new_output_tag, *pss );
				}

				// write this to score-file if applicable
				if ( silent_score_file_ ) {
					silent_score_file_ -> write_silent_struct( *pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
				}
				if ( option[ OptionKeys::out::pdb ] ) fold_pose.dump_pdb( std::string(option[ OptionKeys::out::path::path ]())  + "/" + new_output_tag + ".pdb");
				outsfd.add_structure( pss );
				if ( run_time && option[ OptionKeys::abinitio::no_write_failures ]() ) {  //if we omit decoys we want to count from write-event to write-
					run_time->reset(); //reset clock of TimeEvaluator
				}

			}

			jobdist.dump_silent( outsfd ); // does the same thing as: outsfd.write_all( filename ); //cool bulk-writing makes clusters happy
		} //bProcessDecoy ( false if multi-fast_relax )

		// clean up
		time_t pdb_end_time = time(NULL);
		tr.Info << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (pdb_end_time - pdb_start_time) << " seconds." << std::endl;
		abinitio_protocol.get_checkpoints().clear_checkpoints();
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
/// @detail do fast relax on multiple structures that have been visited during abinitio-protocol.
/// MIKE: please give more documentation to this
bool AbrelaxApplication::multi_fast_relax(
	Protocol& ,
	core::scoring::ScoreFunctionOP ,
	jobdist::PlainSilentFileJobDistributor ,
	int& ,
	jobdist::BasicJobOP&
)
{
	std::cerr << "multi_fast_relax stubbed out" << std::endl;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail full-atom relax of decoys
/// uses either ClassicRelax or FastRelax protocols.
void AbrelaxApplication::relax( pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag ) {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	//bail out if no relax is done
	if ( !option[ OptionKeys::abinitio::relax ]() &&
			!option[ OptionKeys::abinitio::fastrelax ]() ) return;

	// run relax if applicable

	//add cyclic peptide constraints during relax stages if and only if the cyclic peptide is specified
	if ( option[ OptionKeys::abinitio::cyclic_peptide ]() ) {
		protocols::relax::cyclize_pose( pose );
	}

	// remove constraints if option is set
	ConstraintSetOP orig_cst = NULL;
	if ( option[ constraints::no_cst_in_relax ] ) {
		orig_cst = pose.constraint_set()->clone();
		pose.constraint_set( NULL );
	}

	add_fa_constraints_from_cmdline( pose, *scorefxn );

	//fpd  move the detection to _before_ fullatom
	//if ( option[ OptionKeys::abinitio::detect_disulfide_before_relax ] ) {
	// pose.conformation().detect_disulfides();
	//}

	// Support deprecated option
	if ( option[ basic::options::OptionKeys::abinitio::fastrelax ]() ) {
		tr.Error << "WARNING: Using DEPRECATED OPTION fastrelax ! " << std::endl;
		option[ basic::options::OptionKeys::abinitio::relax ].def( true );
		option[ basic::options::OptionKeys::relax::fast ].def( true );
	}
	if ( option[ OptionKeys::abinitio::relax ]() ) {
		if ( !abrelax_checkpoints_.recover_checkpoint( pose, tag, "relax", true, true) ) {
			relax::relax_pose( pose, scorefxn, tag );
			abrelax_checkpoints_.checkpoint( pose, tag, "relax", true ); //since relax_protocol throws away its checkpoints right here
		}
		abrelax_checkpoints_.debug( tag, "relax", (*scorefxn)( pose ) );
	}

	/// Do a final clean relax without the constraints ?
	if ( option[ OptionKeys::loops::final_clean_fastrelax ]() ) {
		if ( !abrelax_checkpoints_.recover_checkpoint( pose,  tag, "final_fastrelax", true, true) ) {
			pose.constraint_set( NULL );
			relax::FastRelax fast_relax( scorefxn );
			fast_relax.set_current_tag( tag );
			fast_relax.apply( pose );
			abrelax_checkpoints_.checkpoint( pose, tag, "final_fastrelax", true ); //since relax_protocol throws away its checkpoints right here
		}
		abrelax_checkpoints_.debug( tag, "final_fastrelax", (*scorefxn)( pose ) );
	}
} // AbrelaxApplication::relax

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail after setup() run either fold() or rerun()
void AbrelaxApplication::run() {
	using namespace basic::options::OptionKeys;
	setup();

	if ( option [ OptionKeys::abinitio::rerun ] ) {
		do_rerun();
		return;
	}

	if ( option [ OptionKeys::abinitio::jdist_rerun ] ) {
		do_distributed_rerun();
		return;
	}

	// setup pose and abinitio
	ProtocolOP prot_ptr;
	pose::Pose init_pose;
#ifdef BOINC
	std::cerr << "Setting up folding (abrelax) ..." << std::endl;
#endif
	setup_fold( init_pose, prot_ptr ); //init_pose may or may not be fullatom (depending on flag option[ in:file:fullatom] )

#ifdef BOINC
	std::cerr << "Beginning folding (abrelax) ... " << std::endl;
#endif
	fold( init_pose, prot_ptr );
	return;
}

//============= Implementation of PoseEvaluators ======================================
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail
void ShowViolation::apply( pose::Pose &pose, std::string, io::silent::SilentStruct& ) const {
	using namespace basic::options::OptionKeys;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		//don't remove its not active ( on BOINC ) if you don't say -viol !!!
		pose.constraint_set()->show_violations(  std::cout, pose, option[ constraints::viol_level ] );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @detail
void ComputeTotalDistCst::apply( pose::Pose &pose, std::string, io::silent::SilentStruct &pss ) const {
	using namespace basic::options::OptionKeys;
	if ( pose.constraint_set()->has_residue_pair_constraints() ) {
		if ( !constraints_ ) {
			constraints_ = ConstraintSetOP( new ConstraintSet( *pose.constraint_set() ) ); //get rid of MAX_SEQ_SEP
		}
		if ( option[ constraints::compute_total_dist_cst ]() ) {
			pose::Pose my_pose ( pose ); //copy of pose, don't want to change any score terms.
			my_pose.constraint_set( constraints_ );
			core::scoring::ScoreFunction scfxn;
			scfxn.set_weight( core::scoring::atom_pair_constraint, option[ constraints::cst_weight ]() );
			scfxn( my_pose );
			pss.add_energy("total_dist_cst", my_pose.energies().total_energies()[ core::scoring::atom_pair_constraint ] );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PcaEvaluator::apply( pose::Pose& pose, std::string , io::silent::SilentStruct &pss ) const {
	PCA::ProjectionVector proj;
	pca_->eval( pose, proj );
	pss.add_energy ( "pca1", proj[1] );
	pss.add_energy ( "pca2", proj[2] );
}

} //abinitio
} //protocols
