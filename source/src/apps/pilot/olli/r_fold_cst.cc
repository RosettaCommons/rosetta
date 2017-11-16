// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file application to run fold-constraint protocol
/// @brief

// libRosetta headers


#include <core/types.hh>
#include <devel/init.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/viewer/viewers.hh>

using namespace protocols::abinitio;

void* my_main( void * )
{
	AbrelaxApplication my_app;
	my_app.run();
	return NULL;
}

int
main( int argc, char * argv [] )
{
	AbrelaxApplication::register_options();
	//FoldConstraint::register_options()
	devel::init( argc, argv );
	protocols::viewer::viewer_main( my_main );
	return 0;
}


#if 0
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <protocols/abinitio/util.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/conformation/ResidueFactory.hh>
#include <protocols/viewer/visualize.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

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

#include <core/pose/util.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/JumpingFoldConstraints.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/relax_protocols.hh>

#include <protocols/evaluation/PCA.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>


#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>

//#include <protocols/loops/Loops.hh>
//#include <protocols/loops/loops_main.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <utility/vector1.hh>

//numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>

static basic::Tracer tr( "r_fold_cst" );

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;
using namespace jumping;
using namespace evaluation;

#include "helper_code_r_trjconv.cc"

using namespace core;

#include <devel/simple_options/option.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using namespace devel::option;

class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
	void process_decoy( pose::Pose &pose, std::string tag, io::silent::ProteinSilentStruct & );
	void add_constraints( pose::Pose &pose );
	void do_rerun();
	void fold();
	void run();
	void setup();
	void fix_chainbreaks( pose::Pose &pose );
private:
	std::ofstream score_file_;
	io::silent::SilentFileDataOP silent_score_file_;
	pose::PoseOP rmsd_pose_;
	pose::PoseOP native_pose_;
	PCA_OP pca_;
	scoring::ScoreFunctionOP scorefxn_;
	std::string sequence_;
	scoring::constraints::ConstraintSetOP cstset_;
	BaseJumpSetupOP jump_def_;
	MetaPoseEvaluatorOP evaluator_;
};

ThisApplication::ThisApplication() :
	rmsd_pose_( NULL),
	pca_( NULL ),
	cstset_( NULL ),
	evaluator_ ( new MetaPoseEvaluator )
{}


using namespace basic::options;
using namespace basic::options::OptionKeys;


#define OPT(akey)                     \
	basic::options::option.add_relevant( akey )

#define NEW_OPT(akey,help,adef)                  \
	basic::options::option.add( akey , help ).def( adef ); 	\
	OPT( akey )

#define OPT_KEY( type, key )                                    \
	namespace core { namespace options { namespace OptionKeys {		\
				type##OptionKey const key( #key );											\
	} } }

OPT_KEY( Boolean, rerun )
OPT_KEY( Boolean, steal )
OPT_KEY( Boolean, start_extended )
OPT_KEY( File, pca )
OPT_KEY( File, rmsd_target )
OPT_KEY( File, sf )
OPT_KEY( Boolean, viol )
OPT_KEY( Integer, viol_level )
OPT_KEY( String, viol_type )
OPT_KEY( StringVector, tag_selector )
//OPT_KEY( Boolean, integrin_jump )
OPT_KEY( Real, perturb )
OPT_KEY( File, jumps )
OPT_KEY( File, jump_lib )
OPT_KEY( IntegerVector, rmsd_residues )
OPT_KEY( Boolean, bGDT )
OPT_KEY( Real, chainbreak_weight_stage1 )
OPT_KEY( Real, chainbreak_weight_stage2 )
OPT_KEY( Real, chainbreak_weight_stage3 )
OPT_KEY( Real, chainbreak_weight_stage4 )
OPT_KEY( Boolean, early_chainbreak_penalty )
OPT_KEY( Boolean, steal_3mers )
OPT_KEY( Boolean, steal_9mers )
OPT_KEY( Boolean, fix_chainbreak )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::native );
	OPT( in::file::silent ); // input silent file
	OPT( out::file::silent );
	OPT( in::file::frag3 );
	OPT( in::file::frag9 );
	OPT( in::file::fasta );
	OPT( constraints::cst_file );
	// OPT( constraints::cst_weight );
	OPT( out::nstruct );
	NEW_OPT( OptionKeys::rerun, "go through intput structures and evaluate ( pca, rmsd, cst-energy )", false );
	NEW_OPT( steal, "use fragments of native (starting) structure", true);
	NEW_OPT( start_extended, "start from extended structure (instead of native)", false );
	NEW_OPT( pca, "compute PCA projections", "");
	NEW_OPT( rmsd_target, "compute rmsd against this structure", "" );
	NEW_OPT( rmsd_residues, "give start and end residue for rmsd calcul.", -1 );
	NEW_OPT( sf, "filename for score output", "score.fsc" );
	NEW_OPT( viol, "show violations", true );
	NEW_OPT( viol_level, "how much detail for violation output", 1 );
	NEW_OPT( viol_type, "work only on these types of constraints", "");
	NEW_OPT( tag_selector, "work only on these tag(s)","");
	// NEW_OPT( integrin_jump, "use beta-sheet jumps for integrin", false );
	NEW_OPT( perturb, "add some perturbation (gaussian) to phi/psi of native", 0.0);
	NEW_OPT( jumps, "read jump_file", "" );
	NEW_OPT( jump_lib, "read jump_library_file for automatic jumps", "" );
	NEW_OPT( chainbreak_weight_stage1, "the weight on chainbreaks", 1.0 );
	NEW_OPT( chainbreak_weight_stage2, "the weight on chainbreaks", 1.0 );
	NEW_OPT( chainbreak_weight_stage3, "the weight on chainbreaks", 1.0 );
	NEW_OPT( chainbreak_weight_stage4, "the weight on chainbreaks", 1.0 );
	NEW_OPT( early_chainbreak_penalty, "use full chainbreak weight also in stage 1+2", false );
	NEW_OPT( bGDT, "compute gdtmmm", false );
	NEW_OPT( steal_3mers, "use 3mers from native", true );
	NEW_OPT( steal_9mers, "use 9mers from native", true );
	NEW_OPT( fix_chainbreak, "minimize to fix ccd in re-runs", false );
	JumpingFoldConstraints::register_options();
}


// little helper class
class PcaEvaluator : public PoseEvaluator {
public:
	PcaEvaluator ( PCA_OP pca ) : pca_( pca ) {};
	void apply( pose::Pose& pose, std::string tag, io::silent::ProteinSilentStruct &pss ) const;
private:
	PCA_OP pca_;
};

using namespace scoring::constraints;
class ShowViolation : public PoseEvaluator {
public:
	ShowViolation( ) : constraints_( NULL ) {};
	void apply( pose::Pose& pose, std::string tag, io::silent::ProteinSilentStruct &pss ) const {
		if ( pose.constraint_set() ) {
			tr.Info << "compute full constraint energy for " << tag << std::endl;
			pose.constraint_set()->show_violations(  std::cout, pose, option[ viol_level ] );
			if ( !constraints_ ) {
				constraints_ = new ConstraintSet( *pose.constraint_set() ); //get rid of MAX_SEQ_SEP
			}
			pose::Pose my_pose ( pose ); //copy of pose, don't want to change any score terms.
			my_pose.constraint_set( constraints_ );
			my_pose.constraint_set()->show_violations(  std::cout, pose, option[ viol_level ] );
			scoring::ScoreFunction scfxn;
			scfxn.set_weight( scoring::atom_pair_constraint, 1.0 );
			scfxn( my_pose );
			pss.add_energy("total_dist_cst", my_pose.energies().total_energies()[ scoring::atom_pair_constraint ] );
		}
	}
private:
	mutable ConstraintSetOP constraints_;
};

void
PcaEvaluator::apply( pose::Pose& pose, std::string , io::silent::ProteinSilentStruct &pss ) const {
	PCA::ProjectionVector proj;
	pca_->eval( pose, proj );
	// tr.Info << "PCAEvaluator:  " << proj[1] << " " << proj[2] << std::endl;
	pss.add_energy ( "pca1", proj[1] );
	pss.add_energy ( "pca2", proj[2] );
}

void ThisApplication::setup() {
	basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "~/minirosetta_database");

	silent_score_file_ = new io::silent::SilentFileData;
	silent_score_file_-> set_filename( std::string( option[ sf ]()  ) );

	// read native pose
	native_pose_ = new pose::Pose;
	core::import_pose::pose_from_pdb( *native_pose_, option[ in::file::native ]() );
	core::util::switch_to_residue_type_set( *native_pose_, chemical::CENTROID );

	// specify sequence -- from fasta file or native_pose
	if ( option [ in::file::fasta ].user() ) {
		sequence_ = read_fasta ( option[ in::file::fasta ]() );
	} else {
		sequence_ = native_pose_->sequence();
	}

	// set rmsd native
	if ( option[ rmsd_residues ].user() ) {
		Size start = option[ rmsd_residues ]()[ 1 ];
		Size end = option[ rmsd_residues ]()[ 2 ];
		evaluator_->add_evaluation( new RmsdEvaluator( native_pose_, start, end,  "_native", option[ bGDT ]() ) );
	} else {
		evaluator_->add_evaluation( new RmsdEvaluator( native_pose_, "_native", option[ bGDT ]() ) );
	}

	// set rmsd_target
	if ( option[ rmsd_target ].user() ) {
		rmsd_pose_ = new pose::Pose;
		rmsd_pose_->clear();
		core::import_pose::pose_from_pdb( *rmsd_pose_, option[ rmsd_target]() );
		if ( option[ rmsd_residues ].user() ) {
			Size start = option[ rmsd_residues ]()[ 1 ];
			Size end = option[ rmsd_residues ]()[ 2 ];
			evaluator_->add_evaluation( new RmsdEvaluator( rmsd_pose_, start, end, "", option[ bGDT ]() ) );
		} else {
			evaluator_->add_evaluation( new RmsdEvaluator( rmsd_pose_, "", option[ bGDT ]() ) );
		}

	}


	if ( option [ viol ] ) evaluator_->add_evaluation( new ShowViolation );

	// read PCA info
	if ( option[ pca ].user() ) {
		pca_ = new PCA;
		pca_->read_eigvec_file( option[ pca ](), *native_pose_, 2 );
		if ( tr.Trace.visible() ) pca_->show( std::cout );
		evaluator_->add_evaluation( new PcaEvaluator( pca_ ) );
	}

	// set score function
	scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scorefxn_->set_weight( scoring::atom_pair_constraint, option [ constraints::cst_weight ] );

	// setup jumps
	if ( option[ jumps ].user() && option[ jump_lib ].user() ) {
		tr.Error << "you can only define either -jumps or -jump_lib\n";
		utility_exit();
	}

	jump_def_ = NULL;
	if ( option[ jumps ].user() ) {
		JumpSetup *ptr = new JumpSetup;
		ptr->read_file( option[ jumps ] );
		jump_def_ = ptr;
		scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	}

	if ( option[ jump_lib ].user() ) {
		pose::set_ss_from_phipsi( *native_pose_ );
		JumpSelector *ptr = new JumpSelector( native_pose_->secstruct() );
		ptr->read_file( option[ jump_lib ] );
		jump_def_ = ptr;
		scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	}
}

void  ThisApplication::fix_chainbreaks( pose::Pose &pose ) {
	using namespace kinematics;
	FoldTree const &f ( pose.fold_tree() );
	for ( int i = 1; i<=f.num_cutpoint(); i++ ) {
		Size cut = f.cutpoint( i );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cut );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cut+1 );
	}
	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	// allow bb moves
	movemap->set_bb( true );
	protocols::simple_moves::MinMoverOP min_move_ = new protocols::simple_moves::MinMover;
	min_move_->movemap( movemap() );
	min_move_->min_type( "lbfgs_armijo_nonmonotone" );
	//get currently used score_function...
	min_move_->score_function( scorefxn_ );
	min_move_->apply( pose );
}

void ThisApplication::process_decoy( pose::Pose &pose, std::string tag, io::silent::ProteinSilentStruct &pss ) {
	//pose.dump_pdb("test.pdb");

	//fix_chainbreaks ... should probably become part of silent-file reader...
	if ( option[ fix_chainbreak ]() && pose.fold_tree().num_cutpoint() ) fix_chainbreaks( pose );

	// would like to put the following two also in an PoseEvaluator
	// ScoreEvaluator
	// StructureDumper
	( *scorefxn_ )( pose );
	if ( tr.Info.visible() ) scorefxn_->show( std::cout, pose );

	pss.fill_struct( pose, tag );

	evaluator_->apply( pose, tag, pss );

	if ( silent_score_file_ ) {
		silent_score_file_ -> write_silent_struct( pss,  silent_score_file_->filename(), true /* bWriteScoresOnly */ );
	}

}


void ThisApplication::add_constraints( pose::Pose& pose ){
	using namespace core::scoring::constraints;
	if ( ! cstset_ ) {
		//temporary hacd: also register HARMONIC
		core::scoring::constraints::ConstraintIO::get_func_factory().add_type( "HARMONIC", new core::scoring::constraints::HarmonicFunc(0,0) );
		if ( option[ constraints::cst_file ].user() ) {
			cstset_ = ConstraintIO::read( option[ constraints::cst_file ], new ConstraintSet, pose );
		}
	} else {
		pose.constraint_set( cstset_ );
	}
}

void ThisApplication::do_rerun() {
	using namespace core;
	using namespace io::silent;
	using namespace pose;

	//read silent file for input
	SilentFileData sfd;
	sfd.read_file( *(option[ in::file::silent ]().begin()) );

	// add native structure to the list
	sfd.add_structure( new ProteinSilentStruct( *native_pose_, "NATIVE", false ) );  /// !!! VERY BAD IF NOT IDEALIZED !!! */
	// run thru all structures
	for ( SilentFileData::const_iterator it=sfd.begin_const(), eit=sfd.end_const(); it!=eit; ++it ) {
		Pose pose;
		std::string tag = it->decoy_tag();
		if ( option[ tag_selector ].user() == 0 || std::find( option[ tag_selector ]().begin(), option[ tag_selector ]().end(), tag ) != option[ tag_selector ]().end() ) {
			it->fill_pose( pose, *chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ));
			if ( tag == "NATIVE" ) pose = *native_pose_; // replace structure with NATIVE so that we don't suffer from non-idealized stuff
			add_constraints( pose );
			tr.Info << tag << " " ;

			ProteinSilentStruct ss;
			process_decoy( pose, tag, ss );
		}
	}
}

void ThisApplication::fold() {
	// ==========================================================================
	///  --------- fold()-specific setup --------------------------
	// ==========================================================================
	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset3mer->read_fragment_file( option [ in::file::frag3 ]);
	fragset9mer->read_fragment_file( option [ in::file::frag9 ], 25 );

	if ( option [ steal ] || ( option[ steal_3mers ] && option[ steal_9mers ] ) ) {
		if ( option[ steal_3mers ] ) steal_constant_length_frag_set_from_pose( *native_pose_, *fragset3mer );
		if ( option[ steal_9mers ] ) steal_constant_length_frag_set_from_pose( *native_pose_, *fragset9mer );
	};

	pose::Pose extended_pose;
	core::pose::make_pose_from_sequence(  extended_pose, sequence_,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
	);

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose.size(); pos++ ) {
		extended_pose.set_phi( pos, -45 );
		extended_pose.set_psi( pos, -45 );
		extended_pose.set_omega( pos, 180 );
	}

	// requires that the sequences match at the beginning -- > use sequence alignment later
	if ( ! option [ start_extended ] ) {
		tr.Info << " *** use native structure as starting template *** \n";
		// determine length of segment to copy from native
		Size seg_len = std::min(extended_pose.size(), native_pose_->size() );
		fragment::Frame long_frame(1, seg_len);

		//create apropriate length FragData object
		FragData frag; //there should be some kind of factory to do this.
		Size nbb ( 3 ); //3 backbone torsions to steal
		for ( Size pos = 1; pos<= seg_len; pos++ ) {
			frag.add_residue( new BBTorsionSRFD( nbb, native_pose_->secstruct( pos ), 'X' ) );
		};
		// get torsion angles from native pose
		frag.steal( *native_pose_, long_frame);

		// apply native torsions to extended structue
		frag.apply(extended_pose, long_frame);
	}; // if option [ start extended ]

	extended_pose.dump_pdb( "starting_structure.pdb" );
	add_constraints( extended_pose );

	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	// allow bb moves
	movemap->set_bb( true );

	FoldConstraintsOP prot_ptr;
	if ( jump_def_ ) {
		tr.Info << "run JumpingFoldConstraints.....\n";
		JumpingFoldConstraints* pp;
		pp = new JumpingFoldConstraints( fragset3mer, fragset9mer, movemap, jump_def_ );
		if ( native_pose_ ) pp->set_native_pose( *native_pose_ ); //to steal native jumps
		if ( option[ early_chainbreak_penalty ]() ) pp->set_defeat_purpose( true );
		prot_ptr = pp;
	} else {
		tr.Info << "run FoldConstraints.....\n";
		prot_ptr = new FoldConstraints( fragset3mer, fragset9mer, movemap );
	}
	FoldConstraints& abinitio_protocol( *prot_ptr ); // hide the fact that protocol is a pointer


	abinitio_protocol.init( extended_pose );

	if ( evaluator_->size() ) abinitio_protocol.set_evaluator( evaluator_ );

	if ( jump_def_ ) { // move that into JumpingFoldConstraints
		// set low chainbreak energy for stage 1 + 2
		Real chainbreak_score_1 = option[ chainbreak_weight_stage1 ]();
		Real chainbreak_score_2 = option[ chainbreak_weight_stage2 ]();
		Real chainbreak_score_3 = option[ chainbreak_weight_stage3 ]();
		Real chainbreak_score_4 = option[ chainbreak_weight_stage4 ]();

		abinitio_protocol.set_score_weight( scoring::chainbreak, chainbreak_score_1, 1 );
		abinitio_protocol.set_score_weight( scoring::chainbreak, chainbreak_score_2, 2 );

		// set full chainbreak energy for stages 3 + 4
		abinitio_protocol.set_score_weight( scoring::chainbreak, chainbreak_score_3, 3 );
		abinitio_protocol.set_score_weight( scoring::chainbreak, chainbreak_score_4, 4 );
	}

	// ==========================================================================
	// ==========================================================================
	// ==========================================================================


	//     Start running abinitio


	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;
	utility::vector1< BasicJobOP > input_jobs;
	int const nstruct = std::max( 1, option [ out::nstruct ]() );
	BasicJobOP job = new BasicJob("classic_abinitio_relax", "classic_abinitio_relax", nstruct);
	input_jobs.push_back( job );
	PlainSilentFileJobDistributor< BasicJobOP > jobdist( input_jobs );
	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();
	while ( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		std::cout << "Starting " << jobdist.get_output_tag(curr_job, curr_nstruct) << " ...\n";

		// retrieve starting pose
		pose::Pose fold_pose ( extended_pose );

		// perturb phi/psi randomly
		if ( option[ perturb ].user() ) {
			Real sig = option[ perturb ];
			for ( Size pos = 1; pos <= fold_pose.size(); pos++ ) {
				fold_pose.set_phi( pos, fold_pose.phi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_psi( pos, fold_pose.psi( pos ) + numeric::random::gaussian()*sig );
				fold_pose.set_omega( pos, fold_pose.omega( pos ) );
			}
		}

		// run abinitio
		abinitio_protocol.set_output_tag( "debug_"+jobdist.get_output_tag(curr_job, curr_nstruct ) );
		abinitio_protocol.apply( fold_pose );

		if ( tr.Info.visible() ) fold_pose.constraint_set()->show_violations(  std::cout , fold_pose, option [ viol_level ] );

		// analyse result
		io::silent::ProteinSilentStruct pss;
		process_decoy( fold_pose, jobdist.get_output_tag( curr_job, curr_nstruct ), pss );
		tr.Debug << "RMSD in pose is " << fold_pose.energies().total_energies()[ scoring::rms ] << std::endl;

		// write to silent file
		jobdist.dump_silent( curr_nstruct, pss );
		prev_job = curr_job;
		num_structures_processed += 1;

		// clean up
		time_t pdb_end_time = time(NULL);
		std::cout << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (pdb_end_time - pdb_start_time) << " seconds.\n";
		abinitio_protocol.clear_checkpoints();
	}
	jobdist.shutdown();
}

void ThisApplication::run() {
	setup();

	if ( option [ rerun ] ) {
		do_rerun();
		return;
	};
	fold();
	return;

}

void* my_main( void * )
{
	ThisApplication my_app;
	my_app.run();
	return NULL;
}

int
main( int argc, char * argv [] )
{
	try{
		ThisApplication::register_options();
		//FoldConstraint::register_options()
		devel::init( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

#endif
