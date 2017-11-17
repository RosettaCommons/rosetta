// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @details
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/IterativeBase.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>


// Package Headers
#include <protocols/noesy_assign/NoesyModule.hh>
#include <protocols/noesy_assign/NoesyModule.impl.hh>

// to test broker setup-file
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/ConstraintEvaluatorWrapper.hh>
#include <protocols/topology_broker/util.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/util.hh>


#include <core/chemical/ChemicalManager.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#include <core/scoring/dssp/PairingsList.hh>
#endif

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/DecoySetEvaluation.impl.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>


#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/ScoreEvaluator.hh>
#include <protocols/simple_filters/RDC_Evaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>

#include <protocols/cluster/cluster.hh>
#include <protocols/toolbox/Cluster.hh>
#include <protocols/toolbox/Cluster.impl.hh>


#include <core/fragment/SecondaryStructure.hh>
// ObjexxFCL Headers

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/random/random.hh>
#include <utility/DereferenceIterator.hh>

// Option Headers
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
// Third-party Headers
#include <boost/functional/hash.hpp>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

// Utility headers
#include <basic/options/option_macros.hh>
#include <numeric/random/random_permutation.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/noesy_assign/CrossPeakList.impl.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

static basic::Tracer tr( "protocols.iterative" );
static basic::MemTracer mem_tr;


using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_2GRP_KEY( File, iterative, enumerate, broker )
OPT_2GRP_KEY( Integer, iterative, enumerate, Naccept )
OPT_2GRP_KEY( Boolean, iterative, enumerate, skip_half )
//OPT_1GRP_KEY( Real, iterative, turnover_rate )
OPT_1GRP_KEY( Integer, iterative, rmsf_nstruct )
OPT_1GRP_KEY( File, iterative, cen_score )
OPT_1GRP_KEY( File, iterative, cen_score_patch )

OPT_1GRP_KEY( File, iterative, fa_score )
OPT_1GRP_KEY( File, iterative, fa_score_patch )
OPT_1GRP_KEY( Real, iterative, chainbreak_evaluator_weight )
OPT_1GRP_KEY( IntegerVector, iterative, max_nstruct )
OPT_1GRP_KEY( Real, iterative, perturb_resampling )
OPT_1GRP_KEY( Boolean, iterative, mix_frags )
OPT_1GRP_KEY( Real, iterative, min_core_fraction_to_score )
OPT_1GRP_KEY( RealVector, iterative, min_diversity )
OPT_1GRP_KEY( RealVector, iterative, accept_ratio )
OPT_1GRP_KEY( Boolean, iterative, copy_pool_for_convergence_check )
OPT_1GRP_KEY( Real, iterative, safety_hatch_scorecut )
OPT_1GRP_KEY( Boolean, iterative, scored_ss_core )
OPT_1GRP_KEY( File, iterative, force_scored_region )
OPT_1GRP_KEY( Boolean, iterative, cluster )
OPT_1GRP_KEY( File, iterative, fix_core )
OPT_1GRP_KEY( Real, iterative, cenpool_noesy_cst_weight )
OPT_1GRP_KEY( String, iterative, chemicalshift_column )
OPT_1GRP_KEY( Real, iterative, cenpool_chemicalshift_weight )
OPT_1GRP_KEY( Real, iterative, centroid_before_quickrelax_weight )
OPT_1GRP_KEY( Real, iterative, fullatom_after_quickrelax_weight )
OPT_1GRP_KEY( Integer, iterative, limit_decoys_for_noe_assign )
OPT_1GRP_KEY( Boolean, iterative, centroid_quickrelax )
OPT_1GRP_KEY( String, iterative, super_quick_relax_protocol )
OPT_1GRP_KEY( String, iterative, super_quick_relax_patch )
OPT_1GRP_KEY( Integer, iterative, skip_redundant_constraints )
OPT_1GRP_KEY( String, iterative, initial_noe_auto_assign_csts )
OPT_1GRP_KEY( Integer, iterative, delay_noesy_reassign )
OPT_1GRP_KEY( Boolean, iterative, nolazy_noesy_reassign )
OPT_1GRP_KEY( String, iterative, auto_assign_scheme )
OPT_1GRP_KEY( Real, iterative, dcut )
OPT_1GRP_KEY( String, iterative, initial_beta_topology )
OPT_1GRP_KEY( Integer, iterative, recompute_beta_Naccept )
OPT_1GRP_KEY( String, iterative, flags_fullatom )
OPT_1GRP_KEY( Boolean, iterative, force_topology_resampling )
OPT_2GRP_KEY( Boolean, iterative, normalize, sampling )
OPT_1GRP_KEY( Boolean, iterative, update_noesy_filter_cst )
OPT_1GRP_KEY( Boolean, iterative, never_update_noesy_filter_cst )

OPT_1GRP_KEY( Boolean, iterative, split_autoNOE_restraints )
OPT_1GRP_KEY( Real, iterative, randomize_elimination_candidates )

OPT_1GRP_KEY( File, iterative, refine_auto_noe_flags )
OPT_1GRP_KEY( FileVector, iterative, staged_auto_noe_flags )

std::string const NOESY_CST_FILE_NAME("noe_auto_assign.cst");

bool protocols::abinitio::IterativeBase::options_registered_( false );


//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::abinitio::IterativeBase::register_options() {
	if ( !options_registered_ ) {
		Parent::register_options();
		noesy_assign::NoesyModule::register_options();
		OPT( in::file::silent );
		OPT( cluster::limit_cluster_size           );
		OPT( cluster::limit_clusters               );
		OPT( cluster::limit_total_structures       );

		NEW_OPT( iterative::enumerate::Naccept, "use enumerated pairings until Naccept decoys were added to archive", 5000 );
		NEW_OPT( iterative::enumerate::broker, "broker-file for enumerated_pairings", "" );
		NEW_OPT( iterative::enumerate::skip_half, "run half of the batches without enumerated pairings -- even before Naccept is reached", false );
		//  NEW_OPT( iterative::turnover_rate, "exchange of X percent of archive before new batch is started", 0.1 );
		NEW_OPT( iterative::accept_ratio, "switch to new stage if accept_ratio is lower than", 0.2 );
		NEW_OPT( iterative::rmsf_nstruct, "how many structures of pool used for computations of cores", 30 );
		NEW_OPT( iterative::cen_score, "energy function for centroid pool", "score3" );
		NEW_OPT( iterative::cen_score_patch, "patch of centroi_pool energy function", "NOPATCH" );
		NEW_OPT( iterative::fa_score, "energy function for centroid pool", "score12_full" );
		NEW_OPT( iterative::fa_score_patch, "patch of centroi_pool energy function", "NOPATCH" );
		NEW_OPT( iterative::max_nstruct, "give maximum numbers of structures generated for a given stage before switch -- 0 for infinite, -1 for skip stage", 0);

		NEW_OPT( iterative::mix_frags, "mix frags with original fragset", false );
		NEW_OPT( iterative::perturb_resampling, "perturb resample_stage2 start structures by this amount", 0.0 );
		NEW_OPT( iterative::min_core_fraction_to_score, "use core (2, 3, 4A) for discriminiation if more than X% of residues", 2.0); //>1 == off
		NEW_OPT( iterative::min_diversity, "don't accept structures that are closer than X to one of the pool", 1.5 );
		NEW_OPT( iterative::copy_pool_for_convergence_check, "make a copy of the pool for each batch that uses the convergence check", true );
		NEW_OPT( iterative::chainbreak_evaluator_weight, "weight for selection in pool...", 1.0 );

		NEW_OPT( iterative::safety_hatch_scorecut, "in CEN2FULLATOM state use structures for innput that are below scorecut in individual batches", 0.1);
		NEW_OPT( iterative::scored_ss_core, "selecte structures based on score restricted to the selected residues based on secondary structure", false );
		NEW_OPT( iterative::force_scored_region, "as scored_ss_core, but residues for scoring are provided manually", "" );

		NEW_OPT( iterative::cluster, "cluster archive with min_diversity_limit before each batch creation", false );

		NEW_OPT( iterative::fix_core, "RIGID file for RigidChunkClaimer to fix residues in stage 1-3", "" );

		NEW_OPT( iterative::cenpool_noesy_cst_weight, "weight to apply to centroid pool for noesy-autoassigned constraints", 1);
		NEW_OPT( iterative::cenpool_chemicalshift_weight, "weight to apply to chemical shifts in centroid pool rescoring", 5 );
		NEW_OPT( iterative::chemicalshift_column, "column name of the ChemicalShiftEvaluator used for chemical shift rescoring -- allows to have inactive shifts in score", "chem_shift" );
		NEW_OPT( iterative::super_quick_relax_protocol, "provide a sequence file for super quick relax ", "none" );
		NEW_OPT( iterative::super_quick_relax_patch, "provide a patch file for use during super quick relax ", "none" );
		NEW_OPT( iterative::centroid_before_quickrelax_weight, "add the centroid score of the decoy before it went into the <super_quick_relax> with this weight " , 1.0 );
		NEW_OPT( iterative::fullatom_after_quickrelax_weight, "add the centroid score of the decoy before it went into the <super_quick_relax> with this weight " , 0.0 );
		NEW_OPT( iterative::limit_decoys_for_noe_assign, "only use X first decoys for NOE assigment... ", 30 );
		NEW_OPT( iterative::centroid_quickrelax, "run a quick relax on centroid structures... ", false );
		NEW_OPT( iterative::skip_redundant_constraints, "skip constranits that have similar constraints within 1=same residue 2=neighbours (0=inactive)", 0 );

		NEW_OPT( iterative::delay_noesy_reassign, "start reassigning NOE restraints after X structures have been generated", 500 );
		NEW_OPT( iterative::nolazy_noesy_reassign, "do not wait with re-assigning until batch has been generated but do it immediately", false );
		NEW_OPT( iterative::initial_noe_auto_assign_csts, "initial nopool assignment of NOESY data.... ", "../initial_assignment/noe_auto_assign.cst" );
		NEW_OPT( iterative::auto_assign_scheme,"select CONST, ADAPT1, ... ","CONST");
		NEW_OPT( iterative::dcut,"in ADAPT1 what dcut should be chosen",7);
		NEW_OPT( iterative::initial_beta_topology,"start with this file as beta.top in stage3","" );
		NEW_OPT( iterative::force_topology_resampling,"if strand-fraction is low topology sampling is usually skipped. Override this with this flags", false );
		NEW_OPT( iterative::recompute_beta_Naccept, "recompute beta-topology after minimum of Naccept structures -- if no initial_beta_topology always recompute", 2000 );
		NEW_OPT( iterative::flags_fullatom, "point to flag-file to read flags for fullatom-refinement and loop-closing specify e.g., as ../flags_fullatom ","");
		NEW_OPT( iterative::normalize::sampling,
			"dynamically determined score-variation will be used for patching the sampling stage scores", false );
		NEW_OPT( iterative::update_noesy_filter_cst, "update noesy filter restraints (used for pool selection) as often as possible", false );
		NEW_OPT( iterative::never_update_noesy_filter_cst, "never update noesy filter restraints (used for pool selection) not even once for fullatom stage", false );
		NEW_OPT( iterative::refine_auto_noe_flags, "file with extra flags to modify NOESY-Module behavior in refinement phase (stage5+ )", "" );
		NEW_OPT8( iterative::staged_auto_noe_flags, "files with extra flags for each RASREC stage, say NONE for no option file", "NONE", "NONE", "NONE", "NONE", "NONE", "NONE", "NONE", "NONE" );
		NEW_OPT( iterative::split_autoNOE_restraints, "split the generated restraints into HI and MED accuracy class and use restraint combination only on the 2nd class", false );
		NEW_OPT( iterative::randomize_elimination_candidates, "put the elimination candidates into extra restraint file and only randomly activate them", 0 );


		options_registered_ = true;
	}
}

#define OBSOLETE(key)              \
	if ( option[ key ].user() ) {									\
		tr.Warning << "Option "<< #key<< " is deprecated!" << std::endl; \
	}

void warn_obsolete_flags() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// if ( option[ iterative::copy_pool_for_convergence_check ].user() ) {
	//  tr.Warning << "Option -iterative::copy_pool_for_convergence_check is deprecated!"  << std::endl;
	// }
	OBSOLETE(iterative::copy_pool_for_convergence_check);
	OBSOLETE(iterative::fix_core);
	OBSOLETE(iterative::enumerate::Naccept);
	OBSOLETE(iterative::min_core_fraction_to_score);
	OBSOLETE(iterative::force_scored_region);
	OBSOLETE(iterative::scored_ss_core);

	OBSOLETE( iterative::perturb_resampling);

	OBSOLETE( iterative::chainbreak_evaluator_weight);

	//OBSOLETE( iterative::safety_hatch_scorecut);
	OBSOLETE( iterative::cluster);


	//OBSOLETE( iterative::cenpool_chemicalshift_weight);
	//OBSOLETE( iterative::chemicalshift_column);
	//OBSOLETE( iterative::super_quick_relax_protocol);
	//OBSOLETE( iterative::centroid_before_quickrelax_weight);
	//OBSOLETE( iterative::fullatom_after_quickrelax_weight);
	//OBSOLETE( iterative::limit_decoys_for_noe_assign);
	//OBSOLETE( iterative::centroid_quickrelax);

	//noesy
	//OBSOLETE( iterative::skip_redundant_constraints);
	//OBSOLETE( iterative::cenpool_noesy_cst_weight);
	//OBSOLETE( iterative::delay_noesy_reassign);
	//OBSOLETE( iterative::initial_noe_auto_assign_csts);
	//OBSOLETE( iterative::auto_assign_scheme);
	//OBSOLETE( iterative::dcut);

	//OBSOLETE( iterative::initial_beta_topology);
	//OBSOLETE( iterative::force_topology_resampling);
	//OBSOLETE( iterative::recompute_beta_Naccept);

	//OBSOLETE( iterative::flags_fullatom);
	//OBSOLETE( iterative::normalize::sampling);

}

namespace protocols {
namespace abinitio {


using jd2::archive::EXCN_Archive;


void fix_chainbreak_patch( scoring::ScoreFunctionOP scorefxn, std::string const& patch_name ) {
	//constraints are done via patches
	bool chainbreaks_in_patch( scorefxn->get_weight( scoring::linear_chainbreak ) > 0.0 );
	core::Real extra_chainbreak_weight( basic::options::option[ basic::options::OptionKeys::jumps::increase_chainbreak ] );
	if ( !chainbreaks_in_patch ) { //if chainbreak is not patched .. set it
		core::Real lin_wgt( 4.0/3.0 * extra_chainbreak_weight );
		core::Real overlap_wgt( 1.0 * extra_chainbreak_weight );
		core::Real quad_wgt( 1.0 * extra_chainbreak_weight );
		tr.Info << "no chainbreaks specified in  " << patch_name << ",..."
			<< " set chainbreaks to " << lin_wgt << " and "
			<< overlap_wgt <<" for linear- and overlap-chainbreak, respectively" << std::endl;
		scorefxn->set_weight( scoring::chainbreak, quad_wgt );
		scorefxn->set_weight( scoring::linear_chainbreak, lin_wgt );
		scorefxn->set_weight( scoring::overlap_chainbreak, overlap_wgt );
	}
}

bool decide_on_beta_jumping_from_frags() {
	/// determine if we have enough strands to do jumping...
	bool bDoBetaJumping(true);
	using namespace core::fragment;
	FragSetOP frags_s = FragmentIO( option[ OptionKeys::abinitio::number_3mer_frags ]() ).read_data( option[ in::file::frag3 ] );
	core::fragment::SecondaryStructure ss_def( *frags_s, true /*no JustUseCentralResidue */ );
	Size ct_E( 0 );
	Size ct_H( 0 );
	for ( Size i=1; i<=ss_def.total_residue(); ++i ) {
		if ( ss_def.strand_fraction( i ) > 0.7 ) ++ct_E;
		if ( ss_def.helix_fraction( i ) > 0.7 ) ++ct_H;
	}
	mem_tr << "fragments loaded for ss-prediction" << std::endl;
	//less than 15% strands and less than 40 residues in total and we skip jumping...
	bDoBetaJumping = true;
	if ( ct_E < 0.1*ct_H && ct_E < 40 ) {
		tr.Info << "skip beta-jumping since mostly alpha-helical protein" << std::endl;
		tr.Info << "only " << (ct_E*1.0)/(ss_def.total_residue()*1.0)*100 << "% of residues display significant strand-character" << std::endl;
		tr.Info << "whereas " << (ct_H*1.0)/(ss_def.total_residue()*1.0)*100 << "% of residues display significant helix-character" << std::endl;
		bDoBetaJumping = false;
	}
	if ( !bDoBetaJumping && option[ OptionKeys::iterative::force_topology_resampling ]() ) {
		tr.Info << "force-topology-resampling despite low beta-sheet content due to flag -force_topology_resampling" << std::endl;
		bDoBetaJumping = true;
	}
	return bDoBetaJumping;
}


IterativeBase::IterativeBase(std::string name_in )
: Parent(),

	stage_( ENUMERATION ),
	finish_stage_( LAST_CENTROID_START ),
	first_batch_this_stage_ ( 1 ),
	first_fullatom_batch_( 100000 ),
	//  min_structures_for_accept_statistics_( 500 ),

	bEnumeratedLastTime_( false ),

	last_accepted_decoys_in_idle_( 0 ),

	reference_pose_( /* NULL */ ),

	cen_score_( option[ iterative::cen_score ]() ),
	cen_score_patch_( option[ iterative::cen_score_patch ]() ),

	fa_score_( option[ iterative::fa_score ]() ),
	fa_score_patch_( option[ iterative::fa_score_patch ]() ),
	noesy_assign_float_cycle_( 1.0 ), //changed OCT 20th 2010 ... start only in generation 3 of STAGE2_RESAMPLE with cyana-cycle 2
	first_noesy_cst_file_( "n/a" ),
	first_noesy_fa_cst_file_ ("n/a" ),
	current_noesy_sampling_file_( "n/a" ),
	bCombineNoesyCst_( true ),
	super_quick_relax_of_centroids_( option[ iterative::centroid_quickrelax ]() ),
	use_dynamic_weights_for_sampling_( option[ iterative::normalize::sampling ]() ),
	delay_noesy_reassign_( option[ iterative::delay_noesy_reassign ]() ),
	rdc_data_( /* NULL */ ),
	cst_data_( /* NULL */ ),
	cst_fa_data_( /* NULL */ ),
	vanilla_options_( basic::options::option )
{
	never_switched_noe_filter_=true;
	warn_obsolete_flags();
	basic::show_time_on_cerr = true;

	using namespace ObjexxFCL;
	mem_tr << "IterativeBase CStor-Start" << std::endl;
	//changes for debug mode
	if ( option[ run::test_cycles ] || option[ run::dry_run ] ) {
		delay_noesy_reassign_ = 50;
		//  min_structures_for_accept_statistics_ = 20;
	}

	//name is e.g., centroid_pool or fullatom_pool
	set_name( name_in );

	// --- setup stage-steering parameters
	// maximum of sampled structures per stage
	max_nstruct_list_ = option[ iterative::max_nstruct ]();
	if ( max_nstruct_list_.size() != ( FINISHED-1 ) ) {
		throw CREATE_EXCEPTION(EXCN_Archive, "wrong number "+string_of( max_nstruct_list_.size() )
			+" of values for max_nstruct -- needs exactly "+string_of( FINISHED - 1 )+" values.");
	}
	max_nstruct_list_.push_back( 0 ); //for FINISH_STAGE

	// min cluster radius for diversification
	min_diversity_list_ = option[ iterative::min_diversity ]();
	if ( min_diversity_list_.size() != ( FINISHED-1 ) ) {
		throw CREATE_EXCEPTION(EXCN_Archive, "wrong number "+string_of( min_diversity_list_.size() )
			+" of values for min_diversity -- needs exactly "+string_of( FINISHED - 1 )+" values.");
	}
	min_diversity_list_.push_back( 0 ); //for FINISH_STAGE

	// accept ratio to end stage
	target_accept_ratio_ = option[ iterative::accept_ratio ]();
	if ( target_accept_ratio_.size() != ( FINISHED-1 ) ) {
		throw CREATE_EXCEPTION(EXCN_Archive, "wrong number "+string_of( target_accept_ratio_.size() )
			+" of values for accept_ratio -- needs exactly "+string_of( FINISHED - 1 )+" values.");
	}
	target_accept_ratio_.push_back( 1 ); //for FINISH_STAGE

	/// -- setup native pose
	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		reference_pose_ = native_pose;
		mem_tr << "setup native pose" << std::endl;
	}

	/// --- setup sequence
	if ( option[ in::file::fasta ].user() ) {
		target_sequence_ = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();
		if ( reference_pose_ ) {
			if ( reference_pose_->sequence()!=target_sequence_ ) {
				std::string msg("target sequence and sequence of native pose -in:file:native are inconsistent! ");
				msg+="\nNATIVE: "+reference_pose_->sequence();
				msg+="\nFASTA : "+target_sequence_;
				utility_exit_with_message(msg);
			}
		}
	} else {
		if ( reference_pose_ ) {
			target_sequence_ = reference_pose_->sequence();
		} else {
			throw CREATE_EXCEPTION(EXCN_Archive, "need either fasta-sequence (-in:file:fasta) or native structure ( -in:file:native )");
		}
	}

	if ( tr.Trace.visible() ) {
		if ( reference_pose_ ) reference_pose_->dump_pdb("reference_pose_in_IterativeBase.pdb");
		tr.Trace << "target_sequence_:\n " << target_sequence_ << std::endl;
	}
}

void IterativeBase::initialize() {
	Parent::initialize();
	// --- setup scorefxn
	core::scoring::ScoreFunctionOP scorefxn =
		core::scoring::ScoreFunctionFactory::create_score_function( cen_score(), cen_score_patch() );
	tr.Info << "create Archive Scorefunction with: "
		<<  option[ iterative::cen_score]() << " "
		<<  option[ iterative::cen_score_patch ]() << std::endl;

	mem_tr << "setup cen-scorefxn" << std::endl;

	//manually set the chainbreaks if user has forgotten to specify them ...
	fix_chainbreak_patch( scorefxn, cen_score_patch() );

	set_scorefxn( scorefxn );
	set_overall_cstfilter_weight( scorefxn->get_weight( scoring::atom_pair_constraint ) );
	set_weight( "atom_pair_constraint", 0 ); //this is now done via FILTER mechanism of ConstraintClaimer only !
	scorefxn->set_weight( scoring::atom_pair_constraint, 0 );

	// --- setup pool-evaluation
	if ( evaluate_local() ) {
		set_weight( "score", 1.0 );
		//setup constraint-evaluation for filter-cst from the ConstraintClaimers ( -broker:setup )
		/// @brief set scorefxn used for evaluation
		setup_filter_cst( overall_cstfilter_weight() );
		//  set_weight( "prefa_clean_score3", option[ iterative::centroid_before_quickrelax_weight ]() );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::RDC_Evaluator("rdc") ) ), scorefxn->get_weight( scoring::rdc ) );
		scorefxn->set_weight( scoring::rdc, 0 );
	} else {
		set_weight( "score", 0.0 ); //don't use score that comes back --- but the score_final thing
		set_weight( "score_final", 1.0 );
		set_weight( "rdc_filt", scorefxn->get_weight( scoring::rdc ) );
		scorefxn->set_weight( scoring::rdc, 0 );
	}
	set_scorefxn( scorefxn );
	//will setup autoNOE module if cmd-line options activated -- might set super_quick_relax_of_centroids_ to true.
	setup_autoNOE();

	// setup chemical shift rescoring -- this is done on worker side, here we pick up the column and give it a weight
	if ( option[ iterative::cenpool_chemicalshift_weight ].user() ) {
		chemshift_column_ = option[ iterative::chemicalshift_column ]();
		set_weight( option[ iterative::chemicalshift_column ](), option[ iterative::cenpool_chemicalshift_weight ]() );
		super_quick_relax_of_centroids_=true;
	}

	bDoBetaJumping_ = decide_on_beta_jumping_from_frags();
	if ( !bDoBetaJumping_ ) {
		stage_ = PURE_TOPO_RESAMPLING;
	}

	//if we use super-quick-relax mode setup appropriate score-evaluation
	//and put a 0-Scorefunction into the Archive (as the scores are computed on worker side)
	//keep this at end of constructor since other setup-routines can change this boolean value
	if ( super_quick_relax_of_centroids_ ) {
		set_weight( "prefa_clean_score3", option[ iterative::centroid_before_quickrelax_weight ]() );
		set_weight( "score_fa", option[ iterative::fullatom_after_quickrelax_weight ]() );

		// --- setup scorefxn a NULL scorefxn since we cannot evaluate a centroid score anymore (could do that by changing residue_type_set...)
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		set_scorefxn( scorefxn );
		set_weight( "score", 0 );
	}

	//start stage-end controller
	test_for_stage_end();
	mem_tr << "IterativeBase CStor-End" << std::endl;
}

IterativeBase::~IterativeBase() = default;

/// @details ready for new batch .... if queue is empty batch will be generated any way, but otherwise we only generate if this yields true.
///  logic here: new batch at beginning, but only if we are in startup phase ( not a reload of a full archive )
///              otherwise make new batch if sufficiently many structures have been accepted since last batch
// bool IterativeBase::ready_for_batch() const {
//  if ( !bStartedInitBatch_ && decoys().size() <= nstruct() ) return true;
//  if ( decoys().size() >= nstruct() && accepts_since_last_batch() > nstruct_accept_before_rerun_ && proposed_since_last_batch() > min_structures_for_accept_statistics_ ) return true;
//  return false;
// }

void IterativeBase::idle() {
	if ( last_accepted_decoys_in_idle_ < ( accepts_since_last_batch() + total_accepts() - decoys().size()/2 ) ) {
		basic::show_time( tr,  "compute_cores/rescore: start..." );
		last_accepted_decoys_in_idle_ = accepts_since_last_batch() + total_accepts();
		compute_cores();
		rescore();
		basic::show_time( tr,  "compute_cores/rescore: done..." );
		//if we take out the compute_cores() we should still trigger a rescore() to get the new score-variations...
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//if we are in remote-evaluation mode we want to switch the noesy-filter restraints once after stage3
	if ( !evaluate_local()
			&& option[ iterative::update_noesy_filter_cst ]()
			&& stage_ >= STAGE2_RESAMPLING
			&& noesy_assign::NoesyModule::cmdline_options_activated()
			&& never_switched_noe_filter_
			&& decoys().begin() != decoys().end()
			) {
		basic::show_time( tr,  "non-local scoring replace noesy_filter constraints: start..." );
		never_switched_noe_filter_=false;
		replace_noesy_filter_constraints();
		basic::show_time( tr,  "non-local scoring replace noesy_filter constraints: done! ");
	}
}
///  ----------------- stage control ----------------------
// ---> outsource to extra class ?  --- might be reused by other protocols
/// @brief batch is expired ?
bool IterativeBase::still_interested( jd2::archive::Batch const& batch ) const {
	return Parent::still_interested( batch ) && batch.id() >= first_batch_this_stage_;
}

//after reading new structures we test for energy-saturation
void IterativeBase::read_structures(
	core::io::silent::SilentFileData& sfd,
	core::io::silent::SilentFileData& alternative_decoys,
	jd2::archive::Batch const& batch
) {
	basic::show_time( tr,  "read structures into "+name()+"..." );
	Parent::read_structures( sfd, alternative_decoys, batch );
	test_for_stage_end();
	basic::show_time( tr,  "done reading into "+name() );
}

/// @brief are we ready to switch to next stage ?
void IterativeBase::test_for_stage_end() {
	//switch to next stage ? want to have some significance to this ratio --- hence at least 1000 proposals
	tr.Info << "current accept ratio: " << current_acceptance_ratio() << " this is "
		<< ( current_acceptance_ratio() < target_accept_ratio() ? "" : "not" )
		<< " lower than " << target_accept_ratio()
		<< "\n"
		<< proposed_since_last_batch() << " proposed decoys since last batch "
		<< std::endl;

	int last_stage_N( 0 );
	if ( stage_ > 1 ) last_stage_N = max_nstruct_list_[ stage_ - 1 ];
	if ( max_nstruct_list_[ stage_ ] && ( (int) total_proposed() - last_stage_N ) > max_nstruct_list_[ stage_ ] ) {
		tr.Info << "maximum number of " << max_nstruct_list_[ stage_ ]
			<< " decoys for stage " << stage_ << " is reached...switching!" << std::endl;
		increment_stage();
	}

	if ( current_acceptance_ratio() < target_accept_ratio() ) {
		increment_stage();
	}
}

/// @brief got to next stage
void IterativeBase::increment_stage() {
	if ( stage_ >= finish_stage_ ) return;
	save_to_file( "_stage" + ObjexxFCL::string_of( stage_ ) );

	//save current number of structures:
	max_nstruct_list_[ stage_ ] = total_proposed();  //used in test_for_stage_end to no how-many new structures since last stage

	//switch to next possible stage
	stage_ = IterationStage( 1 + (int) stage_ );
	while ( (( max_nstruct_list_[ stage_ ] < 0 ) || ( !bDoBetaJumping_ && stage_ == NOESY_PHASEII_TOPO )) && stage_ < finish_stage_ ) {
		stage_ = IterationStage( 1 + (int) stage_ );
	}
	if ( !manager_ptr() ) {
		tr.Info << "start with stage " << stage_ << std::endl;
		return;
	}
	tr.Info << "manager_ptr() " << manager_ptr() << " logical " << !manager_ptr() << std::endl;
	tr.Info << "switched to stage: " << stage_ << std::endl;

	//safe status and reset counters
	first_batch_this_stage_ = manager().last_batch_id()+1;
	save_to_file();
	reset_accept_counter();
	manager().cancel_batches_previous_to( manager().last_batch_id(), true /*allow reading of decoys already in production*/ );

	//if we are not doing enumeration anymore we don't need fragments.
	using namespace core::fragment;
	if ( stage_ >= PURE_TOPO_RESAMPLING ) FragmentIO().clean_frag_cache();

}

/// ------------------------ end stage control
/// @detail rescore and sort archive
void  IterativeBase::rescore() {
	Parent::rescore();
	if ( hedge_archive_ ) {
		hedge_archive_->set_evaluators( evaluators(), weights() );
		hedge_archive_->set_weights( weights() );
		hedge_archive_->set_scorefxn( scorefxn().clone() );
		hedge_archive_->rescore();
	}
}


void IterativeBase::collect_hedge_structures( core::io::silent::SilentStructOP evaluated_decoy, jd2::archive::Batch const& batch ) {
	// if ( !evaluate_local() ) return;
	if ( !hedge_archive_ ) {
		hedge_archive_ = HedgeArchiveOP( new HedgeArchive( name()+"_hedge" ) );
		hedge_archive_->initialize();
		hedge_archive_->set_evaluators( evaluators(), weights() );
	}
	hedge_archive_->add_evaluated_structure( evaluated_decoy, nullptr, batch );
}

/// @brief overload to check for pool_convergence data in incoming decoys
bool IterativeBase::add_structure(
	core::io::silent::SilentStructOP new_decoy,
	core::io::silent::SilentStructOP alternative_decoy,
	jd2::archive::Batch const& batch
) {

	core::io::silent::SilentStructOP evaluated_decoy = evaluate_silent_struct( new_decoy );
	if ( decoys().size() == 0 ) {
		first_batch_this_stage_ = batch.id();
	}

	if ( stage_ < CEN2FULLATOM || ( stage_ == CEN2FULLATOM && batch.id() < first_batch_this_stage_ ) ) {
		select_score( evaluated_decoy ); //make sure that _archive_select_score_ is there
		collect_hedge_structures( evaluated_decoy, batch );
	}
	//comes without pool-convergence ? nothing to do
	if ( !new_decoy->has_energy( "pool_converged_tag" ) || min_diversity_list_[ stage() ] == 0 ) {
		return Parent::add_evaluated_structure( evaluated_decoy, alternative_decoy, batch );
	}

	//okay, let's look at the rmsd to closest structure
	runtime_assert( new_decoy->has_energy( "pool_converged_rmsd" ) );
	Real const rmsd_to_pool( new_decoy->get_energy( "pool_converged_rmsd" ) );

	if ( rmsd_to_pool > min_diversity_list_[ stage() ] ) { //structure is sufficiently different -- add via score
		return Parent::add_evaluated_structure( evaluated_decoy, alternative_decoy, batch );
	} else { //structure is close in RMSD to one of archive
		std::string const tag( new_decoy->get_string_value( "pool_converged_tag" ) );

		//find the pool-structure that is redundant with new decoy
		SilentStructs::iterator it;
		for ( it = decoys().begin(); it != decoys().end(); ++it ) {
			if ( (*it)->decoy_tag() == tag ) break;
		}

		if ( it == decoys().end() ) {
			//can't find tag ... (might be that we are close to a centroid structure and now add to fullatom pool ).
			//                    might be that we have swapped away the original structure
			return Parent::add_evaluated_structure( evaluated_decoy, alternative_decoy, batch );
		}

		//improved score ?
		if ( it != decoys().end() ) {
			core::Real const new_score=select_score( evaluated_decoy );
			core::Real const old_score=select_score( *it );
			core::Real const check_score=(*it)->get_energy( "_archive_select_score_" );
			if ( new_score < old_score ) {
				core::Real const min_score=select_score( decoys().front() );
				core::Real const max_score=select_score( decoys().back() );
				core::Real const delta_max=max_score-min_score;
				core::Real const delta_score=new_score-old_score;
				tr.Debug << "swap " << evaluated_decoy->decoy_tag() << " for " << (*it)->decoy_tag()
					<< " new-score " << new_score << " check: " << check_score
					<< " score-improvement: " << delta_score
					<< " at spread of " << delta_max
					<< " ratio " << delta_score/delta_max << std::endl;
				if ( delta_score/delta_max < -0.05 ) {
					erase_decoy( tag );
					return Parent::add_evaluated_structure( evaluated_decoy, alternative_decoy, batch );
				} else {
					tr.Debug << "swap declined because score-improvement below 5% of total spread" << std::endl;
					return false;
				}
			} // new_score < old_score
		}//improve score block

		tr.Trace << "decoy " << evaluated_decoy->decoy_tag() << " with original tag " << evaluated_decoy->get_comment( "tag_in_file" )
			<< " declined because of min_diversity: rmsd is " << rmsd_to_pool
			<< " limit: " << min_diversity_list_[ stage() ]
			<< std::endl;

		return false; //declined ... either tag not found or score not improved...
	}
	return false; //should never get here
}

/// @details generate new batch...
/// type of batch depends on stage_. we switch to next stage based on some convergence criteria:
/// right now it is how many decoys were accepted from last batch.. if this number drops sufficiently ---> next stage...
///    (maybe need to put a safeguard in here: ratio small but at least XXX decoys proposed since last batch... )
///
core::Size IterativeBase::generate_batch( jd2::archive::Batch& batch, core::Size repeat_id ) {

	if ( repeat_id == 0 ) return repeat_id;
	mem_tr << "IterativeBase::generate_batch " << stage_ << " " << batch.batch() << std::endl;
	tr.Info << "\ngenerate batch from " <<name() << " " << batch.batch() << std::endl;

	//want intermediate structures from abinitio runs
	batch.set_intermediate_structs();

	// --- run some of the gen_X methods to generate the type of run we want
	gen_noe_assignments( batch );

	// first 2 stages: enumerate pairings
	if ( (int) stage_ < (int) PURE_TOPO_RESAMPLING ) gen_enumerate_pairings( batch );

	// beta-sheet-topologies
	if ( stage_ == TOPO_RESAMPLING || stage_ == PURE_TOPO_RESAMPLING || stage_ == NOESY_PHASEII_TOPO ) gen_resample_topologies( batch );

	// reuse fragments and restart from stage2 structures
	if ( stage_ == STAGE2_RESAMPLING || stage_ == NOESY_PHASEII_S2_RESAMPLING ) {
		gen_resample_stage2( batch );
	}

	if ( stage_ >= STAGE2_RESAMPLING ) {
		gen_resample_fragments( batch );
	}

	bool result_is_fullatom = false;
	// close loops - fullatom relax
	if ( stage_ == CEN2FULLATOM ) {
		if ( repeat_id == 1  ) {
			gen_cen2fullatom( batch );
			if ( option[ iterative::safety_hatch_scorecut ].user() ) repeat_id = 3; //we subtract 1 from repeat_id at end of sub-routine and thus return 2.
		} else if ( repeat_id == 2 ) {
			gen_cen2fullatom_non_pool_decoys( batch );
			--repeat_id;
		}
		// gen_resample_fragments( batch );
		batch.set_intermediate_structs( false ); //otherwise the intermediate (centroid) structures will be scored by score_13_envhb
		result_is_fullatom = true;
	}
	mem_tr.Debug << "before evaluation output" << std::endl;
	//finalize batch
	gen_evaluation_output( batch, result_is_fullatom );
	gen_dynamic_patches( batch );
	mem_tr.Debug << "evaluation output" << std::endl;

	basic::show_time( tr,  "finalize batch..." );
	test_broker_settings( batch );
	return repeat_id-1;
}

void IterativeBase::generate_batch() {
	//OBSOLET cluster();
	basic::show_time( tr,  "generate_batch" );
	//initialize batch
	mem_tr << "IterativeBase::start_new_batch " << std::endl;
	Size repeat=1;
	while ( repeat ) {
		jd2::archive::Batch& batch( manager().start_new_batch() );
		if ( basic::options::option[ basic::options::OptionKeys::iterative::nolazy_noesy_reassign ]() ) {
			reassign_noesy_data( batch );
		}
		repeat = generate_batch( batch, repeat );
		manager().finalize_batch( batch );
		basic::show_time( tr,  "finalized batch" );
		tr.Info << std::endl;
		// don't want to reset counters too often... if we run out of steam the QUEUE EMPTY pathway will make sure that we do more runs
		mem_tr << "IterativeBase::generated_batch " << std::endl;
		//now it is best time to do this... JobQueue is definitely filled up....
		if ( !basic::options::option[ basic::options::OptionKeys::iterative::nolazy_noesy_reassign ]() ) {
			reassign_noesy_data( batch );
		}
	}
}

/// ============================================================================
/// -----------           methods to make new batches               ------------
/// -----          each method should only "append" to broker and flags    -----
/// ============================================================================


void IterativeBase::do_dynamic_patching( jd2::archive::Batch& batch, utility::io::ozstream& flags, std::string score, utility::options::FileVectorOptionKey const& key ) const {
	if ( ! option[ key ].user() ) return; //if this wasn't patched it doesn't require dynamic balancing...

	core::Real var_score(0); //can't change score so keep this as reference
	if ( get_weight("score")>0.001 ) {
		var_score = score_variation("score");
	} else if ( get_weight("prefa_clean_score3") > 0.001 ) {
		var_score = score_variation("prefa_clean_score3");
	} else {
		tr.Fatal << "Neither 'score' nor 'prefa_clean_score3' has a weight greater than 0.001" << std::endl;
		utility_exit_with_message("Must have either score or prefa_clean_score3.");
	}

	typedef std::map< scoring::ScoreType, std::string > ScoreTypeMap;
	ScoreTypeMap dynamic_scores;
	dynamic_scores[ scoring::rdc ] ="rdc";
	dynamic_scores[ scoring::atom_pair_constraint ] ="atom_pair_constraint";

	//get patched scorefxn
	scoring::ScoreFunctionOP scorefxn  = scoring::ScoreFunctionFactory::create_score_function( score, option[ key ]()[1] );
	std::ostringstream patches;
	bool have_patch( false );
	for ( ScoreTypeMap::const_iterator it = dynamic_scores.begin(); it != dynamic_scores.end(); ++it ) {
		if ( scorefxn->get_weight( it->first ) > 0.00001 ) {
			core::Real var_rel;
			if ( it->second == "atom_pair_constraint" ) {
				core::Real var_sum( 0 ); core::Size ct_cst( 0 );
				WeightMap const& variations( score_variations() );
				for ( auto const & variation : variations ) {
					if ( variation.first.find( "filter_cst" ) != std::string::npos ) {
						var_sum += variation.second;
						++ct_cst;
					}
				}
				runtime_assert( ct_cst );
				var_sum /= ct_cst;
				var_rel=var_sum/var_score;
			} else { //not a atom_pair_constraint
				var_rel=score_variation( it->second );
				var_rel/=var_score;
			}
			core::Real patch = scorefxn->get_weight( it->first )/var_rel;
			patches << it->second << " = " << patch << std::endl;
			have_patch = true;
		} //has positive weight in patched scorefunction
	}
	if ( have_patch ) {
		std::string patch_file ( batch.dir()+"/dynamic_abinitio_"+score+".patch" );
		utility::io::ozstream patch_fd( patch_file );
		patch_fd << patches.str() << std::endl;
		tr.Info << "dynamic patching: write to flag-file " << key.id() << " " << patch_file << std::endl;
		flags << "-" << key.id() << " " << patch_file << std::endl;
		//  std::cerr << "dynamic patching for " << key.id() << std::endl;
	}
}

void IterativeBase::gen_dynamic_patches( jd2::archive::Batch& batch ) {
	if ( !use_dynamic_weights_for_sampling_ ) return;

	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	do_dynamic_patching( batch, flags, "score0",  OptionKeys::abinitio::stage1_patch  );
	do_dynamic_patching( batch, flags, "score1",  OptionKeys::abinitio::stage2_patch  );
	do_dynamic_patching( batch, flags, "score2",  OptionKeys::abinitio::stage3a_patch  );
	do_dynamic_patching( batch, flags, "score5",  OptionKeys::abinitio::stage3b_patch  );
	do_dynamic_patching( batch, flags, "score3",  OptionKeys::abinitio::stage4_patch  );
	do_dynamic_patching( batch, flags, "score12", OptionKeys::score::patch  );
}


void IterativeBase::gen_evaluation_output( jd2::archive::Batch& batch, bool fullatom ) {
	tr.Info << "gen_evaluation_output in " << name() << std::endl;
	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	flags << "-evaluation:jump_nr" << std::endl; //add JumpNrEvaluator

	//add score evaluator
	if ( !fullatom && !super_quick_relax_of_centroids_ ) { //centroid
		// the standard 'score_final' score column
		flags << "-evaluation:extra_score " << option[ iterative::cen_score]() << std::endl
			<< "-evaluation:extra_score_column _final" << std::endl;
		flags << "-evaluation:extra_score_patch " << option[ iterative::cen_score_patch ]() << std::endl;
	} else { //fullatom
		if ( super_quick_relax_of_centroids_ ) {
			//obtain 'empty' score for 'score_final' and put normal fullatom score into 'score_fa'
			flags << "-evaluation:extra_score " << "empty " << option[ iterative::fa_score]() << std::endl
				<< "-evaluation:extra_score_column _final _fa" << std::endl;
			flags << "-evaluation:extra_score_patch " << option[ iterative::fa_score_patch ]() << " NOPATCH " << std::endl;
		} else {
			//score-final standard fullatom mode
			flags << "-evaluation:extra_score " << option[ iterative::fa_score]() << std::endl
				<< "-evaluation:extra_score_column _final" << std::endl;
			flags << "-evaluation:extra_score_patch " << option[ iterative::fa_score_patch ]() << std::endl;
		}
	}
	flags << "-evaluation:built_in_rdc rdc_filt" << std::endl;

	//setup relax sequence protocol for super-quick relax
	if ( super_quick_relax_of_centroids_ ) {
		if ( !fullatom ) {
			if ( option[ iterative::super_quick_relax_protocol ].user() ) {
				flags<< "-relax:sequence_file " <<  option[ iterative::super_quick_relax_protocol ]() << std::endl;
			} else {
				utility::io::ozstream sequence_file( batch.dir()+"/super_quick_relax.txt" );
				sequence_file << "ramp_repack_min 0.02  0.01"      << std::endl;
				sequence_file << "ramp_repack_min 0.250 0.01"      << std::endl;
				sequence_file << "ramp_repack_min 0.550 0.01"      << std::endl;
				sequence_file << "ramp_repack_min 1     0.00001"   << std::endl;
				sequence_file << "accept_to_best"                  << std::endl;
				flags << "-relax:sequence_file " << batch.dir()+"/super_quick_relax.txt" << std::endl;
			}
			flags << "-relax:sequence" << std::endl;
			flags << "-out:user_tag centroid" << std::endl;
			//   utility::io::ozstream patch_file( batch.dir()+"/super_quick_relax.patch" );
			//   patch_file << "atom_pair_constraint = 0.1" << std::endl;
			//   patch_file << "rdc = 0.1" << std::endl;
			if ( option[ iterative::super_quick_relax_patch ].user() ) {
				flags << "-score:patch " << option[ iterative::super_quick_relax_patch ]() << std::endl;
			}
		} else { //!fullatom
			flags << "-out:user_tag fullatom" << std::endl;
		}
	} //end super-quick relax

	//NOESY FILTER Restraints
	if ( noesy_assign::NoesyModule::cmdline_options_activated() ) {
		update_noesy_filter_files( current_noesy_sampling_file_, fullatom );

		flags << "-evaluation::constraints " << ( fullatom ? first_noesy_fa_cst_file_ : first_noesy_cst_file_ ) << std::endl;
		flags << "-evaluation::constraints_column " << ( evaluate_local() ? "noesy_autoassign_cst_worker" : "noesy_autoassign_cst" ) << std::endl;
	}

	//compute pool_convergence_XXX ?
	if ( min_diversity_list_[ stage() ] > 0 ) {
		gen_diversity_pool( batch, fullatom );
	} //setup of min-divergence pool
} //gen_evaluation


void IterativeBase::gen_diversity_pool( jd2::archive::Batch& batch, bool ) {

	utility::io::ozstream flags( batch.flag_file(), std::ios::app );

	//copy pool to '<batch>/pool.in'
	io::silent::SilentFileOptions opts;
	io::silent::SilentFileData sfd( opts );
	for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
		sfd.add_structure( *it ); //not a copy, only adds OP to sfd
	}
	sfd.write_all( batch.dir()+"/pool.in" );

	//empty pool -- write empty file for score-column consistency
	if ( sfd.size() == 0 ) {
		// write a pool anyway so that we have the tags in the returning structures and don't mess
		// up the columns in the pool decoys.out file... not crucial but makes post-analysis with scripts easier.
		utility::io::ozstream empty_pool_file( batch.dir()+"/pool.in" );
		empty_pool_file << "REMARK no structure in pool " << std::endl;
	}

	//setup flags for min-divergence-pool: ----------------

	//common flags
	flags << "-mc:known_structures " << batch.dir() << "/pool.in" << std::endl;
	flags << "-mc:max_rmsd_against_known_structures " << std::max( 0.0, min_diversity_list_[ stage() ] - 0.25 ) << std::endl;

	//exclude residues from RMSD that are flexible anyhow
	if ( basic::options::option[ basic::options::OptionKeys::in::replonly_residues ].user() ) {
		flags << "-mc:excluded_residues_from_rmsd ";
		utility::vector1<Size> replonly_rsd = basic::options::option[ basic::options::OptionKeys::in::replonly_residues ]();
		for ( utility::vector1<Size>::const_iterator it = replonly_rsd.begin(); it != replonly_rsd.end(); ++it ) {
			flags << " "<< *it;
		}
		flags << std::endl;
	}
}
/// @brief in the comp. modelling protocol the topo-resampling stage might also contain a RigidChunkClaimer...
/// provide start-structures for this as -in:file:silent
void IterativeBase::gen_start_structures( jd2::archive::Batch& /*batch*/ ) {
	// OBSOLETE
	//  batch.set_has_silent_in();
	//  io::silent::SilentFileData sfd;
	//  for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
	//   sfd.add_structure( *it ); //not a copy, only adds OP to sfd
	//  }
	//  sfd.write_all( batch.silent_in() );

	//  utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	//  broker << "\nUSE_INPUT_POSE\n" << std::endl;
	//  broker << "CLAIMER RigidChunkClaimer" << std::endl
	//      << "REGION_FILE " << option[ iterative::fix_core ]() << std::endl
	//      << "END_CLAIMER" << std::endl;

	//  //compute nstruct such that we get the usual amount of total structures
	//  batch.nstruct() = std::max( 1, int( 1.0*batch.nstruct() / ( 1.0*decoys().size() ) ) );
}

//creation of the random-ss samplign stage
void IterativeBase::gen_enumerate_pairings( jd2::archive::Batch& batch ) {

	//after Naccept accepted structures stop doing the enumerated pairing.
	//sounds OBSOLETE TO ME: if ( (int) total_accepts() > option[ iterative::enumerate::Naccept ]() ) return;

	//skip every second time (standard setting in flag-lib: Apr 2012)
	if ( option[ iterative::enumerate::skip_half ]() && bEnumeratedLastTime_ ) {
		bEnumeratedLastTime_ = false;
		return;
	}
	bEnumeratedLastTime_ = true;

	tr.Info << "enumerate pairings\n ";
	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	runtime_assert( broker.good() );

	//look for an explicit broker file for enumeration stage...
	if ( option[ iterative::enumerate::broker ].user() ) {
		utility::io::izstream enum_broker( option[ iterative::enumerate::broker ]() );
		if ( !enum_broker.good() ) {
			throw ( CREATE_EXCEPTION(utility::excn::FileNotFound, "-iterative::enumerate::broker: File "
				+std::string( option[ iterative::enumerate::broker ]())+" not found! ") );
		}
		std::string line;
		while ( getline( enum_broker, line ) ) broker << line << std::endl;
	} else {
		//no explicit enumerate file is given -- generate enumeration automatically from ss-content of fragments
		std::string frag_ss_file( batch.dir() + "psipred_ss2.dat" );
		std::string pairings_file( batch.dir() + "pairings_guess.dat" );

		using namespace core::fragment;
		FragSetOP frags_s = FragmentIO(
			option[ OptionKeys::abinitio::number_3mer_frags ]()
			).read_data( option[ in::file::frag3 ] );
		guess_pairings_from_secondary_structure( *frags_s, pairings_file, frag_ss_file );
		broker << "CLAIMER TemplateJumpClaimer\n"
			<< "RANDOM_SHEETS 2\n"
			<< "SS_INFO " << frag_ss_file << "\n"
			<< "PAIRING_FILE " << pairings_file << "\n"
			<< "END_CLAIMER" << std::endl;
	} //else (not enumerate::broker)
}

/// @brief figure out beta-sheet topologies from pooled decoys and run with jumping
void IterativeBase::gen_resample_topologies( jd2::archive::Batch& batch) {
	if ( !bDoBetaJumping_ ) return;
	tr.Info << "resample topologies\n";


	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	broker << "\nCLAIMER TemplateJumpClaimer \n"
		<< "NO_USE_INPUT_POSE\n"
		<< "topol_file "<< batch.dir() << "beta.top\n"
		<< "END_CLAIMER\n\n" << std::endl;
	broker.close();

	//read beta-topology from start file -fix_topol --- use this file for first recompute_beta_Naccept decoys
	if ( option[ iterative::initial_beta_topology ].user() &&
			(int) total_accepts() < option[ iterative::recompute_beta_Naccept ]() ) {
		//copy file to batch.dir() + "beta.top"
		std::string input( option[ iterative::initial_beta_topology ]() );

		utility::io::izstream beta_input( input );
		if ( beta_input.good() ) {
			std::string line;
			utility::io::ozstream beta_batch( batch.dir()+"beta.top" );
			while ( getline( beta_input, line ) ) beta_batch << line << std::endl;
		} else {
			utility_exit_with_message("topology_file "+input+" not found. supply file or remove option -iterative::initial_beta_topology");
		}
	} else { //generate beta-topology from decoys
		PairingStatisticsOP beta_topol = compute_beta_topology();
		utility::io::ozstream file( batch.dir() + "beta.top" );
		file << *beta_topol << std::endl;
	}
}

/// @brief restart runs from stage2-structures that correspond to those in the pool
void IterativeBase::gen_resample_stage2( jd2::archive::Batch& batch ) {
	tr.Info << "resample_stage2 " << std::endl;
	mem_tr << "IterativeBase::gen_resample_stage2 start" << std::endl;
	SilentStructVector start_decoys;

	//collect to decoys() the corresponding stage2 decoys and store in start_decoys
	collect_alternative_decoys( decoys(), "decoys_stage2.out", start_decoys );


	//  typedef std::map< std::string, utility::vector1< std::string > > SourceFiles;
	//  typedef std::map< std::string, utility::vector1< core::io::silent::SilentStructOP > > AlternativeDecoys;

	//  SourceFiles sources;
	//  AlternativeDecoys alternative_decoys;
	//  Size ct_in( 0 );

	//  //to find the stage2 structures collect first all tags for a specific file
	//  for ( const_decoy_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
	//   runtime_assert( (*it)->has_comment( TAG_IN_FILE ) );
	//   std::string tag( (*it)->get_comment( TAG_IN_FILE ) );
	//   utility::file::FileName file( (*it)->get_comment( SOURCE_FILE ) );
	//   std::string stage2_file( file.path()+file.base() + "_stage2." + file.ext() );

	//   //creates map <filename> <list of tags>
	//   sources[ stage2_file ].push_back( tag );
	//   alternative_decoys[ stage2_file ].push_back( (*it) );
	//   ++ct_in;
	//  }

	//  //read selected structures from each file
	//  Size ct_read( 0 );
	//  io::silent::SilentStructOPs start_decoys;
	//  for ( SourceFiles::const_iterator it = sources.begin(); it != sources.end(); ++it ) {
	//   /// it->first is filename, it->second are all tags collected for this file
	//   io::silent::SilentFileData sfd;
	//   try { //read structures
	//    sfd._read_file( it->first, it->second, true /*throw exceptions */ );
	//    if ( sfd.size() != it->second.size() ) {
	//     tr.Warning << "multiple decoys with same tag detected in file " << it->first << std::endl;
	//    }
	//    copy( sfd.begin(), sfd.end(), std::back_inserter( start_decoys ) );
	//    ct_read += sfd.size();
	//   } catch ( utility::excn::EXCN_IO& excn ) { //ERROR
	//    tr.Warning << "Problem reading silent-file " << it->first << " for " << it->second.size() << " structures " << std::endl;
	//    excn.show( tr.Warning );
	//    tr.Warning << std::endl;
	//    tr.Warning << "use the respective structures in the pool as starting structure instead" << std::endl;
	//    copy( alternative_decoys[ it->first ].begin(), alternative_decoys[ it->first ].end(), std::back_inserter( start_decoys ) );
	//    ct_read += alternative_decoys[ it->first ].size();
	//   }
	//  }

	//  tr.Debug << "structures from pool" << ct_in << " structure retrieved from stage2-files "
	//       << ct_read << " start structs: " << start_decoys.size() << std::endl;
	//  if ( start_decoys.size() != decoys().size() ) {
	//   tr.Warning << "why do we have a different number of decoys in pool and start_decoys ? " << std::endl;
	//  }

	///write flags and broker-file
	if ( start_decoys.size() ) {
		numeric::random::random_permutation( start_decoys, numeric::random::rg() ); //permute to get rid of high fluctuation in acceptance rate
		batch.set_has_silent_in();
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		for ( SilentStructVector::const_iterator
				it = start_decoys.begin(); it != start_decoys.end(); ++it ) {
			sfd.add_structure( **it );
		}
		sfd.write_all( batch.silent_in() );

		//make broker file for RESTART calculations
		utility::io::ozstream broker( batch.broker_file(), std::ios::app );
		broker << "\nUSE_INPUT_POSE\n"
			<< "CLAIMER StartStructClaimer\n"
			<< "PERTURB " << option[ iterative::perturb_resampling ] << std::endl
			<< "END_CLAIMER\n\n"
			<< "CLAIMER FragmentJumpClaimer\n"
			<< "END_CLAIMER\n\n" << std::endl;

		//skip stage 1
		utility::io::ozstream flags( batch.flag_file(), std::ios::app );
		flags << "-abinitio:skip_stages 1 " << std::endl;
		flags << "-out:no_nstruct_label" << std::endl;
		flags << "-out:suffix " << "_b" << std::setfill('0') << std::setw(4) << batch.id() << std::endl;
	} else {
		tr.Warning << "no stage2 decoys found ! " << std::endl;
	}

	//compute nstruct such that we get the usual amount of total structures
	batch.nstruct() = std::max( 1, int( 1.0*batch.nstruct() / ( 1.0*start_decoys.size() ) ) );
	mem_tr << "IterativeBase::gen_resample_stage2 end" << std::endl;
}


void IterativeBase::gen_resample_fragments( jd2::archive::Batch& batch ) {
	using namespace core::fragment;
	ConstantLengthFragSet frags_9mer( 9 ); //steal good old 9mers
	ConstantLengthFragSet frags_3mer( 3 ); //steal good old 9mers
	mem_tr << "IterativeBase::gen_resample_fragments start" << std::endl;

	Size ct( 1 );
	Size const max_frags( 500 );
	if ( option[ iterative::mix_frags ]() ) {
		FragSetOP frags_l = FragmentIO(
			option[ OptionKeys::abinitio::number_9mer_frags ]()
			).read_data( option[ in::file::frag9 ] );
		FragSetOP frags_s = FragmentIO(
			option[ OptionKeys::abinitio::number_3mer_frags ]()
			).read_data( option[ in::file::frag3 ] );
		frags_9mer.add( *frags_l );//assuming we really have read 9mers and 3mers.
		frags_3mer.add( *frags_s );
	}
	for ( const_decoy_iterator it = decoys().begin(); it != decoys().end() && ct <= max_frags; ++it, ++ct ) {
		pose::Pose pose;
		//std::string tag;// = it->decoy_tag();
		(*it)->fill_pose( pose );
		steal_constant_length_frag_set_from_pose( pose, frags_9mer );
		steal_constant_length_frag_set_from_pose( pose, frags_3mer );
	}
	std::string file_9mer( batch.dir() + "frags_9mer.dat" ); //using here gzipped files causes weird seg-fault on BG/P
	std::string file_3mer( batch.dir() + "frags_3mer.dat" );
	FragmentIO().write_data( file_9mer, frags_9mer );
	FragmentIO().write_data( file_3mer, frags_3mer );

	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	flags << "-frag3 " << file_3mer << std::endl
		<< "-frag9 " << file_9mer << std::endl
		<< "-abinitio:number_3mer_frags 0" << std::endl
		<< "-abinitio:number_9mer_frags 0" << std::endl;

	mem_tr << "IterativeBase::gen_resample_fragments end" << std::endl;
}


void IterativeBase::update_noesy_filter_files(
	std::string const& current,
	bool fullatom
) {
	current_noesy_sampling_file_ = current;
	if ( first_noesy_cst_file_ == "n/a" ) {
		if ( super_quick_relax_of_centroids_ ) {
			first_noesy_cst_file_ = current_noesy_sampling_file_+".filter";
		} else {
			first_noesy_cst_file_ = current_noesy_sampling_file_+".filter.centroid";
		}
	}
	if ( fullatom && first_noesy_fa_cst_file_ == "n/a" ) {
		if ( option[ iterative::never_update_noesy_filter_cst ]() ) {
			first_noesy_fa_cst_file_ = first_noesy_cst_file_.substr(0,first_noesy_cst_file_.find(".centroid"));
		} else {
			first_noesy_fa_cst_file_ = current_noesy_sampling_file_+".filter";
		}
	}
	tr.Info << "update_noesy_filter_files for " << name() << std::endl;
	tr.Info << "STATUS now: ";
	save_status( tr.Info );
	tr.Info << std::endl;
}

//generate ConstraintClaimer files for NOE restraints
void IterativeBase::gen_noe_assignments( jd2::archive::Batch& batch ) {
	tr.Info << "gen_noe_assignments in " << name() << std::endl;
	if ( !noesy_assign::NoesyModule::cmdline_options_activated() ) return;
	bool bCombine( bCombineNoesyCst_ );
	mem_tr << "IterativeBase::gen_noe_assignments start" << std::endl;

	/* OBSOLETED by update_noesy_filter_files
	//if just started IterativeFullatom we are in a new instance and the first_noesy_cst_file and
	// the current_noesy_sampling_file_ are not set .. go through backward through past batches and take first (i.e., youngest) file.
	if ( batch.id() > 1 && stage() > CEN2FULLATOM && first_noesy_cst_file_ == "n/a"  ) {
	tr.Info << "Figure out current NOESY CST File by going backwards from " << batch.id() << std::endl;
	for ( Size back = batch.id() - 1; back >= 1; --back ) {
	Batch last_batch( back );
	runtime_assert( last_batch.id() < batch.id() );
	current_noesy_sampling_file_ = last_batch.dir()+"/"+NOESY_CST_FILE_NAME;
	if ( utility::file::file_exists( current_noesy_sampling_file_ ) ) break;
	current_noesy_sampling_file_ = option[ iterative::initial_noe_auto_assign_csts ]();
	}
	//  if ( option[ iterative::update_noesy_filter_cst ]() ) {
	first_noesy_cst_file_ = current_noesy_sampling_file_+".filter";
	first_noesy_fa_cst_file_ = current_noesy_sampling_file_+".filter";
	//} else {
	//first_noesy_cst_file_ = option[ iterative::initial_noe_auto_assign_csts ]()+".filter";
	//first_noesy_fa_cst_file_ = option[ iterative::initial_noe_auto_assign_csts ]()+".filter";
	//}
	}
	*/
	//
	std::string cst_file( current_noesy_sampling_file_ );
	std::string best_cst_file( cst_file + ".good");
	std::string cst_centroid_file( cst_file + ".centroid");
	std::string best_cst_centroid_file( best_cst_file + ".centroid");

	//make broker files with FULL-ATOM and CENTROID CLaimers
	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	bool bSplitted( option[ iterative::split_autoNOE_restraints ]() );

	std::string drop_candidates( cst_file + ".drop" );
	std::string drop_candidates_centroid( drop_candidates + ".centroid" );

	using namespace utility::file; //for file_exists
	bool bRandomRestraints( option[ iterative::randomize_elimination_candidates ] > 0.0 && file_exists( drop_candidates ) );


	if ( stage() >= CEN2FULLATOM ) { //don't add fa constraints in early stages... slows down super-quick relax to much...
		if ( bSplitted ) {
			broker << "\nCLAIMER ConstraintClaimer \n"
				<< "CST_FILE " << best_cst_file << "\n"
				<< "NO_CENTROID\n"
				<< "FULLATOM\n"
				<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
			broker << "COMBINE_RATIO 1" << std::endl;
			broker << "END_CLAIMER\n" << std::endl;
		}
		broker << "\nCLAIMER ConstraintClaimer \n"
			<< "CST_FILE " << cst_file << "\n"
			<< "NO_CENTROID\n"
			<< "FULLATOM\n"
			<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
		if ( bCombine || bSplitted ) {
			broker << "COMBINE_RATIO " << option[ basic::options::OptionKeys::constraints::combine ]() << std::endl;
		}
		broker << "END_CLAIMER\n" << std::endl;
		if ( bRandomRestraints ) {
			broker << "\nCLAIMER ConstraintClaimer \n"
				<< "CST_FILE " << drop_candidates << "\n"
				<< "NO_CENTROID\n"
				<< "FULLATOM\n"
				<< "DROP_RANDOM_RATE " << option[ iterative::randomize_elimination_candidates ]() << "\n"
				<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
			if ( bCombine || bSplitted ) {
				broker << "COMBINE_RATIO " << option[ basic::options::OptionKeys::constraints::combine ]() << std::endl;
			}
			broker << "END_CLAIMER\n" << std::endl;
		} //bRandom
	}
	// if ( bCombine ) broker << "COMBINE_RATIO " << basic::options::option[ basic::options::OptionKeys::constraints::combine ]() << "\n";
	if ( bSplitted ) {
		broker << "\nCLAIMER ConstraintClaimer \n"
			<< "CST_FILE " << best_cst_centroid_file << "\n"
			<< "CENTROID\n"
			<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
		broker << "COMBINE_RATIO 1" << std::endl;
		broker << "END_CLAIMER\n" << std::endl;
	}
	broker << "\nCLAIMER ConstraintClaimer \n"
		<< "CST_FILE " << cst_centroid_file << "\n"
		<< "CENTROID\n"
		<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
	if ( bCombine || bSplitted ) {
		broker << "COMBINE_RATIO " << option[ basic::options::OptionKeys::constraints::combine ]() << std::endl;
	}
	broker << "END_CLAIMER\n" << std::endl;

	if ( bRandomRestraints ) {
		broker << "\nCLAIMER ConstraintClaimer \n"
			<< "CST_FILE " << drop_candidates_centroid << "\n"
			<< "CENTROID\n"
			<< "DROP_RANDOM_RATE " << option[ iterative::randomize_elimination_candidates ]() << "\n"
			<< "SKIP_REDUNDANT "<<option[ iterative::skip_redundant_constraints ]() << "\n";
		if ( bCombine || bSplitted ) {
			broker << "COMBINE_RATIO " << option[ basic::options::OptionKeys::constraints::combine ]() << std::endl;
		}
		broker << "END_CLAIMER\n" << std::endl;
	} //bRandom

	mem_tr << "IterativeBase::gen_noe_assignments end" << std::endl;
}


void IterativeBase::gen_cen2fullatom( jd2::archive::Batch& batch ) {
	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	broker << "\nUSE_INPUT_POSE\n"
		<< "CLAIMER StartStructClaimer\n"
		<< "END_CLAIMER\n\n"
		<< "CLAIMER FragmentJumpClaimer\n"
		<< "END_CLAIMER\n\n" << std::endl;


	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	flags << "-abinitio::skip_stages 1 2 3 4" << std::endl;

	add_fullatom_flags( batch );
	io::silent::SilentStructOPs start_decoys;
	std::copy( decoys().begin(), decoys().end(), std::back_inserter( start_decoys ) );

	if ( start_decoys.size() ) {
		batch.set_has_silent_in();
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		for ( core::io::silent::SilentStructOPs::const_iterator
				it = start_decoys.begin(); it != start_decoys.end(); ++it ) {
			sfd.add_structure( **it );
		}
		sfd.write_all( batch.silent_in() );
	}

	batch.nstruct() = std::max( 1, int( 1.0*batch.nstruct() / ( 1.0*start_decoys.size() ) ) );

}


void IterativeBase::collect_hedgeing_decoys_from_batches(
	jd2::archive::Batch const& batch,
	core::io::silent::SilentStructOPs& start_decoys,
	Real score_cut_per_batch
) {
	using namespace core::io::silent;

	//count total decoys
	Size total( 0 );
	for ( auto const & it : manager().batches() ) {
		if ( it.id() >= first_fullatom_batch_ ) break;
		if ( !it.has_silent_in() ) continue;
		total += it.decoys_returned();
		tr.Debug << "harvest old decoys for safety hatch: batch " << it.id() << " " << it.decoys_returned() << " " << total << std::endl;
	}
	basic::show_time( tr,  "generate safety_hatch: counted total decoys" );

	for ( auto const & it : manager().batches() ) {
		Real percentage_per_batch( 1.0*batch.nstruct() / (1.0*total) );
		if ( it.id() >= first_fullatom_batch_ ) break;
		if ( !it.has_silent_in() ) continue; //usually only the resampling decoys are interesting...
		if ( !it.decoys_returned() ) continue; //avoid looking for empty files
		//  it->silent_out();
		basic::show_time( tr,  "generate safety_hatch: access batch "+it.batch() );
		SilentFileOptions opts;
		SilentFileData sfd( opts );
		std::list< std::pair< core::Real, SilentStructOP > > score_cut_decoys;
		Size ct( 0 );
		tr.Debug << "read and score decoys in " << it.silent_out() << "..." << std::endl;
		sfd.read_file( it.silent_out() );
		for ( SilentFileData::iterator sit=sfd.begin(), esit=sfd.end(); sit!=esit; ++sit ) {
			//std::string tag = sit->decoy_tag();
			sit->set_decoy_tag( "harvest_"+batch.batch()+"_"+ObjexxFCL::lead_zero_string_of( ++ct, 6 ) );

			//note this does nothing but return *it, if b_evaluate_incoming_decoys_ is false
			//Solve this by having another Archive that is once and only once filled with decoys from all relevant batches.
			//This can already be done when the batches are being made. the only difference is that the rule for acceptance is different..
			//the archive has no upper-size limit and always takes the bestX% of all structures.
			//it has to be maintained with updates of evaluators and such together with the main-archive.
			//but one could also choose to *not* update evaluators to have a more different (more hedging approach. )
			//this should probably be renamed into hedging !!!
			core::io::silent::SilentStructOP pss = *sit;//cannot afford this ever: evaluate_silent_struct( *sit );
			score_cut_decoys.push_back( std::make_pair( select_score( pss ), pss ) );
		}
		score_cut_decoys.sort();
		tr.Debug << "select " << percentage_per_batch*100 << "% from batch from the lowest scoring " << score_cut_per_batch*100 << "% of structures" << std::endl;
		basic::show_time( tr,  "generate safety_hatch: collected decoys batch "+it.batch() );
		// if we have less structures below score cut than what we want to harvest,.... take them all
		while ( score_cut_per_batch < percentage_per_batch ) {
			Size ind_max( static_cast< Size > ( score_cut_decoys.size()*score_cut_per_batch ) );
			for ( std::list< std::pair< core::Real, core::io::silent::SilentStructOP > >::const_iterator sit = score_cut_decoys.begin();
					sit != score_cut_decoys.end(); ++sit ) {
				start_decoys.push_back( sit->second );
				if ( --ind_max <= 1 ) break;
			}
			percentage_per_batch-=score_cut_per_batch;
		}
		basic::show_time( tr,  "generate safety_hatch: generated start_decoys batch "+it.batch());
		// for the remaining structures we want to harvest, they clearly will be less than what the score-cut yields... choose randomly...
		if ( percentage_per_batch > 0.01 ) {
			Size ind_max( static_cast< Size > ( score_cut_decoys.size()*score_cut_per_batch ) );
			for ( std::list< std::pair< core::Real, core::io::silent::SilentStructOP > >::const_iterator sit = score_cut_decoys.begin();
					sit != score_cut_decoys.end(); ++sit ) {
				if ( numeric::random::rg().uniform() < ( percentage_per_batch / score_cut_per_batch ) ) {
					start_decoys.push_back( sit->second );
				}
				if ( --ind_max <= 1 ) break;
			}
		}
	}
	basic::show_time( tr,  "generate safety_hatch: done collecting");
}

void IterativeBase::gen_cen2fullatom_non_pool_decoys( jd2::archive::Batch& batch ) {
	using namespace core::io::silent;
	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	broker << "\nUSE_INPUT_POSE\n"
		<< "CLAIMER StartStructClaimer\n"
		<< "END_CLAIMER\n\n"
		<< "CLAIMER FragmentJumpClaimer\n"
		<< "END_CLAIMER\n\n" << std::endl;


	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	flags << "-abinitio::skip_stages 1 2 3 4" << std::endl;
	add_fullatom_flags( batch );


	//now go thru batches and select percentage_per_batch structures randomly from the pool created by score_cut_per_batch...
	//use score_final for this ? or can we use EvaluatedArchive methods to get a reasonable score ???
	//ACHTUNG:  need also to know final centroid batch id...
	SilentStructOPs start_decoys;
	Real score_cut_per_batch( option[ OptionKeys::iterative::safety_hatch_scorecut ] );
	if ( !hedge_archive_ ) {
		tr.Info << "collect stage2 decoys from old batches ..." << std::endl;
		collect_hedgeing_decoys_from_batches( batch, start_decoys, score_cut_per_batch );
	} else {
		tr.Info << "collect stage2 decoys from hedge_archive ..." << std::endl;
		hedge_archive_->collect( batch, start_decoys );
	}

	if ( start_decoys.size() ) {
		batch.set_has_silent_in();
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd(opts);
		Size ct( 0 );
		for ( core::io::silent::SilentStructOPs::const_iterator
				it = start_decoys.begin(); it != start_decoys.end(); ++it ) {
			if ( ++ct > batch.nstruct() ) break;
			sfd.add_structure( **it );
		}
		sfd.write_all( batch.silent_in() );
	}
	basic::show_time( tr,  "generate safety_hatch: done writing");
	batch.nstruct() = std::max( 1, int( 1.0*batch.nstruct() / ( 1.0*start_decoys.size() ) ) );

}

/// -----------------------------------------------------------------------------
///         END gen_xxx ( batch ) methods
/// ==============================================================================

void IterativeBase::reassign_noesy_data( jd2::archive::Batch& batch ) {
	if ( !noesy_assign::NoesyModule::cmdline_options_activated() ) return;
	Size const n_decoys( option[ iterative::limit_decoys_for_noe_assign ] );
	if ( batch.id() == 1 || ( total_proposed() < delay_noesy_reassign_ && stage() < CEN2FULLATOM )
			|| ( stage() >= CEN2FULLATOM && total_proposed() < 2*n_decoys ) ) return; //don't do this at very beginning
	//this takes a while... make sure that backedup version of archive is up-to-date
	manager().save_archive();

	mem_tr << "IterativeBase reassign_noesy_data start" << std::endl;

	SilentStructs calibration_decoys;

	std::ostringstream hash_string;
	hash_string << "NO_POOL " << std::endl;

	//get the calibration subset of decoys...
	if ( decoys().size() < n_decoys && decoys().size() < nstruct() ) {
		calibration_decoys = decoys();
		hash_string << "ONLY "<< decoys().size() << std::endl;
	} else {
		Size ct2 ( 1 );
		for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end() && ct2 <= n_decoys; ++ct2, ++it ) {
			calibration_decoys.push_back( *it );
			hash_string << (*it)->decoy_tag() << std::endl; //create hash-string out of decoy tags
		}
	}

	//check if the low-energy decoy-tags have changed via their hash-string
	size_t hash_val( hasher( hash_string.str() ) );
	if ( hash_val == noesy_assign_hash_ && !option[ iterative::nolazy_noesy_reassign ]() ) {
		tr.Info << "do not create new noesy assignment since low-energy structures have not changed..." << std::endl;
		return; // don't do anything
	}
	noesy_assign_hash_ = hash_val;


	//manipulate params if we have some "SCHEME"

	bool use_hardwired_phaseII_params( false );
	option=vanilla_options_;
	std::string const& option_file( option[ iterative::staged_auto_noe_flags ][ stage() ] );
	if ( option_file != "NONE" ) {
		tr.Info << "Loading Auto-NOE options for stage " << stage() << " from @" << option_file << std::endl;
		using utility::file::file_exists;
		if ( !file_exists( option_file ) ) {
			utility_exit_with_message( "couldn't find flag-file "+option_file );
		}
		option.load_options_from_file( option_file );
	} else {
		tr.Info << "No stage specific Auto-NOE option file found for stage" << stage() << std::endl;
		if ( stage() > STAGE2_RESAMPLING ) {
			if ( option[ iterative::refine_auto_noe_flags ].user() ) {
				using utility::file::file_exists;
				if ( !file_exists( option[ iterative::refine_auto_noe_flags ]()  ) ) {
					utility_exit_with_message( "couldn't find flag-file "+std::string(option[ iterative::refine_auto_noe_flags ]()) );
				}

				basic::options::option.load_options_from_file( option[ iterative::refine_auto_noe_flags ] );
			} else {
				use_hardwired_phaseII_params = true;
			}
		}
	}

	noesy_assign::PeakAssignmentParameters::reset(); //re-read cmd-line options
	noesy_assign::PeakAssignmentParameters& params( *noesy_assign::PeakAssignmentParameters::get_nonconst_instance() );

	if ( use_hardwired_phaseII_params ) {
		tr.Info << "Use hardwired phaseII params for NoesyModule " << std::endl;
		params.set_cycle( 7 );
		params.use_local_distviol_ = true;
	}

	std::string scheme( option[ iterative::auto_assign_scheme ]() );
	if ( scheme == "CONST" ) {
		tr.Info << " reassign NOESY data with cmd-line settings " << std::endl;
		bCombineNoesyCst_ = stage() <= STAGE2_RESAMPLING;
	} else if ( scheme == "ADAPT1" ) {
		tr.Info << " ADAPT1 scheme selected for NOESY data: " << std::endl;
		tr.Info << "stage1-3: cycle7 - nodistViol - combine-cst " << std::endl;
		tr.Info << "stage4: cycle7 - [dcut = -iterative:dcut ] - do not combine " << std::endl;
		if ( stage() < STAGE2_RESAMPLING ) {
			params.dcut_ = -1;
			bCombineNoesyCst_ = true;
		} else {
			params.dcut_ = option[ iterative::dcut ]();
			bCombineNoesyCst_ = false;
		}
	} else {
		utility_exit_with_message("unknown auto_assign_scheme for NOESY data: choose ADAPT1 or CONST" + scheme);
	}

	//initialize noesy-module
	if ( !noesy_module_ ) {
		noesy_module_ = protocols::noesy_assign::NoesyModuleOP( new protocols::noesy_assign::NoesyModule( target_sequence_ ) );
	} else {
		noesy_module_->reset();
	}

	//dump params on tracer and in README file
	params.show_on_tracer();
	{
		utility::io::ozstream param_out( batch.dir()+"/README_noe_auto_assign", std::ios::app );
		params.show( param_out );
	}

	//assign corss-peaks using the calibration decoys
	noesy_module_->assign(
		utility::DereferenceIterator< SilentStructs >( calibration_decoys.begin() ),
		utility::DereferenceIterator< SilentStructs >( calibration_decoys.end() )
	);
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//make constraints
	core::pose::Pose aPose;  //make a pose to make constraint-files
	core::pose::make_pose_from_sequence(
		aPose,
		target_sequence_,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ))
	);
	std::string cst_file( batch.dir()+"/"+NOESY_CST_FILE_NAME ); //cst-file names
	std::string cst_centroid_file( cst_file + ".centroid");
	std::string best_cst_file( cst_file + ".good");
	std::string best_cst_centroid_file( best_cst_file + ".centroid");

	if ( reference_pose_ ) {
		noesy_module_->add_dist_viol_to_assignments( *reference_pose_ );
	}
	bool bSplitted( option[ iterative::split_autoNOE_restraints ] );
	bool bRandomRestraints( option[ iterative::randomize_elimination_candidates ] > 0.0 );
	if ( bSplitted ) {
		noesy_module_->generate_constraint_files( aPose, best_cst_file, best_cst_centroid_file, 5, 0, 2 ); //generate constraints
		noesy_module_->generate_constraint_files( aPose, cst_file, cst_centroid_file, 5, 3, 4 ); //generate constraints
	} else {
		noesy_module_->generate_constraint_files( aPose, cst_file, cst_centroid_file, 5, 0, 5, !bRandomRestraints ); //generate constraints
	}
	noesy_module_->generate_constraint_files( aPose, cst_file+".filter", cst_file+".filter.centroid", 2, 0, 5 );

	std::string drop_candidates( cst_file + ".drop" );
	std::string drop_candidates_centroid( drop_candidates + ".centroid" );

	if ( bRandomRestraints ) {
		noesy_module_->generate_constraint_files( aPose, drop_candidates, drop_candidates_centroid, 5, 0, 5, !bRandomRestraints, true );
	}

	//write NOE_out.dat
	noesy_module_->write_assignments( batch.dir()+"/NOE_out.dat" );

	//remember current cst-file
	current_noesy_sampling_file_ = cst_file;

	//switch out filter constraints --- if evaluate-local() otherwise we do this only once ( triggered from idle() )
	if ( evaluate_local() && option[ iterative::update_noesy_filter_cst ]() ) { //change filter restraints if we are local-evaluators
		topology_broker::ConstraintClaimerOP cst( new topology_broker::ConstraintClaimer( current_noesy_sampling_file_+".filter", "noesy_autoassign_cst" ) );
		cst->set_combine_ratio( bCombineNoesyCst_ ? 2 : 1 );
		cst->set_fullatom( true );
		cst->set_centroid( false );
		cst->set_skip_redundant( option[ iterative::skip_redundant_constraints ]() );
		cst->set_filter_weight( get_weight( "noesy_autoassign_cst" )/overall_cstfilter_weight_ );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new topology_broker::ConstraintEvaluatorWrapper( cst->tag(), cst ) ) ), cst->filter_weight()*overall_cstfilter_weight_ );
		rescore(); //rescore now, since we probably have more time now, when later when the decoys are arriving...
	}
	basic::prof_show();
	mem_tr << "IterativeBase reassign_noesy_data end" << std::endl;
}


void IterativeBase::guess_pairings_from_secondary_structure(
	core::fragment::FragSet const& frags,
	std::string const& out_pairings_file,
	std::string const& out_frag_ss_file
) const {
	core::fragment::SecondaryStructure ss_def( frags, true /*no JustUseCentralResidue */ );
	utility::io::ozstream out_ss( out_frag_ss_file );
	ss_def.write_psipred_ss2( out_ss, target_sequence() ); // << ss_def << std::endl;

	loops::Loops strands;
	int last_sheet_pos=-1;
	int new_sheet_start=0;

	//find all sheets longer than 3 residues and put them into "strands"
	for ( int pos=1; pos <= (int) ss_def.total_residue(); pos++ ) {
		if ( ss_def.sheet_fraction( pos ) >= 0.7 ) {
			if ( pos - 1 != last_sheet_pos ) {
				if ( new_sheet_start > 0 && (last_sheet_pos - new_sheet_start >= 3 ) ) {
					strands.push_back( new_sheet_start, last_sheet_pos );
				}
				new_sheet_start = pos;
			}
			last_sheet_pos = pos;
		}
	}
	// last detected strand is still not pushed into "strands"
	if ( new_sheet_start > 0 && (last_sheet_pos - new_sheet_start >= 3 ) ) {
		strands.push_back( new_sheet_start, last_sheet_pos );
	}
	tr.Debug << "secondary structure in current fragset reveals following beta-strands: " << std::endl << strands;

	//create pairings between all possible strand pairings
	utility::io::ozstream os(out_pairings_file);
	for ( Size i1 = 1; i1 <= strands.size(); i1++ ) {
		for ( Size i2 = i1+1; i2 <= strands.size(); i2++ ) {
			for ( Size pos1=strands[i1].start(); pos1<=strands[i1].stop(); pos1++ ) {
				for ( Size pos2=strands[i2].start(); pos2<=strands[i2].stop(); pos2++ ) {
					if ( pos2-pos1 > 25 ) {
						os << pos1 << " " << pos2 << " A 1" << std::endl;
						os << pos1 << " " << pos2 << " A 2" << std::endl;
						os << pos1 << " " << pos2 << " P 1" << std::endl;
						os << pos1 << " " << pos2 << " P 2" << std::endl;
					}
				}
			}
		}
	}
}


//this is called for non-local evaluation to switch-out the filter restraints once after stage3.
//this causes cancelling of batches and other shit, that is the reason why called rarely/once
//is not called if switching is deactivated
void IterativeBase::replace_noesy_filter_constraints() {
	tr.Info << "replacing NOE filter constraints" << std::endl;
	//are we full-atom ?
	pose::Pose pose;
	(*decoys().begin())->fill_pose( pose );
	bool const fullatom( pose.is_fullatom() );
	update_noesy_filter_files( current_noesy_sampling_file_, fullatom );

	// this needs work in EvaluatedArchive to make sure the right constraint files are used...
	if ( noesy_assign::NoesyModule::cmdline_options_activated() ) {
		std::string const cst_file( fullatom ? first_noesy_fa_cst_file_ : first_noesy_cst_file_ );
		Real weight( get_weight( "noesy_autoassign_cst" )/overall_cstfilter_weight_ );
		topology_broker::ConstraintClaimerOP cst( new topology_broker::ConstraintClaimer( cst_file, "noesy_autoassign_cst" ) );
		cst->set_combine_ratio( 2 );
		cst->set_fullatom( true );
		cst->set_centroid( false );
		cst->set_skip_redundant( option[ OptionKeys::iterative::skip_redundant_constraints ]() );
		cst->set_filter_weight( weight );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new topology_broker::ConstraintEvaluatorWrapper( cst->tag(), cst ) ) ), cst->filter_weight()*overall_cstfilter_weight_ );
	}
	rescore_nonlocal_archive();
}

//this is called for non-local evaluation to switch filter-restraints
void IterativeBase::rescore_nonlocal_archive() {
	if ( evaluate_local() ) return;

	//cancel all running batches...
	manager().cancel_batches_previous_to( manager().last_batch_id(), false /*don't allow reading of old decoys...*/ );

	//make sure that archive gets re-evaluated...
	set_evaluate_local( true );//set this temporarily

	core::scoring::ScoreFunctionOP scfxn( nullptr );
	std::string score_name;
	std::string score_patch;

	//check fullatom:
	core::pose::Pose a_pose;
	(*decoys().begin())->fill_pose( a_pose );
	bool const fullatom( a_pose.is_fullatom() );

	if ( !fullatom ) {
		score_name=option[ iterative::cen_score]();
		score_patch=option[ iterative::cen_score_patch ]();
	} else {
		score_name=option[ iterative::fa_score ]();
		score_patch=option[ iterative::fa_score_patch ]();
	}

	if ( score_patch != "NOPATCH" ) {
		scfxn = core::scoring::ScoreFunctionFactory::create_score_function( score_name, score_patch );
	} else {
		scfxn = core::scoring::ScoreFunctionFactory::create_score_function( score_name );
	}

	bool chainbreaks_in_patch( scfxn->get_weight( scoring::linear_chainbreak ) > 0.0 );
	core::Real extra_chainbreak_weight( basic::options::option[ basic::options::OptionKeys::jumps::increase_chainbreak ] );
	if ( !chainbreaks_in_patch ) { //if chainbreak is not patched .. set it
		core::Real lin_wgt( 4.0/3.0 * extra_chainbreak_weight );
		core::Real overlap_wgt( 1.0 * extra_chainbreak_weight );
		core::Real quad_wgt( 1.0 * extra_chainbreak_weight );
		tr.Info << "no chainbreaks specified in  " << score_patch << ",..."
			<< " set chainbreaks to " << lin_wgt << " and "
			<< overlap_wgt <<" for linear- and overlap-chainbreak, respectively" << std::endl;
		scfxn->set_weight( scoring::chainbreak, quad_wgt );
		scfxn->set_weight( scoring::linear_chainbreak, lin_wgt );
		scfxn->set_weight( scoring::overlap_chainbreak, overlap_wgt );
	}

	if ( super_quick_relax_of_centroids_ && stage_ < CEN2FULLATOM ) {
		std::string score_patch=option[ iterative::fa_score_patch ]();
		if ( score_patch != "NOPATCH" ) {
			scfxn = core::scoring::ScoreFunctionFactory::create_score_function( "empty", score_patch );
		} else {
			scfxn = core::scoring::ScoreFunctionFactory::create_score_function( "empty" );
		}
		// cannot compute this from fullatom-pose... keep energy term around in rescore...
		//    std::string cen_score = option[ iterative::cen_score ]();
		//    core::scoring::ScoreFunctionOP cen_scfxn( NULL );
		//    cen_scfxn = core::scoring::ScoreFunctionFactory::create_score_function( cen_score );
		//    add_evaluation( new evaluation::TruncatedScoreEvaluator( "prefa_centroid_score", selection, cen_scfxn, true /*fullname*/ ), 1.0 );

		std::string fa_score = option[ iterative::fa_score ]();
		core::scoring::ScoreFunctionOP fa_scfxn( nullptr );
		fa_scfxn = core::scoring::ScoreFunctionFactory::create_score_function( fa_score );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::ScoreEvaluator( "score_fa", fa_scfxn, true /*fullname*/ ) ) ), option[ iterative::fullatom_after_quickrelax_weight ]() );
	}

	set_scorefxn( scfxn );
	add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::ScoreEvaluator( "_final", scfxn ) ) ), 1.0 );
	rescore();
	save_to_file();
	set_scorefxn( nullptr );
	remove_evaluation("score_final");
	set_weight("score_final",1.0); // need this for further scoring
	//can only get here, if the mode is evaluate_local == false, then its setting to true was only temporarily
	set_evaluate_local( false );//have set this above to true, it was temporarily
}

void IterativeBase::add_fullatom_flags( jd2::archive::Batch& batch ) {
	if ( first_fullatom_batch_ > batch.id() ) first_fullatom_batch_ = batch.id();
	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	//these are many.... should there be some input file ?
	// flags << "@../flags_fullatom" << std::endl;
	if ( !option[ iterative::flags_fullatom ].user() ) {
		utility_exit_with_message( "option -iterative:flags_fullatom is required !!! ");
	}
	flags << "@"<< option[ iterative::flags_fullatom ]() << std::endl;

	if ( option[ constraints::cst_fa_file ].user() ) {
		utility::io::ozstream broker( batch.broker_file(), std::ios::app );
		broker << "\nCLAIMER ConstraintClaimer \n"
			<< "CMD_FLAG\n"
			<< "FULLATOM\n"
			<< "NO_CENTROID\n";
		broker << "END_CLAIMER\n"
			<< std::endl;
	}
	if ( option[ constraints::cst_file ].user() ) {
		utility::io::ozstream broker( batch.broker_file(), std::ios::app );
		broker << "\nCLAIMER ConstraintClaimer \n"
			<< "CMD_FLAG\n"
			//      << "FULLATOM\n"
			//   if ( option[ constraints::cst_fa_file ].user() ) {
			<< "CENTROID \n"
			<< "END_CLAIMER\n"
			<< std::endl;
	}
}

/// ============================================================================
/// -----------           methods to compute stuff from archive     ------------
/// -----------             these may be called from idle()         ------------
/// ============================================================================

PairingStatisticsOP IterativeBase::compute_beta_topology() {
	using namespace ObjexxFCL;
	PairingStatisticsOP beta_topol;

	//use -jumps::max_strand_gap_allowed 10 -jumps:contact_score 0.2
	if ( !option[ jumps::max_strand_gap_allowed ].user() ) {
		option[ jumps::max_strand_gap_allowed ].def( 10 );
	}
	if ( !option[ jumps::contact_score ].user() ) {
		option[ jumps::contact_score ].def( 0.2 );
	}

	beta_topol = PairingStatisticsOP( new PairingStatistics );
	PairingStatistics::ModelFreq model_freq;
	core::Size ct( 1 );
	for ( const_decoy_iterator it = decoys().begin(); it != decoys().end(); ++it, ++ct ) {
		pose::Pose pose;
		std::string tag;// = it->decoy_tag();
		(*it)->fill_pose( pose );

		// get strand pairings
		core::scoring::dssp::StrandPairingSet strand_pairings( pose );

		tag = right_string_of( ct, 4, '0');
		beta_topol->add_topology( strand_pairings, tag );
		model_freq[ tag.substr(0, 4) ] += 1;

	} // for decoys
	// set score terms
	beta_topol->compute_model_weights( model_freq );
	return beta_topol;
}


void get_core( toolbox::DecoySetEvaluation& eval, core::Real cutoff, loops::Loops& rigid ) {
	ObjexxFCL::FArray1D_double weights( eval.n_atoms(), 1.0 );
	eval.superimpose();
	eval.wRMSD( cutoff, 0.00001, weights );

	utility::vector1< Real > result;
	eval.rmsf( result );

	rigid.clear();
	for ( Size i=1; i<=result.size(); ++i ) {
		if ( result[ i ] < cutoff ) rigid.add_loop( loops::Loop(  i, i ), 5 );
	}
	tr.Debug << "make rigid with cutoff " << cutoff << "  " << std::endl << rigid << std::endl;
}


void IterativeBase::compute_cores() {
	tr.Info << "compute converged regions " << std::endl;
	basic::DynamicProfileThis prof_here("ARCHIVE_COMPUTE_CONVERGED_REGIONS");
	core::Size const opt_nstruct( option[ iterative::rmsf_nstruct ]() );
	core::Size const nstruct( std::min( int( opt_nstruct ), int( decoys().size() ) ));

	if ( nstruct==0 ) return;
	toolbox::DecoySetEvaluation eval;
	eval.reserve( nstruct );
	//Size nres( (*decoys().begin())->nres() );
	// Size ct( 1 );
	//for ( const_decoy_iterator iss = decoys().begin(); iss != decoys().end(); ++iss, ++ct ) {
	//  if ( ct > nstruct ) break;
	//  pose::Pose pose;
	//  (*iss)->fill_pose( pose );
	//  eval.push_back( pose );
	//  nres = pose.size();
	// }
	typedef toolbox::OP_const_iterator<SilentStructs::const_iterator, core::io::silent::SilentStructOP > decoy_iterator;
	eval.push_back_CA_xyz_from_silent_file( nstruct, decoy_iterator( decoys().begin() ), decoy_iterator( decoys().end() ), false );

	loops::Loops old_core15 = core15_;
	loops::Loops old_core2 = core2_;
	loops::Loops old_core3 = core3_;
	loops::Loops old_core4 = core4_;

	get_core( eval, 1.5, core15_ );
	get_core( eval, 2, core2_ );
	get_core( eval, 3, core3_ );
	get_core( eval, 4, core4_ );

	if ( old_core15 != core15_ ) add_core_evaluator( core15_, "_core15" );
	if ( old_core2 != core2_ ) add_core_evaluator( core2_, "_core2" );
	if ( old_core3 != core3_ ) add_core_evaluator( core3_, "_core3" );
	if ( old_core4 != core4_ ) add_core_evaluator( core4_, "_core4" );
	tr.Info << "finished computing converged regions" << std::endl;
}

void IterativeBase::add_core_evaluator( loops::Loops const& core, std::string const& core_tag ) {
	utility::vector1< Size> selection;
	core.get_residues( selection );
	if ( reference_pose_ ) add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::SelectRmsdEvaluator( reference_pose_, selection, core_tag ) ) ) );
	core.write_loops_to_file( name()+"/"+core_tag+".gen.rigid", "RIGID" ); //so we have them for other evaluations
}

void IterativeBase::restore_status( std::istream& is ) {
	Parent::restore_status( is );
	hedge_archive_ = HedgeArchiveOP( new HedgeArchive( name()+"_hedge" ) );
	hedge_archive_->restore_from_file();
	int bla; std::string tag;
	is >> tag >> bla;
	runtime_assert( tag == "IterationStage:" );
	stage_ = IterationStage( bla );
	tr.Info << "restored iteration stage: " << stage_ << std::endl;
	is >> tag >> bla;
	runtime_assert( tag == "first_batch_this_stage:");
	first_batch_this_stage_ = bla;
	is >> tag >> bla;
	runtime_assert( tag == "first_fullatom_batch:" );
	first_fullatom_batch_ = bla;
	compute_cores();
	is >> tag;
	if ( is.good() && tag == "SCORED_CORE:" ) {
		//OBSOLETE TAG: read anyway for backward compatibility
		loops::PoseNumberedLoopFileReader reader;
		reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, name()+"STATUS file", false /*no strict checking */ );
		tr.Warning << "found obsolete tag SCORE_CORE in status file" << std::endl;
	}

	if ( basic::options::option[ iterative::never_update_noesy_filter_cst ]() ) {
		//bugfix because of new sequence of events the fullatom pool gets initialized with the wrong filter restraints
		//should get rid of all these different options, as we have now worked out the optimal algorithm
		//but for now, a hack at this point

		//  this will update the fullatom_pool to have the initial sampling restraints rather than an n/a
		update_noesy_filter_files( current_noesy_sampling_file_, true );
	}

	is >> tag;
	if ( is.good() && tag == "NOESY_CYCLE:" ) {
		is >> noesy_assign_float_cycle_;
	}
	is >> tag;
	if ( is.good() && tag == "NOESY_FIRST_CST:" ) {
		is >> first_noesy_cst_file_;
	}
	is >> tag;
	if ( is.good() && tag == "NOESY_FIRST_FA_CST:" ) {
		is >> first_noesy_fa_cst_file_;
	}
	is >> tag;
	if ( is.good() && tag == "NOESY_CURRENT_CST:" ) {
		is >> current_noesy_sampling_file_;
	}
	update_noesy_filter_files( current_noesy_sampling_file_, false );
	if ( evaluate_local() && current_noesy_sampling_file_ != "n/a" ) { //change filter restraints if we are local-evaluators
		topology_broker::ConstraintClaimerOP cst( new topology_broker::ConstraintClaimer( current_noesy_sampling_file_+".filter", "noesy_autoassign_cst" ) );
		cst->set_combine_ratio( bCombineNoesyCst_ ? 2 : 1 );
		cst->set_fullatom( true );
		cst->set_centroid( false );
		cst->set_skip_redundant( option[ iterative::skip_redundant_constraints ]() );
		cst->set_filter_weight( get_weight( "noesy_autoassign_cst" )/overall_cstfilter_weight_ );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new topology_broker::ConstraintEvaluatorWrapper( cst->tag(), cst ) ) ), cst->filter_weight()*overall_cstfilter_weight_ );
	}
	bCombineNoesyCst_ = stage() < STAGE2_RESAMPLING; //will be overwritten in next reassign NOESY... take guess until then...
}


void IterativeBase::save_status( std::ostream& os ) const {
	Parent::save_status( os );
	if ( hedge_archive_ ) hedge_archive_->save_to_file();

	os << "IterationStage: " << stage_;
	os << "   first_batch_this_stage: " << first_batch_this_stage_;
	os << "   first_fullatom_batch: " << first_fullatom_batch_;
	os << std::endl;
	os << "NOESY_CYCLE: " << noesy_assign_float_cycle_ << std::endl;
	os << "NOESY_FIRST_CST: " << first_noesy_cst_file_ << std::endl;
	os << "NOESY_FIRST_FA_CST: " << first_noesy_fa_cst_file_ << std::endl;
	os << "NOESY_CURRENT_CST: " << current_noesy_sampling_file_ << std::endl;

}


void IterativeBase::setup_default_evaluators() {
	Parent::setup_default_evaluators();
	add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::JumpNrEvaluator ) ) );
}


void IterativeBase::cluster() {
	///OBSOLET
	//  using namespace protocols::cluster;
	//  // using namespace basic::options::OptionKeys;
	//  using namespace basic::options::OptionKeys::cluster;
	//  using namespace basic::options;
	//  using namespace toolbox;

	//  //jump out if inactive
	//  if ( !option[ iterative::cluster ]() ) return;
	//  if ( decoys().size() < 50 ) return;
	//   if ( min_diversity_list_[ stage() ] == 0 ) return;

	//  mem_tr << "IterativeBase cluster-start" << std::endl;

	//  SilentStructs kept_decoys;
	//  toolbox::ClusterOptions cluster_opts( false /*don't change tags to c.XXX.NNN */ );
	//  cluster_opts.cluster_radius = min_diversity_list_[ stage() ];
	//  cluster_opts.keep_center = false; /* keep the lowest energy structures -- not interested in the most central structure */

	//  //read CA coords into DecoySetEvaluation
	//   DecoySetEvaluation CA_set;
	//  CA_set.push_back_CA_xyz_from_silent_file( decoys().size(),
	//   utility::DereferenceIterator< SilentStructs >( decoys().begin() ),
	//   utility::DereferenceIterator< SilentStructs >( decoys().end() ),
	//   false /*don't store plain energies */
	//  );

	//  //we have our special score... so need to gather this information
	//  utility::vector1< core::Real > all_energies;
	//  all_energies.reserve( decoys().size() );
	//  for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
	//   all_energies.push_back( select_score( *it ) );
	//  }
	//  CA_set.set_all_energies( all_energies );

	//  //now do the clustering
	//  toolbox::cluster_silent_structs( CA_set,
	//    utility::DereferenceIterator< SilentStructs >( decoys().begin() ),
	//   utility::DereferenceIterator< SilentStructs >( decoys().end() ),
	//   kept_decoys,
	//   cluster_opts
	//  );

	//  //how many were removed ?
	//  Size n_removed = decoys().size() - kept_decoys.size();
	//  tr.Info << "removed " << n_removed << " structures.   " << kept_decoys.size() << " structures remaining after clustering " << std::endl;
	//  decoys()=kept_decoys;
	//  count_removed_structures( n_removed );

	//  //finally...
	//  mem_tr << "IterativeBase cluster-end" << std::endl;
}


/// @detail before we can apply score-fxn we have to add extra data: RDC, NOES, (not supported yet: PCS, ... )
void IterativeBase::score( pose::Pose & pose ) const {
	//to speed up things we cache the RDC data in the archive
	if ( basic::options::option[ basic::options::OptionKeys::in::file::rdc ].user() ) {
		if ( !rdc_data_ ) rdc_data_ = core::scoring::ResidualDipolarCouplingOP( new core::scoring::ResidualDipolarCoupling );
		core::scoring::store_RDC_in_pose( rdc_data_, pose );
	}
	Parent::score( pose );
}


/// Helper functions


void
IterativeBase::test_broker_settings( jd2::archive::Batch const& batch ) {
	tr.Debug << "test broker settings...." << std::endl;
	OptionCollection vanilla_options( option );
	option.load_options_from_file( batch.flag_file() );
	try {
		topology_broker::TopologyBrokerOP topology_broker( new topology_broker::TopologyBroker() );
		topology_broker::add_cmdline_claims( *topology_broker );
		tr.Debug << "setting of broker::setup  " << std::endl;
		utility::vector1< std::string > files( option[ OptionKeys::broker::setup ]() );
		std::copy( files.begin(), files.end(), std::ostream_iterator<std::string>( tr.Debug, " "));
		tr.Debug << std::endl;
	} catch ( utility::excn::Exception &excn ) {  // clean up options and rethrow
		utility_exit_with_message( "[ERROR] problems with broker setup in "+batch.all_broker_files()+" aborting... ");
		tr.Error << "problems with broker setup in " << batch.all_broker_files() << " aborting... " << std::endl;
		// excn.show( tr.Error );
		option = vanilla_options;
		throw ( CREATE_EXCEPTION(EXCN_Archive,  batch.all_broker_files() + " contains errors: " + excn.msg() ) );
	}
	option = vanilla_options;
}


/// @detail load decoys into archive from -archive:input_pool or so
void IterativeBase::init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) {
	//make bogus batch that contains init-file

	//if non-local evaluation we need to add score_final to decoys --- switch temporarily to local evaluation
	bool b_old_eval_state( evaluate_local() );
	if ( !b_old_eval_state ) {
		tr.Debug << "switch to local evaluation for reading of initial pool" << std::endl;
		set_evaluate_local( true );//set this temporarily
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new simple_filters::ScoreEvaluator( "_final", scorefxn_non_const() ) ) ), 1.0 );
	}

	//read decoys and evaluate
	ArchiveBase::init_from_decoy_set( sfd );

	//switch back to non-local evaluation if applicable
	if ( !b_old_eval_state ) {
		remove_evaluation( "score_final" );
		set_weight( "score_final", 1.0 );
		set_evaluate_local( b_old_eval_state );
	}
}


//setup filter-cst from broker-setup:
void IterativeBase::setup_filter_cst( core::Real overall_weight ) {
	using namespace topology_broker;
	tr.Info << "setup filter-cst module in IterativeBase"<<std::endl;
	TopologyBrokerOP topology_broker( new TopologyBroker() );
	add_cmdline_claims( *topology_broker, false );
	tr.Trace << "topology_broker is initiailized with " << topology_broker->num_claimers() << " claimers "<< std::endl;
	core::Size ct( 1 );
	for ( auto it = topology_broker->begin();
			it != topology_broker->end(); ++it, ++ct ) {
		tr.Trace << "found claimer of type " << (*it)->type() << "trying to cast now..."<< std::endl;
		ConstraintClaimerCOP cst_claimer = utility::pointer::dynamic_pointer_cast< ConstraintClaimer const >( *it );
		if ( cst_claimer ) {
			tr.Info << "found cst-claimer with filter_name " << cst_claimer->filter_name()
				<< " filter_weight " << cst_claimer->filter_weight() << std::endl;
			if ( cst_claimer->filter_weight() < 0.01 ) continue;
			std::string name( cst_claimer->filter_name() );
			core::Real weight( cst_claimer->filter_weight() );
			if ( !name.size() ) {
				name = "filter_cst_"+ObjexxFCL::lead_zero_string_of( ct, 2 );
			} else {
				name = "filter_cst_"+name+"_"+ObjexxFCL::lead_zero_string_of( ct, 2 );
			}
			add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new ConstraintEvaluatorWrapper( name, cst_claimer ) ) ), weight*overall_weight );
		}//if cst_claimer
	} // for claimer
}

void IterativeBase::setup_autoNOE() {
	if ( !noesy_assign::NoesyModule::cmdline_options_activated() ) return;

	tr.Debug << "setup autoNOE module in IterativeBase"<<std::endl;
	super_quick_relax_of_centroids_=true;
	current_noesy_sampling_file_ = option[ iterative::initial_noe_auto_assign_csts ]();
	using utility::file::file_exists;
	if ( !file_exists( current_noesy_sampling_file_ ) ) {
		utility_exit_with_message( "couldn't find initial NOESY-autoassign file "+current_noesy_sampling_file_ );
	}

	//another cst set for filtering
	set_weight( "noesy_autoassign_cst", option[ iterative::cenpool_noesy_cst_weight ]()*overall_cstfilter_weight_ );

	//use restraints combination -- this might be switched off later when convergence is better
	bCombineNoesyCst_ = true;

	//add full-atom restraints as Evaluator (since a full-atom structure is sent back after relax)
	if ( evaluate_local() ) {
		tr.Debug << "cool we have local-scoring active"<<std::endl;
		topology_broker::ConstraintClaimerOP cst( new topology_broker::ConstraintClaimer( current_noesy_sampling_file_+".filter", "noesy_autoassign_cst" ) );
		cst->set_combine_ratio( bCombineNoesyCst_ ? 2 : 1 );
		cst->set_fullatom( true );
		cst->set_centroid( false );
		cst->set_filter_weight( option[ iterative::cenpool_noesy_cst_weight ]() );
		add_evaluation( evaluation::PoseEvaluatorCOP( evaluation::PoseEvaluatorOP( new topology_broker::ConstraintEvaluatorWrapper( cst->tag(), cst ) ) ), cst->filter_weight()*overall_cstfilter_weight_ );
	}

	//hash-string of the low30 decoy-tags to quickly determine if a new AutoNOE run should be done.
	std::ostringstream hash_string;
	hash_string << "NO_POOL " << std::endl;
	noesy_assign_hash_ = hasher( hash_string.str() );
}

} //abinitio
} //protocols
