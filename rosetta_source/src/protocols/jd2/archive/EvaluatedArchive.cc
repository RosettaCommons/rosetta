// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu

// Unit headers
#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>

#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/constraints/ConstraintIO.hh>
//#include <core/scoring/constraints/util.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/util.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

//for DebugArchive
// AUTO-REMOVED #include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

//for setup_default_evaluators
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//#include <core/scoring/ResidualDipolarCoupling.hh>
//#include <core/scoring/constraints/ConstraintSet.hh>
#include <utility/vector1.hh>


static basic::Tracer tr("protocols.jd2.Archive");

//OPT_1GRP_KEY( Real, iterative, chainbreak_evaluator_exponent )
//OPT_1GRP_KEY( Boolean, iterative, simulate_bg4_cbtreatment )
OPT_1GRP_KEY( Boolean, iterative, evaluate_only_on_slaves )
OPT_1GRP_KEY( Real, iterative, penalize_initial_decoys )

std::string const SPECIAL_INITIAL_DECOY_PENALTY( "special_initial_decoy_penalty" );

bool protocols::jd2::archive::EvaluatedArchive::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::jd2::archive::EvaluatedArchive::register_options() {
	Parent::register_options();
	if ( !options_registered_ ) {
		//		NEW_OPT( iterative::chainbreak_evaluator_exponent, "exponent for nrjump_weighted_chainbreaks", 1.5 );
		//	NEW_OPT( iterative::simulate_bg4_cbtreatment, "this gives special cb weights", false );
		NEW_OPT( iterative::evaluate_only_on_slaves,"do not re-evaluate decoys when they are read into archvie (e.g. on BlueGene)", false );
		NEW_OPT( iterative::penalize_initial_decoys, "decoys read from input_pool have X extra score", 1000.0 );
		options_registered_ = true;
	}
}


namespace protocols {
namespace jd2 {
namespace archive {
// using namespace basic::options;
// using namespace basic::options::OptionKeys;
using namespace core;
using namespace core::io::silent;

EvaluatedArchive::~EvaluatedArchive() {}


EvaluatedArchive::EvaluatedArchive()
	: scorefxn_( NULL ),
		b_evaluate_incoming_decoys_( !basic::options::option[ basic::options::OptionKeys::iterative::evaluate_only_on_slaves ]() ) ///yields bottleneck on BG
{
	runtime_assert( options_registered_ );
	setup_default_evaluators();
}

EvaluatedArchive::EvaluatedArchive( ArchiveManagerAP ptr )
	: ArchiveBase( ptr ),
		scorefxn_( NULL ),
		b_evaluate_incoming_decoys_( !basic::options::option[ basic::options::OptionKeys::iterative::evaluate_only_on_slaves ]() ) ///yields bottleneck on BG
{
	runtime_assert( options_registered_ );
	setup_default_evaluators();
}


bool EvaluatedArchive::add_structure( core::io::silent::SilentStructOP from_batch ) {
	core::io::silent::SilentStructOP evaluated_decoy = evaluate_silent_struct( from_batch );
	return add_evaluated_structure( evaluated_decoy );
}

bool EvaluatedArchive::add_evaluated_structure( core::io::silent::SilentStructOP evaluated_decoy ) {
	//get score
	Real new_decoy_score ( select_score( evaluated_decoy ) );
	if ( numeric::isnan( new_decoy_score ) ) return false;

	//find position in sorted list to insert
	SilentStructs::iterator iss = decoys().begin();
	while ( iss != decoys().end() && new_decoy_score >= select_score( *iss ) ) ++iss;

	//if we are at the end this decoy has a worse score than all others
	if ( iss != decoys().end() || decoys().size() < nstruct() ) {
		add_structure_at_position( iss, evaluated_decoy );
		return true;
	}
	// here if score was not good enough to be added
	return false;
}

//overloaded to add some tracer output
void EvaluatedArchive::read_structures( core::io::silent::SilentFileData& sfd, Batch const& batch ) {
	tr.Info << "structures are scored with the following weights: " << std::endl;
	for ( WeightMap::const_iterator it=select_weights_.begin(); it != select_weights_.end(); ++it ) {
		std::string const& name( it->first );
		core::Real const& weight( it->second );
		tr.Info << name << " " << weight << std::endl;
	}
	Parent::read_structures( sfd, batch );
}

///@details evaluate decoy... if non-local evaluation just copy silent-struct
core::io::silent::SilentStructOP
EvaluatedArchive::evaluate_silent_struct( core::io::silent::SilentStructOP iss ) const {

	//non-local evalution ? just return input
	if ( !evaluate_local() ) {
		tr.Trace << "don't evaluate local: trust energies saved in " << iss->decoy_tag() << std::endl;
		return iss;
	}

	tr.Trace << "evaluate local for " << iss->decoy_tag() << std::endl;
	core::io::silent::SilentStructOP pss = iss->clone();

	//	pss->clear_energies(); //need to get rid of e.g., _archive_select_score_
		//note: clear_energies should not kill the comments (TAG_IN_FILE, SOURCE_FILE ) should still be present...
	//need to keep prefa_centroid_score...
	using namespace core::io::silent;
	utility::vector1< SilentEnergy > old_energies( pss->energies() );

	//make pose for scoring and evaluation purposes
	PROF_START( basic::ARCHIVE_FILL_POSE );
	pose::Pose pose;
	pss->fill_pose( pose ); //has to reread RDC file for each pose!
	pose.data().clear();
	pose.energies().clear();
	PROF_STOP( basic::ARCHIVE_FILL_POSE );

	//now evaluate the pose
	pss = evaluate_pose( pss, pose );

	//add all old energy terms that have not been re-computed, e.g., prefa_centroid_score
	for ( utility::vector1< SilentEnergy >::const_iterator it = old_energies.begin(); it != old_energies.end(); ++it ) {
		if ( it->name() != "_archive_select_score_" && !pss->has_energy( it->name() )) {
			pss->add_energy( it->name(), it->value(), it->weight() );
		}
	}
	return pss;

}

///@detail evaluate a pose ... store info in iss (and return the same for convenience )
core::io::silent::SilentStructOP
EvaluatedArchive::evaluate_pose( core::io::silent::SilentStructOP iss, core::pose::Pose& pose ) const {
	tr.Trace << "evaluate decoy " << iss->decoy_tag() << std::endl;

	//apply scorefxn
	PROF_START( basic::ARCHIVE_SCORE_POSE );
	score( pose );
	iss->energies_from_pose( pose );
	PROF_STOP( basic::ARCHIVE_SCORE_POSE );

	//apply evaluators
	PROF_START( basic::ARCHIVE_EVALUATORS );
 	for ( EvaluatorMap::const_iterator it=evaluators_.begin(), eit=evaluators_.end();
				it!=eit; ++it ) {
		it->second->apply( pose, iss->decoy_tag(), *iss );
	}
	PROF_STOP( basic::ARCHIVE_EVALUATORS );

	return iss;
}

///@detail compute the score of an evaluated decoy
Real EvaluatedArchive::select_score( SilentStructOP evaluated_decoy ) {

	///are cached scores clean ?
	if ( !scores_are_clean_ && evaluate_local() ) rescore();

	///is there a cached score for this decoy ?
	if ( evaluated_decoy->has_energy( "_archive_select_score_" ) ) {
		return evaluated_decoy->get_energy( "_archive_select_score_" );
	}

	///no cached score: compute score and cache it
	Real sum( 0.0 );

	///iterate over evaluators
	for ( WeightMap::const_iterator it=select_weights_.begin(); it != select_weights_.end(); ++it ) {
		std::string const& name( it->first );
		Real const& weight( it->second );
		if ( weight == 0.0 ) continue;

		if ( name == SPECIAL_INITIAL_DECOY_PENALTY ) {		///special evaluator name--- compute score-contribution for this one
			//crude hack to test things out ...
			if ( evaluated_decoy->get_comment( "source_file" ).find( "batch_000000/" ) != std::string::npos ) {
				sum += weight;
			}
		} else { //normal score-contribution --- column should be present
			if ( !evaluated_decoy->has_energy( name ) ) {
				throw EXCN_Archive( "energy name "+name+" not found in returned decoys -- run with rescoring in archive to avoid this or fix your batches" );
			} // add weighted column-value to final score
			sum += weight * evaluated_decoy->get_energy( name );
		}
	}

	//tracer output with full result
	tr.Trace << "evaluated select_score for " << evaluated_decoy->decoy_tag()
					 << " that was tagged as " << evaluated_decoy->get_comment( "tag_in_file")
					 << " : " << sum << " with " << select_weights_.size() << " evaluators" << std::endl;

	//add energy to cache
	evaluated_decoy->add_energy( "_archive_select_score_", sum );
	return sum;
}


///@detail restore archive and sort
bool EvaluatedArchive::restore_from_file() {
	bool b_have_restored = Parent::restore_from_file();
	sort();
	return b_have_restored;
}

// ---------------------- sort archive -------------------

class SortPredicate {
public:
	SortPredicate( EvaluatedArchive& arc ) : arc_( arc ) {};
	bool operator() (SilentStructOP const& pss1, SilentStructOP const& pss2 ) {
		return arc_.select_score( pss1 ) < arc_.select_score( pss2 );
	}
	EvaluatedArchive& arc_;
};

void EvaluatedArchive::sort() {
	decoys().sort( SortPredicate( *this ) );
}

// --------------------------- end sort ------------------------------


///@detail rescore and sort archive
void EvaluatedArchive::rescore() {
	tr.Debug << "rescore " << name() << " decoys " << std::endl;

	//rescore all decoys
	for ( SilentStructs::iterator iss = decoys().begin(); iss != decoys().end(); ++iss ) {
		*iss = evaluate_silent_struct( *iss );
	}
	scores_are_clean_ = true; //do this first, since SortPredicate asks for the select_score

	sort();
	tr.Debug << "...done rescoring and sorting " << std::endl;
// 	if ( tr.Trace.visible() ) {
// 		for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
// 			tr.Trace << select_score( *it ) << " " << (*it)->decoy_tag() << std::endl;
// 		}
// 	}
}


/* =================== maintenance of evaluators and weights ====================== */

void EvaluatedArchive::add_evaluation( evaluation::PoseEvaluatorOP eval, Real weight ) {
	for ( Size i=1; i<= eval->size(); i++ ) {
		std::string const& column( eval->name( i ) );
		select_weights_[ column ] = weight;
		evaluators_[ column ] = eval;
	}
	scores_are_clean_ = false;
}

void EvaluatedArchive::remove_evaluation( std::string const& name ) {
	std::string const& column( name );

	EvaluatorMap::iterator iter = evaluators_.find( column );
	if ( iter != evaluators_.end() ) 	evaluators_.erase( iter );

	WeightMap::iterator iter2 = select_weights_.find( column );
	if ( iter2 != select_weights_.end() ) select_weights_.erase( iter2 );

	scores_are_clean_ = false;
}

void EvaluatedArchive::set_weight( std::string const& column, core::Real weight ) {
	//	runtime_assert( has_evaluator( column ) ); or part of score!
	select_weights_[ column ] = weight;
}

core::Real EvaluatedArchive::get_weight( std::string const& column ) {
	//	runtime_assert( has_evaluator( column ) ); or part of score!
	if ( has_evaluator( column ) ) {
		return select_weights_[ column ];
	};
	return 0.0;
}

bool EvaluatedArchive::has_evaluator( std::string const& column ) {
	EvaluatorMap::const_iterator iter = evaluators_.find( column );
	return iter != evaluators_.end();
}

void EvaluatedArchive::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ) {
	scorefxn_ = scorefxn;
	scores_are_clean_ = false;
}

core::scoring::ScoreFunction const &
EvaluatedArchive::scorefxn() const {
	runtime_assert( scorefxn_ );
	return *scorefxn_;
}

core::scoring::ScoreFunctionOP
EvaluatedArchive::scorefxn_non_const() {
	runtime_assert( scorefxn_ );
	return scorefxn_;
}

void EvaluatedArchive::setup_default_evaluators() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//using namespace scoring::constraints;
// 	if ( evaluate_local() && option[ constraints::cst_file ].user() ) {
// 		std::string filename( option[ constraints::cst_file ]()[ 1 ] );
// 		evaluation::PoseEvaluatorOP ev_cst ( new evaluation::ConstraintEvaluator( "cmdline", filename ) );
// 		add_evaluation( ev_cst, option[ constraints::cst_weight ] );
// 		//The ConstraintEvaluator creates two columns: cmdline_cst and cmdline_viol...
// 		set_weight( "cmdline_viol", 0.0 );
// 	}
	set_weight( "score", 1.0 );
	if ( option[ OptionKeys::iterative::penalize_initial_decoys ]() > 0.0 ) {
		set_weight( SPECIAL_INITIAL_DECOY_PENALTY , option[ OptionKeys::iterative::penalize_initial_decoys ]() );
	}

	evaluation::MetaPoseEvaluator cmdline_evals;
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(cmdline_evals);
	for ( evaluation::MetaPoseEvaluator::EvaluatorList::const_iterator it = cmdline_evals.evaluators().begin();
				it != cmdline_evals.evaluators().end(); ++it ) {
 		add_evaluation( *it );
	}
}

/* =================== end maintenance of evaluators and weights ====================== */


}//archive
}//jd2
}//protocols
