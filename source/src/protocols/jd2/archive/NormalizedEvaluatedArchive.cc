// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu

// Unit headers
#include <protocols/jd2/archive/NormalizedEvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>
#include <protocols/jd2/archive/VarianceStatisticsArchive.hh>

//#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/constraints/ConstraintIO.hh>
//#include <core/scoring/constraints/util.hh>
//#include <core/scoring/Energies.hh>

//#include <core/pose/Pose.hh>
//#include <basic/datacache/BasicDataCache.hh>

//#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/util.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

//for DebugArchive

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

//for setup_default_evaluators
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.jd2.NormalizedEvaluatedArchive" );

//OPT_1GRP_KEY( Real, iterative, chainbreak_evaluator_exponent )
//OPT_1GRP_KEY( Boolean, iterative, simulate_bg4_cbtreatment )
OPT_2GRP_KEY( Boolean, iterative, normalize, activate )
OPT_2GRP_KEY( Boolean, iterative, normalize, extra_archive  )
OPT_2GRP_KEY( Real, iterative, normalize, keep_adding )
OPT_2GRP_KEY( Integer, iterative, normalize, num_struct)
OPT_2GRP_KEY( Integer, iterative, normalize, start )
OPT_2GRP_KEY( StringVector, iterative, normalize, force_zero )
OPT_2GRP_KEY( Boolean, iterative, normalize, lower_quartile )

std::string const SPECIAL_INITIAL_DECOY_PENALTY( "special_initial_decoy_penalty" );

bool protocols::jd2::archive::NormalizedEvaluatedArchive::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::jd2::archive::NormalizedEvaluatedArchive::register_options() {
	Parent::register_options();
	if ( !options_registered_ ) {
		NEW_OPT( iterative::normalize::activate, "score-variations are determined to normalize scores", false );
		NEW_OPT( iterative::normalize::extra_archive, "determine score variations from extra archive", false );
		NEW_OPT( iterative::normalize::keep_adding, "keep adding X percent of the incoming structures, throw out old ones randomly", 0.0 );
		NEW_OPT( iterative::normalize::num_struct, "number of structures in varaiance archive", 1000 );
		NEW_OPT( iterative::normalize::start, "do not normalize until X structures have been accumulated, if 0 this will be set to normalize:nstruct structures", 0 );
		NEW_OPT2( iterative::normalize::force_zero, "for scores whose name starts with XXX compute variance as 0...Q3 instead of Q1..Q3", "rdc", "filter_cst" );
		NEW_OPT( iterative::normalize::lower_quartile, "use Lo..Q1 instead of Q1..Q3 to determine range", false );
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

NormalizedEvaluatedArchive::~NormalizedEvaluatedArchive() {}


NormalizedEvaluatedArchive::NormalizedEvaluatedArchive() {
	runtime_assert( options_registered_ );
	score_variations_are_clean_ = false;
	init_from_options();
}

NormalizedEvaluatedArchive::NormalizedEvaluatedArchive( ArchiveManagerAP ptr )
: EvaluatedArchive( ptr )
{
	runtime_assert( options_registered_ );
	score_variations_are_clean_ = false;
	init_from_options();
}

void NormalizedEvaluatedArchive::init_from_options() {
	using namespace basic::options;
	insertion_prob_ = option[ OptionKeys::iterative::normalize::keep_adding ]();
	keep_adding_to_statistics_ = insertion_prob_ > 0.0;
	use_variance_archive_ = option[ OptionKeys::iterative::normalize::extra_archive ]();
	nstruct_for_statistics_ = option[ OptionKeys::iterative::normalize::num_struct ]();
	min_decoys_for_statistics_ = option[ OptionKeys::iterative::normalize::start ]();
	if ( !min_decoys_for_statistics_ ) min_decoys_for_statistics_ = nstruct_for_statistics_;
	activated_ = option[ OptionKeys::iterative::normalize::activate ]();
	positive_scores_ = option[ OptionKeys::iterative::normalize::force_zero ]();
	lower_quartile_ = option[ OptionKeys::iterative::normalize::lower_quartile ]();
	if ( option[ OptionKeys::run::test_cycles ] || option[ OptionKeys::run::dry_run ] ) {
		nstruct_for_statistics_ = 10;
		min_decoys_for_statistics_ = 10;
	}
}

void NormalizedEvaluatedArchive::initialize() {
	if ( use_variance_archive_ ) {
		variance_archive_ = VarianceStatisticsArchiveOP( new VarianceStatisticsArchive( name()+"_variance" ) );
		variance_archive_->set_nstruct( nstruct_for_statistics_ );
		variance_archive_->set_insertion_prob( insertion_prob_ );
		variance_archive_->initialize();
	}
}

//completely overwrites the EvaluatedArchive  version. Doesn't call Parent function.
bool NormalizedEvaluatedArchive::add_evaluated_structure(
	core::io::silent::SilentStructOP evaluated_decoy,
	core::io::silent::SilentStructOP alternative_decoy,
	Batch const& batch
) {
	bool added( Parent::add_evaluated_structure( evaluated_decoy, alternative_decoy, batch ) );
	bool added_variance( added );
	if ( variance_archive_ ) {
		added_variance = false;
		if ( keep_adding_to_statistics_ || total_proposed() < nstruct_for_statistics_ ) {
			added_variance = variance_archive_->add_evaluated_structure( evaluated_decoy, alternative_decoy, batch );
		}
	}
	tr.Info << "offered structure that was " << ( added ? "" : "not" ) << " relevant to the variances " << std::endl;
	score_variations_are_clean_ = score_variations_are_clean_ && !added_variance;
	return added;
}

void NormalizedEvaluatedArchive::save_to_file( std::string suffix ) {
	Parent::save_to_file( suffix );
	if ( variance_archive_ ) variance_archive_->save_to_file( suffix );
}

/// @detail restore archive and sort
bool NormalizedEvaluatedArchive::restore_from_file() {
	bool b_have_restored = Parent::restore_from_file();
	if ( use_variance_archive_ ) {
		runtime_assert( variance_archive_ != 0 );
		variance_archive_->restore_from_file();
	}
	score_variations_are_clean_ = false;
	return b_have_restored;
}
// --------------------------- end sort ------------------------------

/// @detail determine variations of the non-zero weighted (select_weight_) scores by taking the difference Q3-Q1 (upper / lower quartil)
bool NormalizedEvaluatedArchive::determine_score_variations() const {
	if ( score_variations_are_clean_ ) return false; //not changed
	score_variations_are_clean_ = true;
	SilentStructs const& my_decoys( variance_archive_ ? variance_archive_->decoys() : decoys() );
	core::Size ndecoys( my_decoys.size() );
	tr.Info << "determine score variations in NormalizedEvaluatedArchive " << name() << "... " << std::endl;
	tr.Info << "use " << ndecoys << " decoys from " << (variance_archive_ ? variance_archive_->name() : name() ) << std::endl;
	utility::vector1<core::Real> values;
	for ( WeightMap::const_iterator it = weights().begin(); it != weights().end(); ++it ) {
		if ( it->first == "special_initial_decoy_penalty" ) continue;
		if ( it->second > 0.01 ) {
			if ( ndecoys >= min_decoys_for_statistics_ && activated_ ) {

				std::string const& name( it->first );
				Size ct( 1 );
				core::Size half( ndecoys / 2 );
				core::Size lowQ( half / 2 );
				core::Size highQ( half + lowQ );
				if ( lower_quartile_ ) {
					highQ = lowQ;
					lowQ = 1 ;
				}
				//score_variations_.clear(); not really needed. should be faster without
				values.resize( ndecoys );
				for ( SilentStructs::const_iterator iss = my_decoys.begin(); iss != my_decoys.end(); ++iss, ++ct ) {
					if ( !(*iss)->has_energy( name ) ) {
						throw EXCN_Archive( "energy name "+name+" not found in returned decoys -- run with rescoring in archive to avoid this or fix your batches" );
					} // add weighted column-value to final score
					values[ct]=(*iss)->get_energy( name );
				}
				runtime_assert( lowQ > 0 && highQ < ndecoys );
				std::sort(values.begin(), values.end());
				if ( is_start_zero_score( it->first ) ) {
					score_variations_[ it->first ] = values[highQ];
					tr.Info << "score variation of " << score_variations_[ it->first ] << " for " << name
						<< " between 0 (forced)"
						<< " and "     << values[highQ] << " at " << highQ << std::endl;
				} else {
					score_variations_[ it->first ] = std::abs( values[highQ]-values[lowQ] );
					tr.Info << "score variation of " << score_variations_[ it->first ] << " for " << name
						<< " between " << values[lowQ]  << " at " << lowQ
						<< " and "     << values[highQ] << " at " << highQ << std::endl;
				}
			} else { //not enough decoys or not activated
				score_variations_[ it->first ] = 1.0;
			}
			//cutoff to avoid division by 0
			if ( score_variations_[ it->first ]< 1e-20 ) {
				score_variations_[ it->first ]= 1e-20;
			}
		} // if weight > 0.01
	} //for select_weights
	return true; //changed
} //determine_score_variations

bool NormalizedEvaluatedArchive::is_start_zero_score( std::string const& str ) const {
	for ( utility::vector1< std::string >::const_iterator it =positive_scores_.begin(); it != positive_scores_.end(); ++it ) {
		if ( str.substr( 0, it->size() ) == *it ) return true;
	}
	return false;
}

/// @detail rescore and sort archive
void NormalizedEvaluatedArchive::rescore() {
	Parent::rescore();
	score_variations_are_clean_ = false;
	if ( variance_archive_ ) {
		variance_archive_->set_evaluators( evaluators(), weights() );
		variance_archive_->set_scorefxn( scorefxn().clone() );
		variance_archive_->rescore();
	}
}

NormalizedEvaluatedArchive::WeightMap const& NormalizedEvaluatedArchive::score_variations() const {
	tr.Info << "ask for score_variations. They are " << (score_variations_are_clean_? "clean" : "not clean") << std::endl;
	if ( !score_variations_are_clean_ ) determine_score_variations();
	return score_variations_;
}

core::Real NormalizedEvaluatedArchive::score_variation( std::string const& column ) const {
	tr.Info << "ask for score_variations. They are " << (score_variations_are_clean_? "clean" : "not clean") << std::endl;
	if ( !score_variations_are_clean_ ) determine_score_variations();
	WeightMap::const_iterator iter = score_variations_.find( column );
	if ( iter != score_variations_.end() ) return iter->second;
	else return 1.0;
}

/* =================== end maintenance of evaluators and weights ====================== */


}//archive
}//jd2
}//protocols
