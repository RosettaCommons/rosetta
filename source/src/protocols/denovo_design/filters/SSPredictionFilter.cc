// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/SSPredictionFilter.cc
/// @brief Filter to determine agreement with SSPrediction for secondary structure prediction
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

// unit headers
#include <protocols/denovo_design/filters/SSPredictionFilter.hh>
#include <protocols/denovo_design/filters/SSPredictionFilterCreator.hh>

// package headers

#include <protocols/denovo_design/components/StructureData.hh>

// project headers
#include <core/io/external/PsiPredInterface.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/ss_prediction/SS_predictor.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <fstream>
#include <cstdlib>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.filters.SSPredictionfilter" );

namespace protocols {
namespace denovo_design {
namespace filters {
		
using namespace core::io::external;
// general constructor
SSPredictionFilter::SSPredictionFilter()
: protocols::filters::Filter( "SSPrediction" ),
	threshold_( 0 ),
	cmd_( "" ),
	blueprint_( /* NULL */ ),
	use_probability_( false ),
	mismatch_probability_(false),
	use_confidence_( false ),
	use_svm_( true ),
	temp_( 0.6 ),
	ss_predictor_( /* NULL */ ),
	psipred_interface_( /* NULL */ )
{}

// value constructor
SSPredictionFilter::SSPredictionFilter( core::Real const threshold,
	std::string const & cmd,
	std::string const & blueprint_filename,
	bool const use_probability,
	bool const mismatch_probability )
: protocols::filters::Filter( "SSPrediction" ),
	threshold_( threshold ),
	cmd_( cmd ),
	blueprint_( protocols::jd2::parser::BluePrintOP( new protocols::jd2::parser::BluePrint( blueprint_filename ) ) ),
	use_probability_( use_probability ),
	mismatch_probability_ (mismatch_probability),
	use_confidence_( false ),
	use_svm_( true ),
	temp_( 0.6 ),
	ss_predictor_( /* NULL */ ),
	psipred_interface_( PsiPredInterfaceOP( new PsiPredInterface( cmd ) ) )
{}

SSPredictionFilter::~SSPredictionFilter()
{}

// virtual functions that need to be overridden
bool
SSPredictionFilter::apply( core::pose::Pose const & pose ) const {
	core::Real value( compute( pose ) );
	if ( value >= threshold_ && ! use_probability_ ) {
		return true;
	} else if ( value <= threshold_ && use_probability_ ) {
		return true;
	} else {
		return false;
	}
}

protocols::filters::FilterOP
SSPredictionFilter::clone() const {
	return protocols::filters::FilterOP( new SSPredictionFilter( *this ) );
}

protocols::filters::FilterOP
SSPredictionFilter::fresh_instance() const {
	return protocols::filters::FilterOP( new SSPredictionFilter() );
}

void SSPredictionFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out << "SSPredictionFilter returns " << compute( pose ) << std::endl;
}

core::Real
SSPredictionFilter::report_sm( core::pose::Pose const & pose ) const {
	return compute( pose );
}

/// @brief computes the weighted boltzmann sum of the passed vector
core::Real
SSPredictionFilter::compute_boltz_sum( utility::vector1< core::Real > const & probabilities ) const
{
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=probabilities.size(); ++i ) {
		sum += exp( -probabilities[i] / temp_ );
	}
	TR << "Boltzmann sum is " << sum << " -- Final value=" << sum/probabilities.size() << std::endl;
	return sum/probabilities.size();
}

/// @brief computes the overall probability P as the product of each residue probability p.
/// Then returns the Nth root of P for N residues.
core::Real
SSPredictionFilter::compute_mismatch_prob( utility::vector1< core::Real > const & probabilities ) const
{
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=probabilities.size(); ++i ) {
		sum += log( probabilities[i]);
		//TR << i << "  " << probabilities[i] << "  " <<log(probabilities[i]) << std::endl;
	}

	TR << "Probability of correct secondary structure is " << exp(sum/probabilities.size()) << "^(n_residues) -- Final value=" << 1-exp(sum/probabilities.size()) << std::endl;
	return 1-exp(sum/probabilities.size());
}

// @brief Calculates the score. if probabilities are used, it compute the boltzmann sum, which will be a number between 0 and 1, where 0 is the best, and 1 is the worst.
core::Real
SSPredictionFilter::compute( core::pose::Pose const & pose ) const {
	std::string wanted_ss;

	// if a blueprint is specified, use it. If not, use Dssp or pose metadata.
	if ( blueprint_ ) {
		wanted_ss = blueprint_->secstruct();
	} else if ( components::StructureData::has_cached_string( pose ) ) {
		components::StructureDataOP sd = components::StructureData::create_from_pose( pose, "" );
		debug_assert( sd );
		wanted_ss = sd->ss();
	} else {
		core::scoring::dssp::Dssp dssp( pose );
		wanted_ss = dssp.get_dssp_secstruct();
	}

	// strip ligands from the ss string -- dssp now includes a character for ligand
	runtime_assert( pose.total_residue() == wanted_ss.size() );
	std::string pruned_ss( "" );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			pruned_ss += wanted_ss[i-1];
		}
	}
	wanted_ss = pruned_ss;

	if ( use_svm_ ) {
		runtime_assert( ss_predictor_ != 0 );
		std::string sequence;
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue( i ).is_protein() ) sequence += pose.residue( i ).name1();
		}
		runtime_assert( sequence.size() == wanted_ss.size() );
		utility::vector1< utility::vector1< core::Real > > ss_pred( ss_predictor_->predict_ss( sequence ) );
		utility::vector1< core::Real > probabilities;
		if ( use_probability_ ) {
			for ( core::Size i=1; i<=wanted_ss.size(); ++i ) {
				probabilities.push_back( protocols::ss_prediction::get_prob( wanted_ss[i-1], ss_pred[i] ) );
			}
			return compute_boltz_sum( probabilities );
		} else {
			core::Size count( 0 );
			// check to see if the prediction matches the desired SS at each prediction
			for ( core::Size i=1; i<=wanted_ss.size(); ++i ) {
				if ( protocols::ss_prediction::get_label( ss_pred[i] ) == wanted_ss[i-1] ) {
					++count;
				}
			}
			return count / sequence.size();
		}
	} else {
		runtime_assert( psipred_interface_ != 0 );
		PsiPredResult const psipred_result = psipred_interface_->run_psipred( pose, wanted_ss );
		if ( use_probability_ ) {
			TR << "Blueprint SS = " << wanted_ss << std::endl;
			if ( use_confidence_ ) {
				runtime_assert( wanted_ss.size() == psipred_result.psipred2_confidence.size() );
				return compute_boltz_sum( generate_prob( psipred_result, wanted_ss ) );
			} else {
				if ( mismatch_probability_ ) {
					runtime_assert( wanted_ss.size() == psipred_result.psipred2_confidence.size() );
					return compute_mismatch_prob( psipred_result.psipred_prob );
				} else {
					TR.Debug << "Wanted SS: " << wanted_ss.size() << " Result SS: " << psipred_result.psipred_prob.size() << std::endl;
					runtime_assert( wanted_ss.size() == psipred_result.psipred_prob.size() );
					return compute_boltz_sum( psipred_result.psipred_prob );
				}
			}
		} else {
			core::Real const count( psipred_result.nres );
			TR << count << " residues had correct secondary structure." << std::endl;
			return count / psipred_result.psipred_prob.size();
		}
	}
}

//parse the rosetta scripts xml
void SSPredictionFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );
	use_probability_ = tag->getOption< bool >( "use_probability", use_probability_ );
	mismatch_probability_ = tag->getOption< bool >( "mismatch_probability", mismatch_probability_ );
	use_confidence_ = tag->getOption< bool >( "use_confidence", use_confidence_ );
	temp_ = tag->getOption< core::Real >( "temperature", temp_ );
	use_svm_ = tag->getOption< bool >( "use_svm", use_svm_ );

	cmd_ = tag->getOption< std::string >( "cmd", "" );
	if ( cmd_ == "" && ! use_svm_ ) {
		utility_exit_with_message("The SSPrediction Filter requires the psipred executable be set with the cmd option in the XML tag if SVM is not being used. Exiting now...");
	}

	// now that cmd is set, create the psipred interface
	if ( use_svm_ ) {
		ss_predictor_ = protocols::ss_prediction::SS_predictorOP( new protocols::ss_prediction::SS_predictor( "HLE" ) );
	} else {
		psipred_interface_ = PsiPredInterfaceOP( new PsiPredInterface( cmd_ ) );
	}

	std::string blueprint_file = tag->getOption< std::string >( "blueprint", "" );
	if ( blueprint_file != "" ) {
		TR << "Dssp-derived secondary structure will be overridden by user-specified blueprint file." << std::endl;
		blueprint_ = protocols::jd2::parser::BluePrintOP( new protocols::jd2::parser::BluePrint( blueprint_file ) );
		if ( ! blueprint_ ) {
			utility_exit_with_message("There was an error getting the blueprint file loaded.");
		}
	}
}

protocols::filters::FilterOP
SSPredictionFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SSPredictionFilter() );
}

std::string
SSPredictionFilterCreator::keyname() const {
	return "SSPrediction";
}


//namespaces
}
}
}
