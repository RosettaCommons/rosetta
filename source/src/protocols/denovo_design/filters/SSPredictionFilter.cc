// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/SSPredictionFilter.cc
/// @brief Filter to determine agreement with SSPrediction for secondary structure prediction
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

// unit headers
#include <protocols/denovo_design/filters/SSPredictionFilter.hh>
#include <protocols/denovo_design/filters/SSPredictionFilterCreator.hh>

// protocol headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// project headers
#include <core/io/external/PsiPredInterface.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/ss_prediction/SS_predictor.hh>

// utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

static basic::Tracer TR( "protocols.denovo_design.filters.SSPredictionfilter" );

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
	psipred_interface_( /* NULL */ ),
	secstruct_(""),
	scratch_dir_("")
{}

// value constructor
SSPredictionFilter::SSPredictionFilter( core::Real const threshold,
	std::string const & cmd,
	std::string const & blueprint_filename,
	bool const use_probability,
	bool const mismatch_probability,
	bool const use_scratch_dir /* = false */ )
: protocols::filters::Filter( "SSPrediction" ),
	threshold_( threshold ),
	cmd_( cmd ),
	blueprint_( utility::pointer::make_shared< protocols::parser::BluePrint >( blueprint_filename ) ),
	use_probability_( use_probability ),
	mismatch_probability_ (mismatch_probability),
	use_confidence_( false ),
	use_svm_( true ),
	temp_( 0.6 ),
	ss_predictor_( /* NULL */ ),
	psipred_interface_( /* NULL */ ),
	secstruct_(""),
	scratch_dir_("")
{
	if ( use_scratch_dir ) set_scratch_dir();
	set_psipred_interface( utility::pointer::make_shared< PsiPredInterface >( cmd, scratch_dir_ ) );
}

SSPredictionFilter::~SSPredictionFilter() = default;

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
	return utility::pointer::make_shared< SSPredictionFilter >( *this );
}

protocols::filters::FilterOP
SSPredictionFilter::fresh_instance() const {
	return utility::pointer::make_shared< SSPredictionFilter >();
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
SSPredictionFilter::compute( core::pose::Pose const & pose ) const
{
	std::string wanted_ss;

	// if a blueprint is specified, use it. If not, use Dssp or pose metadata.
	if ( blueprint_ ) {
		wanted_ss = blueprint_->secstruct();
	} else if ( ! secstruct_.empty() ) {
		if ( secstruct_.length() != pose.size() ) {
			utility_exit_with_message("The secondary structure string provided to SSPredictionFilter is not the same length as the pose");
		}
		wanted_ss = secstruct_;
	} else if ( components::StructureDataFactory::get_instance()->observer_attached( pose ) ) {
		wanted_ss = components::StructureDataFactory::get_instance()->get_from_const_pose( pose ).ss();
	} else {
		core::scoring::dssp::Dssp dssp( pose );
		wanted_ss = dssp.get_dssp_secstruct();
	}

	// strip ligands from the ss string -- dssp now includes a character for ligand
	std::string pruned_ss( "" );
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( pose.residue(i).is_protein() ) {
			pruned_ss += wanted_ss[i-1];
		}
	}
	wanted_ss = pruned_ss;

	if ( use_svm_ ) {
		runtime_assert( ss_predictor_ != nullptr );
		std::string sequence;
		for ( core::Size i=1; i<=pose.size(); ++i ) {
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
			return core::Real(count) / sequence.size();
		}
	} else {
		runtime_assert( psipred_interface_ != nullptr );
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

void SSPredictionFilter::set_secstruct(std::string const secstruct){
	secstruct_ = secstruct ;
}

void
SSPredictionFilter::set_threshold( core::Real const threshold )
{
	threshold_ = threshold;
}

void
SSPredictionFilter::set_use_probability( bool const use_probability )
{
	use_probability_ = use_probability;
}

void
SSPredictionFilter::set_mismatch_probability( bool const mismatch_probability )
{
	mismatch_probability_ = mismatch_probability;
}

void
SSPredictionFilter::set_use_confidence( bool const use_confidence )
{
	use_confidence_ = use_confidence;
}

void
SSPredictionFilter::set_temperature( core::Real const temp )
{
	temp_ = temp;
}

void
SSPredictionFilter::set_use_svm( bool const use_svm )
{
	use_svm_ = use_svm;
}

void SSPredictionFilter::set_scratch_dir() {
	// There is a default value, but I'm 99% sure no one would ever want to use it
	if ( ! basic::options::option[ basic::options::OptionKeys::out::path::scratch ].user() ) {
		utility_exit_with_message("Set SSPrediction Filter use_scratch_dir=\"true\" but didn't set -out:path:scratch !");
	}
	scratch_dir_ = basic::options::option[ basic::options::OptionKeys::out::path::scratch ]();
}

//parse the rosetta scripts xml
void SSPredictionFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
){
	set_threshold( tag->getOption< core::Real >( "threshold", threshold_ ) );
	set_use_probability( tag->getOption< bool >( "use_probability", use_probability_ ) );
	set_mismatch_probability( tag->getOption< bool >( "mismatch_probability", mismatch_probability_ ) );
	set_use_confidence( tag->getOption< bool >( "use_confidence", use_confidence_ ) );
	set_temperature( tag->getOption< core::Real >( "temperature", temp_ ) );
	set_use_svm( tag->getOption< bool >( "use_svm", use_svm_ ) );

	set_cmd( tag->getOption< std::string >( "cmd", "" ) );
	if ( cmd_ == "" && ! use_svm_ ) {
		utility_exit_with_message("The SSPrediction Filter requires the psipred executable be set with the cmd option in the XML tag if SVM is not being used. Exiting now...");
	}

	if ( tag->getOption< bool >( "use_scratch_dir", false ) ) set_scratch_dir();

	// now that cmd is set, create the psipred interface
	if ( use_svm_ ) {
		ss_predictor_ = utility::pointer::make_shared< protocols::ss_prediction::SS_predictor >( "HLE" );
	} else {
		set_psipred_interface( utility::pointer::make_shared< PsiPredInterface >( cmd_, scratch_dir_ ) );
	}

	std::string blueprint_file = tag->getOption< std::string >( "blueprint", "" );
	if ( !blueprint_file.empty() ) {
		set_blueprint_file( blueprint_file );
	}

	std::string secstruct = tag->getOption< std::string >( "secstruct", "" );
	if ( !secstruct.empty() ) {
		if ( !secstruct.empty() and !blueprint_file.empty() ) {
			utility_exit_with_message("Please provide either a blueprint or secondary structure sting, not both.");
		}
		set_secstruct(secstruct);
	}
}

void
SSPredictionFilter::set_cmd( std::string const & cmd )
{
	cmd_ = cmd;
}

void
SSPredictionFilter::set_blueprint( protocols::parser::BluePrintOP blueprint )
{
	blueprint_ = blueprint;
}

void
SSPredictionFilter::set_blueprint_file( std::string const & blueprint_file )
{
	TR << "Dssp-derived secondary structure will be overridden by user-specified blueprint file." << std::endl;
	set_blueprint( utility::pointer::make_shared< protocols::parser::BluePrint >( blueprint_file ) );
	if ( ! blueprint_ ) {
		utility_exit_with_message("There was an error getting the blueprint file loaded.");
	}
}

std::string SSPredictionFilter::name() const {
	return class_name();
}

std::string SSPredictionFilter::class_name() {
	return "SSPrediction";
}

void SSPredictionFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"threshold", xsct_real,
		"If threshold is set and use_probability is true, the filter returns "
		"true only if the calculated value is less than this number. If use_"
		"probability is false, the filter returns true if the calculated "
		"value is greater than this number." )
		+ XMLSchemaAttribute(
		"use_probability", xsct_rosetta_bool,
		"If true, the probability information that psipred calculates will be "
		"used to determine the score. IF false, the filter will return "
		"the percentage of residues that match." )
		+ XMLSchemaAttribute(
		"mismatch_probability", xsct_rosetta_bool,
		"If true AND use_probability is true, the score is determined as the "
		"geometric average of the probability getting a WRONG secondary "
		"structure type at each position. Just as in regular use_probability, "
		"you want to minimize this score." )
		+ XMLSchemaAttribute(
		"use_confidence", xsct_rosetta_bool,
		"XRW TO DO" )
		+ XMLSchemaAttribute(
		"temperature", xsct_real,
		"XRW TO DO" )
		+ XMLSchemaAttribute(
		"use_svm", xsct_rosetta_bool,
		"If set, an SVM will be used to make secondary structure predictions "
		"instead of psipred. This requires downloading some database files. "
		"If false, the psipred executable specified by cmd will be used." )
		+ XMLSchemaAttribute(
		"cmd", xs_string,
		"Full path to runpsipred_single or runpsipred executable. "
		"Must be specified if use_svm=false" )
		+ XMLSchemaAttribute(
		"blueprint", xs_string,
		"If specified, the filter will take desired secondary structure from "
		"a blueprint file, rather from DSSP on the pose." )
		+ XMLSchemaAttribute(
		"secstruct", xs_string,
		"If specified, the filter will take desired secondary structure from "
		"this string, rather from DSSP on the pose. Format is: 'L' for loop, 'H' for helix"
		"and 'E' ")
		+ XMLSchemaAttribute(
		"use_scratch_dir", xsct_rosetta_bool,
		"If specified, use the directory specified with -out:path:scratch for the temporary "
		"files produced by PsiPred. The path must exist.");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Uses the sequence in the pose to generate secondary structure predictions "
		"for each position. Secondary structure predictions are then compared to "
		"the desired secondary structure to determine a score. If use_probability "
		"is true, the score returned is a value between 0 and 1, where 0 is complete "
		"secondary structure agreement, and 1 is no agreement. The following "
		"equation is used to determine the score: sum(i=1;N;e^(-p[i]/T)), where N is "
		"the number of residues, p[i] is the probability of correct secondary "
		"structure at position i, and T is a temperature factor set to 0.6 by "
		"default. If use_probability is false, the filter returns the fraction of "
		"residues that match the desired secondary structure as a number between 0 "
		"and 1. If use_probability is true AND mismatch_probability is true, the "
		"score is the geometric average of the probability of picking the WRONG "
		"secondary structure type at all residue positions. Minimizing this number "
		"will maximize the geometric average of the probability of picking the "
		"CORRECT secondary structure type at all residue positions. This option "
		"should be the most correct method for comparing two sequences to determine "
		"their expected fragment quality based on their predicted secondary "
		"structure probabilities at each residue.",
		attlist );
}

std::string SSPredictionFilterCreator::keyname() const {
	return SSPredictionFilter::class_name();
}

protocols::filters::FilterOP
SSPredictionFilterCreator::create_filter() const {
	return utility::pointer::make_shared< SSPredictionFilter >();
}

void SSPredictionFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SSPredictionFilter::provide_xml_schema( xsd );
}



//namespaces
}
}
}
