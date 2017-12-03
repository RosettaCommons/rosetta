// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/RotamerBoltzmannWeight2.cc
/// @brief Next-generation RotamerBoltzmannWeight filter
/// @author Tom Linsky (tlinsky@uw.edu)

#include <protocols/simple_filters/RotamerBoltzmannWeight2.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight2Creator.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/AlaScan.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/toolbox/EnergyLandscapeEvaluator.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

static basic::Tracer TR( "protocols.simple_filters.RotamerBoltzmannWeight2" );

namespace protocols {
namespace simple_filters {

core::Real const RotamerBoltzmannWeight2::default_temperature_ = 0.8;
core::Real const RotamerBoltzmannWeight2::default_lambda_ = 0.5;

RotamerBoltzmannWeight2::RotamerBoltzmannWeight2():
	protocols::filters::Filter( "RotamerBoltzmannWeight2" ),
	score_type_( MEAN_PROBABILITY ),
	calculator_id_( new_calculator_id() ),
	scorefxn_( core::scoring::get_score_function() ),
	calculator_( new protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator( core::scoring::get_score_function(), default_temperature_ ) )
{
	calculator_->set_lazy( false );
	TR << "Created new calculator id: " << calculator_id_ << std::endl;
}

RotamerBoltzmannWeight2::~RotamerBoltzmannWeight2()
= default;

RotamerBoltzmannWeight2::RotamerBoltzmannWeight2( RotamerBoltzmannWeight2 const & rval ):
	protocols::filters::Filter( rval ),
	score_type_( rval.score_type_ ),
	calculator_id_( new_calculator_id() ),
	calculator_( new protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculator( *rval.calculator_ ) )
{
	TR << "Created new calculator id: " << calculator_id_ << std::endl;
}

void
RotamerBoltzmannWeight2::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	core::Real temperature = default_temperature_ ;
	if ( tag->hasOption( "temperature" ) ) {
		temperature = tag->getOption< core::Real >( "temperature" ) ;
	}
	core::Real lambda = default_lambda_ ;
	if ( tag->hasOption( "lambda" ) ) {
		lambda = tag->getOption< core::Real >( "lambda" ) ;
	}
	std::string const probability_type = tag->getOption< std::string >( "probability_type", "" );
	if ( !probability_type.empty() ) {
		if ( probability_type == "BOLTZMANN_SUM" ) {
			set_energy_landscape_evaluator( protocols::toolbox::EnergyLandscapeEvaluatorOP(
				new protocols::toolbox::RotamerBoltzmannWeightEvaluator( temperature, true ) ) );
		} else if ( probability_type == "PNEAR" ) {
			set_energy_landscape_evaluator( protocols::toolbox::EnergyLandscapeEvaluatorOP(
				new protocols::toolbox::MulliganPNearEvaluator( temperature, lambda ) ) );
		} else {
			std::stringstream msg;
			msg << "RotamerBoltzmannWeight2::parse_my_tag(): invalid probability_type specified: " << probability_type << std::endl;
			msg << "Valid score types are: BOLTZMANN_SUM, PNEAR" << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
	}

	std::string const score_type = tag->getOption< std::string >( "score_type", "" );
	if ( !score_type.empty() ) {
		if ( score_type == "MEAN_PROBABILITY" ) {
			set_score_type( MEAN_PROBABILITY );
		} else if ( score_type == "MAX_PROBABILITY" ) {
			set_score_type( MAX_PROBABILITY );
		} else if ( score_type == "MODIFIED_DDG" ) {
			set_score_type( MODIFIED_DDG );
		} else if ( score_type == "PROBABILITY" ) {
			set_score_type( PROBABILITY );
		} else if ( score_type == "MIN_PROBABILITY" ) {
			set_score_type( MIN_PROBABILITY );
		} else {
			std::stringstream msg;
			msg << "RotamerBoltzmannWeight2::parse_my_tag(): invalid score_type specified: " << score_type << std::endl;
			msg << "Valid score types are: MEAN_PROBABILITY, MAX_PROBABILITY, MODIFIED_DDG, MIN_PROBABILITY, PROBABILITY" << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
	}

	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, data );
	if ( selector ) {
		set_residue_selector( selector );
	}

	core::scoring::ScoreFunctionOP scorefxn = protocols::rosetta_scripts::parse_score_function( tag, data );
	if ( scorefxn ) {
		set_scorefxn( scorefxn );
	}
}

protocols::filters::FilterOP
RotamerBoltzmannWeight2::clone() const
{
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight2( *this ) );
}


protocols::filters::FilterOP
RotamerBoltzmannWeight2::fresh_instance() const
{
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight2 );
}

std::string
RotamerBoltzmannWeight2::get_name() const
{
	return RotamerBoltzmannWeight2::class_name();
}

bool
RotamerBoltzmannWeight2::apply( core::pose::Pose const & ) const
{
	return true;
}

core::Real
RotamerBoltzmannWeight2::report_sm( core::pose::Pose const & pose ) const
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	if ( !calculator_ ) {
		utility_exit_with_message( "RotamerBoltzmannWeight2::report_sm(): calculator is not set!" );
	}

	// create and register calculator if it's not already registered
	register_calculator();

	// retrieve calculator value
	basic::MetricValue< RotamerProbabilities > probabilities;
	pose.metric( calculator_id_, "probabilities", probabilities );

	return compute_score( pose, probabilities.value() );
}

void
RotamerBoltzmannWeight2::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	// create and register calculator if it's not already registered
	register_calculator();

	basic::MetricValue< RotamerProbabilities > probabilities_metricvalue;
	pose.metric( calculator_id_, "probabilities", probabilities_metricvalue );
	RotamerProbabilities const & probabilities = probabilities_metricvalue.value();

	out << "Residue\tProbability\tEnergy_Reduction" << std::endl;
	for ( auto const & probabilitie : probabilities ) {
		core::Real const energy_reduction = compute_modified_residue_energy( probabilitie.second );
		out << pose.residue(probabilitie.first).name() << probabilitie.first << '\t'
			<< probabilitie.second << '\t' << energy_reduction << std::endl;
	}

	out << "Final score: " << compute_score( pose, probabilities ) << std::endl;
}

core::Real
compute_mean_probability( protocols::toolbox::pose_metric_calculators::RotamerProbabilities const & probs )
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	core::Real sum = 0.0;
	for ( auto const & prob : probs ) {
		sum += prob.second;
	}
	return sum / static_cast< core::Real >( probs.size() );
}

core::Real
compute_probability( protocols::toolbox::pose_metric_calculators::RotamerProbabilities const & probs )
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	core::Real total_prob = 1.0;
	for ( auto const & prob : probs ) {
		total_prob = total_prob * prob.second;
	}
	return total_prob ;
}

core::Real
compute_max_probability( protocols::toolbox::pose_metric_calculators::RotamerProbabilities const & probs )
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	if ( probs.empty() ) return 0.0;

	core::Real max = probs.begin()->second;
	for ( auto rp=(++probs.begin()); rp!=probs.end(); ++rp ) {
		if ( rp->second > max ) max = rp->second;
	}
	return max;
}

core::Real
compute_min_probability( protocols::toolbox::pose_metric_calculators::RotamerProbabilities const & probs )
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	if ( probs.empty() ) return 0.0;

	core::Real min = probs.begin()->second;
	for ( auto rp=(++probs.begin()); rp!=probs.end(); ++rp ) {
		if ( rp->second < min ) min = rp->second;
	}
	return min;
}

core::Real
RotamerBoltzmannWeight2::compute_modified_residue_energy(
	toolbox::pose_metric_calculators::RotamerProbability const & probability ) const
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbability;
	static core::Real const energy_reduction_factor = 0.5;
	static RotamerProbability const threshold_prob = 1.0;
	if ( probability <= threshold_prob ) {
		return -1 * default_temperature_ * log( probability ) * energy_reduction_factor;
	} else {
		return 0.0;
	}
}

core::Real
RotamerBoltzmannWeight2::compute_modified_ddg(
	core::pose::Pose const & pose,
	toolbox::pose_metric_calculators::RotamerProbabilities const & probs ) const
{
	using protocols::toolbox::pose_metric_calculators::RotamerProbability;
	using protocols::toolbox::pose_metric_calculators::RotamerProbabilities;

	core::Real const ddg_in = compute_ddg( pose );
	core::Real modified_ddG = ddg_in;
	for ( auto const & prob : probs ) {
		core::Size const res( prob.first );
		core::Real const energy = compute_modified_residue_energy( prob.second );
		modified_ddG += energy;
		TR.Debug << pose.residue( res ).name3() << pose.pdb_info()->number( res ) << '\t'
			<< '\t' << prob.second << '\t' << energy << std::endl;
	}
	TR.Debug << "ddG before, after modification: " << ddg_in << ", " << modified_ddG << std::endl;
	return modified_ddG;
}

core::Real
RotamerBoltzmannWeight2::compute_ddg( core::pose::Pose const & pose ) const
{
	static int const jump = 1;
	static core::Size const nrepeats = 3;
	protocols::simple_filters::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn_->clone(), jump, nrepeats );
	ddg_filter.repack( true );

	return ddg_filter.compute( pose );
}

core::Real
RotamerBoltzmannWeight2::compute_score(
	core::pose::Pose const & pose,
	protocols::toolbox::pose_metric_calculators::RotamerProbabilities const & probabilities ) const
{
	if ( score_type_ == MEAN_PROBABILITY ) {
		return -compute_mean_probability( probabilities );
	} else if ( score_type_ == MAX_PROBABILITY ) {
		return -compute_max_probability( probabilities );
	} else if ( score_type_ == MODIFIED_DDG ) {
		return compute_modified_ddg( pose, probabilities );
	} else if ( score_type_ == MIN_PROBABILITY ) {
		return -compute_min_probability( probabilities );
	} else if ( score_type_ == PROBABILITY ) {
		return -compute_probability( probabilities );
	} else {
		utility_exit_with_message( "Unknown score type" );
	}

	return static_cast< core::Real >(probabilities.size());
}

void
RotamerBoltzmannWeight2::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	calculator_->set_residue_selector( selector );
	unregister_calculator();
}

void
RotamerBoltzmannWeight2::set_energy_landscape_evaluator( protocols::toolbox::EnergyLandscapeEvaluatorCOP evaluator )
{
	calculator_->set_energy_landscape_evaluator( evaluator );
	unregister_calculator();
}

void
RotamerBoltzmannWeight2::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	calculator_->set_scorefxn( scorefxn->clone() );
	scorefxn_ = scorefxn->clone();
	unregister_calculator();
}

void
RotamerBoltzmannWeight2::set_score_type( ScoreType const scoretype )
{
	score_type_ = scoretype;
}

void
RotamerBoltzmannWeight2::register_calculator() const
{
	// register calculator if it doesn't already exist
	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calculator_id_ ) ) {
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( calculator_id_, calculator_ );
		TR << "Registered " << calculator_id_ << std::endl;
	}
}

void
RotamerBoltzmannWeight2::unregister_calculator()
{
	// un-register calculator if it exists
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calculator_id_ ) ) {
		core::pose::metrics::CalculatorFactory::Instance().remove_calculator( calculator_id_ );
		TR << "Unregistered " << calculator_id_ << std::endl;
		calculator_id_ = new_calculator_id();
	}
}

std::string
RotamerBoltzmannWeight2::new_calculator_id()
{
	core::Size const id = IdManager< core::Size >::get_instance()->register_new_id();
	std::stringstream calc_id;
	calc_id << RotamerBoltzmannWeight2::class_name() << "_" << id;
	return calc_id.str();
}

std::string const &
RotamerBoltzmannWeight2::calculator_id() const
{
	return calculator_id_;
}

/////////////// Creator ///////////////

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP RotamerBoltzmannWeight2Creator::create_filter() const
// XRW TEMP {
// XRW TEMP  return protocols::filters::FilterOP( new RotamerBoltzmannWeight2 );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RotamerBoltzmannWeight2Creator::keyname() const
// XRW TEMP {
// XRW TEMP  return RotamerBoltzmannWeight2::class_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RotamerBoltzmannWeight2::class_name()
// XRW TEMP {
// XRW TEMP  return "RotamerBoltzmannWeight2";
// XRW TEMP }

std::string RotamerBoltzmannWeight2::name() const {
	return class_name();
}

std::string RotamerBoltzmannWeight2::class_name() {
	return "RotamerBoltzmannWeight2";
}

void RotamerBoltzmannWeight2::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "probability_type", xs_string, "How probabilities are calculated: BOLTZMANN_SUM -- A simple boltzmann sum defined by 1/SUM( e^( -(score - score_0)/temp ) ), where score_0 is the score with the rotamer in the input pose. The score of the input pose is included in this sum, and therefore it can range from (0, 1]. All rotamers not in the input pose are essentially considered equally \"far\" from the input. PNEAR -- Probability of being near the input state, using a Gaussian to define nearness in sidechain heavy atom RMSD to the input state as proposed by Vikram Mulligan in the Baker Lab. The lambda parameter defines the width of the Gaussian. For a lambda of 0.5, a rotamer with sidechain RMSD of 0.5 will be scored as as half \"near\" and half \"far\" from the input state. PNEAR is defined as NUMERATOR/DENOMINATOR, where NUMERATOR = SUM( e^(-rms^2/lambda^2)*e(-(score - score_0)/temp) ), and DENOMINATOR = SUM( e^( -(score - score_0)/temp ).","BOLTZMANN_SUM")
		+ XMLSchemaAttribute::attribute_w_default( "score_type", xs_string, " The method used to combine the rotamer probabilities into the single number reported by the filter. Available score types are: MEAN_PROBABILITY -- Returns the mean of all computed probabilities. MAX_PROBABILITY -- Returns the maximum of all computed probabilities. MIN_PROBABILITY -- Returns the minimum of all computed probabilities. PROBABILITY -- Returns the product of all computed probabilities. MODIFIED_DDG -- Returns a ddG value weighted by the rotamer probabilities as described in Fleishman et al. (2011) Protein Sci. 20:753.", "MODIFIED_DDG")
		+ XMLSchemaAttribute::attribute_w_default( "lambda", xs_string, "The \"lambda\" value to be used in computing \"nearness\" to the input state during rotamer probability calculation. This is only used if \"PNEAR\" is the probability method.", "0.5")
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xs_string, "The temperature to be used in computing rotamer probabilities", "0.8");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "XRW TO DO");

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string RotamerBoltzmannWeight2Creator::keyname() const {
	return RotamerBoltzmannWeight2::class_name();
}

protocols::filters::FilterOP
RotamerBoltzmannWeight2Creator::create_filter() const {
	return protocols::filters::FilterOP( new RotamerBoltzmannWeight2 );
}

void RotamerBoltzmannWeight2Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotamerBoltzmannWeight2::provide_xml_schema( xsd );
}


//////////////// Calculator id management /////////////////////

template< class T >
IdManager< T >::IdManager():
	utility::SingletonBase< IdManager< T > >(),
	used_ids_()
{}

template< class T >
T const &
IdManager< T >::register_new_id()
{
	T test_id = 0;
	while ( true ) {
		if ( ! id_exists( test_id ) ) {
			return *(used_ids_.insert( test_id ).first);
		}
		++test_id;
	}
	utility_exit_with_message("Cannot register a new id!");
	return *(used_ids_.end());
}

template< class T >
bool
IdManager< T >::id_exists( T const & id ) const
{
	return ( used_ids_.find( id ) != used_ids_.end() );
}

template< class T >
void
IdManager< T >::register_id( T const & id )
{
	if ( id_exists( id ) ) {
		std::stringstream msg;
		msg << "IdManager::register_id(): Attempting to register an id (" << id << ") that already exists!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	used_ids_.insert( id );
}

template< class T >
void
IdManager< T >::unregister_id( T const & id )
{
	std::stringstream msg;
	msg << "IdManager::register_id(): Attempting to unregister an id (" << id << ") that doesn't exist!" << std::endl;
	utility_exit_with_message( msg.str() );
}

} //protocols
} //simple_filters

